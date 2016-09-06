#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -o Intergenic.sh.stdout
#PBS -e Intergenic.sh.stderr
#PBS -V
#PBS -N Intergenic-split-sense

cd /gs/project/eav-760-aa/RNA-CKII-Spt6

rm -rf output/intergenic
mkdir -p output/intergenic

# Load bedtools.
module add mugqic/samtools/1.3
module add mugqic/bedtools/2.22.1

# UCSC-obtained data has "chr" in front of chromosome number. Remvoe that.
sed -e 's/^chr//' input/UCSC_SGD_Genes.bed > output/intergenic/UCSC_SGD_Genes.nochr.bed

# Define a list of all intergenic regions
bedtools subtract -a input/All.bed -b output/intergenic/UCSC_SGD_Genes.nochr.bed > output/intergenic/intergenic.bed

# Split and annotate all intergenic regions based on their proximity to the closest gene.
perl scripts/annotate_split.pl output/intergenic/UCSC_SGD_Genes.nochr.bed 1000 < output/intergenic/intergenic.bed > output/intergenic/intergenic.split.bed


for type in promoter downstream
do
    grep -e "$type.*+$" output/intergenic/intergenic.split.bed > $type.plus.bed
    grep -e "$type.*-$" output/intergenic/intergenic.split.bed > $type.minus.bed
done

grep -e "intergenic" output/intergenic/intergenic.split.bed > intergenic.bed

# Split the bam files
for lib in CK2 WT Spt6-SA
do
    # Get the bam file's name.
    BAMFILE=output/bwa/$lib.sorted.bam
     
    for type in promoter downstream
    do
        for direction in sense antisense
        do
            plus_samtools_arg="-f"
            minus_samtools_arg="-F"
            if [ $direction = "antisense" ]
            then
                plus_samtools_arg="-F"
                minus_samtools_arg="-f"        
            fi
            
            prefix=output/intergenic/$lib.$type.$direction
            
            cat $lib.$type.$direction.generate_bam.sh <<EOF
samtools view -b -h $plus_samtools_arg 16 -L $type.plus.bed $BAMFILE > $prefix.1.bam
samtools view -b -h $minus_samtools_arg 16 -L $type.minus.bed $BAMFILE > $prefix.2.bam
samtools merge $prefix.bam $prefix.1.bam $prefix.2.bam
samtools sort -o $prefix.sorted.bam $prefix.bam
samtools index $prefix.sorted.bam
rm $prefix.1.bam $prefix.2.bam $prefix.bam        
EOF

            
         done
    done
    
    cat $lib.intergenic.generate_bam.sh <<EOF
samtools view -b -h -L intergenic.bed $BAMFILE > output/intergenic/$lib.intergenic.bam
samtools sort -o output/intergenic/$lib.intergenic.sorted.bam output/intergenic/$lib.intergenic.bam
samtools index output/intergenic/$lib.intergenic.sorted.bam
EOF

done

module load mugqic/python/2.7.8
for lib in CK2 WT Spt6-SA
do
    for type in promoter downstream
    do
        grep -e "$type" output/intergenic/intergenic.split.bed > $type.bed
    
        for direction in sense antisense
        do
            prefix=output/intergenic/$lib.$type.$direction
            BAMFILE=$prefix.sorted.bam
            BIGWIG=$prefix.sorted.bigwig
                
            cat $lib.$type.$direction.generate_matrix.sh <<EOF
./bin/bamCoverage --skipNonCoveredRegions --binSize 10 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
./bin/computeMatrix reference-point --nanAfterEnd --numberOfProcessors 8 --regionsFileName $type.bed --scoreFileName $BIGWIG --outFileName $prefix.sorted.matrix.gz
EOF

        done
    done
    
    # And finally, do intergenic.
    BAMFILE=output/intergenic/$lib.intergenic.sorted.bam
    BIGWIG=output/intergenic/$lib.intergenic.sorted.bigwig
    
    cat $lib.$type.$direction.generate_matrix.sh <<EOF
./bin/bamCoverage --skipNonCoveredRegions --binSize 10 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
./bin/computeMatrix scale-regions --numberOfProcessors 8 --regionsFileName intergenic.bed --scoreFileName $BIGWIG --outFileName $prefix.sorted.matrix.gz
EOF

done

