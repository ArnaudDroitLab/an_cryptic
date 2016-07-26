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

# UCSC-obtained data has "chr" in front of chromosome number. Remove that.
sed -e 's/^chr//' input/UCSC_SGD_Genes.bed > output/intergenic/UCSC_SGD_Genes.nochr.bed

# Define a list of all intergenic regions
bedtools subtract -a input/All.bed -b output/intergenic/UCSC_SGD_Genes.nochr.bed > output/intergenic/intergenic.bed

# Split and annotate all intergenic regions based on their proximity to the closest gene.
perl scripts/annotate_split.pl output/intergenic/UCSC_SGD_Genes.nochr.bed 1000 < output/intergenic/intergenic.bed > output/intergenic/intergenic.split.bed

# Split whole annotation into promoter, downstream and intergenic files.
for type in promoter downstream
do
    grep -e "$type.*+$" output/intergenic/intergenic.split.bed > $type.plus.bed
    grep -e "$type.*-$" output/intergenic/intergenic.split.bed > $type.minus.bed
done
grep -e "intergenic" output/intergenic/intergenic.split.bed > intergenic.bed

# Declare a function for submitting jobs using preset parameters.
RAP_ID=eav-760-aa
function submit_job {
    if [ "$2" = "" ]
    then
        dependFlag=""
    else
        dependFlag="-W depend=afterok:$2"
    fi
    
    qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d `pwd` \
         -j oe -o $1.stdout \
         -e $1.stderr \
         -N $1 \
         -l walltime=1:00:0 -q metaq -l nodes=1:ppn=8 -l pmem=1700m \
         $dependFlag $1 | grep "[0-9]"
}

# Declare a named array to keep job IDs to keep track of dependencies.
declare -A job_ids

# Split the bam files
for lib in CK2 WT Spt6-SA CK2-2 WT-2 Spt6-SA-2
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
            
            cat > $lib.$type.$direction.generate_bam.sh <<EOF
module add mugqic/samtools/1.3
samtools view -b -h $plus_samtools_arg 16 -L $type.plus.bed $BAMFILE > $prefix.1.bam
samtools view -b -h $minus_samtools_arg 16 -L $type.minus.bed $BAMFILE > $prefix.2.bam
samtools merge $prefix.sorted.bam $prefix.1.bam $prefix.2.bam
samtools index $prefix.sorted.bam
rm $prefix.1.bam $prefix.2.bam   
EOF
            job_ids[$lib.$type.$direction]=$(submit_job $lib.$type.$direction.generate_bam.sh )
        done
    done
    
    cat > $lib.intergenic.generate_bam.sh <<EOF
module add mugqic/samtools/1.3
samtools view -b -h -L intergenic.bed $BAMFILE > output/intergenic/$lib.intergenic.sorted.bam
samtools index output/intergenic/$lib.intergenic.sorted.bam
EOF
    job_ids[$lib.intergenic]=$(submit_job $lib.intergenic.generate_bam.sh)
done


for lib in CK2 WT Spt6-SA CK2-2 WT-2 Spt6-SA-2
do
    for type in promoter downstream
    do
        grep -e "$type" output/intergenic/intergenic.split.bed > $type.bed
    
        for direction in sense antisense
        do
            prefix=output/intergenic/$lib.$type.$direction
            BAMFILE=$prefix.sorted.bam
            BIGWIG=$prefix.sorted.bigwig
                
            cat > $lib.$type.$direction.generate_matrix.sh <<EOF
module load mugqic/python/2.7.8            
./bin/bamCoverage --skipNonCoveredRegions --binSize 10 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
./bin/computeMatrix reference-point --nanAfterEnd --numberOfProcessors 8 --regionsFileName $type.bed --scoreFileName $BIGWIG --outFileName $prefix.sorted.matrix.gz
EOF
            submit_job $lib.$type.$direction.generate_matrix.sh ${job_ids[$lib.$type.$direction]}
        done
    done
    
    # And finally, do intergenic.
    BAMFILE=output/intergenic/$lib.intergenic.sorted.bam
    BIGWIG=output/intergenic/$lib.intergenic.sorted.bigwig
    
    cat > $lib.intergenic.generate_matrix.sh <<EOF
module load mugqic/python/2.7.8    
./bin/bamCoverage --skipNonCoveredRegions --binSize 10 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
./bin/computeMatrix scale-regions --numberOfProcessors 8 --regionsFileName intergenic.bed --scoreFileName $BIGWIG --outFileName $prefix.sorted.matrix.gz
EOF
    submit_job $lib.intergenic.generate_matrix.sh ${job_ids[$lib.intergenic]}
done

