#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -o Metagene.sh.stdout
#PBS -e Metagene.sh.stderr
#PBS -V
#PBS -N Metagene-split-sense

# Extract sense and antisense exons.
cd "/gs/project/eav-760-aa/RNA-CKII-Spt6"

grep -P "\t\+\t" < input/exons.bed > input/sense.exons.bed
grep -P "\t-\t" < input/exons.bed > input/antisense.exons.bed

module load mugqic_dev/samtools/1.3

# Split bam files in two: one for sense, and one for antisense. 
for lib in CK2 WT Spt6-SA
do
    BAMFILE=output/bwa/$lib.sorted.bam
    samtools view -b -h -f 16 -L input/sense.exons.bed $BAMFILE > output/bwa/$lib.sense1.bam
    samtools view -b -h -F 16 -L input/antisense.exons.bed $BAMFILE > output/bwa/$lib.sense2.bam
    samtools merge output/bwa/$lib.sense.bam output/bwa/$lib.sense1.bam output/bwa/$lib.sense2.bam
    samtools sort -o "output/bwa/"$lib"_sense_sorted.bam" output/bwa/$lib.sense.bam
    samtools index "output/bwa/"$lib"_sense_sorted.bam"
    rm output/bwa/$lib.sense1.bam output/bwa/$lib.sense2.bam output/bwa/$lib.sense.bam
    
    samtools view -h -f 16 -L input/antisense.exons.bed $BAMFILE > output/bwa/$lib.antisense1.bam
    samtools view -h -F 16 -L input/sense.exons.bed $BAMFILE > output/bwa/$lib.antisense2.bam
    samtools merge output/bwa/$lib.antisense.bam output/bwa/$lib.antisense1.bam output/bwa/$lib.antisense2.bam
    samtools sort -o "output/bwa/"$lib"_antisense_sorted.bam" output/bwa/$lib.antisense.bam
    samtools index "output/bwa/"$lib"_antisense_sorted.bam"
    rm output/bwa/$lib.antisense.bam output/bwa/$lib.antisense1.bam output/bwa/$lib.antisense2.bam
done

module load mugqic_dev/python/2.7.8
for lib in CK2 WT Spt6-SA
do
    for direction in sense antisense
    do
        BAMFILE="output/bwa/"$lib"_"$direction"_sorted.bam"
        BIGWIG="output/bwa/"$lib"_"$direction"_sorted.bigwig"
        
        #./bin/bamCoverage --skipNonCoveredRegions --binSize 20 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
        ./bin/computeMatrix scale-regions --numberOfProcessors 8 --regionsFileName input/exons.cut.bed --scoreFileName $BIGWIG --outFileName "output/bwa/"$lib"_"$direction"_sorted.matrix.gz"
    done
done
