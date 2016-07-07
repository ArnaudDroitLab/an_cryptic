#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -V

# USAGE: ReadCount5-3.sh bam_path output_dir result_label
BAMFILE=$1
outdir=$2
lib=$3


# Extract sense and antisense exons.
cd "/gs/project/eav-760-aa/RNA-CKII-Spt6"
mkdir -p $outdir
grep -P "\t\+\t" < input/exons.bed > input/sense.exons.bed
grep -P "\t-\t" < input/exons.bed > input/antisense.exons.bed

module load mugqic_dev/samtools/1.3

# Split bam files in two: one for sense, and one for antisense.
function senseSpecific {
    forwardBed=$1
    reverseBed=$2
    label=$3

    samtools view -b -h -f 16 -L $forwardBed $BAMFILE > $outdir/$lib.$label.1.bam
    samtools view -b -h -F 16 -L $reverseBed $BAMFILE > $outdir/$lib.$label.2.bam
    samtools merge $outdir/$lib.$label.bam $outdir/$lib.$label.1.bam $outdir/$lib.$label.2.bam
    samtools index $outdir/$lib.$label.bam
    rm $outdir/$lib.$label.1.bam $outdir/$lib.$label.2.bam
}

# Sense reads
#senseSpecific input/sense.exons.bed input/antisense.exons.bed sense
#senseSpecific input/antisense.exons.bed input/sense.exons.bed antisense

module load mugqic/python/2.7.8

echo "TEST!"

for direction in sense antisense
do
    BAMFILE=$outdir/$lib.$direction.bam
    BIGWIG=$outdir/$lib.$direction.bw
    
    echo $BAMFILE
    echo $BIGWIG
    
    echo "./bin/bamCoverage --skipNonCoveredRegions --binSize 20 --numberOfProcessors 4 -b $BAMFILE -o $BIGWIG"
    ./bin/bamCoverage --skipNonCoveredRegions --binSize 20 --numberOfProcessors 1 -b $BAMFILE -o $BIGWIG
    echo "./bin/computeMatrix scale-regions --numberOfProcessors 4 --regionsFileName input/exons.cut.bed --scoreFileName $BIGWIG --outFileName $outdir/$lib.$direction.sorted.matrix.gz"
    ./bin/computeMatrix scale-regions --numberOfProcessors 1 --regionsFileName input/exons.cut.bed --scoreFileName $BIGWIG --outFileName $outdir/$lib.$direction.sorted.matrix.gz
done

