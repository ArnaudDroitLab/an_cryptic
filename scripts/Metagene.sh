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
    samtools merge $outdir/$lib.sense.bam $outdir/$lib.$label.1.bam $outdir/$lib.$label.2.bam
    samtools index $outdir/$lib.sense.bam
    rm $outdir/$lib.$label.1.bam $outdir/$lib.$label.2.bam
}

# Sense reads
senseSpecific input/sense.exons.bed input/antisense.exons.bed sense
senseSpecific input/antisense.exons.bed input/sense.exons.bed antisense

for direction in sense antisense
do
    BAMFILE=$outdir/$lib.$direction.bam
    BIGWIG=$outdir/$lib.$direction.bw
    
    ./bin/bamCoverage --skipNonCoveredRegions --binSize 20 --numberOfProcessors 8 -b $BAMFILE -o $BIGWIG
    ./bin/computeMatrix scale-regions --numberOfProcessors 8 --regionsFileName input/exons.cut.bed --scoreFileName $BIGWIG --outFileName $outdir/$lib.$direction.sorted.matrix.gz
done

