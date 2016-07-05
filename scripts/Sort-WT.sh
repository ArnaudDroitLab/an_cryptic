#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -o Sort-WT.sh.stdout
#PBS -e Sort-WT.sh.stderr
#PBS -V
#PBS -N sort-WT


module add mugqic_dev/samtools/1.3
cd /gs/project/eav-760-aa/RNA-CKII-Spt6

samtools sort -@ 4 -o output/bwa/WT.sorted.bam output/bwa/WT.bam
samtools index output/bwa/WT.sorted.bam