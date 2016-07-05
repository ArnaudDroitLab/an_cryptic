#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -o Sort-Spt6-SA.sh.stdout
#PBS -e Sort-Spt6-SA.sh.stderr
#PBS -V
#PBS -N sort-Spt6-SA


module add mugqic_dev/samtools/1.3
cd /gs/project/eav-760-aa/RNA-CKII-Spt6

samtools sort -@ 4 -o output/bwa/Spt6-SA.sorted.bam output/bwa/Spt6-SA.bam
samtools index output/bwa/Spt6-SA.sorted.bam