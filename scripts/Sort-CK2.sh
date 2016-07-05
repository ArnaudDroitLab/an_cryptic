#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=2:00:00
#PBS -A eav-760-aa
#PBS -o Sort-CK2.sh.stdout
#PBS -e Sort-CK2.sh.stderr
#PBS -V
#PBS -N sort-CK2


module add mugqic_dev/samtools/1.3
cd /gs/project/eav-760-aa/RNA-CKII-Spt6

samtools sort -@ 4 -o output/bwa/CK2.sorted.bam output/bwa/CK2.bam
samtools index output/bwa/CK2.sorted.bam