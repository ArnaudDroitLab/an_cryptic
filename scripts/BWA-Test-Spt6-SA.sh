#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -A eav-760-aa
#PBS -o Spt6-SA-BWA-test.sh.stdout
#PBS -e Spt6-SA-BWA-test.sh.stderr
#PBS -V
#PBS -N Spt6-SA-BWA-test.sh

module add mugqic/bwa/0.7.10 
module add mugqic_dev/samtools/1.3

cd '/gs/project/eav-760-aa/RNA-CKII-Spt6'

indexPath=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa
inputFile=input/HI.2998.007.Index_8.Spt6-SA_R1.fastq.gz
libraryName=Spt6-SA

mkdir -p output/bwa
bwa mem  -t 10 $indexPath $inputFile | samtools view -h -b - > output/bwa/$libraryName.bam
