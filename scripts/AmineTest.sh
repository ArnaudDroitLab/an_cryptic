#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:00:00
#PBS -A eav-760-aa
#PBS -o AmineTest.sh.stdout
#PBS -e AmineTest.sh.stderr
#PBS -V
#PBS -N sense-antisense-introns

module add mugqic_dev/samtools/1.3

processLibrary () { 
  library=$1
  alignment=output/alignment/$library/$library.sorted.mdup.bam

  while read p; do
    chr=`echo "$p" | cut -f 1`
    start=`echo "$p" | cut -f 2`
    end=`echo "$p" | cut -f 3`
    name=`echo "$p" | cut -f 4`
    strand=`echo "$p" | cut -f 6`

    forward=`samtools view -f 16 -c $alignment $chr:$start-$end`
    reverse=`samtools view -F 16 -c $alignment $chr:$start-$end`
      
    if [ $strand = "+" ]
    then
      sense=$forward
      antisense=$reverse
    else
      sense=$reverse
      antisense=$forward
    fi
      
    echo -e "$name\t$sense\t$antisense" >> results.full.$library
  done <exons.bed
}

cd /gs/project/eav-760-aa/RNA-CKII-Spt6
rm results.full.*

processLibrary "Spt6-SA" &
processLibrary "CK2" &
processLibrary "WT" &

wait
