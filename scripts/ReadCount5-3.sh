#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=4:00:00
#PBS -V


# USAGE: ReadCount5-3.sh bam_path output_dir result_label
inputBAM=$1
outdir=$2
label=$3

module add mugqic_dev/samtools/1.3


cd /gs/project/eav-760-aa/RNA-CKII-Spt6


mkdir -p $outdir
echo -e "Gene\tChr\tStart\tEnd\tStrand\tSense-All\tAntisense-All\tSense-5'\tAntisense-5'\tSense-3'\tAntisense-3'" >> $outdir/$label.results

while read p; do
    chr=`echo "$p" | cut -f 1`
    start=`echo "$p" | cut -f 2`
    end=`echo "$p" | cut -f 3`
    name=`echo "$p" | cut -f 4`
    strand=`echo "$p" | cut -f 6`
    
    
    # Count reads for the full transcript.
    forward=`samtools view -f 16 -c $inputBAM $chr:$start-$end`
    reverse=`samtools view -F 16 -c $inputBAM $chr:$start-$end`
      
    if [ $strand = "+" ]
    then
        sense=$forward
        antisense=$reverse
    else
        sense=$reverse
        antisense=$forward
    fi
    
    # Count reads in the 5' and 3' ends.
    let tenpercent=(end-start)/10
    
    if [ $strand = "+" ]
    then
        let start5=start
        let end5=start+tenpercent
        
        let start3=end-tenpercent
        let end3=end
    else
        let start3=start
        let end3=start+tenpercent
        
        let start5=end-tenpercent
        let end5=end
    fi
    
    forward5=`samtools view -f 16 -c $inputBAM $chr:$start5-$end5`
    reverse5=`samtools view -F 16 -c $inputBAM $chr:$start5-$end5`
  
    forward3=`samtools view -f 16 -c $inputBAM $chr:$start3-$end3`
    reverse3=`samtools view -F 16 -c $inputBAM $chr:$start3-$end3`
    
    if [ $strand = "+" ]
    then
        sense5=$forward5
        antisense5=$reverse5
        
        sense3=$forward3
        antisense3=$reverse3
    else
        sense5=$reverse5
        antisense5=$forward5
        
        sense3=$reverse3
        antisense3=$forward3
    fi
      
    echo -e "$name\t$chr\t$start\t$end\t$strand\t$sense\t$antisense\t$sense5\t$antisense5\t$sense3\t$antisense3" >> $outdir/$label.results
done <input/exons.bed

