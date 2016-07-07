#!/bin/bash
for libraryPath in output-new/alignment/*/
do
  echo $i
  libraryName=`basename $libraryPath`
  echo "scripts/ReadCount5-3.sh $libraryPath/$libraryName.sorted.mdup.bam output-new/5-3Counts $libraryName" |
  qsub -N $libraryName.5-3 -o $libraryName.5-3.stdout -e $libraryName.5-3.stderr -A "eav-760-aa" -d `pwd`
done
