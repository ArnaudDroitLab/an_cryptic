#!/bin/bash
for libraryPath in output-new/alignment/*/
do
  echo $i
  libraryName=`basename $libraryPath`
  echo "./scripts/Metagene.sh $libraryPath/$libraryName.sorted.mdup.bam output-new/Metagene $libraryName" |
  qsub -N $libraryName.metagene -o $libraryName.metagene.stdout -e $libraryName.metagene.stderr -A "eav-760-aa"
done
