#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq PBSScheduler Job Submission Bash script
# Version: 2.1.1
# Created on: 2015-12-06T01:32:36
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 3 jobs
#   merge_trimmomatic_stats: 1 job
#   star: 8 jobs
#   picard_merge_sam_files: 0 job... skipping
#   picard_sort_sam: 3 jobs
#   picard_mark_duplicates: 3 jobs
#   picard_rna_metrics: 3 jobs
#   estimate_ribosomal_rna: 3 jobs
#   rnaseqc: 2 jobs
#   wiggle: 12 jobs
#   raw_counts: 3 jobs
#   raw_counts_metrics: 4 jobs
#   cufflinks: 3 jobs
#   cuffmerge: 1 job
#   cuffquant: 3 jobs
#   cuffdiff: 3 jobs
#   cuffnorm: 1 job
#   fpkm_correlation_matrix: 2 jobs
#   gq_seq_utils_exploratory_analysis_rnaseq: 3 jobs
#   differential_expression: 1 job
#   differential_expression_goseq: 4 jobs
#   TOTAL: 66 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/project/eav-760-aa/RNA-CKII-Spt6/output
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt6-SA
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt6-SA.70dcaa4f63713f783abced560fff8fbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt6-SA.70dcaa4f63713f783abced560fff8fbd.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.32 && \
mkdir -p trim/Spt6-SA && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/RNA-CKII-Spt6/input/HI.2998.007.Index_8.Spt6-SA_R1.fastq.gz \
  trim/Spt6-SA/Spt6-SA.trim.single.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/Spt6-SA/Spt6-SA.trim.log
trimmomatic.Spt6-SA.70dcaa4f63713f783abced560fff8fbd.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.CK2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.CK2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.CK2.150ac345c84c7030e7ead68311ff46aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.CK2.150ac345c84c7030e7ead68311ff46aa.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.32 && \
mkdir -p trim/CK2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/RNA-CKII-Spt6/input/HI.2998.007.Index_9.CK2_R1.fastq.gz \
  trim/CK2/CK2.trim.single.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/CK2/CK2.trim.log
trimmomatic.CK2.150ac345c84c7030e7ead68311ff46aa.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.WT
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.WT
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.WT.68c97d2b15825961f254ecbc3a8534d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.WT.68c97d2b15825961f254ecbc3a8534d1.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.32 && \
mkdir -p trim/WT && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/RNA-CKII-Spt6/input/HI.2998.007.Index_3.WT_R1.fastq.gz \
  trim/WT/WT.trim.single.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/WT/WT.trim.log
trimmomatic.WT.68c97d2b15825961f254ecbc3a8534d1.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.bfda13232a1cfe9cef52443c254e33e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.bfda13232a1cfe9cef52443c254e33e9.mugqic.done'
module load mugqic/pandoc/1.13.1 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt6-SA/Spt6-SA.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt6-SA	Spt6-SA	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/CK2/CK2.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/CK2	CK2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/WT/WT.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/WT	WT	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=32 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.bfda13232a1cfe9cef52443c254e33e9.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: star
#-------------------------------------------------------------------------------
STEP=star
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: star_1_JOB_ID: star_align.1.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.Spt6-SA
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/star/star_align.1.Spt6-SA.e9d37a891417e3270536aa71a566e0dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.Spt6-SA.e9d37a891417e3270536aa71a566e0dd.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment_1stPass/Spt6-SA/Spt6-SA && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/star_index/Ensembl77.sjdbOverhang99 \
  --readFilesIn \
    trim/Spt6-SA/Spt6-SA.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/Spt6-SA/Spt6-SA/ \
  --outSAMattrRGline ID:"Spt6-SA" 	PL:"ILLUMINA" 			SM:"Spt6-SA" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.Spt6-SA.e9d37a891417e3270536aa71a566e0dd.mugqic.done
)
star_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_2_JOB_ID: star_align.1.CK2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.CK2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/star/star_align.1.CK2.a5620bf5d240d9c2a4f1f83afb9d06c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.CK2.a5620bf5d240d9c2a4f1f83afb9d06c5.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment_1stPass/CK2/CK2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/star_index/Ensembl77.sjdbOverhang99 \
  --readFilesIn \
    trim/CK2/CK2.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/CK2/CK2/ \
  --outSAMattrRGline ID:"CK2" 	PL:"ILLUMINA" 			SM:"CK2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.CK2.a5620bf5d240d9c2a4f1f83afb9d06c5.mugqic.done
)
star_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_3_JOB_ID: star_align.1.WT
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.WT
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/star/star_align.1.WT.46add2a83e3b3bf93906e298c9b73a57.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.WT.46add2a83e3b3bf93906e298c9b73a57.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment_1stPass/WT/WT && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/star_index/Ensembl77.sjdbOverhang99 \
  --readFilesIn \
    trim/WT/WT.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/WT/WT/ \
  --outSAMattrRGline ID:"WT" 	PL:"ILLUMINA" 			SM:"WT" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.WT.46add2a83e3b3bf93906e298c9b73a57.mugqic.done
)
star_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_4_JOB_ID: star_index.AllSamples
#-------------------------------------------------------------------------------
JOB_NAME=star_index.AllSamples
JOB_DEPENDENCIES=$star_1_JOB_ID:$star_2_JOB_ID:$star_3_JOB_ID
JOB_DONE=job_output/star/star_index.AllSamples.be22b25e2039b86fa0d70da08d9d0bd5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_index.AllSamples.be22b25e2039b86fa0d70da08d9d0bd5.mugqic.done'
module load mugqic/star/2.4.0f1 && \
cat \
  alignment_1stPass/Spt6-SA/Spt6-SA/SJ.out.tab \
  alignment_1stPass/CK2/CK2/SJ.out.tab \
  alignment_1stPass/WT/WT/SJ.out.tab | \
awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p reference.Merged && \
STAR --runMode genomeGenerate \
  --genomeDir reference.Merged \
  --genomeFastaFiles /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  --runThreadN 12 \
  --limitGenomeGenerateRAM 50000000000 \
  --sjdbFileChrStartEnd alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 4000000000 \
  --sjdbOverhang 99
star_index.AllSamples.be22b25e2039b86fa0d70da08d9d0bd5.mugqic.done
)
star_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_5_JOB_ID: star_align.2.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.Spt6-SA
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$star_4_JOB_ID
JOB_DONE=job_output/star/star_align.2.Spt6-SA.e7325abd1311885a96e23290bed27efe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.Spt6-SA.e7325abd1311885a96e23290bed27efe.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment/Spt6-SA/Spt6-SA && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/Spt6-SA/Spt6-SA.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/Spt6-SA/Spt6-SA/ \
  --outSAMattrRGline ID:"Spt6-SA" 	PL:"ILLUMINA" 			SM:"Spt6-SA" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f Spt6-SA/Aligned.sortedByCoord.out.bam alignment/Spt6-SA/Spt6-SA.sorted.bam
star_align.2.Spt6-SA.e7325abd1311885a96e23290bed27efe.mugqic.done
)
star_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_6_JOB_ID: star_align.2.CK2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.CK2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID:$star_4_JOB_ID
JOB_DONE=job_output/star/star_align.2.CK2.164bcc5ddd2181fa7df9cf2118a6ffd2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.CK2.164bcc5ddd2181fa7df9cf2118a6ffd2.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment/CK2/CK2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/CK2/CK2.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/CK2/CK2/ \
  --outSAMattrRGline ID:"CK2" 	PL:"ILLUMINA" 			SM:"CK2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f CK2/Aligned.sortedByCoord.out.bam alignment/CK2/CK2.sorted.bam
star_align.2.CK2.164bcc5ddd2181fa7df9cf2118a6ffd2.mugqic.done
)
star_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_7_JOB_ID: star_align.2.WT
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.WT
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID:$star_4_JOB_ID
JOB_DONE=job_output/star/star_align.2.WT.c4b852c1d783e1e97384dda7f75842db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.WT.c4b852c1d783e1e97384dda7f75842db.mugqic.done'
module load mugqic/star/2.4.0f1 && \
mkdir -p alignment/WT/WT && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/WT/WT.trim.single.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/WT/WT/ \
  --outSAMattrRGline ID:"WT" 	PL:"ILLUMINA" 			SM:"WT" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f WT/Aligned.sortedByCoord.out.bam alignment/WT/WT.sorted.bam
star_align.2.WT.c4b852c1d783e1e97384dda7f75842db.mugqic.done
)
star_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_8_JOB_ID: star_report
#-------------------------------------------------------------------------------
JOB_NAME=star_report
JOB_DEPENDENCIES=$star_5_JOB_ID:$star_6_JOB_ID:$star_7_JOB_ID
JOB_DONE=job_output/star/star_report.e23cea8170fa6b18e1b6a65226273205.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_report.e23cea8170fa6b18e1b6a65226273205.mugqic.done'
module load mugqic/pandoc/1.13.1 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Saccharomyces_cerevisiae" \
  --variable assembly="R64-1-1" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.e23cea8170fa6b18e1b6a65226273205.mugqic.done
)
star_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_sort_sam
#-------------------------------------------------------------------------------
STEP=picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_1_JOB_ID: picard_sort_sam.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.Spt6-SA
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.Spt6-SA.667a76e520db96bf36d0694ccc84e92c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.Spt6-SA.667a76e520db96bf36d0694ccc84e92c.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.bam \
  OUTPUT=alignment/Spt6-SA/Spt6-SA.QueryNameSorted.bam \
  SORT_ORDER=queryname \
  MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.Spt6-SA.667a76e520db96bf36d0694ccc84e92c.mugqic.done
)
picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_2_JOB_ID: picard_sort_sam.CK2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.CK2
JOB_DEPENDENCIES=$star_6_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.CK2.c2848511972dc1b4c1b3961cbcfabce7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.CK2.c2848511972dc1b4c1b3961cbcfabce7.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/CK2/CK2.sorted.bam \
  OUTPUT=alignment/CK2/CK2.QueryNameSorted.bam \
  SORT_ORDER=queryname \
  MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.CK2.c2848511972dc1b4c1b3961cbcfabce7.mugqic.done
)
picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_3_JOB_ID: picard_sort_sam.WT
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.WT
JOB_DEPENDENCIES=$star_7_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.WT.b857d91299510eef613bfc3ee5e2588c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.WT.b857d91299510eef613bfc3ee5e2588c.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/WT/WT.sorted.bam \
  OUTPUT=alignment/WT/WT.QueryNameSorted.bam \
  SORT_ORDER=queryname \
  MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.WT.b857d91299510eef613bfc3ee5e2588c.mugqic.done
)
picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt6-SA
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt6-SA.6257dfffaf7a619e5ab0c1d6acf535e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt6-SA.6257dfffaf7a619e5ab0c1d6acf535e9.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.bam \
  OUTPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  METRICS_FILE=alignment/Spt6-SA/Spt6-SA.sorted.mdup.metrics \
  MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.Spt6-SA.6257dfffaf7a619e5ab0c1d6acf535e9.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.CK2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.CK2
JOB_DEPENDENCIES=$star_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.CK2.7988dc9f2d1f42bd54c0059c4aeb331a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.CK2.7988dc9f2d1f42bd54c0059c4aeb331a.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/CK2/CK2.sorted.bam \
  OUTPUT=alignment/CK2/CK2.sorted.mdup.bam \
  METRICS_FILE=alignment/CK2/CK2.sorted.mdup.metrics \
  MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.CK2.7988dc9f2d1f42bd54c0059c4aeb331a.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.WT
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.WT
JOB_DEPENDENCIES=$star_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.WT.6936ea3914a7aff567e1ad22a47df3d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.WT.6936ea3914a7aff567e1ad22a47df3d6.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/WT/WT.sorted.bam \
  OUTPUT=alignment/WT/WT.sorted.mdup.bam \
  METRICS_FILE=alignment/WT/WT.sorted.mdup.metrics \
  MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.WT.6936ea3914a7aff567e1ad22a47df3d6.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_rna_metrics
#-------------------------------------------------------------------------------
STEP=picard_rna_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_1_JOB_ID: picard_rna_metrics.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.Spt6-SA
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.Spt6-SA.6e3edfe308bb1f5d306738965a082cd0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.Spt6-SA.6e3edfe308bb1f5d306738965a082cd0.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p metrics/Spt6-SA && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectMultipleMetrics.jar \
  PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=/gs/scratch/$USER \
  REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  OUTPUT=metrics/Spt6-SA/Spt6-SA \
  MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectRnaSeqMetrics.jar \
  VALIDATION_STRINGENCY=SILENT  \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  OUTPUT=metrics/Spt6-SA/Spt6-SA \
  REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.ref_flat.tsv \
  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
  MINIMUM_LENGTH=200 \
  REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.Spt6-SA.6e3edfe308bb1f5d306738965a082cd0.mugqic.done
)
picard_rna_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_2_JOB_ID: picard_rna_metrics.CK2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.CK2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.CK2.d6284635fcfdf5956888753b264d6ce0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.CK2.d6284635fcfdf5956888753b264d6ce0.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p metrics/CK2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectMultipleMetrics.jar \
  PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=/gs/scratch/$USER \
  REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  INPUT=alignment/CK2/CK2.sorted.mdup.bam \
  OUTPUT=metrics/CK2/CK2 \
  MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectRnaSeqMetrics.jar \
  VALIDATION_STRINGENCY=SILENT  \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/CK2/CK2.sorted.mdup.bam \
  OUTPUT=metrics/CK2/CK2 \
  REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.ref_flat.tsv \
  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
  MINIMUM_LENGTH=200 \
  REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.CK2.d6284635fcfdf5956888753b264d6ce0.mugqic.done
)
picard_rna_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_3_JOB_ID: picard_rna_metrics.WT
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.WT
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.WT.144bcdd739b023f26dbdbd116e779fbe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.WT.144bcdd739b023f26dbdbd116e779fbe.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p metrics/WT && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectMultipleMetrics.jar \
  PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=/gs/scratch/$USER \
  REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  INPUT=alignment/WT/WT.sorted.mdup.bam \
  OUTPUT=metrics/WT/WT \
  MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/CollectRnaSeqMetrics.jar \
  VALIDATION_STRINGENCY=SILENT  \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/WT/WT.sorted.mdup.bam \
  OUTPUT=metrics/WT/WT \
  REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.ref_flat.tsv \
  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
  MINIMUM_LENGTH=200 \
  REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.WT.144bcdd739b023f26dbdbd116e779fbe.mugqic.done
)
picard_rna_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: estimate_ribosomal_rna
#-------------------------------------------------------------------------------
STEP=estimate_ribosomal_rna
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_1_JOB_ID: bwa_mem_rRNA.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.Spt6-SA
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.Spt6-SA.ee2c0d72215153a8ec02f4a47bb4d107.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.Spt6-SA.ee2c0d72215153a8ec02f4a47bb4d107.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bvatools/1.5 mugqic/bwa/0.7.10 mugqic/picard/1.123 mugqic/mugqic_tools/2.1.0 mugqic/python/2.7.8 && \
mkdir -p alignment/Spt6-SA/Spt6-SA metrics/Spt6-SA/Spt6-SA && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/Spt6-SA/Spt6-SA/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:Spt6-SA	SM:Spt6-SA	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/rrna_bwa_index/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.rrna.fa \
  /dev/stdin | \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=/dev/stdin \
  OUTPUT=metrics/Spt6-SA/Spt6-SA/Spt6-SArRNA.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/Spt6-SA/Spt6-SA/Spt6-SArRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  -o metrics/Spt6-SA/Spt6-SA/Spt6-SArRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.Spt6-SA.ee2c0d72215153a8ec02f4a47bb4d107.mugqic.done
)
estimate_ribosomal_rna_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_2_JOB_ID: bwa_mem_rRNA.CK2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.CK2
JOB_DEPENDENCIES=$star_6_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.CK2.831d862d95409d0ac39695e4a9e7ae58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.CK2.831d862d95409d0ac39695e4a9e7ae58.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bvatools/1.5 mugqic/bwa/0.7.10 mugqic/picard/1.123 mugqic/mugqic_tools/2.1.0 mugqic/python/2.7.8 && \
mkdir -p alignment/CK2/CK2 metrics/CK2/CK2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/CK2/CK2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:CK2	SM:CK2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/rrna_bwa_index/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.rrna.fa \
  /dev/stdin | \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=/dev/stdin \
  OUTPUT=metrics/CK2/CK2/CK2rRNA.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/CK2/CK2/CK2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  -o metrics/CK2/CK2/CK2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.CK2.831d862d95409d0ac39695e4a9e7ae58.mugqic.done
)
estimate_ribosomal_rna_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_3_JOB_ID: bwa_mem_rRNA.WT
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.WT
JOB_DEPENDENCIES=$star_7_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.WT.a81f5a9201d6b24c20633821569f6e6f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.WT.a81f5a9201d6b24c20633821569f6e6f.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bvatools/1.5 mugqic/bwa/0.7.10 mugqic/picard/1.123 mugqic/mugqic_tools/2.1.0 mugqic/python/2.7.8 && \
mkdir -p alignment/WT/WT metrics/WT/WT && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/WT/WT/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:WT	SM:WT	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/rrna_bwa_index/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.rrna.fa \
  /dev/stdin | \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=/dev/stdin \
  OUTPUT=metrics/WT/WT/WTrRNA.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/WT/WT/WTrRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  -o metrics/WT/WT/WTrRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.WT.a81f5a9201d6b24c20633821569f6e6f.mugqic.done
)
estimate_ribosomal_rna_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: rnaseqc
#-------------------------------------------------------------------------------
STEP=rnaseqc
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: rnaseqc_1_JOB_ID: rnaseqc
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc.16dbe8800f9b952b9114cf342998fe5f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc.16dbe8800f9b952b9114cf342998fe5f.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.10 mugqic/rnaseqc/1.1.7 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
Spt6-SA	alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam	RNAseq
CK2	alignment/CK2/CK2.sorted.mdup.bam	RNAseq
WT	alignment/WT/WT.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $RNASEQC_JAR \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.transcript_id.gtf \
  -ttype 2 \
  -singleEnd && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep
rnaseqc.16dbe8800f9b952b9114cf342998fe5f.mugqic.done
)
rnaseqc_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$rnaseqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: rnaseqc_2_JOB_ID: rnaseqc_report
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc_report.1f34e6d475e7a6029ec66e88d3767d1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc_report.1f34e6d475e7a6029ec66e88d3767d1b.mugqic.done'
module load mugqic/python/2.7.8 mugqic/pandoc/1.13.1 && \
mkdir -p report && \
cp metrics/rnaseqRep.zip report/reportRNAseqQC.zip && \
python -c 'import csv; csv_in = csv.DictReader(open("metrics/rnaseqRep/metrics.tsv"), delimiter="	")
print "	".join(["Sample", "Aligned Reads", "Alternative Alignments", "%", "rRNA Reads", "Coverage", "Exonic Rate", "Genes"])
print "\n".join(["	".join([
    line["Sample"],
    line["Mapped"],
    line["Alternative Aligments"],
    str(float(line["Alternative Aligments"]) / float(line["Mapped"]) * 100),
    line["rRNA"],
    line["Mean Per Base Cov."],
    line["Exonic Rate"],
    line["Genes Detected"]
]) for line in csv_in])' \
  > report/trimAlignmentTable.tsv.tmp && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F"	" 'FNR==NR{raw_reads[$1]=$2; surviving_reads[$1]=$3; surviving_pct[$1]=$4; next}{OFS="	"; if ($2=="Aligned Reads"){surviving_pct[$1]="%"; aligned_pct="%"; rrna_pct="%"} else {aligned_pct=($2 / surviving_reads[$1] * 100); rrna_pct=($5 / surviving_reads[$1] * 100)}; printf $1"	"raw_reads[$1]"	"surviving_reads[$1]"	"surviving_pct[$1]"	"$2"	"aligned_pct"	"$3"	"$4"	"$5"	"rrna_pct; for (i = 6; i<= NF; i++) {printf "	"$i}; print ""}' \
  metrics/trimSampleTable.tsv \
  report/trimAlignmentTable.tsv.tmp \
  > report/trimAlignmentTable.tsv
else
  cp report/trimAlignmentTable.tsv.tmp report/trimAlignmentTable.tsv
fi && \
rm report/trimAlignmentTable.tsv.tmp && \
trim_alignment_table_md=`if [[ -f metrics/trimSampleTable.tsv ]] ; then cut -f1-13 report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8), sprintf("%\47d", $9), sprintf("%.1f", $10), sprintf("%.2f", $11), sprintf("%.2f", $12), sprintf("%\47d", $13)}}' ; else cat report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.2f", $6), sprintf("%.2f", $7), $8}}' ; fi`
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.rnaseqc.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.rnaseqc.md \
  --variable trim_alignment_table="$trim_alignment_table_md" \
  --to markdown \
  > report/RnaSeq.rnaseqc.md
rnaseqc_report.1f34e6d475e7a6029ec66e88d3767d1b.mugqic.done
)
rnaseqc_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$rnaseqc_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: wiggle
#-------------------------------------------------------------------------------
STEP=wiggle
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: wiggle_1_JOB_ID: wiggle.Spt6-SA.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.Spt6-SA.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.Spt6-SA.forward_strandspec.f7720af70e0ae164e5230f9542297944.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.Spt6-SA.forward_strandspec.f7720af70e0ae164e5230f9542297944.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
samtools view -bh -F 256 -f 81 \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  > alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  > alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.forward.bam \
  OUTPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.forward.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.forward.bam alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.forward.bam
wiggle.Spt6-SA.forward_strandspec.f7720af70e0ae164e5230f9542297944.mugqic.done
)
wiggle_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_2_JOB_ID: wiggle.Spt6-SA.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.Spt6-SA.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.Spt6-SA.reverse_strandspec.f5e4b858a4210b4df2bf7da55fd05e83.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.Spt6-SA.reverse_strandspec.f5e4b858a4210b4df2bf7da55fd05e83.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p tracks/Spt6-SA tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  > alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  > alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.reverse.bam \
  OUTPUT=alignment/Spt6-SA/Spt6-SA.sorted.mdup.reverse.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp1.reverse.bam alignment/Spt6-SA/Spt6-SA.sorted.mdup.tmp2.reverse.bam
wiggle.Spt6-SA.reverse_strandspec.f5e4b858a4210b4df2bf7da55fd05e83.mugqic.done
)
wiggle_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_3_JOB_ID: wiggle.Spt6-SA.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.Spt6-SA.forward
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.Spt6-SA.forward.13a6360dfdbb39ef7ffa87abee89bb52.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.Spt6-SA.forward.13a6360dfdbb39ef7ffa87abee89bb52.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/Spt6-SA tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/Spt6-SA/Spt6-SA.forward.bedGraph && \
bedGraphToBigWig \
  tracks/Spt6-SA/Spt6-SA.forward.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/Spt6-SA.forward.bw
wiggle.Spt6-SA.forward.13a6360dfdbb39ef7ffa87abee89bb52.mugqic.done
)
wiggle_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_4_JOB_ID: wiggle.Spt6-SA.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.Spt6-SA.reverse
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.Spt6-SA.reverse.a4e5ccbd3b62fc6ed127a9c5bf48030c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.Spt6-SA.reverse.a4e5ccbd3b62fc6ed127a9c5bf48030c.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/Spt6-SA tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/Spt6-SA/Spt6-SA.reverse.bedGraph && \
bedGraphToBigWig \
  tracks/Spt6-SA/Spt6-SA.reverse.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/Spt6-SA.reverse.bw
wiggle.Spt6-SA.reverse.a4e5ccbd3b62fc6ed127a9c5bf48030c.mugqic.done
)
wiggle_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_5_JOB_ID: wiggle.CK2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.CK2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.CK2.forward_strandspec.b1ff2da1fe5b070692bfaba264fbcc70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.CK2.forward_strandspec.b1ff2da1fe5b070692bfaba264fbcc70.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
samtools view -bh -F 256 -f 81 \
  alignment/CK2/CK2.sorted.mdup.bam \
  > alignment/CK2/CK2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/CK2/CK2.sorted.mdup.bam \
  > alignment/CK2/CK2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/CK2/CK2.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/CK2/CK2.sorted.mdup.tmp2.forward.bam \
  OUTPUT=alignment/CK2/CK2.sorted.mdup.forward.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/CK2/CK2.sorted.mdup.tmp1.forward.bam alignment/CK2/CK2.sorted.mdup.tmp2.forward.bam
wiggle.CK2.forward_strandspec.b1ff2da1fe5b070692bfaba264fbcc70.mugqic.done
)
wiggle_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_6_JOB_ID: wiggle.CK2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.CK2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.CK2.reverse_strandspec.c11a8b6ed1665ddb1dd59ce6d4c73b6b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.CK2.reverse_strandspec.c11a8b6ed1665ddb1dd59ce6d4c73b6b.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p tracks/CK2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/CK2/CK2.sorted.mdup.bam \
  > alignment/CK2/CK2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/CK2/CK2.sorted.mdup.bam \
  > alignment/CK2/CK2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/CK2/CK2.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/CK2/CK2.sorted.mdup.tmp2.reverse.bam \
  OUTPUT=alignment/CK2/CK2.sorted.mdup.reverse.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/CK2/CK2.sorted.mdup.tmp1.reverse.bam alignment/CK2/CK2.sorted.mdup.tmp2.reverse.bam
wiggle.CK2.reverse_strandspec.c11a8b6ed1665ddb1dd59ce6d4c73b6b.mugqic.done
)
wiggle_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_7_JOB_ID: wiggle.CK2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.CK2.forward
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.CK2.forward.d43ced24dca3ff399cf8daba547769f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.CK2.forward.d43ced24dca3ff399cf8daba547769f1.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/CK2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/CK2/CK2.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/CK2/CK2.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/CK2/CK2.forward.bedGraph && \
bedGraphToBigWig \
  tracks/CK2/CK2.forward.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/CK2.forward.bw
wiggle.CK2.forward.d43ced24dca3ff399cf8daba547769f1.mugqic.done
)
wiggle_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_8_JOB_ID: wiggle.CK2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.CK2.reverse
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.CK2.reverse.2681af88ce8d1aa18a8202f393df26bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.CK2.reverse.2681af88ce8d1aa18a8202f393df26bc.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/CK2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/CK2/CK2.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/CK2/CK2.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/CK2/CK2.reverse.bedGraph && \
bedGraphToBigWig \
  tracks/CK2/CK2.reverse.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/CK2.reverse.bw
wiggle.CK2.reverse.2681af88ce8d1aa18a8202f393df26bc.mugqic.done
)
wiggle_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_9_JOB_ID: wiggle.WT.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.WT.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.WT.forward_strandspec.d678bfec10da040f869c07f15d2dd833.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.WT.forward_strandspec.d678bfec10da040f869c07f15d2dd833.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
samtools view -bh -F 256 -f 81 \
  alignment/WT/WT.sorted.mdup.bam \
  > alignment/WT/WT.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/WT/WT.sorted.mdup.bam \
  > alignment/WT/WT.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/WT/WT.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/WT/WT.sorted.mdup.tmp2.forward.bam \
  OUTPUT=alignment/WT/WT.sorted.mdup.forward.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/WT/WT.sorted.mdup.tmp1.forward.bam alignment/WT/WT.sorted.mdup.tmp2.forward.bam
wiggle.WT.forward_strandspec.d678bfec10da040f869c07f15d2dd833.mugqic.done
)
wiggle_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_10_JOB_ID: wiggle.WT.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.WT.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.WT.reverse_strandspec.1f80d924968e3c625a60c650899f8ed2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.WT.reverse_strandspec.1f80d924968e3c625a60c650899f8ed2.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p tracks/WT tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/WT/WT.sorted.mdup.bam \
  > alignment/WT/WT.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/WT/WT.sorted.mdup.bam \
  > alignment/WT/WT.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/gs/scratch/$USER \
  INPUT=alignment/WT/WT.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/WT/WT.sorted.mdup.tmp2.reverse.bam \
  OUTPUT=alignment/WT/WT.sorted.mdup.reverse.bam \
  MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/WT/WT.sorted.mdup.tmp1.reverse.bam alignment/WT/WT.sorted.mdup.tmp2.reverse.bam
wiggle.WT.reverse_strandspec.1f80d924968e3c625a60c650899f8ed2.mugqic.done
)
wiggle_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_11_JOB_ID: wiggle.WT.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.WT.forward
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.WT.forward.b28250dd101fc59f3589ec072a753a0c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.WT.forward.b28250dd101fc59f3589ec072a753a0c.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/WT tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/WT/WT.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/WT/WT.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/WT/WT.forward.bedGraph && \
bedGraphToBigWig \
  tracks/WT/WT.forward.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/WT.forward.bw
wiggle.WT.forward.b28250dd101fc59f3589ec072a753a0c.mugqic.done
)
wiggle_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_12_JOB_ID: wiggle.WT.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.WT.reverse
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.WT.reverse.1d3bd21784471b6d5e4ddbae95eb509b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.WT.reverse.1d3bd21784471b6d5e4ddbae95eb509b.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.22.1 mugqic/ucsc/20140212 && \
mkdir -p tracks/WT tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/WT/WT.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/WT/WT.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  > tracks/WT/WT.reverse.bedGraph && \
bedGraphToBigWig \
  tracks/WT/WT.reverse.bedGraph \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa.fai \
  tracks/bigWig/WT.reverse.bw
wiggle.WT.reverse.1d3bd21784471b6d5e4ddbae95eb509b.mugqic.done
)
wiggle_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------
STEP=raw_counts
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_1_JOB_ID: htseq_count.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.Spt6-SA
JOB_DEPENDENCIES=$picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.Spt6-SA.298ad5115aa7769a77c1c28f81020ee3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.Spt6-SA.298ad5115aa7769a77c1c28f81020ee3.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/python/2.7.8 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/Spt6-SA/Spt6-SA.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  > raw_counts/Spt6-SA.readcounts.csv
htseq_count.Spt6-SA.298ad5115aa7769a77c1c28f81020ee3.mugqic.done
)
raw_counts_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.CK2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.CK2
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.CK2.3c5b79beb8dc2d7cbb7286ba492dff6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.CK2.3c5b79beb8dc2d7cbb7286ba492dff6c.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/python/2.7.8 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/CK2/CK2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  > raw_counts/CK2.readcounts.csv
htseq_count.CK2.3c5b79beb8dc2d7cbb7286ba492dff6c.mugqic.done
)
raw_counts_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_3_JOB_ID: htseq_count.WT
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.WT
JOB_DEPENDENCIES=$picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.WT.458fd13debfe66e593d87c101142b922.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.WT.458fd13debfe66e593d87c101142b922.mugqic.done'
module load mugqic/samtools/0.1.19 mugqic/python/2.7.8 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/WT/WT.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  > raw_counts/WT.readcounts.csv
htseq_count.WT.458fd13debfe66e593d87c101142b922.mugqic.done
)
raw_counts_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: raw_counts_metrics
#-------------------------------------------------------------------------------
STEP=raw_counts_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_1_JOB_ID: metrics.matrix
#-------------------------------------------------------------------------------
JOB_NAME=metrics.matrix
JOB_DEPENDENCIES=$raw_counts_1_JOB_ID:$raw_counts_2_JOB_ID:$raw_counts_3_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.9029dbff174b2d1daab82f48b08f8d07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.matrix.9029dbff174b2d1daab82f48b08f8d07.mugqic.done'
module load mugqic/mugqic_tools/2.1.0 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/Spt6-SA.readcounts.csv \
  raw_counts/CK2.readcounts.csv \
  raw_counts/WT.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | cut -d. -f1) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt
metrics.matrix.9029dbff174b2d1daab82f48b08f8d07.mugqic.done
)
raw_counts_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_2_JOB_ID: metrics.wigzip
#-------------------------------------------------------------------------------
JOB_NAME=metrics.wigzip
JOB_DEPENDENCIES=$wiggle_3_JOB_ID:$wiggle_4_JOB_ID:$wiggle_7_JOB_ID:$wiggle_8_JOB_ID:$wiggle_11_JOB_ID:$wiggle_12_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done'
zip -r tracks.zip tracks/bigWig
metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
)
raw_counts_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_3_JOB_ID: rpkm_saturation
#-------------------------------------------------------------------------------
JOB_NAME=rpkm_saturation
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/rpkm_saturation.334d93e95d6ad72272c1de29a52db9a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rpkm_saturation.334d93e95d6ad72272c1de29a52db9a5.mugqic.done'
module load mugqic/R_Bioconductor/3.1.2_3.0 mugqic/mugqic_tools/2.1.0 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  11 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation
rpkm_saturation.334d93e95d6ad72272c1de29a52db9a5.mugqic.done
)
raw_counts_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q lm -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_4_JOB_ID: raw_count_metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=raw_count_metrics_report
JOB_DEPENDENCIES=$raw_counts_metrics_2_JOB_ID:$raw_counts_metrics_3_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/raw_count_metrics_report.456b8a0d760b92e99a4de2b5cc6e4d0d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'raw_count_metrics_report.456b8a0d760b92e99a4de2b5cc6e4d0d.mugqic.done'
module load mugqic/pandoc/1.13.1 && \
mkdir -p report && \
cp metrics/rnaseqRep/corrMatrixSpearman.txt report/corrMatrixSpearman.tsv && \
cp tracks.zip report/ && \
cp metrics/saturation.zip report/ && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.raw_counts_metrics.md \
  --variable corr_matrix_spearman_table="`head -16 report/corrMatrixSpearman.tsv | cut -f-16| awk -F"	" '{OFS="	"; if (NR==1) {$0="Vs"$0; print; gsub(/[^	]/, "-"); print} else {printf $1; for (i=2; i<=NF; i++) {printf "	"sprintf("%.2f", $i)}; print ""}}' | sed 's/	/|/g'`" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.raw_counts_metrics.md \
  > report/RnaSeq.raw_counts_metrics.md
raw_count_metrics_report.456b8a0d760b92e99a4de2b5cc6e4d0d.mugqic.done
)
raw_counts_metrics_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cufflinks
#-------------------------------------------------------------------------------
STEP=cufflinks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cufflinks_1_JOB_ID: cufflinks.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.Spt6-SA
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.Spt6-SA.328b653d1f7facee163c03d96fef0ac4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.Spt6-SA.328b653d1f7facee163c03d96fef0ac4.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/Spt6-SA && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/Spt6-SA \
  --num-threads 12 \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam
cufflinks.Spt6-SA.328b653d1f7facee163c03d96fef0ac4.mugqic.done
)
cufflinks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_2_JOB_ID: cufflinks.CK2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.CK2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.CK2.2cdc5639b6773a12cb9ab1a84d30d680.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.CK2.2cdc5639b6773a12cb9ab1a84d30d680.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/CK2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/CK2 \
  --num-threads 12 \
  alignment/CK2/CK2.sorted.mdup.bam
cufflinks.CK2.2cdc5639b6773a12cb9ab1a84d30d680.mugqic.done
)
cufflinks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_3_JOB_ID: cufflinks.WT
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.WT
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.WT.bcfbe06951d3d0e0445614b518a47947.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.WT.bcfbe06951d3d0e0445614b518a47947.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/WT && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/WT \
  --num-threads 12 \
  alignment/WT/WT.sorted.mdup.bam
cufflinks.WT.bcfbe06951d3d0e0445614b518a47947.mugqic.done
)
cufflinks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffmerge
#-------------------------------------------------------------------------------
STEP=cuffmerge
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffmerge_1_JOB_ID: cuffmerge
#-------------------------------------------------------------------------------
JOB_NAME=cuffmerge
JOB_DEPENDENCIES=$cufflinks_1_JOB_ID:$cufflinks_2_JOB_ID:$cufflinks_3_JOB_ID
JOB_DONE=job_output/cuffmerge/cuffmerge.99c865caa1eede84c5ac75f2dd48a43b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffmerge.99c865caa1eede84c5ac75f2dd48a43b.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/AllSamples && \
`cat > cufflinks/cuffmerge.samples.txt << END
cufflinks/Spt6-SA/transcripts.gtf
cufflinks/CK2/transcripts.gtf
cufflinks/WT/transcripts.gtf
END
  
` && \
cuffmerge  \
  --ref-gtf /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf \
  --ref-sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  -o cufflinks/AllSamples \
  --num-threads 12 \
  cufflinks/cuffmerge.samples.txt
cuffmerge.99c865caa1eede84c5ac75f2dd48a43b.mugqic.done
)
cuffmerge_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffmerge_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffquant
#-------------------------------------------------------------------------------
STEP=cuffquant
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffquant_1_JOB_ID: cuffquant.Spt6-SA
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.Spt6-SA
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.Spt6-SA.26c2a07fb91bccaf2d45cb46ebb5944e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.Spt6-SA.26c2a07fb91bccaf2d45cb46ebb5944e.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/Spt6-SA && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/Spt6-SA \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/Spt6-SA/Spt6-SA.sorted.mdup.bam
cuffquant.Spt6-SA.26c2a07fb91bccaf2d45cb46ebb5944e.mugqic.done
)
cuffquant_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_2_JOB_ID: cuffquant.CK2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.CK2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.CK2.2372d1e003ad2558c545967381fd8c9f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.CK2.2372d1e003ad2558c545967381fd8c9f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/CK2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/CK2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/CK2/CK2.sorted.mdup.bam
cuffquant.CK2.2372d1e003ad2558c545967381fd8c9f.mugqic.done
)
cuffquant_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_3_JOB_ID: cuffquant.WT
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.WT
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.WT.43f445dd35884734c06f2a1ccda9f297.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.WT.43f445dd35884734c06f2a1ccda9f297.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/WT && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/WT \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/WT/WT.sorted.mdup.bam
cuffquant.WT.43f445dd35884734c06f2a1ccda9f297.mugqic.done
)
cuffquant_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffdiff
#-------------------------------------------------------------------------------
STEP=cuffdiff
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffdiff_1_JOB_ID: cuffdiff.Spt6-WT
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Spt6-WT
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Spt6-WT.f57d12bd935a57840d884d9fa12e907a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Spt6-WT.f57d12bd935a57840d884d9fa12e907a.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Spt6-WT && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Spt6-WT \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/WT/abundances.cxb \
  cufflinks/Spt6-SA/abundances.cxb
cuffdiff.Spt6-WT.f57d12bd935a57840d884d9fa12e907a.mugqic.done
)
cuffdiff_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_2_JOB_ID: cuffdiff.CK2-WT
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.CK2-WT
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.CK2-WT.d396cb84d388696bd1556dd4b13f942b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.CK2-WT.d396cb84d388696bd1556dd4b13f942b.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/CK2-WT && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/CK2-WT \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/WT/abundances.cxb \
  cufflinks/CK2/abundances.cxb
cuffdiff.CK2-WT.d396cb84d388696bd1556dd4b13f942b.mugqic.done
)
cuffdiff_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_3_JOB_ID: cuffdiff.Spt6-CK2
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Spt6-CK2
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Spt6-CK2.c4c1c62d7e84c0f27e688765505e0ea1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Spt6-CK2.c4c1c62d7e84c0f27e688765505e0ea1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Spt6-CK2 && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/Saccharomyces_cerevisiae.R64-1-1.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Spt6-CK2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/CK2/abundances.cxb \
  cufflinks/Spt6-SA/abundances.cxb
cuffdiff.Spt6-CK2.c4c1c62d7e84c0f27e688765505e0ea1.mugqic.done
)
cuffdiff_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffnorm
#-------------------------------------------------------------------------------
STEP=cuffnorm
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffnorm_1_JOB_ID: cuffnorm
#-------------------------------------------------------------------------------
JOB_NAME=cuffnorm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffnorm/cuffnorm.7bf37cc3b8a05f1ea72cbef62015c1b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm.7bf37cc3b8a05f1ea72cbef62015c1b1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffnorm && \
cuffnorm -q  \
  --library-type fr-firststrand \
  --output-dir cuffnorm \
  --num-threads 12 \
  --labels Spt6-SA,CK2,WT \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/Spt6-SA/abundances.cxb \
  cufflinks/CK2/abundances.cxb \
  cufflinks/WT/abundances.cxb
cuffnorm.7bf37cc3b8a05f1ea72cbef62015c1b1.mugqic.done
)
cuffnorm_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffnorm_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: fpkm_correlation_matrix
#-------------------------------------------------------------------------------
STEP=fpkm_correlation_matrix
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: fpkm_correlation_matrix_1_JOB_ID: fpkm_correlation_matrix_transcript
#-------------------------------------------------------------------------------
JOB_NAME=fpkm_correlation_matrix_transcript
JOB_DEPENDENCIES=$cuffnorm_1_JOB_ID
JOB_DONE=job_output/fpkm_correlation_matrix/fpkm_correlation_matrix_transcript.4ab8e237f03fba4552cb540d82d3b810.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'fpkm_correlation_matrix_transcript.4ab8e237f03fba4552cb540d82d3b810.mugqic.done'
module load mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p metrics && \
R --no-save --no-restore <<-EOF
dataFile=read.table("cuffnorm/isoforms.fpkm_table",header=T,check.names=F)
fpkm=cbind(dataFile[,2:ncol(dataFile)])
corTable=cor(log2(fpkm+0.1))
corTableOut=rbind(c('Vs.',colnames(corTable)),cbind(rownames(corTable),round(corTable,3)))
write.table(corTableOut,file="metrics/transcripts_fpkm_correlation_matrix.tsv",col.names=F,row.names=F,sep="	",quote=F)
print("done.")

EOF
fpkm_correlation_matrix_transcript.4ab8e237f03fba4552cb540d82d3b810.mugqic.done
)
fpkm_correlation_matrix_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$fpkm_correlation_matrix_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: fpkm_correlation_matrix_2_JOB_ID: fpkm_correlation_matrix_gene
#-------------------------------------------------------------------------------
JOB_NAME=fpkm_correlation_matrix_gene
JOB_DEPENDENCIES=$cuffnorm_1_JOB_ID
JOB_DONE=job_output/fpkm_correlation_matrix/fpkm_correlation_matrix_gene.99b6f4694190d208bddff755684057f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'fpkm_correlation_matrix_gene.99b6f4694190d208bddff755684057f9.mugqic.done'
module load mugqic/R_Bioconductor/3.1.2_3.0 && \
R --no-save --no-restore <<-EOF
dataFile=read.table("cuffnorm/genes.fpkm_table",header=T,check.names=F)
fpkm=cbind(dataFile[,2:ncol(dataFile)])
corTable=cor(log2(fpkm+0.1))
corTableOut=rbind(c('Vs.',colnames(corTable)),cbind(rownames(corTable),round(corTable,3)))
write.table(corTableOut,file="metrics/gene_fpkm_correlation_matrix.tsv",col.names=F,row.names=F,sep="	",quote=F)
print("done.")

EOF
fpkm_correlation_matrix_gene.99b6f4694190d208bddff755684057f9.mugqic.done
)
fpkm_correlation_matrix_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$fpkm_correlation_matrix_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: gq_seq_utils_exploratory_analysis_rnaseq
#-------------------------------------------------------------------------------
STEP=gq_seq_utils_exploratory_analysis_rnaseq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID:$cuffnorm_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/gq_seq_utils_exploratory_analysis_rnaseq.08bf3d37c7f34563166cd6dd30feb973.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'gq_seq_utils_exploratory_analysis_rnaseq.08bf3d37c7f34563166cd6dd30feb973.mugqic.done'
module load mugqic/R_Bioconductor/3.1.2_3.0 mugqic/mugqic_R_packages/1.0.1 && \
mkdir -p exploratory && \
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryAnalysisRNAseq(htseq.counts.path="DGE/rawCountMatrix.csv", cuffnorm.fpkms.dir="cuffnorm", genes.path="/cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.genes.tsv", output.dir="exploratory")
desc = readRDS(file.path("exploratory","index.RData"))
write.table(desc,file=file.path("exploratory","index.tsv"),sep='	',quote=F,col.names=T,row.names=F)
print("done.")

EOF
gq_seq_utils_exploratory_analysis_rnaseq.08bf3d37c7f34563166cd6dd30feb973.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=00:30:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_report
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_report
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/gq_seq_utils_exploratory_analysis_rnaseq_report.257a1030c1bc8d70845b805daa071e68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'gq_seq_utils_exploratory_analysis_rnaseq_report.257a1030c1bc8d70845b805daa071e68.mugqic.done'
mkdir -p report && \
cp -r exploratory/ report/ && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.md \
  report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.md && \
cut -f3-4 exploratory/index.tsv | sed '1d' | perl -pe 's/^([^	]*)	(.*)$/* [\1](\2)/' \
  >> report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.md
gq_seq_utils_exploratory_analysis_rnaseq_report.257a1030c1bc8d70845b805daa071e68.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID: cuffnorm_report
#-------------------------------------------------------------------------------
JOB_NAME=cuffnorm_report
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/cuffnorm_report.7181548353ee5f5cefe324c01650d400.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm_report.7181548353ee5f5cefe324c01650d400.mugqic.done'
mkdir -p report && \
zip -r report/cuffAnalysis.zip cufflinks/ cuffdiff/ cuffnorm/ && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.cuffnorm.md \
  report/RnaSeq.cuffnorm.md
cuffnorm_report.7181548353ee5f5cefe324c01650d400.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: differential_expression
#-------------------------------------------------------------------------------
STEP=differential_expression
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: differential_expression_1_JOB_ID: differential_expression
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/differential_expression/differential_expression.ef0184f11fab0a0e29bfc034dedcf5a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression.ef0184f11fab0a0e29bfc034dedcf5a9.mugqic.done'
module load mugqic/mugqic_tools/2.1.0 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p DGE && \
Rscript $R_TOOLS/edger.R \
  -d ../input/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE && \
Rscript $R_TOOLS/deseq.R \
  -d ../input/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE
differential_expression.ef0184f11fab0a0e29bfc034dedcf5a9.mugqic.done
)
differential_expression_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: differential_expression_goseq
#-------------------------------------------------------------------------------
STEP=differential_expression_goseq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: differential_expression_goseq_1_JOB_ID: differential_expression_goseq.dge.Spt6-WT
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq.dge.Spt6-WT
JOB_DEPENDENCIES=$differential_expression_1_JOB_ID
JOB_DONE=job_output/differential_expression_goseq/differential_expression_goseq.dge.Spt6-WT.c0520169ea91728b4a7cbff922d2a608.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression_goseq.dge.Spt6-WT.c0520169ea91728b4a7cbff922d2a608.mugqic.done'
module load mugqic/mugqic_tools/2.1.0 mugqic/R_Bioconductor/3.1.2_3.0 && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.genes.length.tsv \
  -G /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.GO.tsv \
  -d DGE/Spt6-WT/dge_results.csv \
  -c 1,6 \
  -o DGE/Spt6-WT/gene_ontology_results.csv
differential_expression_goseq.dge.Spt6-WT.c0520169ea91728b4a7cbff922d2a608.mugqic.done
)
differential_expression_goseq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_goseq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: differential_expression_goseq_2_JOB_ID: differential_expression_goseq.dge.CK2-WT
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq.dge.CK2-WT
JOB_DEPENDENCIES=$differential_expression_1_JOB_ID
JOB_DONE=job_output/differential_expression_goseq/differential_expression_goseq.dge.CK2-WT.82ba8b39009b8808744aa0bb37d5b70c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression_goseq.dge.CK2-WT.82ba8b39009b8808744aa0bb37d5b70c.mugqic.done'
module load mugqic/mugqic_tools/2.1.0 mugqic/R_Bioconductor/3.1.2_3.0 && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.genes.length.tsv \
  -G /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.GO.tsv \
  -d DGE/CK2-WT/dge_results.csv \
  -c 1,6 \
  -o DGE/CK2-WT/gene_ontology_results.csv
differential_expression_goseq.dge.CK2-WT.82ba8b39009b8808744aa0bb37d5b70c.mugqic.done
)
differential_expression_goseq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_goseq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: differential_expression_goseq_3_JOB_ID: differential_expression_goseq.dge.Spt6-CK2
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq.dge.Spt6-CK2
JOB_DEPENDENCIES=$differential_expression_1_JOB_ID
JOB_DONE=job_output/differential_expression_goseq/differential_expression_goseq.dge.Spt6-CK2.962084481d80b0ab4a62e7fa3bf9f640.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression_goseq.dge.Spt6-CK2.962084481d80b0ab4a62e7fa3bf9f640.mugqic.done'
module load mugqic/mugqic_tools/2.1.0 mugqic/R_Bioconductor/3.1.2_3.0 && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.genes.length.tsv \
  -G /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/annotations/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.GO.tsv \
  -d DGE/Spt6-CK2/dge_results.csv \
  -c 1,6 \
  -o DGE/Spt6-CK2/gene_ontology_results.csv
differential_expression_goseq.dge.Spt6-CK2.962084481d80b0ab4a62e7fa3bf9f640.mugqic.done
)
differential_expression_goseq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_goseq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: differential_expression_goseq_4_JOB_ID: differential_expression_goseq_report
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression_goseq_report
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID:$cuffdiff_1_JOB_ID:$cuffdiff_2_JOB_ID:$cuffdiff_3_JOB_ID:$differential_expression_1_JOB_ID:$differential_expression_goseq_1_JOB_ID:$differential_expression_goseq_2_JOB_ID:$differential_expression_goseq_3_JOB_ID
JOB_DONE=job_output/differential_expression_goseq/differential_expression_goseq_report.297d0c9ec1418472106e85c75f7f8243.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression_goseq_report.297d0c9ec1418472106e85c75f7f8243.mugqic.done'
module load mugqic/pandoc/1.13.1 && \
mkdir -p report && \
cp /gs/project/eav-760-aa/RNA-CKII-Spt6/input/design.txt report/design.tsv && \
cp DGE/rawCountMatrix.csv report/ && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.differential_expression.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.1/bfx/report/RnaSeq.differential_expression.md \
  --variable design_table="`head -7 report/design.tsv | cut -f-8 | awk -F"	" '{OFS="	"; if (NR==1) {print; gsub(/[^	]/, "-")} print}' | sed 's/	/|/g'`" \
  --variable raw_count_matrix_table="`head -7 report/rawCountMatrix.csv | cut -f-8 | awk -F"	" '{OFS="	"; if (NR==1) {print; gsub(/[^	]/, "-")} print}' | sed 's/	/|/g'`" \
  --to markdown \
  > report/RnaSeq.differential_expression.md && \
for contrast in Spt6-WT CK2-WT Spt6-CK2
do
  mkdir -p report/DiffExp/$contrast/
  echo -e "\n#### $contrast Results\n" >> report/RnaSeq.differential_expression.md
  cp DGE/$contrast/dge_results.csv report/DiffExp/$contrast/${contrast}_Genes_DE_results.tsv
  echo -e "\nTable: Differential Gene Expression Results (**partial table**; [download full table](DiffExp/$contrast/${contrast}_Genes_DE_results.tsv))\n" >> report/RnaSeq.differential_expression.md
  head -7 report/DiffExp/$contrast/${contrast}_Genes_DE_results.tsv | cut -f-8 | sed '2i ---	---	---	---	---	---	---	---' | sed 's/	/|/g' >> report/RnaSeq.differential_expression.md
  sed '1s/^tracking_id/test_id/' cuffdiff/$contrast/isoforms.fpkm_tracking | awk -F"	" 'FNR==NR{line[$1]=$0; next}{OFS="	"; print line[$1], $0}' - cuffdiff/$contrast/isoform_exp.diff | python -c 'import csv,sys; rows_in = csv.DictReader(sys.stdin, delimiter="	"); rows_out = csv.DictWriter(sys.stdout, fieldnames=["test_id", "gene_id", "tss_id","nearest_ref_id","class_code","gene","locus","length","log2(fold_change)","test_stat","p_value","q_value"], delimiter="	", extrasaction="ignore"); rows_out.writeheader(); rows_out.writerows(rows_in)' > report/DiffExp/$contrast/${contrast}_Transcripts_DE_results.tsv
  echo -e "\n---\n\nTable: Differential Transcript Expression Results (**partial table**; [download full table](DiffExp/$contrast/${contrast}_Transcripts_DE_results.tsv))\n" >> report/RnaSeq.differential_expression.md
  head -7 report/DiffExp/$contrast/${contrast}_Transcripts_DE_results.tsv | cut -f-8 | sed '2i ---	---	---	---	---	---	---	---' | sed 's/	/|/g' >> report/RnaSeq.differential_expression.md
  if [ `wc -l DGE/$contrast/gene_ontology_results.csv | cut -f1 -d\ ` -gt 1 ]
  then
    cp DGE/$contrast/gene_ontology_results.csv report/DiffExp/$contrast/${contrast}_Genes_GO_results.tsv
    echo -e "\n---\n\nTable: GO Results of the Differentially Expressed Genes (**partial table**; [download full table](DiffExp/${contrast}/${contrast}_Genes_GO_results.tsv))\n" >> report/RnaSeq.differential_expression.md
    head -7 report/DiffExp/${contrast}/${contrast}_Genes_GO_results.tsv | cut -f-8 | sed '2i ---	---	---	---	---	---	---	---' | sed 's/	/|/g' >> report/RnaSeq.differential_expression.md
  else
    echo -e "\nNo FDR adjusted GO enrichment was significant (p-value too high) based on the differentially expressed gene results for this design.\n" >> report/RnaSeq.differential_expression.md
  fi
done
differential_expression_goseq_report.297d0c9ec1418472106e85c75f7f8243.mugqic.done
)
differential_expression_goseq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_goseq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n04&ip=10.241.129.14&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates,picard_rna_metrics,estimate_ribosomal_rna,rnaseqc,wiggle,raw_counts,raw_counts_metrics,cufflinks,cuffmerge,cuffquant,cuffdiff,cuffnorm,fpkm_correlation_matrix,gq_seq_utils_exploratory_analysis_rnaseq,differential_expression,differential_expression_goseq&samples=3" --quiet --output-document=/dev/null

