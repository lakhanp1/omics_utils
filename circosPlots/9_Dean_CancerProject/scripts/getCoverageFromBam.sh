#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N cov 
#$ -j y
#$ -A cparsania
#$ -q fhs.q
#$ -pe orte 2
#$ -m bes
#$ -M chirag.parsania@gmail.com

. /opt/gridengine/default/common/settings.sh
. /etc/profile.d/modules.sh

module add samtools/1.2

recal_bam=(`ls -a ../*/*_recal_reads.bam`)

for (( i=0; i<${#recal_bam[@]}; i++ ));

do
echo ${recal_bam[$i]}
/share/apps/samtools/1.2/bin/samtools mpileup ${recal_bam[$i]} | awk '{print $4}' | perl ~/Projects/cancer_project_20151119/20151119/hiseq_160114_D00691_0040_AC7G8RANXX/finishedSamples/mean_coverage.pl
echo ""
done