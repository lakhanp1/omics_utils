#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N bwa_preprocess
#$ -j y
#$ -A cparsania
#$ -q fhs.q
#$ -pe orte 48 
#$ -m bes
#$ -M yb57653@umac.mo

. /opt/gridengine/default/common/settings.sh
. /etc/profile.d/modules.sh

SCRIPT_HOME=`pwd`

module add java/jdk1.8.0_25
module add gatk/3.3

#java -jar /share/apps/gatk/gcc/bin/GenomeAnalysisTK.jar -h

fastqPath="/home/hiseq/160114_D00691_0040_AC7G8RANXX/Data/Intensities/BaseCalls/MiaoKai"
genome="/home/cparsania/Database/human/b37/genomes/human_g1k_v37_renamed.fasta"

r1=(`ls -a ${fastqPath}/*L004_R1*`)
r2=(`ls -a ${fastqPath}/*L004_R2*`)

for (( i=0; i<${#r1[@]}; i++ ));

do

fileName=`basename ${r1[$i]} .fastq.gz`
date=`echo ${fileName} | cut -d'_' -f1`
sample=`echo ${fileName} | cut -d'_' -f1`
mkdir `pwd`/${sample}
outDir=`pwd`/${sample}
id=`echo ${sample}`

echo "**************************************************************************************************$sample bwa started**************************************************************************************************"

/share/apps/blastall/2.2.23/bin/bwa  mem -t $NSLOTS \
${genome} \
${r1[i]} \
${r2[i]} \
-M -R "@RG\tID:${id}\tSM:${sample}\tPL:illumina\tLB:patient_data\tPU:unit1" > ${outDir}/${sample}.sam


echo "**************************************************************************************************bwa Finished **************************************************************************************************"
echo "\n\n"
echo "**************************************************************************************************Sorting Started **************************************************************************************************"

java -Xmx500g -Djava.io.tmpdir=`pwd`  -jar /share/apps/gatk/src/picard-tools-1.119/SortSam.jar \
TMP_DIR=`pwd`/tmp \
INPUT=${outDir}/${sample}.sam \
OUTPUT=${outDir}/${sample}_sorted.bam \
SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=5000000



echo "**************************************************************************************************Sorting Finished **************************************************************************************************"
echo "\n\n"
echo "**************************************************************************************************Mark Dulpicate Started **************************************************************************************************"

java -Djava.io.tmpdir=`pwd`/tmp  \
-jar /share/apps/gatk/src/picard-tools-1.119/MarkDuplicates.jar \
TMP_DIR=`pwd`/tmp \
INPUT=${outDir}/${sample}_sorted.bam \
OUTPUT=${outDir}/${sample}_dedup_reads.bam \
ASSUME_SORTED=true \
METRICS_FILE=${outDir}/${sample}_metrics.txt

echo "**************************************************************************************************Mark Duplicate Finished **************************************************************************************************"
echo "\n\n"
echo "**************************************************************************************************Get Alignment Stats **************************************************************************************************"  
 
java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/picard-tools-1.119/CollectAlignmentSummaryMetrics.jar \
INPUT=${outDir}/${sample}_sorted.bam \
REFERENCE_SEQUENCE=${genome} \
OUTPUT=${outDir}/${sample}_stats.txt

done

    
echo "************************************************************************************************** Alignment Stats Finished **************************************************************************************************"
echo "************************************************************************************************** Alignment Stats Finished **************************************************************************************************"
echo "************************************************************************************************** Alignment Stats Finished **************************************************************************************************"   