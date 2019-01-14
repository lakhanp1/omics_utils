#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N realign_to_bqsr
#$ -j y
#$ -A cparsania
#$ -q fhs.q
#$ -pe orte 36 
#$ -m bes
#$ -M cparsania@umac.mo

. /opt/gridengine/default/common/settings.sh
. /etc/profile.d/modules.sh

SCRIPT_HOME=`pwd`

module add java/jdk1.8.0_25
module add gatk/3.3

#java -jar /share/apps/gatk/gcc/bin/GenomeAnalysisTK.jar -h

dedupBamFiles=(`find . -name "*_dedup_reads.bam"`)
ref="/home/cparsania/Database/human/b37/genomes/human_g1k_v37_renamed.fasta"

for (( i=0; i<${#dedupBamFiles[@]}; i++ ))

do
 
echo "**************************************************************************************************Create Interval For Realignment**************************************************************************************************"    

java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-nt $NSLOTS \
-T RealignerTargetCreator \
-I ${dedupBamFiles[i]} \
-R ${ref} \
-o ${dedupBamFiles[i]}.intervals \
-known:myvcf,VCF /home/cparsania/Database/human/b37/VCFs/Mills_and_1000G_gold_standard.indels.b37.vcf \
  
  
echo "**************************************************************************************************Interval For Realignment Created Successfully **************************************************************************************************"    
echo "**************************************************************************************************Interval For Realignment Created Successfully **************************************************************************************************"    

echo "**************************************************************************************************Indel Realignemt Near Interval Startred **************************************************************************************************" 
  
  
java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T IndelRealigner \
-I ${dedupBamFiles[i]} \
-R ${ref} \
-targetIntervals ${dedupBamFiles[i]}.intervals \
-known:myvcf,VCF /home/cparsania/Database/human/b37/VCFs/Mills_and_1000G_gold_standard.indels.b37.vcf \
-o ${dedupBamFiles[i]}_realigned.bam

echo "**************************************************************************************************Indel Realignemt Near Interval Finished **************************************************************************************************" 
echo "**************************************************************************************************Indel Realignemt Near Interval Finished **************************************************************************************************" 

echo "**************************************************************************************************Analyze patterns of covariation in the sequence dataset*****************************************************************"

java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R ${ref} \
-I ${dedupBamFiles[i]}_realigned.bam \
-knownSites /home/cparsania/Database/human/b37/VCFs/dbsnp_138.b37.vcf \
-knownSites /home/cparsania/Database/human/b37/VCFs/Mills_and_1000G_gold_standard.indels.b37.vcf \
-o ${dedupBamFiles[i]}_recal_data.table

echo "**************************************************************************************************Analyze patterns of covariation in the sequence dataset Finished*****************************************************************"
echo "**************************************************************************************************Analyze patterns of covariation in the sequence dataset Finished*****************************************************************"

echo "**************************************************************************************************Do a second pass to analyze covariation remaining after recalibration*****************************************************************"
 
java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R ${ref} \
-I ${dedupBamFiles[i]}_realigned.bam \
-knownSites /home/cparsania/Database/human/b37/VCFs/dbsnp_138.b37.vcf \
-knownSites /home/cparsania/Database/human/b37/VCFs/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR ${dedupBamFiles[i]}_recal_data.table \
-o ${dedupBamFiles[i]}_post_recal_data.table

echo "**************************************************************************************************Do a second pass to analyze covariation remaining after recalibration Finished*****************************************************************"
echo "**************************************************************************************************Do a second pass to analyze covariation remaining after recalibration Finished*****************************************************************"

echo "*************************************************************************************************Generate before/after plots*****************************************************************"
 
java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R ${ref} \
-before ${dedupBamFiles[i]}_recal_data.table \
-after ${dedupBamFiles[i]}_post_recal_data.table \
-plots recalibration_plots.pdf
 
 echo "*************************************************************************************************Apply the recalibration to your sequence data*****************************************************************"
 
java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T PrintReads \
-R ${ref} \
-I ${dedupBamFiles[i]}_realigned.bam \
-BQSR ${dedupBamFiles[i]}_recal_data.table \
-o ${dedupBamFiles[i]}_recal_reads.bam
 
 
done

                                                                                                                                                                                                                                                                                                            157,1         80%
 echo "*************************************************************************************************All Finished *****************************************************************"




