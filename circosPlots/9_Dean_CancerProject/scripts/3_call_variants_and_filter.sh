#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N mt_variant_discovery
#$ -j y
#$ -A cparsania
#$ -q urgent.q
#$ -pe orte 48
#$ -m bes
#$ -M cparsania@umac.mo

. /opt/gridengine/default/common/settings.sh
. /etc/profile.d/modules.sh

SCRIPT_HOME=`pwd`

module add java/jdk1.8.0_25
module add gatk/3.3

ref=/home/cparsania/Database/human/b37/genomes/human_g1k_v37_renamed.fasta
hapmap=/home/cparsania/Database/human/b37/VCFs/hapmap_3.3.b37.vcf
omni=/home/cparsania/Database/human/b37/VCFs/1000G_omni2.5.b37.vcf
thousnedGenome=/home/cparsania/Database/human/b37/VCFs/1000G_phase1.snps.high_confidence.b37.vcf
dbsnp=/home/cparsania/Database/human/b37/VCFs/dbsnp_138.b37.vcf
millsIndel=/home/cparsania/Database/human/b37/VCFs/Mills_and_1000G_gold_standard.indels.b37.vcf

recalReadsBamFiles=(`find . -name "*recal_reads.bam"`)

for (( i=0; i<${#recalReadsBamFiles[@]}; i++ ))

do

echo "************************************************************************************************** Call variants in your sequence data*************************************************************************************************"        

java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${ref} \
-I ${recalReadsBamFiles[i]} \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
-o ${recalReadsBamFiles[i]}_raw_variants.vcf


echo "**************************************************************************************************  Build the SNP recalibration model*************************************************************************************************"    

java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref} \
-input ${recalReadsBamFiles[i]}_raw_variants.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
-resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousnedGenome} \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile ${recalReadsBamFiles[i]}_recalibrate_SNP.recal \
-tranchesFile ${recalReadsBamFiles[i]}_recalibrate_SNP.tranches \
-rscriptFile ${recalReadsBamFiles[i]}_recalibrate_SNP_plots.R


echo "**************************************************************************************************  Apply the desired level of recalibration to the SNPs in the call set *************************************************************************************************" 


java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref} \
-input ${recalReadsBamFiles[i]}_raw_variants.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile ${recalReadsBamFiles[i]}_recalibrate_SNP.recal \
-tranchesFile ${recalReadsBamFiles[i]}_recalibrate_SNP.tranches \
-o ${recalReadsBamFiles[i]}_recalibrated_snps_raw_indels.vcf


echo "**************************************************************************************************Build the Indel recalibration model**************************************************************************************************"

java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref} \
-input ${recalReadsBamFiles[i]}_recalibrated_snps_raw_indels.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 ${millsIndel} \
-an QD \
-an DP \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile ${recalReadsBamFiles[i]}_recalibrate_SNP.recal \
-tranchesFile ${recalReadsBamFiles[i]}_recalibrate_SNP.tranches \
-rscriptFile ${recalReadsBamFiles[i]}_recalibrate_INDEL_plots.R


echo "**************************************************************************************************Apply the desired level of recalibration to the Indels in the call set**************************************************************************************************"


java -Xmx50g -Xms50g -Djava.io.tmpdir=`pwd`/tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref} \
-input ${recalReadsBamFiles[i]}_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile ${recalReadsBamFiles[i]}_recalibrate_SNP.recal \
-tranchesFile ${recalReadsBamFiles[i]}_recalibrate_SNP.tranches \
-o ${recalReadsBamFiles[i]}recalibrated_variants.vcf

done

echo "**************************************************************************************************All Finished **************************************************************************************************"