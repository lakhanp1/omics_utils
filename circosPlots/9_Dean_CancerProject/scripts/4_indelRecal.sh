#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N mt_variant_discovery
#$ -j y
#$ -A lakhanp
#$ -q fhs.q
#$ -pe orte 6
#$ -m bes
#$ -M lakhanp@umac.mo

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

SAMPLE=""

echo "********** Build the SNP recalibration model **********"    


java -Xmx50g -Xms50g -Djava.io.tmpdir=../tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref} \
-input ${SAMPLE}_dedup_reads.bam_recal_reads.bam_raw_variants.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
-resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${thousnedGenome} \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.recal \
-tranchesFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.tranches \
-rscriptFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP_plots.R


echo "********** Apply the desired level of recalibration to the SNPs in the call set **********" 


java -Xmx50g -Xms50g -Djava.io.tmpdir=../tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref} \
-input ${SAMPLE}_dedup_reads.bam_recal_reads.bam_raw_variants.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.recal \
-tranchesFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.tranches \
-o ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrated_snps_raw_indels.vcf


echo "********** Build the Indel recalibration model **********"

java -Xmx50g -Xms50g -Djava.io.tmpdir=../tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref} \
-input ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrated_snps_raw_indels.vcf \
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
-recalFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.recal \
-tranchesFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.tranches \
-rscriptFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_INDEL_plots.R


echo "********** Apply the desired level of recalibration to the Indels in the call set **********"


java -Xmx50g -Xms50g -Djava.io.tmpdir=../tmp -jar /share/apps/gatk/src/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref} \
-input ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.recal \
-tranchesFile ${SAMPLE}_dedup_reads.bam_recal_reads.bam_recalibrate_SNP.tranches \
-o ${SAMPLE}_dedup_reads.bam_recal_reads.bamrecalibrated_variants.vcf

echo "********** Done **********"
