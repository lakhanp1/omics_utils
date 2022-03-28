#!/bin/bash

## RNAseq pipeline script: hisat2 maping -> samtools indexing -> stringtie

set -e
set -u
set -o pipefail

##----------------------------------------------------------------
## argument parsing
usage="
Usage: bash RNAseq_process.sh -i /path/to/HiSAT2/index -g /path/to/annotation.gtf -c sampleInfo.tab
-c, --conf        FILE: Tabular file with three columns: <sampleId> <R1.fastq.gz> <R2.fastq.gz>
-g, --gtf         FILE: GTF annotation file
-i, --index       CHR: HiSat index prefix
-h, --help        This small usage guide
"

## read command line arguments
PARSED_OPT=$(getopt -o hi:g:c: --long "help,index:,gtf:,conf:"  -- "$@")

if [ $? != "0" ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi
if [ $# -eq 0 ] ; then printf "Error: No arguments\n${usage}" >&2 ; exit 1 ; fi

eval set -- "$PARSED_OPT"

while true ; do
	case "$1" in
	-h|--help)
		printf "${usage}" >&2; exit 1;;
	-i|--index)
		index=$2; shift 2;;
    -g|--gtf)
    gtf=$2; shift 2;;
	-c|--conf)
		file_config=$2; shift 2;;
	--) shift ; break ;;
	*) echo "Internal error!" >&2; exit 1;;
	esac
done


if [ -z "$index" ]; then
	printf "Error: Missing --index argument\n" 1>&2
	exit 1
fi

if [ -z "$gtf" ]; then
	printf "Error: Missing --gtf argument\n" 1>&2
	exit 1
fi

if [ -z "$file_config" ]; then
	printf "Error: Missing --conf argument\n" 1>&2
	exit 1
fi

##----------------------------------------------------------------
## check if config file exists
if [ ! -f "${file_config}" ]; then
	printf "Error: reference config file does not exist...\n" 1>&2
	exit 1
fi

## check if GTF file exists
if [ ! -f "${gtf}" ]; then
	printf "Error: GTF file does not exist...\n" 1>&2
	exit 1
fi

#the error_exit function will terminate the script if any of the command fails
function error_exit
{
	finishTime=`date "+%T %Y/%m/%d"`
	if [ $1 != "0" ]; then
		printf "Error: Failed at $finishTime\n" 1>&2
		exit 1
	else
		printf "Done... $finishTime\n\n" 1>&2
	fi
}

export -f error_exit

function process_start
{
	startTime=`date "+%T %Y/%m/%d"`
	printf "Started at $startTime: $1\n" 1>&2
}

export -f process_start

## check if tools are installed
for tool in hisat2 stringtie samtools
do
	printf "Checking installer for $tool: "
	which $tool
	error_exit $?
done

##----------------------------------------------------------------

while IFS=$'\t' read -r sampleId read1 read2 ; do

	printf "Processing sample: sampleId:%s**read1:%s**read2:%s**\n" $sampleId $read1 $read2

	## check for non empty string
	if [ -z "$read1" ]; then
		printf "Error: Provide valid R1 file name: $read1...\n" 1>&2
		exit 1
	fi

	## check for non empty string
	if [ -z "$read2" ]; then
		printf "Error: Provide valid R2 file name: *$read2*...\n" 1>&2
		exit 1
	fi

	## check if FASTQ files are present
	for fqFile in $(echo ${read1} ${read2} | tr "," "\n")
	do
		if [ ! -f "${fqFile}" ]; then
			printf "Error: File not found: %s ...\n" ${fqFile} 1>&2
			exit 1
		fi
	done

	outDir=$sampleId
	[ ! -d ${outDir} ] && mkdir ${outDir}

	## align using HiSAT2
	process_start hisat2
	hisat2 -p 4 --summary-file ${outDir}/hisat.summary  -x ${index} -1 ${read1} -2 ${read2} | \
	samtools view -bS - | \
	samtools sort  -O bam -o ${outDir}/${sampleId}_hisat2.bam

	error_exit $?

	## mapping stats
	process_start samtools_index
	samtools index ${outDir}/${sampleId}_hisat2.bam
	error_exit $?
	samtools flagstat ${outDir}/${sampleId}_hisat2.bam > ${outDir}/alignment.stats

	## run StringTie
	process_start stringtie
	stringtie ${outDir}/${sampleId}_hisat2.bam -p 4 -e -B -G ${gtf} -o ${outDir}/stringTie_${sampleId}/${sampleId}.gtf
	error_exit $?

	printf "Sample $sampleId done\n"

done < ${file_config}


##----------------------------------------------------------------


