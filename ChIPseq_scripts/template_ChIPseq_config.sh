#!/bin/bash

## This script prints the reference file variable for the organism of interest
## it uses 'yq' to read the YAML config file

usage="Usage: bash template_ChIP_process.sh -o orgID -c reference_config.yaml
  -o, --org         CHR: Organism ID
  -c, --conf        FILE: reference data YAML config file
  --polII           FLAG: Whether to run polII signal calculation script
  -h, --help        This small usage guide
"

## read command line arguments
PARSED_OPT=$(getopt -o ho:c: --long "help,org:,conf:,polII"  -- "$@")

if [ $? != "0" ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi
if [ $# -eq 0 ] ; then printf "Error: No arguments\n${usage}" >&2 ; exit 1 ; fi

eval set -- "$PARSED_OPT"

POLII=false

while true ; do
	case "$1" in
	-h|--help)
		printf "${usage}" >&2; exit 1;;
	-o|--org)
		orgId=$2; shift 2;;
	-c|--conf)
		file_config=$2; shift 2;;
	--polII)
		POLII=true; shift;;
	--) shift ; break ;;
	*) echo "Internal error!" >&2; exit 1;;
	esac
done

if [ -z "$orgId" ]; then
	printf "Error: Missing --org argument\n" 1>&2
	exit 1
fi

if [ -z "$file_config" ]; then
	printf "Error: Missing --conf argument\n" 1>&2
	exit 1
fi

## check if config file exists
if [ ! -f "${file_config}" ]; then
	printf "Error: reference config file does not exist...\n" 1>&2
	exit 1
fi

## get reference variables for organism of interest
chrSize=$(yq eval .${orgId}.chrSize ${file_config})
file_polIIFeatures=$(yq eval .${orgId}.polII_cds ${file_config})


if [ $chrSize == "null" ]; then
	printf "Error: ${orgId}.chrSize not found...\n" 1>&2
	exit 1
fi


if [ $POLII == "true" ] && [ $file_polIIFeatures == "null" ]; then
	printf "Error: ${orgId}.polII_cds not found...\n" 1>&2
	exit 1
fi


## print the variables
printf "##reference config\n"

## basic files for ChIPseq data processing
printf "conf_orgId=\"${orgId}\"\n"
printf "conf_chrSize=\"${chrSize}\"\n"

## BED file for polII signal calculation using Miao' script
if [ $POLII == "true" ]; then
	printf "conf_polIIFeatures=\"${file_polIIFeatures}\"\n"
fi

printf "##\n\n"



