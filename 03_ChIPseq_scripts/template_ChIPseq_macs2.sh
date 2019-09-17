
##macs2:START
conf_macs2_g=$(yq r /home/lakhanp/database/reference_genomes.yaml ${conf_orgId}.macs2_g)

if [ -f "${control}" ]
then
	macs2_ctrl="-c ${control}"
	macs2_outDir="macs2_withCtrl"
	macs2_name='withCtrl'
else
	macs2_ctrl=' '
	macs2_outDir="macs2_withoutCtrl"
	macs2_name='withoutCtrl'
fi

if [ "${peakType}" = "broad" ]
then
	macs2_broad=" --broad"
	macs2_outDir+="_broad"
else
	macs2_broad=" "
	macs2_outDir+="_narrow"
fi

##macs2 narrowPeak calling: 
macs2 callpeak -t SAMPLE_ID_bt2.bam --name SAMPLE_ID.${macs2_name} ${macs2_ctrl} --outdir ${macs2_outDir} -g ${conf_macs2_g} --nomodel --extsize 200 -B --SPMR ${macs2_broad}
error_exit $?

##generate fold-enrichment track
macs2 bdgcmp -t "${macs2_outDir}/SAMPLE_ID.${macs2_name}_treat_pileup.bdg" -c "${macs2_outDir}/SAMPLE_ID.${macs2_name}_control_lambda.bdg" -o "${macs2_outDir}/SAMPLE_ID.${macs2_name}.FE.bdg" -m FE

##bedSort and bedGraph to bigWig conversion
bedSort "${macs2_outDir}/SAMPLE_ID.${macs2_name}.FE.bdg" "${macs2_outDir}/SAMPLE_ID.${macs2_name}.FE.bdg"
bedGraphToBigWig "${macs2_outDir}/SAMPLE_ID.${macs2_name}.FE.bdg" ${conf_chrSize} "${macs2_outDir}/SAMPLE_ID.${macs2_name}.FE.bw"
error_exit $?
rm ${macs2_outDir}/SAMPLE_ID*.bdg

##macs2:END
