## mapping stats
samtools index SAMPLE_ID_bt2.bam
samtools flagstat SAMPLE_ID_bt2.bam > alignment.stats

mappedReads=`grep -P ' 0 mapped \(' alignment.stats | grep -P -o '^\d+'`
scale=`perl -e "printf('%.3f', 1000000/$mappedReads)"`

##macs2 pileup with 200bp extension
macs2 pileup --extsize 200 -i SAMPLE_ID_bt2.bam -o SAMPLE_ID_pileup.bdg -f BAM
error_exit $?

##normalize
printf "Normalizing SAMPLE_ID_pileup.bdg with factor %s\n" $scale
macs2 bdgopt -i SAMPLE_ID_pileup.bdg -m multiply -p $scale -o temp_normalized.bdg
error_exit $?

##Remove the first line
sed -n '2,$p' temp_normalized.bdg > SAMPLE_ID_normalized.bdg
error_exit $?
rm temp_normalized.bdg

##bedSort and bedGraph to bigWig conversion
bedSort SAMPLE_ID_normalized.bdg SAMPLE_ID_normalized.bdg
bedGraphToBigWig SAMPLE_ID_normalized.bdg ${conf_chrSize} SAMPLE_ID_normalized.bw
error_exit $?

##Pol-II expression value
if [ -f "${conf_polIIFeatures}" ]; then
	perl $CODE_DIR/miaoScripts/zqWinSGR-v2.pl -feature_file ${conf_polIIFeatures} -socre_file SAMPLE_ID_normalized.bdg -chrom_column 1 -start_column 2 -end_column 3  -direction_column 6 -bin_count 1 -output_folder $PWD -outout_name SAMPLE_ID_polii_expr.tab
	error_exit $?
fi

gzip *.bdg

