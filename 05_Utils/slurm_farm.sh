#!/bin/bash

if [ $# -ne 2 ] ; then
printf "Error: incorrect arguments
Usage: slurm_farm.sh sample_name.list mpi-short\n\n" >&2
exit 1
fi

userName=`whoami`
file_samples="$1"
partition="$2"
partString=$(echo $partition | cut -c1-8)

if [ ! -f "$file_samples" ]; then
printf "Error: $file_samples does not exists...\n" 1>&2
exit 1
fi

qjobs=$(qstat -a -n -1 -u ${userName} | grep "${partString} .* Q " | wc -l)
printf "qjobs = $qjobs\n";


## submit jobs or wait for other jobs to finish
for sampleId in `cat ${file_samples}`
do
	
	while true; do
		qjobs=$(qstat -a -n -1 -u ${userName} | grep "${partString} .* Q " | wc -l)
		
		if [ "$qjobs" -lt 3 ]; then
			cd ${sampleId}
			printf "Submitting job for sample $sampleId...\n"
			sbatch -J $sampleId -p $partition generalJob.sh
			echo ""
			cd ..
			break
		else
			sleep 60
		fi
	done

done

