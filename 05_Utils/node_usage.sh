#/bin/bash

header="\n %-20s %-15.15s  %-15s %-10.10s %5s %5s\n"
format=" %-20s %-15.15s  %-15s %-10.10s %5d %5d\n"

printf "$header" "user" "name" "node" "partition" "cpu" "mem"
printf "%0.s=" {1..80}
printf "\n"

qstat -a -n -1 | grep 'FHS_.* R ' | sed -r 's/\s+/\t/' | cut -f1 | while read -r jobId
do
user=$(scontrol show -d job $jobId | grep -oP 'UserId=[^ ]+' | sed -r 's/\w+=//')
name=$(scontrol show -d job $jobId | grep -oP 'JobName=[^ ]+$' | sed -r 's/\w+=//')
node=$(scontrol show -d job $jobId | grep -oP ' NodeList=[^ ]+$' | sed -r 's/\s*\w+=//')
partition=$(scontrol show -d job $jobId | grep -oP ' Partition=[^ ]+' | sed -r 's/\s*\w+=//')
cpu=$(scontrol show -d job $jobId | grep -oP 'MinCPUsNode=\d+' | sed -r 's/\w+=//')
mem=$(scontrol show -d job $jobId | grep -oP 'MinMemoryNode=\d+' | sed -r 's/\w+=//')
printf "$format" $user $name $node $partition $cpu $mem
done


