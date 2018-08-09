#!/bin/bash

RESULTS_PATH="$1"
#RESULTS_PATH="/home/magda/Dane"  ON THE LAPTOP

##lanes separately
#for l in "$RESULTS_PATH"/Powroty/panel/fastq/MARCEL*R1_001.fastq;
#   do path=$(echo $l | cut -d/ -f9); 
#   #do path=$(echo $l | cut -d/ -f8); ON THE LAPTOP
#   sample=$(echo $path | cut -d_ -f1-3);
#   echo $sample;

   ##must remove beginning of log because logs from two runs were contatenated
   #mv "$RESULTS_PATH"/Powroty/panel/diagnose_sequencing/results/"$sample".log "$RESULTS_PATH"/Powroty/panel/diagnose_sequencing/results/"$sample".log.old;
   #tail -n +112063 "$RESULTS_PATH"/Powroty/panel/diagnose_sequencing/results/"$sample".log.old > "$RESULTS_PATH"/Powroty/panel/diagnose_sequencing/results/"$sample".log;
   #rm "$RESULTS_PATH"/Powroty/panel/diagnose_sequencing/results/"$sample".log.old;

#   python extract_results.py $sample 
#   done

#merged lanes

for l in "$RESULTS_PATH"/*.log
   do path=$(echo ${l##*/});
   sample=$(echo $path | cut -d. -f1);
   echo -e "\n$sample";
   python extract_results.py $sample $RESULTS_PATH
   done
