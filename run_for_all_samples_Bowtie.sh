#!/bin/bash

FASTQ_PATH="$1"
RESULTS_PATH="$2"

for l in "$FASTQ_PATH"/*_L001_R1_001.fastq;
   do path=$(echo ${l##*/});
   sample=$(echo $path | cut -d_ -f1-2);
   echo -e "\n$sample";
   bash run_diagnostic_pipeline_with_Bowtie.sh $sample $FASTQ_PATH $RESULTS_PATH &> "$RESULTS_PATH$sample".log 
   done


