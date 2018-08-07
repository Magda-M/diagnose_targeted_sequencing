#!/bin/bash

for l in /mnt/chr11/Data/magda/Powroty/panel/fastq/MARCEL*_L001_R1_001.fastq;
   do path=$(echo $l | cut -d/ -f9); 
   sample=$(echo $path | cut -d_ -f1-2);
   echo -e "\n$sample";
   bash run_diagnostic_pipeline.sh $sample &> /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/minlen50/"$sample".log 
   done


#MARCEL4_S14_L001
#MARCEL2_S12_L002_R1_001.fastq
