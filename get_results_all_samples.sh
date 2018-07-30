#!/bin/bash

for l in /mnt/chr11/Data/magda/Powroty/panel/fastq/MARCEL*R1_001.fastq;
   do path=$(echo $l | cut -d/ -f9); 
   sample=$(echo $path | cut -d_ -f1-3);
   echo $sample;
   #must remove beginning of log because logs from two runs were contatenated
   mv /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"$sample".log /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"$sample".log.old;
   tail -n +112063 /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"$sample".log.old > /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"$sample".log;
   rm /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"$sample".log.old;
   python extract_results.py $sample 
   done