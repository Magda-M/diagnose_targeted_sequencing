#!/bin/bash

for l in /mnt/chr11/Data/magda/Powroty/panel/fastq/MARCEL*_001.fastq;
   do sample=$(echo $l | cut -d_R -f1);
   echo $sample;
   done


#MARCEL4_S14_L001
#MARCEL2_S12_L002_R1_001.fastq