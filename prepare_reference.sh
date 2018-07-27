#!/bin/bash

#index reference for bwa
/home/magda/bwa/bwa index /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa

#index reference for GATK
java -Xmx4g -Xms4g -jar /home/magda/picard.jar CreateSequenceDictionary R=/mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa O=/mnt/chr11/Data/magda/Powroty/panel/ref/hg38.dict
samtools faidx /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa