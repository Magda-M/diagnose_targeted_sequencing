#!/bin/bash

RESULTS_PATH="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/minlen50/"
PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"

for l in /mnt/chr11/Data/magda/Powroty/panel/fastq/MARCEL*_L001_R1_001.fastq;
   do path=$(echo $l | cut -d/ -f9); 
   SAMPLE_NAME=$(echo $path | cut -d_ -f1-2);
   echo -e "\n$SAMPLE_NAME";

   echo -e "\n####Calculate depth of coverage"
   ##DepthOfCoverage is not present in GATK4 - GATK3 must be used
   java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_primary_target_coverage \
   -L "$PRIMARY_TARGET" \
   -ct 1 -ct 10 -ct 20 >> /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/minlen50/"$SAMPLE_NAME".log 2>&1

   java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_capture_target_coverage \
   -L "$CAPTURE_TARGET" \
   -ct 1 -ct 10 -ct 20 >> /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/minlen50/"$SAMPLE_NAME".log 2>&1

   done 