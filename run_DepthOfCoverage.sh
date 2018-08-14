#!/bin/bash

PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"

FASTQ_PATH="$1"
RESULTS_PATH="$2"

for l in "$FASTQ_PATH"/*_L001_R1_001.fastq;
   do path=$(echo ${l##*/});
   SAMPLE_NAME=$(echo $path | cut -d_ -f1-2);
   echo -e "\n$SAMPLE_NAME";


   echo -e "\n####Calculate depth of coverage"
   ##DepthOfCoverage is not present in GATK4 - GATK3 must be used
   java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /home/magda//Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_primary_target_coverage \
   -L "$PRIMARY_TARGET" \
   -ct 1 -ct 10 -ct 20 >> "$RESULTS_PATH$SAMPLE_NAME".log 2>&1

   java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /home/magda//Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_capture_target_coverage \
   -L "$CAPTURE_TARGET" \
   -ct 1 -ct 10 -ct 20 >> "$RESULTS_PATH$SAMPLE_NAME".log 2>&1

   done 