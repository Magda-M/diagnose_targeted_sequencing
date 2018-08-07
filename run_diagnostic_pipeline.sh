#!/bin/bash

#Based on Evaluate NimbleGen SeqCap EZ Target Enrichment Data Technical Note from Roche (version March 2014)
#and code obtained from Bartosz Wojtas

SAMPLE_NAME="$1"
FASTQ_PATH="/mnt/chr11/Data/magda/Powroty/panel/fastq/"
RESULTS_PATH="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/minlen50/"
PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"
PADDED_CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded.bed"

echo "\n####Process sample:" "$SAMPLE_NAME"

#echo "\n####Unpack"
#gunzip "$FASTQ_PATH$SAMPLE_NAME"_L001_R1_001.fastq.gz
#gunzip "$FASTQ_PATH$SAMPLE_NAME"_L001_R2_001.fastq.gz
#gunzip "$FASTQ_PATH$SAMPLE_NAME"_L002_R1_001.fastq.gz
#gunzip "$FASTQ_PATH$SAMPLE_NAME"_L002_R2_001.fastq.gz

echo "\n####Trimming bad reads and getting read of adapters - lane 1"
java -Xms4g -Xmx4g -jar /home/magda/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -phred33 \
   "$FASTQ_PATH$SAMPLE_NAME"_L001_R1_001.fastq \
   "$FASTQ_PATH$SAMPLE_NAME"_L001_R2_001.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R1_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R1_unpaired_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R2_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R2_unpaired_minlen50.fastq \
   ILLUMINACLIP:/home/magda/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50

echo "\n####Trimming bad reads and getting read of adapters - lane 2"
java -Xms4g -Xmx4g -jar /home/magda/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -phred33 \
   "$FASTQ_PATH$SAMPLE_NAME"_L002_R1_001.fastq \
   "$FASTQ_PATH$SAMPLE_NAME"_L002_R2_001.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R1_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R1_unpaired_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R2_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R2_unpaired_minlen50.fastq \
   ILLUMINACLIP:/home/magda/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50

echo "\n####Aligning reads - lane 1"
/home/magda/bwa/bwa mem -R '@RG\tID:'"$SAMPLE_NAME"'_L001\tPL:illumina\tLB:NIMBLEGEN\tSM:'"$SAMPLE_NAME" \
   /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa -t 40 \
   -M "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R1_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L001_R2_trimmed_minlen50.fastq | \
   samtools view -Sb - > "$RESULTS_PATH$SAMPLE_NAME"_L001_unsorted.bam

echo "\n####Aligning reads - lane 2"
/home/magda/bwa/bwa mem -R '@RG\tID:'"$SAMPLE_NAME"'_L002\tPL:illumina\tLB:NIMBLEGEN\tSM:'"$SAMPLE_NAME" \
   /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa -t 40 \
   -M "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R1_trimmed_minlen50.fastq \
   "$FASTQ_PATH"/minlen50/"$SAMPLE_NAME"_L002_R2_trimmed_minlen50.fastq | \
   samtools view -Sb - > "$RESULTS_PATH$SAMPLE_NAME"_L002_unsorted.bam

echo "\n####Merge bam for lanes"
samtools merge "$RESULTS_PATH$SAMPLE_NAME"_unsorted.bam "$RESULTS_PATH$SAMPLE_NAME"_L001_unsorted.bam "$RESULTS_PATH$SAMPLE_NAME"_L002_unsorted.bam

echo "\n####Sort bam file"
samtools sort "$RESULTS_PATH$SAMPLE_NAME"_unsorted.bam -o "$RESULTS_PATH$SAMPLE_NAME"_sorted.bam

echo "\n####Mark duplicates"
java -Xmx10g -Xms10g -jar /home/magda/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="$RESULTS_PATH$SAMPLE_NAME"_sorted.bam OUTPUT="$RESULTS_PATH$SAMPLE_NAME".bam METRICS_FILE="$RESULTS_PATH$SAMPLE_NAME"_picard_markduplicates_metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true

echo "\n####Index bam file"
samtools index  "$RESULTS_PATH$SAMPLE_NAME".bam

echo "\n####Extract basic mapping metrics"
samtools flagstat "$RESULTS_PATH$SAMPLE_NAME".bam > "$RESULTS_PATH$SAMPLE_NAME"_samtools_flagstat_metrics.txt

echo "\n####Estimate Insert Size Distribution"
java -Xmx10g -jar /home/magda/picard.jar CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE="$RESULTS_PATH$SAMPLE_NAME"_picard_insert_size_plot.pdf INPUT="$RESULTS_PATH$SAMPLE_NAME".bam OUTPUT="$RESULTS_PATH$SAMPLE_NAME"_picard_insert_size_metrics.txt

echo "\n####Count on target rate"
bedtools intersect -bed -abam "$RESULTS_PATH$SAMPLE_NAME".bam -b "$PRIMARY_TARGET" | wc -l > "$RESULTS_PATH$SAMPLE_NAME"_primary_target_count
bedtools intersect -bed -abam  "$RESULTS_PATH$SAMPLE_NAME".bam -b "$CAPTURE_TARGET" | wc -l > "$RESULTS_PATH$SAMPLE_NAME"_capture_target_count
bedtools intersect -bed -abam "$RESULTS_PATH$SAMPLE_NAME".bam -b "$PADDED_CAPTURE_TARGET" | wc -l > "$RESULTS_PATH$SAMPLE_NAME"_padded_target_count

echo "\n####Calculate depth of coverage"
##DepthOfCoverage is not present in GATK4 - GATK3 must be used
java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_primary_target_coverage \
   -L "$PRIMARY_TARGET" -ct 1 -ct 10 -ct 20 -nt 10

java -Xmx10g -Xms10g -jar /home/magda/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R /mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa \
   -I "$RESULTS_PATH$SAMPLE_NAME".bam \
   -o "$RESULTS_PATH$SAMPLE_NAME"_gatk_capture_target_coverage \
   -L "$CAPTURE_TARGET" -ct 1 -ct 10 -ct 20 -nt 10

echo "\n###Create Picard Interval Lists"
samtools view -H "$RESULTS_PATH$SAMPLE_NAME".bam > "$RESULTS_PATH$SAMPLE_NAME"_bam_header.txt

echo "\n###Create a Picard Target Interval List Body"
cat "$PRIMARY_TARGET" | awk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > "$RESULTS_PATH$SAMPLE_NAME"_target.body.txt
cat "$CAPTURE_TARGET" | awk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > "$RESULTS_PATH$SAMPLE_NAME"_bait.body.txt

echo "\n###Concatenate to create a picard bait interval list"
cat "$RESULTS_PATH$SAMPLE_NAME"_bam_header.txt "$RESULTS_PATH$SAMPLE_NAME"_bait.body.txt > "$RESULTS_PATH$SAMPLE_NAME"_bait_intervals.txt
cat "$RESULTS_PATH$SAMPLE_NAME"_bam_header.txt "$RESULTS_PATH$SAMPLE_NAME"_target.body.txt > "$RESULTS_PATH$SAMPLE_NAME"_target_intervals.txt

echo "\n###Hybrid Selection analysis"
java -Xmx10g -Xms10g -jar /home/magda/picard.jar CollectHsMetrics  BAIT_INTERVALS="$RESULTS_PATH$SAMPLE_NAME"_bait_intervals.txt TARGET_INTERVALS="$RESULTS_PATH$SAMPLE_NAME"_target_intervals.txt INPUT="$RESULTS_PATH$SAMPLE_NAME".bam OUTPUT="$RESULTS_PATH$SAMPLE_NAME"_picard_hs_metrics.txt METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=/mnt/chr11/Data/magda/Powroty/panel/ref/hg38.fa VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
