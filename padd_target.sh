#!/bin/bash

bedtools slop -i /mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed -b 100 -g /mnt/chr11/Data/magda/Powroty/hg38.chrom.sizes > /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded_unmerged_unsorted.bed

bedtools sort -i /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded_unmerged_unsorted.bed > /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded_unmerged_sorted.bed

bedtools merge -i /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded_unmerged_sorted.bed > /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded.bed