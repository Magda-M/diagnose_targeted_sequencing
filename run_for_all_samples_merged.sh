#!/bin/bash

samples_array=( "MARCEL1" "MARCEL2" "MARCEL3" "MARCEL4" )
for sample in "${samples_array[@]}"
   do
   echo $sample;
   bash run_diagnostic_pipeline.sh $sample &> /mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/merged/"$sample".log;
   done


