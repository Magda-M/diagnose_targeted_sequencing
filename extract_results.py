import sys
import pandas as pd
import numpy as np


FASTQ_PATH="/mnt/chr11/Data/magda/Powroty/panel/fastq/"
RESULTS_PATH="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"
PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"
PADDED_CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded.bed"

RESULTS_CSV='../results/summary.csv'

def add_column_and_value(results, existing_columns, new_sample_series, column_name, value):
	if column_name not in existing_columns:
		results = results.assign( **{column_name : pd.Series( len(results) * [np.nan], index = results.index)})
		existing_columns.append(column_name)
	new_sample_series[column_name] = value	

def extract_results(sample_name):
	try:
		results = pd.read_csv(RESULTS_CSV, sep='\t')
	except:
		results = pd.DataFrame({'sample_name' : []})
		
	existing_columns = list(results.columns)

	LOG_FILE='../results/' + sample_name + '.log'

	#create new row for the sample
	values = [sample_name]
	values.extend((len(results.columns) - 1) * [np.nan])
	new_sample_series = pd.Series(values, index=existing_columns)

	#% of reads kept after trimming
	with open(LOG_FILE, 'r') as log:
		trimmomatic_line = log.read().split('Both Surviving:')[1].split('Forward Only Surviving:')[0]
		both_surviving_num = trimmomatic_line.split(' ')[1]
		both_surviving_percent = trimmomatic_line.split('(')[1].split('%)')[0]
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving num", both_surviving_num)
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving perc", both_surviving_percent)

	#duplicates metrics
	with open('../results/' + sample_name + '_picard_markduplicates_metrics.txt', 'r') as markduplicates:
		markduplicates_results_line = markduplicates.readlines()[7]
		perc_duplication = float(markduplicates_results_line.split('\t')[-2]) * 100
	add_column_and_value(results, existing_columns, new_sample_series, "Percent duplication", perc_duplication)	

	#mapping metrics from flagstat
	with open('../results/' + sample_name + '_samtools_flagstat_metrics.txt', 'r') as flagstat:
		flagstat_lines = flagstat.readlines()
		flagstat_results = {}
		for line in flagstat_lines:
			field_name = line.split(' ')[3].strip()
			if field_name == "secondary":
				flagstat_results['secondary_mapping'] = line.split(' ')[0]
			if field_name == "duplicates":
				flagstat_results['duplicates'] = line.split(' ')[0]
			if field_name == "mapped":
				flagstat_results['mapped'] = line.split(' ')[0]
				flagstat_results['perc_mapped'] = line.split('(')[1].split('%')[0]
		flagstat_results['total_qc_passed_reads'] = flagstat_lines[0].split(' ')[0]
	add_column_and_value(results, existing_columns, new_sample_series, "Flagstat total qc passed reads", flagstat_results['total_qc_passed_reads'])	
	add_column_and_value(results, existing_columns, new_sample_series, "Flagstat mapped", flagstat_results['mapped'])	
	add_column_and_value(results, existing_columns, new_sample_series, "Flagstat perc mapped", flagstat_results['perc_mapped'])
	add_column_and_value(results, existing_columns, new_sample_series, "Flagstat secondary mapping", flagstat_results['secondary_mapping'])
	add_column_and_value(results, existing_columns, new_sample_series, "Flagstat duplicates", flagstat_results['duplicates'])		

	#insert size metrics
	with open('../results/' + sample_name + '_picard_insert_size_metrics.txt', 'r') as picard_insert_size:
		metrics = picard_insert_size.readlines()[7].strip().split('\t')
		median = metrics[0]
		mean = metrics[5]
	add_column_and_value(results, existing_columns, new_sample_series, "Insert size median", median)
	add_column_and_value(results, existing_columns, new_sample_series, "Insert size mean", mean)

	#on target rate
	primary_target_count = open('../results/' + sample_name + '_primary_target_count', 'r').read().strip()
	capture_target_count = open('../results/' + sample_name + '_capture_target_count', 'r').read().strip()
	padded_target_count = open('../results/' + sample_name + '_padded_target_count', 'r').read().strip()

	add_column_and_value(results, existing_columns, new_sample_series, "Perc on primary target", float(primary_target_count)/float(flagstat_results['mapped']))
	add_column_and_value(results, existing_columns, new_sample_series, "Perc on capture target", float(capture_target_count)/float(flagstat_results['mapped']))
	add_column_and_value(results, existing_columns, new_sample_series, "Perc on padded target", float(padded_target_count)/float(flagstat_results['mapped']))

	results = results.append(new_sample_series, ignore_index=True)
	results = results.reindex_axis(existing_columns, axis=1)
	results.to_csv(RESULTS_CSV, sep='\t', index=False)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python %s sample_name" % sys.argv[0])
		sys.exit(1)
	sample = sys.argv[1]
	extract_results(sample)