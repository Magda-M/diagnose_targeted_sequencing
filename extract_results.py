import sys, os
import pandas as pd
import numpy as np
import datetime

PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"
PADDED_CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded.bed"


def add_column_and_value(results, existing_columns, new_sample_series, column_name, value):
	if column_name not in existing_columns:
		results = results.assign( **{column_name : pd.Series( len(results) * [np.nan], index = results.index)})
		existing_columns.append(column_name)
	new_sample_series[column_name] = value	

def extract_results(sample_name, results_path, results_csv):
	try:
		results = pd.read_csv(results_csv, sep='\t')
	except:
		results = pd.DataFrame({'sample_name' : []})
		
	existing_columns = list(results.columns)

	LOG_FILE = os.path.join(results_path, sample_name + '.log')

	#create new row for the sample
	values = [sample_name]
	values.extend((len(results.columns) - 1) * [np.nan])
	new_sample_series = pd.Series(values, index=existing_columns)

	#% of reads kept after trimming
	with open(LOG_FILE, 'r') as log:
		#for lane 1
		trimmomatic_line1 = log.read().split('Both Surviving:')[1].split('Forward Only Surviving:')[0]
		both_surviving_num1 = trimmomatic_line1.split(' ')[1]
		both_surviving_percent1 = trimmomatic_line1.split('(')[1].split('%)')[0]
		#for lane 2
		trimmomatic_line2 = log.read().split('Both Surviving:')[2].split('Forward Only Surviving:')[0]
		both_surviving_num2 = trimmomatic_line2.split(' ')[1]
		both_surviving_percent2 = trimmomatic_line2.split('(')[1].split('%)')[0]

		both_surviving_num_total = int(both_surviving_num1) + int(both_surviving_num2)

	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving num lane 1", both_surviving_num1)
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving perc lane 1", both_surviving_percent1)
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving num lane 2", both_surviving_num2)
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving perc lane 2", both_surviving_percent2)
	add_column_and_value(results, existing_columns, new_sample_series, "Trimmomatic both surviving num total", both_surviving_num_total)

	#duplicates metrics
	with open(os.path.join(results_path, sample_name + '_picard_markduplicates_metrics.txt'), 'r') as markduplicates:
		markduplicates_results_line = markduplicates.readlines()[7]
		perc_duplication = float(markduplicates_results_line.split('\t')[-2]) * 100
	add_column_and_value(results, existing_columns, new_sample_series, "Percent duplication", perc_duplication)	

	#mapping metrics from flagstat
	with open(os.path.join(results_path, sample_name + '_samtools_flagstat_metrics.txt'), 'r') as flagstat:
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
	with open(os.path.join(results_path, sample_name + '_picard_insert_size_metrics.txt'), 'r') as picard_insert_size:
		metrics = picard_insert_size.readlines()[7].strip().split('\t')
		median = metrics[0]
		mean = metrics[5]
	add_column_and_value(results, existing_columns, new_sample_series, "Insert size median", median)
	add_column_and_value(results, existing_columns, new_sample_series, "Insert size mean", mean)

	#on target rate
	primary_target_count = open(os.path.join(results_path, sample_name + '_primary_target_count'), 'r').read().strip()
	capture_target_count = open(os.path.join(results_path, sample_name + '_capture_target_count'), 'r').read().strip()
	padded_target_count = open(os.path.join(results_path, sample_name + '_padded_target_count'), 'r').read().strip()

	add_column_and_value(results, existing_columns, new_sample_series, "Perc on primary target", float(primary_target_count)/float(flagstat_results['mapped']))
	add_column_and_value(results, existing_columns, new_sample_series, "Perc on capture target", float(capture_target_count)/float(flagstat_results['mapped']))
	add_column_and_value(results, existing_columns, new_sample_series, "Perc on padded target", float(padded_target_count)/float(flagstat_results['mapped']))

    #coverage on primary and capture target
	for target_type in ['primary', 'capture']:
		gatk_coverage = pd.read_csv(os.path.join(results_path, sample_name + '_gatk_%s_target_coverage.sample_summary' % target_type), sep='\t')
		gatk_mean = gatk_coverage.iloc[0]['mean']
		gatk_3rd = gatk_coverage.iloc[0]['granular_third_quartile']
		gatk_median = gatk_coverage.iloc[0]['granular_median']
		gatk_1st = gatk_coverage.iloc[0]['granular_first_quartile']
		gatk_1 = gatk_coverage.iloc[0]['%_bases_above_1']
		gatk_10 = gatk_coverage.iloc[0]['%_bases_above_10']
		gatk_20 = gatk_coverage.iloc[0]['%_bases_above_20']

		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target coverage mean" % target_type, gatk_mean)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target coverage 3rd quart" % target_type, gatk_3rd)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target coverage median" % target_type, gatk_median)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target coverage 1st quart" % target_type, gatk_1st)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target perc bases > 1" % target_type, gatk_1)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target perc bases > 10" % target_type, gatk_10)
		add_column_and_value(results, existing_columns, new_sample_series, "GATK %s target perc bases > 20" % target_type, gatk_20)

	#HsMetrics
	picard_hsmetrics = pd.read_csv(os.path.join(results_path, sample_name + '_picard_hs_metrics.txt'),
								  skiprows=[0,1,2,3,4,5], nrows=1, sep='\t')
	for col in list(picard_hsmetrics.columns):
		metric_value = picard_hsmetrics.iloc[0][col]
		add_column_and_value(results, existing_columns, new_sample_series, col, metric_value)
			

	results = results.append(new_sample_series, ignore_index=True)
	results = results.reindex_axis(existing_columns, axis=1)
	results.to_csv(results_csv, sep='\t', index=False)

	#transpose final table
	results_updated = pd.read_csv(results_csv, sep='\t')
	results_csv_T = results_csv.replace('_horizontal', '')
	results_updated.T.to_csv(results_csv_T, sep='\t', header=False)

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python %s sample_name results_path" % sys.argv[0])
		sys.exit(1)
	sample = sys.argv[1]
	results_path = sys.argv[2]

	today_date = datetime.datetime.today()
	date = today_date.strftime('%Y%m%d')
	results_csv = os.path.join(results_path, 'summary_horizontal_%s.csv' % date)

	extract_results(sample, results_path, results_csv)