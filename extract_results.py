import sys
import pandas as pd
import numpy as np


FASTQ_PATH="/mnt/chr11/Data/magda/Powroty/panel/fastq/"
RESULTS_PATH="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/"
PRIMARY_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_primary_targets.bed"
CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/Symfonia_v2_capture_targets.bed"
PADDED_CAPTURE_TARGET="/mnt/chr11/Data/magda/Powroty/panel/diagnose_sequencing/results/Symfonia_v2_capture_targets_padded.bed"

RESULTS_CSV='../results/summary.csv'

def extract_results(sample_name):
	try:
		results = pd.read_csv(RESULTS_CSV, sep='\t')
	except:
		results = pd.DataFrame({'sample_name' : []})
		
	existing_columns = list(results.columns)

	sample_name="TEST"
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

	if "Trimmomatic_both_surviving_num" not in existing_columns:
		results = results.assign( Trimmomatic_both_surviving_num = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Trimmomatic_both_surviving_num')
	new_sample_series["Trimmomatic_both_surviving_num"] = both_surviving_num	

	if "Trimmomatic_both_surviving_perc" not in existing_columns:
		results = results.assign( Trimmomatic_both_surviving_perc = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Trimmomatic_both_surviving_perc')
	new_sample_series["Trimmomatic_both_surviving_perc"] = both_surviving_percent

	#duplicates metrics
	with open('../results/' + sample_name + '_picard_markduplicates_metrics.txt', 'r') as markduplicates:
		markduplicates_results_line = markduplicates.readlines()[7]
		perc_duplication = float(markduplicates_results_line.split('\t')[-2]) * 100
		
	if "Percent_duplication" not in existing_columns:
		results = results.assign( Percent_duplication = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Percent_duplication')
	new_sample_series["Percent_duplication"] = perc_duplication

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
		
	if "Flagstat_total_qc_passed_reads" not in existing_columns:
		results = results.assign( Flagstat_total_qc_passed_reads = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Flagstat_total_qc_passed_reads')
	new_sample_series["Flagstat_total_qc_passed_reads"] = flagstat_results['total_qc_passed_reads']

	if "Flagstat_mapped" not in existing_columns:
		results = results.assign( Flagstat_mapped = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Flagstat_mapped')
	new_sample_series["Flagstat_mapped"] = flagstat_results['mapped']

	if "Flagstat_perc_mapped" not in existing_columns:
		results = results.assign( Flagstat_perc_mapped = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Flagstat_perc_mapped')
	new_sample_series["Flagstat_perc_mapped"] = flagstat_results['perc_mapped']

	if "Flagstat_secondary_mapping" not in existing_columns:
		results = results.assign( Flagstat_secondary_mapping = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Flagstat_secondary_mapping')
	new_sample_series["Flagstat_secondary_mapping"] = flagstat_results['secondary_mapping']

	if "Flagstat_duplicates" not in existing_columns:
		results = results.assign( Flagstat_duplicates = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Flagstat_duplicates')
	new_sample_series["Flagstat_duplicates"] = flagstat_results['duplicates']
		

	#insert size metrics
	with open('../results/' + sample_name + '_picard_insert_size_metrics.txt', 'r') as picard_insert_size:
		metrics = picard_insert_size.readlines()[7].strip().split('\t')
		median = metrics[0]
		mean = metrics[5]

	if "Insert_size_median" not in existing_columns:
		results = results.assign( Insert_size_median = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Insert_size_median')
	new_sample_series["Insert_size_median"] = median

	if "Insert_size_mean" not in existing_columns:
		results = results.assign( Insert_size_mean = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Insert_size_mean')
	new_sample_series["Insert_size_mean"] = mean

	#on target rate
	primary_target_count = open('../results/' + sample_name + '_primary_target_count', 'r').read().strip()
	capture_target_count = open('../results/' + sample_name + '_capture_target_count', 'r').read().strip()
	padded_target_count = open('../results/' + sample_name + '_padded_target_count', 'r').read().strip()

	if "Perc_on_primary_target" not in existing_columns:
		results = results.assign( Perc_on_primary_target = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Perc_on_primary_target')
	new_sample_series["Perc_on_primary_target"] = float(primary_target_count)/float(flagstat_results['mapped'])

	if "Perc_on_capture_target" not in existing_columns:
		results = results.assign( Perc_on_capture_target = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Perc_on_capture_target')
	new_sample_series["Perc_on_capture_target"] = float(capture_target_count)/float(flagstat_results['mapped'])

	if "Perc_on_padded_target" not in existing_columns:
		results = results.assign( Perc_on_padded_target = pd.Series( len(results) * [np.nan], index = results.index))
		existing_columns.append('Perc_on_padded_target')
	new_sample_series["Perc_on_padded_target"] = float(padded_target_count)/float(flagstat_results['mapped'])

	results = results.append(new_sample_series, ignore_index=True)
	results.to_csv(RESULTS_CSV, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python %s sample_name" % sys.argv[0])
        sys.exit(1)
    sample = sys.argv[1]
    extract_results(sample)