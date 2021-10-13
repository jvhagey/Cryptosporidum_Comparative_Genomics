#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2021

from numpy import *
import glob
import os
import pandas as pd
from argparse import ArgumentParser
import numpy as np
from Bio import SeqIO

# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CEND = '\033[0m'
# allowing more output
pd.set_option("max_columns", 20)
pd.set_option("max_rows", None)


def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Consensus_to_matrix.py", description=""".""")
	parser.add_argument("-d", "--dir", dest="directory", action="store", required=True, help="Path for a directory where all fasta files are found. Expected to have an extention of")
	parser.add_argument("-p", "--percent", dest="percent_to_drop", action="store", required=False, help="The cutoff percent of samples with N at any given locus. Columns with a percent of Ns greater than or equal than the cutoff will be dropped. \
	Example: '-p 0.9' will drop any locus position that has greater than or equal to 90% Ns. Thus, this is a high threshold as you will only drop locus positions with 90% Ns at any given location ")
	args = parser.parse_args()
	return args

def make_df(fasta_file):
	column_names = []
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		for i in range(len(seq_record.seq)):
			chromosome = seq_record.id.replace("CpIOWA-ATCC_", "")
			name = chromosome + "." + str(i)
			column_names.append(name)
	snp_df = pd.DataFrame()
	return column_names, snp_df

def clean_list(fasta_files):
	"""This removes samples to be process"""
	Bad_list = ['41897_Cattle-masked-seq.fasta', '39640_Human-masked-seq.fasta', '42129_Cattle-masked-seq.fasta', '40457_Human-masked-seq.fasta', '40236_Calf-masked-seq.fasta', '41573_Cow-masked-seq.fasta', '42134_Cattle-masked-seq.fasta', '41891_Cattle-masked-seq.fasta', '41892_Cattle-masked-seq.fasta']
	for sample in Bad_list:
		if sample in fasta_files:
			fasta_files.remove(sample)
	print(CYEL + "Running the following samples {}.\n".format(fasta_files))
	return fasta_files

def add_fasta_to_dataframe():
	fasta_files = glob.glob("*-masked-seq.fasta")
	fasta_files = clean_list(fasta_files)
	position_names, snp_df = make_df(fasta_files[0])
	file_names = []
	count = 0
	for fasta_file in fasta_files:
		seq_list = []
		file_name = fasta_file.replace("-masked-seq.fasta", "")
		file_names.append(file_name)
		for seq_record in SeqIO.parse(fasta_file, "fasta"):
			seq_list.append(str(seq_record.seq))
		cat_seq = ''.join(seq_list)
		snp_df = pd.concat([snp_df, pd.DataFrame(list(cat_seq))], ignore_index=True, axis = 1)
		count = count + 1
		print(CYEL + "Finished {} samples.".format(count) + CEND)
	snp_df.index = position_names
	snp_df.columns = file_names
	return snp_df.transpose()

def df_clean(snp_df, percent_to_drop):
	#remove duplicate columns
	nunique = snp_df.nunique(axis=0)
	cols_to_drop = nunique[nunique == 1].index
	snp_df = snp_df.drop(cols_to_drop, axis=1)
	#drop columns with greater than 10% N's
	values = snp_df.apply(pd.value_counts).fillna(0)
	count_to_drop = round(snp_df.shape[0]*0.1)
	filtered = pd.DataFrame(values.iloc[3,]>=count_to_drop)
	#get list of column names to drop
	col_list = filtered.index[filtered["N"] == True].tolist()
	snp_df_filtered = snp_df.drop(col_list, axis = 1)
	print(CYEL+"After filtering dataframe the dimentions of the dataframe are {}.".format(snp_df_filtered.shape)+CEND)
	#convert to values
	values = {'N':0,'A':1,'T':2,'C':3,'G':4}
	snp_df_num = snp_df.applymap(lambda x: values.get(x) if x in values else x)
	snp_df_num_90 = snp_df_filtered.applymap(lambda x: values.get(x) if x in values else x)
	print(CYEL+"Finished making snp_df_num dataframe."+CEND)
	snp_df_num.to_csv('snp_df_num.csv', sep='\t', index=True)
	snp_df_num_90.to_csv('snp_df_num_90.csv', sep='\t', index=True)
	snp_df.to_csv('snp_df.csv', sep='\t', index=True)
	snp_df_filtered.to_csv('snp_df_90.csv', sep='\t', index=True)
	return snp_df_filtered, snp_df

def create_metadata(snp_df_num):
	metadata = pd.DataFrame()
	location = pd.DataFrame()
	columns = snp_df_num.columns
	chr_list = []
	pos_list = []
	for column in columns:
		chr_list.append(column.split(".")[0])
		pos_list.append(column.split(".")[1])
	location["CHR"] = chr_list
	location["POS"] = pos_list
	# collect row names in list
	metadata["Samples"] = list(snp_df_num.index.values)
	metadata['Host'] = np.where(metadata.Samples.str.contains("Human", case=False), "Human",
                   np.where(metadata.Samples.str.contains("Cow", case=False), "Cow",
                   np.where(metadata.Samples.str.contains("Calf", case=False), "Cow",
                   np.where(metadata.Samples.str.contains("Cattle", case=False), "Cow", "Unknown"))))
	metadata['Specific_Host'] = np.where(metadata.Samples.str.contains("Human", case=False), "Human",
	                   np.where(metadata.Samples.str.contains("Cow", case=False), "Cow",
	                   np.where(metadata.Samples.str.contains("Calf", case=False), "Calf",
	                   np.where(metadata.Samples.str.contains("Cattle", case=False), "Cattle", "Unknown"))))
	metadata.to_csv('metadata.csv', sep='\t', index=False)
	location.to_csv('location.csv', sep='\t', index=False)

def main():
	args = parse_cmdline()
	os.chdir(args.directory)
	snp_df = add_fasta_to_dataframe()
	snp_df_filtered, snp_df = df_clean(snp_df, args.percent_to_drop)
	#create_metadata(snp_df_num)

if __name__ == '__main__':
	main()
