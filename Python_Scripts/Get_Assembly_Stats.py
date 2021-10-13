#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2020
## Given the outfile (.out) from kneaddata returns read and assembly stats for graphing in R

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import glob, os
import re
import sys

def parse_cmdline():
	'''Parse command-line arguments for script.'''
	parser = ArgumentParser(prog="Get_Sequence_Files.py", description="""Given sample_ID.txt file (or will run through all directories in input directory if not file is given) and input directory it will produce Assembly_stats.csv file that contains 
     in the output directory with Genome_Fraction, N50, Total_Alignment_Length, Duplication_Ratio, Number_of_Contigs, Total_Assembly_Length,
     Largest_Contig, Missasemblies'in the file. Ex: python Get_Assembly_Stats.py -i /scicomp/home/qpk9/Assembly/quast_results/ -o /scicomp/home/qpk9/Assembly/quast_results/ -f sample_IDs.txt -c 20 99 1000 1.5""")
	parser.add_argument("-f", "--file", dest="infile", action="store", default=None, required=False, help="A file with sample_IDs on each line. These will be the same IDs as the folder names that are the output of quast. If no file is given (default) script will run through all directories in input directory.")
	parser.add_argument("-o", "--outdir", dest="outdir",action="store", default=None, required=True, help="Full path for output 'Assembly_stats.csv' to be written to.")
	parser.add_argument("-i", "--indir", type=str, dest="indir", action="store", required=True, default=None, help='Full path to where you can find quast results folders to loop through (required)')
	parser.add_argument("-c", "--cutoffs", type=float, dest="cutoffs", action="store", nargs="+", required=False, default=None, 
						help='A comma separated list of cutoffs for stats you want warnings for. The numbers should be in this order Missasemblies (30), Genome Fraction (95) percent, '
						'Number of_Contigs (2000) and Duplication_Ratio (1.2). Default values in (). Currently, you need to just let it run with default values or set all of them.')
	args = parser.parse_args()
	return args


# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CPUR = '\033[95m'
CEND = '\033[0m'


def enter_folder(infile, indir, outdir):
	""" This function checks to see if an infile was given and changes to the correct directory to loop through assembly stats. """
	df_new = pd.DataFrame(columns=['SeqID', 'Genome_Fraction', 'N50', 'Total_Alignment_Length', 'Duplication_Ratio',
				 'Number_of_Contigs', 'Total_Assembly_Length', 'Largest_Contig', 'Missasemblies'])  # make a empty dataframe
	if infile is not None: 
		df = pd.read_csv(infile, sep='\t', index_col=0, header=None)
		for i in range(df.shape[0]):
			sample = str(df.iloc[i].name)  # folder name that will be looped through code.
			os.chdir(indir + sample)
			df_new = get_assembly_stats(df, indir, outdir, df_new, sample)
	else:
		dirlist = os.listdir(indir) # get files and directories in the input directory
		for folder in dirlist:
			os.chdir(indir + folder)
			#sample=os.path.split(os.path.split(os.path.realpath(__file__))[0])[1]
			sample=os.path.split(os.path.abspath(os.getcwd()))[1]
			df_new = get_assembly_stats(None, indir, outdir, df_new, sample)
	os.chdir(outdir)
	df_new.to_csv('Assembly_stats.csv', sep='\t', index=True)
	return df_new
		
def get_assembly_stats(df, indir, outdir, df_new, sample):
	""" This function enters each quast folder, parses the reports.tsv file into a 'Assembly_stats.csv' file that
	contains the assembly stats. This will be useful information for graphing later. """
	try:
		df_results = pd.read_csv('report.tsv', sep='\t')
		df_results.columns = ['Assembly', 'scaffolds'] # The second column can be called contigs or scaffolds so this is just to make sure its general enough to run.
		n50 = int(df_results.loc[df_results['Assembly'] == 'N50', 'scaffolds'].item())
		genome_fraction = float(df_results.loc[df_results['Assembly'] == 'Genome fraction (%)','scaffolds'].item())
		total_alignment_length = int(df_results.loc[df_results['Assembly'] == 'Total aligned length', 'scaffolds'].item())
		contig_num = int(df_results.loc[df_results['Assembly'] == '# contigs', 'scaffolds'].item())
		total_assembly_length = int(df_results.loc[df_results['Assembly'] == 'Total length', 'scaffolds'].item())
		largest_contig = int(df_results.loc[df_results['Assembly'] == 'Largest contig', 'scaffolds'].item())
		missasemblies = int(df_results.loc[df_results['Assembly'] == '# misassemblies', 'scaffolds'].item())
		Duplication_ratio = float(df_results.loc[df_results['Assembly'] == 'Duplication ratio', 'scaffolds'].item())
		new_row = {'SeqID': sample, 'Genome_Fraction': genome_fraction, 'Duplication_Ratio': Duplication_ratio, 'N50': n50, 'Total_Alignment_Length': total_alignment_length, 'Number_of_Contigs': contig_num, 'Total_Assembly_Length': total_assembly_length, 'Largest_Contig': largest_contig, 'Missasemblies': missasemblies}
		df_new = df_new.append(new_row, ignore_index=True)
	except FileNotFoundError as err:
		print(CRED + ' Error: File not found for sample {}'.format(str(df.iloc[i].name)) , err)
		print('Just FYI you are currently in the directory {}'.format(os.getcwd())+ CEND)
	except:
		print(CRED + " Well this is embarrassing I wasn't able to incorporate the assembly statistics from sample {"
					 "}. There was an unexpected error: {}{}, check line {}.".format(df.iloc[i].name, sys.exc_info()[0],sys.exc_info()[1], sys.exc_info()[2].tb_lineno) + CEND)
		pass
	return df_new


def check_assembly_stats(df, outdir, missass_num, GF_num, contig_num, DR_num):
	""" This function will take the Assembly stats file and check for any samples that might need to be looked at. In \
	its current form it looks for high missassemblies and low genome coverage etc. """
	pd.set_option('display.max_rows', None, 'display.max_columns', None)
	df = pd.read_csv(outdir +'Assembly_stats.csv', sep='\t')  # Read in file that was just written
	dict_to_check = {'Missasemblies': missass_num, 'Genome_Fraction':GF_num, 'Number_of_Contigs': contig_num, 'Duplication_Ratio': float(DR_num)}
	for key,item in dict_to_check.items():
		if key != "Genome_Fraction":
			if (df[key] >= item).any():
				bad_samples = (df.loc[df[key] >= item])['SeqID'].unique()
				list_of_stats = (df.loc[df[key] >= item])[key].unique()
				print(CYEL + '\nLooks like sample number(s) {} has/have {} {}, respectively. Probably a good idea to have a closer ' \
				'look.\n'.format(bad_samples, list_of_stats, key.lower()) + CEND)
			else:
				print(CPUR + 'Yay! There are no samples with {} greater than or equal to {}!\n'.format(key.lower().replace("_"," "), item) + CEND)
		else:
			item = float(item)
			if (df[key] <= item).any():
				bad_sample = (df.loc[df[key] <= item])['SeqID'].unique()
				statistic = (df.loc[df[key] <= item])[key].unique()
				print(CYEL + '\nLooks like sample number(s) {} has/have a {} of {}%. Probably a good idea to have a closer ' \
				'look.\n'.format(bad_sample, key.lower(), statistic) + CEND)
			else:
				print(CPUR + 'Yay! There are no samples with {} less than or equal to {}!\n'.format(key.lower().replace("_"," "), item) + CEND)


def main():
	args = parse_cmdline()
	#os.chdir(args.indir)
	#df = get_read_stats()
	if args.cutoffs is not None:
		tup = tuple(args.cutoffs)  # Make tuple
		missass_num, GF_num, contig_num, DR_num = tup # Unpacking
	else:
		tup = (10, 95, 2000, 1.1)  # Make tuple
		missass_num, GF_num, contig_num, DR_num = tup # Unpacking
	df_assembly_stats = enter_folder(args.infile, args.indir, args.outdir)
	check_assembly_stats(df_assembly_stats, args.outdir, missass_num, GF_num, contig_num, DR_num)


if __name__ == '__main__':
	main()
