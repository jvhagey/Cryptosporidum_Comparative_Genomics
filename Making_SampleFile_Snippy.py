#!/usr/bin/env/ python

# Jill Hagey
# CDC
# jvhagey@gmail.com
# https://github.com/jvhagey
# 2021
# This script was written to aid in the running of multiple samples with snippy

# importing packages
import os, re
import glob
import pandas as pd
from argparse import ArgumentParser

def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Making_SampleFile_Snippy.py", description="""This script will gather the files in 
	the indir (expects the extension _Kclean_R1.fastq) and create tab delimited files that can be used by snippy when running the snippy-multi command.""")
	parser.add_argument("-i", "--indir", dest="indir", action="store", default=None, required=True,
						help="Path for the folder where you can find the fastq files that were cleaned by Kraken and "
							 "concatenated (expects the extension _Kclean_R1.fastq).")
	parser.add_argument("-o", "--outdir", dest="outdir", action="store", default=None, required=True,
						help="Path for the folder you want the Cow_files.tab and Human_files.tab to be written to.")
	args = parser.parse_args()
	return args

def check_for_files(outdir):
	'''Deletes old version of files so it doesn't keep appending in the make_file function.'''
	try:
		os.remove(outdir + "Human_files.tab")
		os.remove(outdir + "Cow_files.tab")
		os.remove(outdir + "All_files.tab")
	except OSError:
		pass

def get_files(indir):
	'''Collects the files need'''
	os.chdir(indir)  # change directory files with extention that will be looped through code
	human_output_files = glob.glob('*_Human_Kclean_R1.fastq')  # grabbing files to loop through
	types = ('*_Cattle_Kclean_R1.fastq', '*_Cow_Kclean_R1.fastq', '*_Calf_Kclean_R1.fastq')  # the tuple of file types
	cow_output_files = []
	for files in types:
		cow_output_files.extend(glob.glob(files))
	return human_output_files, cow_output_files

def make_file(outdir, indir, human_output_files, cow_output_files):
	'''Creating the files for Snippy'''
	for file in human_output_files:
		fname = file.replace("_Kclean_R1.fastq","")
		R1 = indir + "/" + fname + "_Kclean_R1.fastq"
		R2 = indir + "/" + fname + "_Kclean_R2.fastq"
		line = fname + "\t" + R1 + "\t" + R2 + "\n"
		with open(outdir + "/" + "Human_files.tab", "a") as f:
			f.write(line)
		with open(outdir + "/" + "All_files.tab", "a") as f:
			f.write(line)
	for file in cow_output_files:
		fname = file.replace("_Kclean_R1.fastq","")
		R1 = indir + "/" + fname + "_Kclean_R1.fastq"
		R2 = indir + "/" + fname + "_Kclean_R2.fastq"
		line = fname + "\t" + R1 + "\t" + R2 + "\n"
		with open(outdir + "/" + "Cow_files.tab", "a") as f:
			f.write(line)
		with open(outdir + "/" + "All_files.tab", "a") as f:
			f.write(line)

def main():
	args = parse_cmdline()
	#indir = "Z:/Old/Alignment/Kraken_Cleaned_NoPhiX/ConCat"
	#outdir = "Z:/Old/Alignment/Snippy/Reads"
	check_for_files(args.outdir)
	human_output_files, cow_output_files = get_files(args.indir)
	make_file(args.outdir, args.indir, human_output_files, cow_output_files)

if __name__ == '__main__':
	main()