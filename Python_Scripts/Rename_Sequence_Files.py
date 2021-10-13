#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2020
## Given Sample_Read_Locations_2020_08_24.txt and input directory rename files with host name in them.


from argparse import ArgumentParser
import pandas as pd
import shutil
import glob,os
import re
import sys


def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Get_Sequence_Files.py", description="""Takes a txt file of seqIDs on each line and creates a .sh file that will
    create links to each seq file in your current working directory. Files will have the same name they originally did""")
	parser.add_argument("-f", "--file", dest="infile", action="store", default=None, required=True, help="""Full path of
    Sample_Read_Locations_2020_08_24.txt file. Include the file name in the path. Needs to be tab delimited (required)""")
	parser.add_argument("-o", "--outdir", dest="outdir",action="store", default=None, required=True, help="Full path for output to be written to.")
	parser.add_argument("-p", dest="projectpath",action="store", default=None, required=True, help="Full path to Project_dataset.txt file")
	parser.add_argument("-i", "--indir", type=str, dest="indir", action="store", required=True, default=None,
						help='Directory where you can find .fastq.gz files. (required)')
	args = parser.parse_args()
	return args


# set colors for warnings so they are noticable
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CPUR = '\033[95m'
CEND = '\033[0m'


def find_paths(seqIDfile, infile):
	seqID = re.search("^[^_]*", seqIDfile).group(0)
	df = pd.read_csv(infile, sep='\t', header=0)  # Has sequence ID and path information
	df_sm = df[df["SampleID"] == seqID]  # return dataframe with only match SeqID in the sampleID row.
	pd.DataFrame(df_sm.to_csv("SeqPaths.txt", header=True, sep='\t', mode='w', index=None))  # Write to a new file.
	path_list = list(df_sm["Path and File Name"])
	return path_list

def get_instrument(path):
	if "HiSeq" in path[0]:
		instrument = "HiSeq"
	elif "MiSeq" in path[0]:
		instrument = "MiSeq"
	return instrument


def get_host(seqID, Project_dataset):
	df = pd.read_csv(Project_dataset, sep='\t', header=0, dtype={'SampleName': str, 'Host': str}, low_memory=False)  # Has sequence ID and path information
	host = df.loc[df['SampleName'] == str(seqID), 'Host'].iloc[0].capitalize()
	return host

def get_file_name(indir, outdir, infile, Project_dataset):
	# files with extention that will be looped through code
	#os.chdir(indir)
	#os.mkdir(outdir) # You need torun this when you aren't running snakemake
	count = 0
	fastq_files = glob.glob(indir + "/"+ "*.fastq.gz")
	for file in fastq_files:
		file_name = file.replace(indir +"/", "") # Strip the folder off the file name
		if file_name.startswith("sample"):
			file_name = re.sub("sample\d+","",file_name)
			num = re.findall(r"_\d+_",os.path.split(os.path.split(os.path.realpath(file))[0])[1])
			clean_num = num[0].replace("_","")
			file_name = clean_num + file_name
			seqID = clean_num
			end_name = file_name.replace(seqID, "")  #  Capture the rest of the file name
			end_name = re.findall(r"_L\d+_R\d_001.fastq.gz", end_name) #  Capture the rest of the file name
			end_name = end_name[0].replace("_001","")
		else:
			seqID = re.search("^[^_]*", file_name).group(0)  # Capture the Sequence ID at the beginning
			end_name = file_name.replace(seqID, "")  #  Capture the rest of the file name
		host = get_host(seqID, Project_dataset)
		path_list = find_paths(file_name, infile)
		machine = get_instrument(path_list)
		#length = get_seq_length(file)
		os.rename(file, outdir + "/" + seqID + "_" + host + end_name)
		#os.rename(file, seqID + "_" + machine + "_" + host + "_" + end_name)
		count = count + 1
	print(CYEL + "\n There was {} files renamed and can be found in {}\n".format(count, os.listdir(indir)) + CEND)
	sys.stdout.write(CYEL + "\n There was {} files renamed and can be found in {}\n".format(count, os.listdir(indir)) + CEND)
	#shutil.rmtree(indir)

def shorted_file_name(indir):
	for fileName in os.listdir(indir):
		if fileName.endswith(".fastq.gz"):
			new_name = re.sub(r"(^[0-9]+_).*(L00[1-3]_R[1-2]).*(_001)", r'\1\2', fileName)
			os.rename(indir+"/"+fileName, indir+"/"+new_name)
	if not any(fileName.endswith('.fastq.gz') for fileName in os.listdir(indir)):
		print(CRED + "\n" + " There were no '.fastq.gz' files the found in the folder " + indir + " :/\n" + CEND)
		sys.stderr.write(CRED + "\n" + " There were no '.fastq.gz' files the found in the folder " + indir + " :/\n" + CEND)

def main():
	args = parse_cmdline()
	shorted_file_name(args.indir)
	get_file_name(args.indir, args.outdir, args.infile, args.projectpath)
	#get_file_name("C:/Users/qpk9/OneDrive - CDC/Cryto_GWAS", "Sample_Read_Locations_2020_08_24.txt", "C:/Users/qpk9/OneDrive - CDC/Cryto_GWAS/new_name")


if __name__ == '__main__':
	main()
