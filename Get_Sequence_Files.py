#!/usr/bin/env python

## Jill Hagey
## CDC
## qpk9@cdc.gov
## https://github.com/jvhagey/
## 2020
## Given Sample_Read_Locations_2020_11_02.txt and file with SeqIDs it will generate .sh file for creating links to the directory is run in.


from argparse import ArgumentParser
import pandas as pd
import os

# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CEND = '\033[0m'

def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="Get_Sequence_Files.py", description="""Takes a txt file of seqIDs on each line and creates a .sh file that will 
    create links to each seq file in your current working directory. Files will have the same name they originally did""")
	parser.add_argument("-i", "--infile", dest="infile", action="store", default=None, required=True, help="""Full location of 
    Sample_Read_Locations_2020_11_02.txt file. Include the file name in the path. Needs to be tab delimited (required)""")
	parser.add_argument("-o", "--outdir", dest="outdir",action="store", default=None, required=True, help="Directory for where you want the sequences to end up when CreateRawReadLinks.sh is run.")
	parser.add_argument("-f", "--seqID-file", type=str, dest="seqIDfile", action="store", required=True, default=None,
						help='A file with seqID on each line (required)')
	args = parser.parse_args()
	return args


def find_paths(seqIDfile, infile):
	"""Gets paths for select sample IDs"""
	with open(seqIDfile) as file:
		seqID = file.read().splitlines()
	df = pd.read_csv(infile, sep='\t', header=0)  # Has sequence ID and path information
	df_sm = df[df["SampleID"].isin(seqID)]  # return dataframe with only match SeqID in the sampleID row.
	pd.DataFrame(df_sm.to_csv("SeqPaths.txt", header=True, sep='\t', mode='w', index=None))  # Write to a new file.
	path_list = list(df_sm["Best_Sequencing_Run"])
	return path_list


def write_new(path_list, outdir):
	"""Gets paths for select sample IDs"""
	count = 0
	for path in path_list:
		for filename in os.listdir(path):
			if filename.endswith(".fastq.gz"):
				with open("CreateRawReadLinks.sh", "a") as f:
					f.write("ln -s " + os.path.join(path, filename) + " " + os.path.join(outdir, filename) + "\n")
					count = count + 1
	print(CYEL + "\n" + "A total of " + str(count) + " paths were written for linking were written to CreateRawReadLinks.sh" + CEND + "\n")


def main():
	args = parse_cmdline()
	path_list = find_paths(args.seqIDfile, args.infile)
	write_new(path_list, args.outdir)


if __name__ == '__main__':
	main()
