#! /bin/bash -l
## -- begin embedded SGE options --
# save the standard output text to this file instead of the default jobID.o file
#$ -o /scicomp/home/qpk9/Alignment_Eval/Error_Out_Files/Simulated_Seq_Mi250_Alignment2.out
#
# save the standard error text to this file instead of the default jobID.e file
#$ -e /scicomp/home/qpk9/Alignment_Eval/Error_Out_Files/Simulated_Seq_Mi250_Alignment2.err
#
# Rename the job to be this string instead of the default which is the name of the script
#$ -N Simulated_Seq_Mi250_Alignment2
# 
# Requesting shared memory across 6 cpus
#$ -pe smp 6
#
# Requesting 20G of Memory for the job
#$ -l h_vmem=20G
#
# Refer all file reference to work the current working directory which is the directory from which the script was qsubbed
#$ -cwd
#
# Always run in the default queue
#$ -q all.q
#
# Email me updates
#$ -M qpk9@cdc.gov
#
## -- end embedded SGE options --
#

echo I am running on node: node63.hpc.biotech.cdc.gov
mkdir -p /scicomp/home/qpk9/Alignment_Eval/Error_Out_Files
mkdir -p /scicomp/home/qpk9/Alignment_Eval
cd /scicomp/home/qpk9/Reference_Genomes
#
## Running Bowtie2 for alignment
#module load bowtie2/2.3.5.1
## Running local alignment
#time bowtie2 --local --threads 12 -x CparvumIOWAII_Index -1 /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq -2 /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq -S /scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250.local.sam
## Running end-to-end alignment
#time bowtie2 --end-to-end --threads 12 -x CparvumIOWAII_Index -1 /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq -2 /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq -S /scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250.E2E.sam
##
## Running BWA for Alignment
module load bwa/0.7.17
#cd /scicomp/home/qpk9/Reference_Genomes/BWA_Index
## Running BWA for alignment
#time bwa mem -t 12 CparvumIOWAII_Index /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq > /scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250.bwa.sam
##
## Running BBMap for Alignment
#module load BBMap/38.84
#time bbmap.sh ref=/scicomp/home/qpk9/Reference_Genomes/CparvumIOWAII.fasta in=/scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq in2=/scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq out=/scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250.bbmap.sam outu=/scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250_unaligned_reads.bbmap.fa
## Running Stampy for Alignment
## To speed up mapping, use BWA (recommended) and multithreading:
## Running BWA for alignment
module load samtools/1.10
## Find the alignments in suffix array (SA) coordinates of the input reads. The resulting file contains ALL the alignments found by BWA.
time bwa aln -t12 /scicomp/home/qpk9/Reference_Genomes/BWA_Index/CparvumIOWAII_Index /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq > Simulated_Seq_Mi250_1.sai
time bwa aln -t12 /scicomp/home/qpk9/Reference_Genomes/BWA_Index/CparvumIOWAII_Index /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq > Simulated_Seq_Mi250_2.sai
## Converting SA coordinates to chromosomal coordinates. Generate alignments in the SAM format given paired-end reads, which is piped into samtools to make the bam file. Repetitive read pairs will be placed randomly.
time bwa sampe /scicomp/home/qpk9/Reference_Genomes/BWA_Index/CparvumIOWAII_Index Simulated_Seq_Mi250_1.sai Simulated_Seq_Mi250_2.sai /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq | samtools view -Sb > Simulated_Seq_Mi250.bwa.bam
time python2 /scicomp/home/qpk9/bin/stampy-1.0.32/stampy.py --genome /scicomp/home/qpk9/Reference_Genomes/Stampy_Index/CparvumIOWAII_Index --hash /scicomp/home/qpk9/Reference_Genomes/Stampy_Index/CparvumIOWAII_Index --threads 12 --bamkeepgoodreads --map Simulated_Seq_Mi250.bwa.bam
time python2 /scicomp/home/qpk9/bin/stampy-1.0.32/stampy.py --threads 12
module load minimap2/2.17
# Aligning with Minimaps2
time minimap2 -ax sr /scicomp/home/qpk9/Reference_Genomes/MiniMap_Index/CparvumIOWAII_Index.mmi /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_1.fq /scicomp/home/qpk9/Reference_Genomes/Simulated_Seq/Simulated_Seq_Mi250_2.fq > /scicomp/home/qpk9/Alignment_Eval/Simulated_Seq_Mi250.minimap2.sam
