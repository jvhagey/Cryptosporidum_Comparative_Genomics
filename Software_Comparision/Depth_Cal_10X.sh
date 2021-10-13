for f in /$PWD/Alignment_Eval/Results_10X/*.sam

do
	fname=$(basename $f .sam)
	STAMPYDir=/$PWD/bin/stampy-1.0.32
	ErrorOutDir=/$PWD/Alignment_Eval/Error_Out_Files
	indir=/$PWD/Alignment_Eval/Results_10X
	outdir=/$PWD/Alignment_Eval/Depth_10X
	Refdir=/$PWD/Reference_Genomes
	echo $fname
	echo "#! /bin/bash -l" > $fname.depth.sh
	echo "## -- begin embedded SGE options --" >> $fname.depth.sh
	echo "# save the standard output text to this file instead of the default jobID.o file" >> $fname.depth.sh
	echo "#$ -o $ErrorOutDir/${fname}_depth.out" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# save the standard error text to this file instead of the default jobID.e file" >> $fname.depth.sh
	echo "#$ -e $ErrorOutDir/${fname}_depth.err" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# Rename the job to be this string instead of the default which is the name of the script" >> $fname.depth.sh
	echo "#$ -N ${fname}_depth" >> $fname.depth.sh
	echo "# " >> $fname.depth.sh
	echo "# Requesting shared memory across 3 cpus" >> $fname.depth.sh
	echo "#$ -pe smp 3" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# Requesting 20G of Memory for the job" >> $fname.depth.sh
	echo "#$ -l h_vmem=20G" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# Refer all file reference to work the current working directory which is the directory from which the script was qsubbed" >> $fname.depth.sh
	echo "#$ -cwd" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# Always run in the default queue" >> $fname.depth.sh
	echo "#$ -q all.q" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "# Email me updates" >> $fname.depth.sh
	echo "#$ -M qpk9@cdc.gov" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "## -- end embedded SGE options --" >> $fname.depth.sh
	echo "#" >> $fname.depth.sh
	echo "" >> $fname.depth.sh
	echo "echo "I am running on node:" `hostname`" >> $fname.depth.sh
	echo "mkdir -p $ErrorOutDir" >> $fname.depth.sh
	echo "mkdir -p $outdir" >> $fname.depth.sh
	echo "cd $outdir" >> $fname.depth.sh
	echo "## Running BWA for alignment" >> $fname.depth.sh
	echo "module load samtools/1.10" >> $fname.depth.sh
	echo "## Converting from sam to bam file" >> $fname.depth.sh
	echo "time samtools view -bS $indir/${fname}.sam -o $indir/${fname}.bam" >> $fname.depth.sh
	echo "## Sorting bam file" >> $fname.depth.sh
	echo "time samtools sort -@ 6 $indir/${fname}.bam -o $indir/${fname}.sorted.bam" >> $fname.depth.sh
	echo "## Calculated depth coverage" >> $fname.depth.sh
	echo "time samtools depth -a -o $outdir/${fname}.alldepth.txt $indir/${fname}.sorted.bam" >> $fname.depth.sh
	echo "time samtools depth -o $outdir/${fname}.depth.txt $indir/${fname}.sorted.bam" >> $fname.depth.sh
	echo "## Running flagstat to get the stat for the bam file alignment as a .tsv file" >> $fname.depth.sh
	echo "time samtools flagstat --threads 10 -O tsv $indir/${fname}.sorted.bam > $outdir/${fname}.tsv" >> $fname.depth.sh
	echo "time samtools stat --threads 10 --reference $Refdir/CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta $indir/${fname}.bam > $outdir/${fname}.SamTools.stats.txt" >> $fname.depth.sh
	echo "## Running Picard" >> $fname.depth.sh
	echo "module load picard/2.23.0" >> $fname.depth.sh
	echo "time picard CollectAlignmentSummaryMetrics R=$Refdir/CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta I=$indir/${fname}.sorted.bam O=$outdir/${fname}.Picard.stats.txt" >> $fname.depth.sh
	echo "## Running BedTools to calculate genome coverage" >> $fname.depth.sh
	echo "module load BEDTools/2.27.1" >> $fname.depth.sh
	echo "time bedtools genomecov -bga -split -ibam $indir/${fname}.sorted.bam -g $Refdir/CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta > $outdir/${fname}.perSite_depthCoverage.txt" >> $fname.depth.sh
	echo "The Mean read depth is:" >> $fname.depth.sh
	echo "samtools depth -a $indir/${fname}.sorted.bam | awk '{c++;s+=$3}END{print s/c}'" >> $fname.depth.sh
	echo "The breadth of coverage is:" >> $fname.depth.sh
	echo "samtools depth -a $indir/${fname}.sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'" >> $fname.depth.sh
	echo "The % of reads mapped to the reference is:" >> $fname.depth.sh
	echo "samtools flagstat $indir/${fname}.sorted.bam | awk -F \"[(|%]\" 'NR == 3 {print $2}'">> $fname.depth.sh
qsub $fname.depth.sh
done
