for f in /scicomp/home-pure/qpk9/Alignment_Eval/Simulated_Seq_10X/*.log

do
	fname=$(basename $f _10X.log)
	STAMPYDir=/scicomp/home-pure/qpk9/bin/stampy-1.0.32
	ErrorOutDir=/scicomp/home-pure/qpk9/Alignment_Eval/Error_Out_Files
	indir=/scicomp/home-pure/qpk9/Alignment_Eval/Simulated_Seq_10X
	SAI=/scicomp/home-pure/qpk9/Alignment_Eval/SAI_10X
	outdir=/scicomp/home-pure/qpk9/Alignment_Eval/Results_10X
	Refdir=/scicomp/home-pure/qpk9/Reference_Genomes
	Depth=/scicomp/home-pure/qpk9/Alignment_Eval/Depth_10X
	echo $fname
	echo "#! /bin/bash -l" > $fname.Alignment10X.sh
	echo "## -- begin embedded SGE options --" >> $fname.Alignment10X.sh
	echo "# save the standard output text to this file instead of the default jobID.o file" >> $fname.Alignment10X.sh
	echo "#$ -o $ErrorOutDir/${fname}_Alignment_ATCC_10X.out" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# save the standard error text to this file instead of the default jobID.e file" >> $fname.Alignment10X.sh
	echo "#$ -e $ErrorOutDir/${fname}_Alignment_ATCC_10X.err" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# Rename the job to be this string instead of the default which is the name of the script" >> $fname.Alignment10X.sh
	echo "#$ -N ${fname}_Alignment" >> $fname.Alignment10X.sh
	echo "# " >> $fname.Alignment10X.sh
	echo "# Requesting shared memory across 6 cpus" >> $fname.Alignment10X.sh
	echo "#$ -pe smp 6" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# Requesting 20G of Memory for the job" >> $fname.Alignment10X.sh
	echo "#$ -l h_vmem=20G" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# Refer all file reference to work the current working directory which is the directory from which the script was qsubbed" >> $fname.Alignment10X.sh
	echo "#$ -cwd" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# Always run in the default queue" >> $fname.Alignment10X.sh
	echo "#$ -q all.q" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "# Email me updates" >> $fname.Alignment10X.sh
	echo "#$ -M qpk9@cdc.gov" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "## -- end embedded SGE options --" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "" >> $fname.Alignment10X.sh
	echo "echo "I am running on node:" `hostname`" >> $fname.Alignment10X.sh
	echo "mkdir -p $ErrorOutDir" >> $fname.Alignment10X.sh
	echo "mkdir -p $Depth" >> $fname.Alignment10X.sh
	echo "mkdir -p $SAI" >> $fname.Alignment10X.sh
	echo "mkdir -p $outdir" >> $fname.Alignment10X.sh
	echo "cd $Refdir" >> $fname.Alignment10X.sh
	echo "#" >> $fname.Alignment10X.sh
	echo "echo '## Running Bowtie2 for alignment'" >> $fname.Alignment10X.sh
	echo "module load bowtie2/2.3.5.1" >> $fname.Alignment10X.sh
	echo "echo '## Running local alignment'" >> $fname.Alignment10X.sh
	echo "time bowtie2 --local --threads 12 -x CparvumIOWA-ATCC -1 $indir/${fname}_10X_1.fq -2 $indir/${fname}_10X_2.fq -S $outdir/${fname}_10X.local.sam" >> $fname.Alignment10X.sh
	echo "echo '## Running end-to-end alignment'" >> $fname.Alignment10X.sh
	echo "time bowtie2 --end-to-end --threads 12 -x CparvumIOWA-ATCC -1 $indir/${fname}_10X_1.fq -2 $indir/${fname}_10X_2.fq -S $outdir/${fname}_10X.E2E.sam" >> $fname.Alignment10X.sh
	echo "##" >> $fname.Alignment10X.sh
	echo "echo '## Running BWA for Alignment'" >> $fname.Alignment10X.sh
	echo "module load bwa/0.7.17" >> $fname.Alignment10X.sh
	echo "cd $Refdir/BWA_Index_ATCC" >> $fname.Alignment10X.sh
	echo "echo '## Running BWA for alignment'" >> $fname.Alignment10X.sh
	echo "time bwa mem -t 12 CparvumIOWA-ATCC $indir/${fname}_10X_1.fq $indir/${fname}_10X_2.fq > $outdir/${fname}_10X.bwa.sam" >> $fname.Alignment10X.sh
	echo "##" >> $fname.Alignment10X.sh
	echo "echo '## Running BBMap for Alignment'" >> $fname.Alignment10X.sh
	echo "module load BBMap/38.84" >> $fname.Alignment10X.sh
	echo "cd $Refdir/BBMap_Index_ATCC" >> $fname.Alignment10X.sh
	echo "time bbmap.sh -Xmx10g threads=12 in=$indir/${fname}_10X_1.fq in2=$indir/${fname}_10X_2.fq out=$outdir/${fname}_10X.bbmap.sam outu=$outdir/${fname}_unaligned_reads.bbmap.fa refstats=$Depth/${fname}_10X.bbmapStats.txt" >> $fname.Alignment10X.sh
	echo "echo '## Running Stampy for Alignment'" >> $fname.Alignment10X.sh
	echo "## To speed up mapping, use BWA (recommended) and multithreading:" >> $fname.Alignment10X.sh
	echo "echo '## Running BWA for alignment'" >> $fname.Alignment10X.sh
	echo "module load samtools/1.10" >> $fname.Alignment10X.sh
	echo "cd $outdir" >> $fname.Alignment10X.sh
	echo "echo '## Find the alignments in suffix array (SA) coordinates of the input reads. The resulting file contains ALL the alignments found by BWA.'" >> $fname.Alignment10X.sh
	echo "time bwa aln -t12 $Refdir/BWA_Index_ATCC/CparvumIOWA-ATCC $indir/${fname}_10X_1.fq > $SAI/${fname}_10X_1.sai" >> $fname.Alignment10X.sh
	echo "time bwa aln -t12 $Refdir/BWA_Index_ATCC/CparvumIOWA-ATCC $indir/${fname}_10X_2.fq > $SAI/${fname}_10X_2.sai" >> $fname.Alignment10X.sh
	echo "echo '## Converting SA coordinates to chromosomal coordinates. Generate alignments in the SAM format given paired-end reads, which is piped into samtools to make the bam file. Repetitive read pairs will be placed randomly.'" >> $fname.Alignment10X.sh
	echo "time bwa sampe $Refdir/BWA_Index_ATCC/CparvumIOWA-ATCC $SAI/${fname}_10X_1.sai $SAI/${fname}_10X_2.sai $indir/${fname}_10X_1.fq $indir/${fname}_10X_2.fq | samtools view -Sb > $outdir/${fname}_10X.bwa2Stampy.bam" >> $fname.Alignment10X.sh
	echo "time python2 $STAMPYDir/stampy.py --genome $Refdir/Stampy_Index/CparvumIOWA-ATCC_Index --hash $Refdir/Stampy_Index/CparvumIOWA-ATCC_Index --threads 12 --bamkeepgoodreads --map $outdir/${fname}_10X.bwa2Stampy.bam" >> $fname.Alignment10X.sh
	echo "echo '## Convert bam to sam'" >> $fname.Alignment10X.sh
	echo "time samtools view -h -o $outdir/${fname}_10X.bwa2Stampy.sam $outdir/${fname}_10X.bwa2Stampy.bam" >> $fname.Alignment10X.sh
	echo "module load minimap2/2.17" >> $fname.Alignment10X.sh
	echo "echo '## Aligning with Minimaps2'" >> $fname.Alignment10X.sh
	echo "time minimap2 -ax sr $Refdir/MiniMap_Index_ATCC/CparvumIOWA-ATCC.mmi $indir/${fname}_10X_1.fq $indir/${fname}_10X_2.fq > $outdir/${fname}_10X.minimap2.sam" >> $fname.Alignment10X.sh
qsub $fname.Alignment10X.sh
done
