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
	echo "#! /bin/bash -l" > $fname.BBMaps10X.sh
	echo "## -- begin embedded SGE options --" >> $fname.BBMaps10X.sh
	echo "# save the standard output text to this file instead of the default jobID.o file" >> $fname.BBMaps10X.sh
	echo "#$ -o $ErrorOutDir/${fname}_BBMaps_ATCC_10X.out" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# save the standard error text to this file instead of the default jobID.e file" >> $fname.BBMaps10X.sh
	echo "#$ -e $ErrorOutDir/${fname}_BBMaps_ATCC_10X.err" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# Rename the job to be this string instead of the default which is the name of the script" >> $fname.BBMaps10X.sh
	echo "#$ -N ${fname}_BBMaps" >> $fname.BBMaps10X.sh
	echo "# " >> $fname.BBMaps10X.sh
	echo "# Requesting shared memory across 6 cpus" >> $fname.BBMaps10X.sh
	echo "#$ -pe smp 6" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# Requesting 20G of Memory for the job" >> $fname.BBMaps10X.sh
	echo "#$ -l h_vmem=20G" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# Refer all file reference to work the current working directory which is the directory from which the script was qsubbed" >> $fname.BBMaps10X.sh
	echo "#$ -cwd" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# Always run in the default queue" >> $fname.BBMaps10X.sh
	echo "#$ -q all.q" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "# Email me updates" >> $fname.BBMaps10X.sh
	echo "#$ -M qpk9@cdc.gov" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "## -- end embedded SGE options --" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "" >> $fname.BBMaps10X.sh
	echo "echo "I am running on node:" `hostname`" >> $fname.BBMaps10X.sh
	echo "mkdir -p $ErrorOutDir" >> $fname.BBMaps10X.sh
	echo "mkdir -p $Depth" >> $fname.BBMaps10X.sh
	echo "mkdir -p $SAI" >> $fname.BBMaps10X.sh
	echo "mkdir -p $outdir" >> $fname.BBMaps10X.sh
	echo "cd $Refdir" >> $fname.BBMaps10X.sh
	echo "#" >> $fname.BBMaps10X.sh
	echo "## Running BBMap for Alignment" >> $fname.BBMaps10X.sh
	echo "module load BBMap/38.84" >> $fname.BBMaps10X.sh
	echo "cd $Refdir/BBMap_Index_ATCC" >> $fname.BBMaps10X.sh
	echo "time bbmap.sh -Xmx10g threads=12 in=$indir/${fname}_10X_1.fq in2=$indir/${fname}_10X_2.fq out=$outdir/${fname}_10X.bbmap.sam outu=$outdir/${fname}_10X_unaligned_reads.bbmap.fa refstats=$Depth/${fname}_10X.bbmapStats.txt" >> $fname.BBMaps10X.sh
qsub $fname.BBMaps10X.sh
done
