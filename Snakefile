####################################################
#Some python stuff
# importing packages
#import pandas as pd
import os

# set colors for warnings so they are seen
CRED = '\033[91m' + '\nWarning:'
CYEL = '\033[93m'
CEND = '\033[0m'

#pd.set_option("max_columns", None)
#pd.set_option("max_rows", None)
#pd.set_option('display.width', 1000)

####################################################

#Setting local rules to not be run on an HPC
localrules: create_links, renaming_files, checkpoint, assembly_stats, checkpoint2, generate_file_for_snippy, depth_check, snippy_core, copy, rename, combine2, make_bad_list_file

shell.executable('/usr/bin/bash')
shell.prefix('source ~/.bash_profile; ')

configfile: '/scicomp/home-pure/qpk9/GWAS_Study/config.json'
workdir: '/scicomp/home-pure/qpk9/SM_Test'

OLD_NAME_SAMPLES, = glob_wildcards(config['old_data'] + '{filename}.fastq.gz')
NEW_FILE_NAME, = glob_wildcards(config['data'] + '{sample}_R1.fastq.gz')
NL_FILE_NAME, LANE = glob_wildcards(config['data'] + '{nolane_sample}_L{lane}_R1.fastq.gz')
NL_FILE_NAME_LESS, = glob_wildcards(config['workdir'] + 'QC/Concat/{sample}_Kclean_R1.fastq')


########################################################

#getting list of cattle vs Human_Core
NL_FILE_NAME_LESS
Human_samples = []
Cow_samples = []
matches = ["Cattle", "Cow", "Calf"]
for name in NL_FILE_NAME_LESS:
	if "Human" in name:
		Human_samples.append(name)
	if any(x in name for x in matches):
		Cow_samples.append(name)

#cleaning lists
Bad_list = ['41897_Cattle', '39640_Human', '42129_Cattle', '40457_Human', '40236_Calf', '41573_Cow', '42134_Cattle', '41891_Cattle', '41892_Cattle']
#making copy
New_all = NL_FILE_NAME_LESS
for sample in Bad_list:
	if sample in Human_samples:
		Human_samples.remove(sample)
	if sample in Cow_samples:
		Cow_samples.remove(sample)
	if sample in New_all:
		New_all.remove(sample)
########################################################

#When two rules can produce the same output file, snakemake cannot decide which one to use without additional guidance (an AmbiguousRuleException will be thrown).
#Here we add a rule order to help snakemake to decide which rule to use since multiple ones can create the same output file!
ruleorder: concat_L001 > concat_L002 > concat_L003

# Note that this snakemake assumes you have already created bowtie indexes for cow, human and phiX genomes.
rule all:
	input:
		# raw
		expand(config['old_data'] + '{filename}.fastq.gz', filename=OLD_NAME_SAMPLES),
		expand(config['old_data'] + '{filename}.fastq.gz', filename=OLD_NAME_SAMPLES),
		# Raw renamed
		expand(config['data'] + '{sample}_R1.fastq.gz', sample=NEW_FILE_NAME),
		expand(config['data'] + '{sample}_R2.fastq.gz', sample=NEW_FILE_NAME),
		# QC files
		expand(config['workdir'] + 'QC/Kraken_Cleaned/{sample}_Kclean_R1.fastq', sample=NEW_FILE_NAME),
		expand(config['workdir'] + 'QC/Kraken_Cleaned/{sample}_Kclean_R2.fastq', sample=NEW_FILE_NAME),
		expand(config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R1.fastq', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R2.fastq', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'QC/Kraken/{sample}.kraken', sample=NEW_FILE_NAME),
		expand(config['workdir'] + 'QC/Kraken/{sample}_report.txt', sample=NEW_FILE_NAME),
		# Assembly files
		expand(config['workdir'] + 'Assembly/Skesa/{nolane_sample}_contigs_skesa.fa', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Assembly/quast_results/{nolane_sample}/report.tsv', nolane_sample=NL_FILE_NAME),
		config['workdir'] + 'Assembly/Assembly_stats.csv',
		config['workdir'] + 'Assembly/Problem_Assemblies.txt',
		## Making Tree
		config['workdir'] + 'Variant_Call_Test/Tree/Skesa_Contigs.fasta',
		# Alignment files
		expand(config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.sam', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.bam', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.sorted.bam', nolane_sample=NL_FILE_NAME),
		# Variant call test files
		config['workdir'] + 'Variant_Call_Test/Consensus/Depth_stats.csv',
		config['workdir'] + 'Variant_Call_Test/Consensus/combined-masked-seq-RN.fasta',
		expand(config['workdir'] + 'Variant_Call_Test/Freebayes/{nolane_sample}-freebayes-SNP-differences.vcf', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz', nolane_sample=NL_FILE_NAME),
		## Merge vcf files
		config['workdir'] + 'Variant_Call_Test/VCF/all_merge.vcf.gz',
		config['workdir'] + 'Variant_Call_Test/VCF/human_merge.vcf.gz',
		config['workdir'] + 'Variant_Call_Test/VCF/cow_merge.vcf.gz',
		## Get masked sequences
		expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-per-site-depth-coverage.txt', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq.fasta', nolane_sample=NL_FILE_NAME),
		expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-genes.fasta', nolane_sample=NL_FILE_NAME),
		config['workdir'] + 'Variant_Call_Test/Consensus/depth_check_output.txt',
		# Making snp-dist matrix
		config['workdir'] + 'Variant_Call_Test/Consensus/checkpoint3.txt',
		config['workdir'] + 'Variant_Call_Test/distances.tab',
		# Snippy - Variant calling files
		multiext(config['workdir'] + 'VariantCalls/Snippy/All_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf"),
		multiext(config['workdir'] + 'VariantCalls/Snippy/Human_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf"),
		multiext(config['workdir'] + 'VariantCalls/Snippy/Cow_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf")

#You first need to run CreateRawReadLinks.sh
#sh CreateRawReadLinks.sh
#python Get_Sequence_Files.py -i Sample_Read_Locations_2020-11-02.txt -f Sample_IDs.txt -o Test_Sequences/

# Create links for Raw Sequences
rule create_links:
	input:
		script = config['workdir'] + 'createRawReadLinks_n132.sh'
	params:
		old_data = config['old_data'],
		data = config['data']
	output:
		files = expand(config['old_data']+'{filename}.fastq.gz', filename=OLD_NAME_SAMPLES)
	shell:
		'''
		echo 'I am running on node:' `hostname`
		set +u
		mkdir -p {params.old_data}
		mkdir -p {params.data}
		cd {params.old_data}
		sh {input}
		set -u
		'''

# Renaming files to it has the host in it.
#You will mostly likely have to force this rule to run and then run the rest of the workflow:
#time snakemake -j 15 --profile Comparative_Genomics -R renaming_files
rule renaming_files:
	input:
		expand(config['old_data']+'{filename}.fastq.gz', filename=OLD_NAME_SAMPLES),
		sam_locations = config['workdir'] +'Sample_Read_Locations_2020_08_24.txt',
		project_dataset = config['workdir'] +'Project_dataset.txt',
		seq_dir_old = config['old_data'],
		seq_dir = config['data']
	params:
		scripts = config['scripts']
	output:
		expand(config['data']+'{short_sample}_R1.fastq.gz', short_sample=NEW_FILE_NAME),
		expand(config['data']+'{short_sample}_R2.fastq.gz', short_sample=NEW_FILE_NAME)
	message: '''--- Renaming files to shorter name.'''
	shell: 'python {params.scripts}/Rename_Sequence_Files.py -o {input.seq_dir} -f {input.sam_locations} -i {input.seq_dir_old} -p {input.project_dataset}'

#print(expand(config['data']+'{short_sample}_R1.fastq.gz', short_sample=NEW_FILE_NAME))

rule bbduk:
	input:
		file_1 = config['data']+'{sample}_R1.fastq.gz',
		file_2 = config['data']+'{sample}_R2.fastq.gz'
	output:
		file_nophix_1 = config['workdir'] + 'QC/BBDuK/{sample}_nophix_R1.fastq.gz',
		file_nophix_2 = config['workdir'] + 'QC/BBDuK/{sample}_nophix_R2.fastq.gz',
		file_adap_1 = config['workdir'] + 'QC/BBDuK/{sample}_noadp_R1.fastq.gz',
		file_adap_2 = config['workdir'] + 'QC/BBDuK/{sample}_noadp_R2.fastq.gz',
		file_clean_1 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R1.fastq',
		file_clean_2 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R2.fastq'
	params:
		wkdir = config['workdir'],
		outdir = config['workdir'] + 'QC/BBDuK/Stats',
		bbdukdir = config['workdir'] + 'QC/BBDuK',
		ErrorOut = config['ErrorOut'] + 'BBDuK',
		adap = '/apps/x86_64/bbmap/38.84/resources/adapters.fa',
		phix ='/apps/x86_64/bbmap/38.84/resources/phix174_ill.ref.fa.gz',
		fastqc = config['workdir'] + 'QC/BBDuK/fastqc'
	message: '''--- Running BBDuk to remove adapters and trim for quality ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p /scicomp/home-pure/qpk9/SM_Test/QC
		mkdir -p {params.bbdukdir}
		mkdir -p {params.outdir}
		mkdir -p {params.fastqc}
		module load BBMap/38.84
		time bbduk.sh -Xmx1g in={input.file_1} in2={input.file_2} out={output.file_adap_1} out2={output.file_adap_2} ref={params.adap} ktrim=r k=23 mink=11 hdist=1 refstats={params.outdir}/{wildcards.sample}_RefStats.txt tpe tbo
		time bbduk.sh in={output.file_adap_1} in2={output.file_adap_2} out={output.file_nophix_1} out2={output.file_nophix_2} outm={params.bbdukdir}/{wildcards.sample}_matched.fq ref={params.phix} k=31 hdist=1 stats={params.outdir}/{wildcards.sample}_PhiX_Stats.txt
		time bbduk.sh -Xmx1g ftl=10 ftr2=10 qtrim=r trimq=30 minlength=50 in={output.file_nophix_1} in2={output.file_nophix_2} out={output.file_clean_1} out2={output.file_clean_2} stats={params.outdir}/{wildcards.sample}_Stats.txt
		module load fastqc/0.11.5
		time fastqc {output.file_clean_1} -o {params.fastqc}
		time fastqc {output.file_clean_2} -o {params.fastqc}
		'''

rule kraken:
	input:
		file_1 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R1.fastq',
		file_2 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R2.fastq'
	params:
		DBdir = '/scicomp/home-pure/qpk9/Kraken_DB/Updated_DB/',
		outdir = config['workdir'] + 'QC/Kraken',
		ErrorOut = config['ErrorOut'] + 'Kraken',
	output:
		kraken = config['workdir'] + 'QC/Kraken/{sample}.kraken',
		report = config['workdir'] + 'QC/Kraken/{sample}_report.txt'
	message: '''--- Running kraken to check for contamination ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		module load kraken/2.0.8
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		time kraken2 --use-names --threads 10 --db {params.DBdir} --report {output.report} --paired {input.file_1} {input.file_2} --output {output.kraken}
		'''

rule clean_kraken:
	input:
		file_1 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R1.fastq',
		file_2 = config['workdir'] + 'QC/BBDuK/{sample}_clean_R2.fastq',
		kraken = config['workdir'] + 'QC/Kraken/{sample}.kraken',
		report = config['workdir'] + 'QC/Kraken/{sample}_report.txt'
	params:
		outdir = config['workdir'] + 'QC/Kraken_Cleaned',
		ErrorOut = config['ErrorOut'] + 'Kraken_Cleaned'
	output:
		k_file_1 = config['workdir'] + 'QC/Kraken_Cleaned/{sample}_Kclean_R1.fastq',
		k_file_2 = config['workdir'] + 'QC/Kraken_Cleaned/{sample}_Kclean_R2.fastq'
	message: '''--- Running extract_kraken_reads.py to get fastq files without contamination. ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		time python /scicomp/home-pure/qpk9/bin/KrakenTools-master/extract_kraken_reads.py --include-children --fastq-output --taxid 5806 -s {input.file_1} -s2 {input.file_2} -o {output.k_file_1} -o2 {output.k_file_2} --report {input.report} -k {input.kraken}
		'''

rule checkpoint:
	input:
		expand(config['workdir'] + 'QC/Kraken_Cleaned/{sample}_Kclean_R{read}.fastq', sample=NEW_FILE_NAME, read=['1','2'])
	output:
		config['workdir']+ 'checkpoint.txt'
	shell:
		'''
		echo {input} files are finished > checkpoint.txt
		'''

rule concat_L001:
	input:
		checkpoint = config['workdir']+ 'checkpoint.txt',
		k_file_1 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L001_Kclean_R1.fastq',
		k_file_2 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L001_Kclean_R2.fastq'
	params:
		outdir = config['workdir'] + 'QC/Concat',
		kdir = config['workdir'] + 'QC/Kraken_Cleaned',
		ErrorOut = config['ErrorOut'] + 'Concat'
	output:
		con_file_1 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R1.fastq',
		con_file_2 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R2.fastq'
	message: '''--- Concating Sequences ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		cd /scicomp/home-pure/qpk9/SM_Test/QC/Kraken_Cleaned
		if [ -f {params.kdir}/{wildcards.nolane_sample}_L001_Kclean_R1.fastq ]
		then
			cat {input.k_file_1} > {output.con_file_1}
			cat {input.k_file_1} > {output.con_file_2}
		fi
		'''
## Notes for rule below:
#if file doesn't exists (statement tests if there is something in it) and has something in it ("greater than zero") then touch creates file file with nothing in it.
#if the file does have something then nothing happens.
# || is an OR opporator and only occurs if the first statement is false.

rule concat_L002:
	input:
		checkpoint = config['workdir']+ 'checkpoint.txt',
		k_file_1 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L002_Kclean_R1.fastq',
		k_file_2 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L002_Kclean_R2.fastq'
	params:
		outdir = config['workdir'] + 'QC/Concat',
		kdir = config['workdir'] + 'QC/Kraken_Cleaned',
		ErrorOut = config['ErrorOut'] + 'Concat'
	output:
		con_file_1 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R1.fastq',
		con_file_2 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R2.fastq'
	message: '''--- Concating Sequences ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		cd /scicomp/home-pure/qpk9/SM_Test/QC/Kraken_Cleaned
		if [ -f {params.kdir}/{wildcards.nolane_sample}_L002_Kclean_R1.fastq ] #if file exists
		then
			[[ -s {output.con_file_1} ]] || touch {output.con_file_1}
			[[ -s {output.con_file_2} ]] || touch {output.con_file_2}
			cat {params.kdir}/{wildcards.nolane_sample}_L002_Kclean_R1.fastq >> {output.con_file_1}
			cat {params.kdir}/{wildcards.nolane_sample}_L002_Kclean_R2.fastq >> {output.con_file_2}
		fi
		'''

rule concat_L003:
	input:
		checkpoint = config['workdir']+ 'checkpoint.txt',
		k_file_1 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L003_Kclean_R1.fastq',
		k_file_2 = config['workdir'] + 'QC/Kraken_Cleaned/{nolane_sample}_L003_Kclean_R2.fastq'
	params:
		outdir = config['workdir'] + 'QC/Concat',
		kdir = config['workdir'] + 'QC/Kraken_Cleaned',
		ErrorOut = config['ErrorOut'] + 'Concat'
	output:
		con_file_1 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R1.fastq',
		con_file_2 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R2.fastq'
	message: '''--- Concating Sequences ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		cd /scicomp/home-pure/qpk9/SM_Test/QC/Kraken_Cleaned
		if [ -f {params.kdir}/{wildcards.nolane_sample}_L003_Kclean_R1.fastq ]
		then
			[[ -s {output.con_file_1} ]] || touch {output.con_file_1}
			[[ -s {output.con_file_2} ]] || touch {output.con_file_2}
			cat {params.kdir}/{wildcards.nolane_sample}_L003_Kclean_R1.fastq >> {output.con_file_1}
			cat {params.kdir}/{wildcards.nolane_sample}_L003_Kclean_R2.fastq >> {output.con_file_2}
		fi
		'''

##Assembly and QC checks
rule skesa:
	input:
		con_file_1 = config['workdir'] + 'QC/Concat/{sample}_Kclean_R1.fastq',
		con_file_2 = config['workdir'] + 'QC/Concat/{sample}_Kclean_R2.fastq'
	params:
		outdir = config['workdir'] + 'Assembly/Skesa',
		ErrorOut = config['ErrorOut'] + 'Skesa'
	output:
		config['workdir'] + 'Assembly/Skesa/{sample}_contigs_skesa.fa'
	message: '''--- Assemblying with Skesa ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		module load Skesa/2.3.0
		time skesa --fastq {input.con_file_1},{input.con_file_2} --cores 10 --memory 50 --contigs_out {output}
		'''

rule quast:
	input:
		config['workdir'] + 'Assembly/Skesa/{sample}_contigs_skesa.fa'
	params:
		outdir = config['workdir'] + 'Assembly/quast_results',
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta',
		ref_gff = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC.gff',
		ErrorOut = config['ErrorOut'] + 'Quast'
	output:
		config['workdir'] + 'Assembly/quast_results/{sample}/report.tsv',
		directory(config['workdir'] + 'Assembly/quast_results/{sample}')
	message: '''--- Running Quast to Evaluate Assemblies ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		module load quast/5.0
		time quast.py --eukaryote -o {params.outdir}/{wildcards.sample}/ -r {params.ref} {input} --features {params.ref_gff}
		'''

rule assembly_stats:
	input:
		expand(config['workdir'] + 'Assembly/quast_results/{sample}/report.tsv', sample=NL_FILE_NAME)
	params:
		outdir = config['workdir'] + 'Assembly/',
		indir = config['workdir'] + 'Assembly/quast_results/',
		scripts = config['scripts']
	output:
		config['workdir'] + 'Assembly/Assembly_stats.csv',
		probs = config['workdir'] + 'Assembly/Problem_Assemblies.txt'
	message: '''--- Checking Assemblies with Get_Assembly_Stats.py ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		time python {params.scripts}/Get_Assembly_Stats.py -i {params.indir} -o {params.outdir} > {output.probs}
		'''

rule checkpoint2:
	input:
		expand(config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R{read}.fastq', nolane_sample=NL_FILE_NAME, read=['1','2'])
	output:
		config['workdir'] + 'checkpoint2.txt'
	shell:
		'''
		echo {input} files are finished > checkpoint2.txt
		'''

rule make_bad_list_file:
	input: config['workdir'] + 'checkpoint2.txt'
	params:
		file_name = 'bad_list_file.txt',
		tree_dir = config['workdir'] + 'Variant_Call_Test/Tree',
	output:
		output_file = config['workdir'] + 'Variant_Call_Test/Tree/bad_list_file.txt',
	run:
		tree_dir = " ".join({params.tree_dir}) #make str from set
		output_file = " ".join({output.output_file})
		if not os.path.exists(tree_dir):
			os.makedirs(tree_dir)
		#some python to make my bad_list_file for the next rule:
		Bad_list = ['41897_Cattle', '39640_Human', '42129_Cattle', '40457_Human', '40236_Calf', '41573_Cow', '42134_Cattle', '41891_Cattle', '41892_Cattle']
		os.chdir(tree_dir)
		with open(output_file, 'w') as f:
			for item in Bad_list:
				f.write("%s\n" % item)

# Making Tree
rule ksnp:
	input:
		expand(config['workdir'] + 'Assembly/Skesa/{sample}_contigs_skesa.fa', sample=New_all),
		bad_list_file = config['workdir'] + 'Variant_Call_Test/Tree/bad_list_file.txt',
		outputdir = config['workdir'] + 'Variant_Call_Test/Tree',
	params:
		assembly_dir = config['workdir'] + 'Assembly',
		ErrorOut = config['ErrorOut'] + 'Tree',
		ref = "CryptoDB-50_CparvumIOWA-ATCC_Genome",
		ref_seq = config['ref'] + "CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta",
		bad_list_file_loc = '/scicomp/home-pure/qpk9/GWAS_Study/bad_list_file.txt',
	output:
		kSNP3_list = config['workdir'] + 'Variant_Call_Test/Tree/kSNP3_in_list_cleaned',
		contigs = config['workdir'] + 'Variant_Call_Test/Tree/Skesa_Contigs.fasta'
	message: '''--- Using de-novo assembly for Tree building ---'''
	shell:
		'''
		mkdir -p {params.ErrorOut}
		cd {params.assembly_dir}
		module load kSNP/3.021
		if [ ! -f {input.outputdir}/kSNP3_in_list ] #if file doesn't exists
		then
			MakeKSNP3infile Skesa {input.outputdir}/kSNP3_in_list A
			printf "{params.ref_seq}\t{params.ref}" >> {input.outputdir}/kSNP3_in_list
		fi
		if [ ! -f {output.kSNP3_list} ] #if file doesn't exists
		then
			grep -v -F -f {input.bad_list_file} {input.outputdir}/kSNP3_in_list > {output.kSNP3_list}
		fi
		if [ ! -f {output.contigs} ] #if file doesn't exists
		then
			#This will through an error because MakeFasta will prompt you to say continue or exit after first running it.
			MakeFasta {output.kSNP3_list} {output.contigs}
			#Kchooser {output.contigs}
		fi
		kSNP3 -k 21 -in {output.kSNP3_list} -CPU 15 -core -outdir {input.outputdir}
		'''

#Alignment for calculatng depth coverage at each position
##Aligning and mapping to look for areas of low coverage
rule bowtie2:
	input:
		con_file_1 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R1.fastq',
		con_file_2 = config['workdir'] + 'QC/Concat/{nolane_sample}_Kclean_R2.fastq',
		ref = config['ref']
	params:
		outdir = config['workdir'] + 'Alignment',
		ErrorOut = config['ErrorOut'] + 'Bowtie2'
	output:
		sam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.sam'
	message: '''--- Running Bowtie2 to align sequences to reference genome ---'''
	shell:
		'''
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		mkdir -p {params.outdir}/Sam_Bam/
		module load bowtie2/2.3.5.1
		cd {input.ref}
		time bowtie2 --end-to-end --threads 12 -x CparvumIOWA-ATCC -1 {input.con_file_1} -2 {input.con_file_2} -S {output.sam}
		'''

rule samtools:
	input:
		sam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.sam',
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta'
	params:
		dir = config['workdir'] + 'Alignment/Sam_Bam',
		ErrorOut = config['ErrorOut'] + 'Sam_Bam'
	output:
		bam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.bam',
		sorted_bam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.sorted.bam'
	message: '''--- Running Samtools to convert bam file from bowtie2 to sam file, then indexing and sorting ---'''
	shell:
		'''
		module load samtools/1.10
		echo 'I am running on node:' `hostname`
		mkdir -p {params.ErrorOut}
		mkdir -p {params.dir}
		cd {params.dir}
		## Converting from sam to bam format
		time samtools view -bS --reference {input.ref} -o {output.bam} {input.sam}
		## Running sort to speed up the indexing and viewing later
		time samtools sort {output.bam} -o {output.sorted_bam}
		## Indexing of bam file
		time samtools index {output.sorted_bam}
		'''

## Variant calling SM style
rule freebayes:
	input:
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta',
		sorted_bam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.sorted.bam'
	params:
		outdir = "/scicomp/home-pure/qpk9/SM_Test/Variant_Call_Test/Freebayes",
		ErrorOut = config['ErrorOut'] + 'Freebayes'
	output:
		vcf = config['workdir'] + 'Variant_Call_Test/Freebayes/{nolane_sample}-freebayes-SNP-differences.vcf'
	message: ''' --- Running Freebayes to get variant calls ---- '''
	shell:
		'''
		source activate freebayes
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		freebayes --min-base-quality 20 --ploidy 1 --min-coverage 20 --min-alternate-fraction 0.75 --use-mapping-quality -f {input.ref} {input.sorted_bam} -v {output.vcf}
		'''

rule remove_indels:
	input:
		vcf = config['workdir'] + 'Variant_Call_Test/Freebayes/{nolane_sample}-freebayes-SNP-differences.vcf'
	params:
		no_indels_file_basename = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-freebayes-SNP-differences-NO-INDELS',
		outdir = config['workdir'] + 'Variant_Call_Test/VCF',
		ErrorOut = config['ErrorOut'] + 'VCF'
	output:
		no_indels_file = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-freebayes-SNP-differences-NO-INDELS.recode.vcf'
	message: ''' --- Running VCF Tools to remove indels --- '''
	shell:
		'''
		module load vcftools/0.1.17
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		vcftools --vcf {input.vcf} --remove-indels --recode --recode-INFO-all --out {params.no_indels_file_basename}
		'''

rule filter_vcf:
	input:
		vcf = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-freebayes-SNP-differences-NO-INDELS.recode.vcf'
	output:
		filtered_vcf = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf'
	message: ''' --- Running VCF filter to drop sites based on quality --- '''
	shell:
		'''
		module load vcflib/2019.12.31
		vcffilter -f "QUAL > 1" {input.vcf} > {output.filtered_vcf}
		'''

rule zip_vcf:
	input:
		filtered_vcf = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf'
	output:
		zipped_vcf = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz'
	message: ''' --- Zipping vcf --- '''
	shell:
		'''
		module load htslib/1.10
		bgzip {input.filtered_vcf}
		tabix --force --preset vcf {output.zipped_vcf}
		'''

rule vcf_merge:
	input:
		all_vcf_files = expand(config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz', nolane_sample=New_all),
		human_vcf_files = expand(config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz', nolane_sample=Human_samples),
		cow_vcf_files = expand(config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz', nolane_sample=Cow_samples)
	output:
		all_vcf_merge = config['workdir'] + 'Variant_Call_Test/VCF/all_merge.vcf.gz',
		human_vcf_merge = config['workdir'] + 'Variant_Call_Test/VCF/human_merge.vcf.gz',
		cow_vcf_merge = config['workdir'] + 'Variant_Call_Test/VCF/cow_merge.vcf.gz',
	message: ''' --- Running VCF Tools to merge vcf files --- '''
	shell:
		'''
		module load vcftools/0.1.17
		vcf-merge --remove-duplicates {input.all_vcf_files} | bgzip -c > {output.all_vcf_merge}
		vcf-merge --remove-duplicates {input.human_vcf_files} | bgzip -c > {output.human_vcf_merge}
		vcf-merge --remove-duplicates {input.cow_vcf_files} | bgzip -c > {output.cow_vcf_merge}
		'''

rule bcftools:
	input:
	params:
		vcf_file_dir = config['workdir'] + 'Variant_Call_Test/VCF'
	output:
	shell:
		'''
		# make file of gz names
		if [ ! -f {params.vcf_file_dir}/vcffiltered_files ] #if file exists
		then
			cd {params.vcf_file_dir}
			for f in *vcffilter.vcf.gz; do echo $f > vcffiltered_files ; done
		fi
		module load bcftools/1.10.2
		bcftools stats -S vcffiltered_files
		'''

rule consensus:
	input:
		zipped_vcf = config['workdir'] + 'Variant_Call_Test/VCF/{nolane_sample}-vcffilter.vcf.gz',
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta'
	params:
		scripts = config['scripts'],
		outdir = config['workdir'] + 'Variant_Call_Test/Consensus/'
	output:
		consensus_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}.fasta',
		one_line_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-one-line.fasta'
	message: ''' --- Creating consensus fasta file --- '''
	shell:
		"""
		module load vcftools/0.1.17
		mkdir -p /scicomp/home-pure/qpk9/SM_Test/Error_Out_Files/Consensus
		mkdir -p {params.outdir}
		cat {input.ref} | vcf-consensus {input.zipped_vcf} > {output.consensus_fasta}
		# Convert 60 line character fasta file to one line fasta file
		python {params.scripts}/UnwrapFasta.py -f {output.consensus_fasta} -o {output.one_line_fasta}
		"""

rule bedtools:
	input:
		sorted_bam = config['workdir'] + 'Alignment/Sam_Bam/{nolane_sample}.E2E.sorted.bam',
		one_line_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-one-line.fasta'
	params:
		outdir = "/scicomp/home-pure/qpk9/SM_Test/Variant_Call_Test/Consensus",
		scripts = config['scripts'],
		ref_gff = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC.gff',
		ErrorOut = config['ErrorOut'] + 'BEDTools'
	output:
		coverage = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-per-site-depth-coverage.txt',
		masked_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq.fasta',
		genes_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-genes.fasta'
	message: '''--- Running bedtools to get coverage at each position of the genome, mask low coverage sites and identify gene locations. '''
	shell:
		'''
		module load BEDTools/2.27.1
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		bedtools genomecov -bga -split -ibam {input.sorted_bam} > {output.coverage}
		# Use python script to mask sites below 20x depth coverage
		python {params.scripts}/maskSites.py -d 20 -f {input.one_line_fasta} -c {output.coverage} -o {output.masked_fasta}  ## check masking site scripts
		# Use Bedtools to identify gene sequences and their start/end coordinates
		bedtools getfasta -fi {output.masked_fasta} -bed {params.ref_gff} -name -fo > {output.genes_fasta}
		'''

rule depth_check:
	input:
		expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-per-site-depth-coverage.txt', nolane_sample=NL_FILE_NAME)
	params:
		depth_files_dir = config['workdir'] + 'Variant_Call_Test/Consensus/',
		scripts_dir = config['scripts']
	output:
		config['workdir'] + 'Variant_Call_Test/Consensus/Depth_stats.csv',
		output_file = config['workdir'] + 'Variant_Call_Test/Consensus/depth_check_output.txt'
	shell:
		'''
		python {params.scripts_dir}/CheckDepth.py -d 20 -c {params.depth_files_dir} > {output.output_file}
		'''

rule consensus_matrix:
	input:
		input_files = expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq.fasta', nolane_sample=New_all)
	params:
		input_dir = config['workdir'] + 'Variant_Call_Test/Consensus/',
		scripts = config['scripts']
	output:
		loc = config['workdir'] + 'Variant_Call_Test/Consensus/location.csv',
		meta = config['workdir'] + 'Variant_Call_Test/Consensus/metadata.csv',
		snp_df_num = config['workdir'] + 'Variant_Call_Test/Consensus/snp_df_num.csv',
		snp_df = config['workdir'] + 'Variant_Call_Test/Consensus/snp_df.csv'
	shell:
		'''
		python {params.scripts}/Consensus_to_matrix.py -d {params.input_dir}
		'''

#making snp-dist matrix
rule copy:
	input:
		fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq.fasta'
	output:
		fasta_renamed = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq-RN.fasta'
	shell:
		'''
		cp {input.fasta} {output.fasta_renamed}
		'''

rule rename:
	input:
		fastas = expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq-RN.fasta', nolane_sample=New_all)
	params:
		outputdir = config['workdir'] + 'Variant_Call_Test/Consensus/',
		ErrorOut = config['ErrorOut'] + 'Tree'
	output:
		combined_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/checkpoint3.txt'
	shell:
		'''
		mkdir -p {params.outputdir}
		mkdir -p {params.ErrorOut}
		cd {params.outputdir}
		#add the file name to the sequence ID
		for f in *-masked-seq-RN.fasta; do fname=$(basename $f -masked-seq-RN.fasta); sed -i "s/Chr[0-9]/${{fname}}/" $f; done
		echo "Finished Rename" > checkpoint3.txt
		'''

rule combine:
	input:
		fasta = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq-RN.fasta'
	output:
		fasta_oneline = config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq-RN-oneline.fasta'
	shell:
		'''
		#combine into one line
		cat {input.fasta} | sed -e '1!{{/^>.*/d;}}' | sed ':a;N;$!ba;s/\\n//2g' | sed -e '$a\\' > {output.fasta_oneline}
		'''

rule combine2:
	input:
		all_fastas = expand(config['workdir'] + 'Variant_Call_Test/Consensus/{nolane_sample}-masked-seq-RN-oneline.fasta', nolane_sample=New_all)
	output:
		combined_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/combined-masked-seq-RN.fasta'
	shell:
		'''
		#create multi-line fasta with all sequences from sample
		cat {input.all_fastas} > {output.combined_fasta}
		'''

rule snp_dist_consensus:
	input:
		combined_fasta = config['workdir'] + 'Variant_Call_Test/Consensus/combined-masked-seq-RN.fasta'
	output:
		distances = config['workdir'] + 'Variant_Call_Test/distances.tab'
	shell:
		'''
		module load snp-dists/0.7
		snp-dists {input.combined_fasta} > {output.distances}
		'''

##Variant calling with Snippy
rule generate_file_for_snippy:
	input:
		checkpoint = config['workdir' ] + 'checkpoint2.txt',
	params:
		indir = config['workdir'] + 'QC/Concat',
		parentoutdir = config['workdir'] + 'VariantCalls',
		outdir = config['workdir'] + 'VariantCalls/Snippy',
		scripts = config['scripts']
	output:
		outfile_1 = config['workdir'] + 'VariantCalls/Snippy/Human_files.tab',
		outfile_2 = config['workdir'] + 'VariantCalls/Snippy/Cow_files.tab',
		outfile_3 = config['workdir'] + 'VariantCalls/Snippy/All_files.tab'
	shell:
		'''
		mkdir -p {params.parentoutdir}
		mkdir -p {params.outdir}
		time python {params.scripts}/Making_SampleFile_Snippy.py -i {params.indir} -o {params.outdir}
		'''

rule snippy:
	input:
		file_1 = config['workdir'] + 'QC/Concat/{sample}_Kclean_R1.fastq',
		file_2 = config['workdir'] + 'QC/Concat/{sample}_Kclean_R2.fastq',
	params:
		outdir = config['workdir'] + 'VariantCalls/Snippy',
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta',
		ErrorOut = config['ErrorOut'] + 'Snippy'
	output:
		multiext(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps', ".aligned.fa", ".bam", ".bam.bai", ".bed", ".consensus.fa", ".consensus.subs.fa", ".csv", ".filt.vcf", ".gff", ".html", ".log", ".raw.vcf", ".subs.vcf", ".tab", ".txt", ".vcf", ".vcf.gz", ".vcf.gz.csi")
	shell:
		'''
		mkdir -p {params.ErrorOut}
		mkdir -p {params.outdir}
		cd {params.outdir}
		conda activate Snippy
		time snippy --outdir {wildcards.sample} --mincov 20 --basequal 20 --R1 {input.file_1} --R2 {input.file_2} --ref {params.ref} --cpus 16 --force
		'''

rule snippy_core:
	input:
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.aligned.fa', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.bam', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.bam.bai', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.bed', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.consensus.fa', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.consensus.subs.fa', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.csv', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.filt.vcf', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.gff', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.html', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.log', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.raw.vcf', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.subs.vcf', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.tab', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.txt', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.vcf', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.vcf.gz', sample=New_all),
		expand(config['workdir'] + 'VariantCalls/Snippy/{sample}/snps.vcf.gz.csi', sample=New_all),
	params:
		outdir = config['workdir'] + 'VariantCalls/Snippy',
		all_files = New_all,
		Cow_files = Cow_samples,
		Human_files = Human_samples,
		ref = config['ref'] + 'CryptoDB-50_CparvumIOWA-ATCC_Genome.fasta'
	output:
		multiext(config['workdir'] + 'VariantCalls/Snippy/All_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf"),
		multiext(config['workdir'] + 'VariantCalls/Snippy/Human_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf"),
		multiext(config['workdir'] + 'VariantCalls/Snippy/Cow_core', ".ref.fa", ".tab", ".aln", ".txt", ".full.aln", ".vcf")
	shell:
		'''
		cd {params.outdir}
		conda activate Snippy
		time snippy-core --ref '{params.ref}' {params.all_files} --prefix=All_core
		time snippy-core --ref '{params.ref}' {params.Human_files} --prefix=Human_core
		time snippy-core --ref '{params.ref}' {params.Cow_files} --prefix=Cow_core
		'''

rule snp_dist:
	input:
		alignment_full = config['workdir'] + 'VariantCalls/Snippy/All_core.full.aln',
		alignment_human = config['workdir'] + 'VariantCalls/Snippy/Human_core.full.aln',
		alignment_cow = config['workdir'] + 'VariantCalls/Snippy/Cow_core.full.aln'
	output:
		distances_full = config['workdir'] + 'VariantCalls/distances_all.tab',
		distances_human = config['workdir'] + 'VariantCalls/distances_human.tab',
		distances_cow = config['workdir'] + 'VariantCalls/distances_cow.tab'
	shell:
		'''
		module load snp-dists/0.7
		snp-dists {input.alignment_full} > {output.distances_full}
		snp-dists {input.alignment_human} > {output.distances_human}
		snp-dists {input.alignment_cow} > {output.distances_cow}
		'''
