
""" Config = name conventions and paths to dependencies

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

"""

class Config:
	dirtree = '''
	.
	├── python
	└── task_dna
	'''
	fastqc = "/home/ivana/third/FastQC/fastqc"
	unzip  = "/usr/bin/unzip"

	cutadapt = "/usr/local/bin/cutadapt"
	paired_reads_rootnames = ["AH_S1_L001", "CH_S2_L001"]
	read_labels = ["R1", "R2"]

	# TODO how do you make this systematic
	# this is from here https://www.biostars.org/p/371399/
	adapter_seq = {"Illumina Universal Adapter": "AGATCGGAAGAG"}

	bwa = "/home/ivana/third/bwa-0.7.17/bwa"
	reference_fasta = "/storage/databases/ucsc/goldenpath/hg19/hg19.fa"
	samtools = "/home/ivana/third/samtools-1.11/samtools"
	bcftools = "/home/ivana/third/bcftools-1.11/bcftools"

	basic_uniprot = "/storage/databases/uniprot/uniprot_basic_info.tsv"
