params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fq"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fq"

workflow{
	FastQC
}

process FastQC{
	output:
	stdout

	script:
	"""
	fastqc $params.reads1
	fastqc $params.reads2
	"""
}


