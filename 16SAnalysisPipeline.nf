params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fq"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fq"

workflow{
	FastQC | view
}

process FastQC{
	output:
	stdout

	script:
	"""
	fastqc $params.reads1 -o file(params.reads1).parent
	fastqc $params.reads2 -o file(params.reads2).parent
	"""
}


