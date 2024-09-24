params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fastq.gz"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fastq.gz"
params.threads = 4

workflow {
    FastQC | view
}

process FastQC {
    output:
    stdout

    script:
    def output_dir1 = file(params.reads1).parent
    def output_dir2 = file(params.reads2).parent

    """
	echo ** FASTAQC **
	echo fastqc $params.reads1 -o $output_dir1
	echo fastqc $params.reads2 -o $output_dir2
    fastqc --threads $params.threads $params.reads1 -o $output_dir1
    fastqc --threads $params.threads $params.reads2 -o $output_dir2
    """
}