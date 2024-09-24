params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fq"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fq"

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
	echo fastqc $params.reads1 -o $output_dir1
    fastqc $params.reads1 -o $output_dir1
    fastqc $params.reads2 -o $output_dir2
    """
}