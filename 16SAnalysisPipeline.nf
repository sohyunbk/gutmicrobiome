params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fq"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fq"

workflow {
    FastQC | view
}

process FastQC {
    input:
    path reads1 = file(params.reads1)
    path reads2 = file(params.reads2)

    script:
    def output_dir1 = reads1.parent
    def output_dir2 = reads2.parent

    """
    fastqc $reads1 -o $output_dir1
    fastqc $reads2 -o $output_dir2
    """
}