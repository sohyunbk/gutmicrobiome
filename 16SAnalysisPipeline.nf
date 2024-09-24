params.reads1 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_1.fq"
params.reads2 = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_2.fq"

workflow{
	FastQC | view
}

process FastQC{
    input:
    path reads1 from params.reads1
    path reads2 from params.reads2

    output:
    path "*.html" into fastqc_results
    path "*.zip" into fastqc_results

    script:
    """
    fastqc $reads1
    fastqc $reads2
    """
    
    publishDir reads1.parent, mode: 'copy'
}


