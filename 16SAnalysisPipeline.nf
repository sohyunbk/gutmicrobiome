params.dir = "/scratch/sb14489/10.Metagenome/"
params.reads1 = "Leaf0.5_Re1_1.fastq.gz"
params.reads2 = "Leaf0.5_Re1_2.fastq.gz"
params.threads = 4

workflow {
    FastQC | view
}

process FastQC {
    output:
    stdout

    script:
    """
	echo ** FASTAQC **
    fastqc --threads $params.threads  $params.dir/1.RawData/$params.reads1 -o $params.dir/1.RawData/ 
    fastqc --threads $params.threads  $params.dir/1.RawData/$params.reads2 -o $params.dir/1.RawData/
    """
}

process Trimmomatic {
    output:
    def output_dir = file(params.reads1).parent  
    def reads1_base = file(params.reads1).replace('.fastq.gz', '') 
    def reads2_base = file(params.reads1).replace('.fastq.gz', '')

	script:
	"""
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        $params.reads1 $params.reads2 \\
        trimmed_${reads1_base}_paired.fq.gz trimmed_${reads1_base}_unpaired.fq.gz \\
        trimmed_${reads2_base}_paired.fq.gz trimmed_${reads2_base}_unpaired.fq.gz \\
        ILLUMINACLIP:/path/to/adapters.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
