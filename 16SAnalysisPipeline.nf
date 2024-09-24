params.dir = "/scratch/sb14489/10.Metagenome/"
params.reads1 = "Leaf0.5_Re1_1.fastq.gz"
params.reads2 = "Leaf0.5_Re1_2.fastq.gz"
params.threads = 4

workflow {
    FastQC | Trimmomatic | view
}

process FastQC {
    output:
    stdout

    script:
    """
	fastqc --threads $params.threads  $params.dir/1.RawData/$params.reads1 -o $params.dir/1.RawData/ 
    fastqc --threads $params.threads  $params.dir/1.RawData/$params.reads2 -o $params.dir/1.RawData/
    """
}

process Trimmomatic {
    output:
    stdout

	script:
	def re1Name = $params.reads1.replace(".fastq.gz","")
	def re2Name = $params.reads2.replace(".fastq.gz","")
	"""
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        $params.dir/1.RawData/$params.reads1 $params.dir/1.RawData/$params.reads2 \\
        ${re1Name}_trimmed_paired.fq.gz ${re1Name}_trimmed_unpaired.fq.gz \\
        ${re2Name}_trimmed_paired.fq.gz ${re2Name}_trimmed_unpaired.fq.gz \\
        ILLUMINACLIP:/path/to/adapters.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
