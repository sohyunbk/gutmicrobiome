params.dir = "/scratch/sb14489/10.Metagenome/"
params.reads1 = "Leaf0.5_Re1_1.fastq.gz"
params.reads2 = "Leaf0.5_Re1_2.fastq.gz"
params.threads = 4

workflow {
    FastQC | Trimmomatic | MergeReads | view
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
	input:
    stdin

    output:
    stdout

	script:
	def re1Name = params.reads1.replace(".fastq.gz", "")
	def re2Name = params.reads2.replace(".fastq.gz","")
	"""
	mkdir -p $params.dir/2.Trimmomatic/
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        $params.dir/1.RawData/$params.reads1 $params.dir/1.RawData/$params.reads2 \\
        $params.dir/2.Trimmomatic/${re1Name}_trimmed_paired.fq.gz $params.dir/2.Trimmomatic/${re1Name}_trimmed_unpaired.fq.gz \\
        $params.dir/2.Trimmomatic/${re2Name}_trimmed_paired.fq.gz $params.dir/2.Trimmomatic/${re2Name}_trimmed_unpaired.fq.gz \\
        ILLUMINACLIP:/path/to/adapters.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process MergeReads {
	input:
    stdin

    output:
    stdout

	script:
	def samplename = params.reads1.replace("_1.fastq.gz", "")
	"""
	mkdir -p $params.dir/3.Pear/
	pear -f $params.dir/2.Trimmomatic/${re1Name}_trimmed_paired.fq.gz \\
	 -r $params.dir/2.Trimmomatic/${re2Name}_trimmed_paired.fq.gz \\
	  -j $params.threads \\
	  -o $params.dir/3.Pear/${samplename} > $params.dir/3.Pear/${samplename}.log
 
	"""

}