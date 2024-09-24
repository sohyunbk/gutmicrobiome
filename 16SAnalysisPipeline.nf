params.dir = "/scratch/sb14489/10.Metagenome/"
params.samplename = "Leaf0.5_Re1"
params.threads = 4

workflow {
    Channel
        .fromFilePairs([params.dir + "/1.RawData/" + params.samplename +"_1.fastq.gz", params.dir + "/1.RawData/" + params.samplename +"_2.fastq.gz"])
        .set { raw_reads }

    raw_reads | FastQC | Trimmomatic | view
    }

process FastQC {
    input:
    tuple val(reads1), val(reads2)

    output:
    tuple val(reads1), val(reads2)

    script:
    """
    fastqc --threads $params.threads $reads1 -o $params.dir/1.RawData/
    fastqc --threads $params.threads $reads2 -o $params.dir/1.RawData/
    """
}

process Trimmomatic {
    input:
    tuple val(reads1), val(reads2)

    output:
    tuple val(reads1), val(reads2)

    script:
    """
	echo "The value of reads1 is: $reads1"
    mkdir -p $params.dir/2.Trimmomatic/
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        $params.dir/1.RawData/$params.reads1 $params.dir/1.RawData/$params.reads2 \\
        $params.dir/2.Trimmomatic/${params.samplename}_1_trimmed_paired.fq.gz $params.dir/2.Trimmomatic/${params.samplename}_1_trimmed_unpaired.fq.gz \\
        $params.dir/2.Trimmomatic/${params.samplename}_2_trimmed_paired.fq.gz $params.dir/2.Trimmomatic/${params.samplename}_2_trimmed_unpaired.fq.gz \\
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process MergeReads {
    input:
    tuple val(reads1), val(reads2)

    output:
    stdout

    script:
    def samplename = params.reads1.replace("_1.fastq.gz", "")
    """
    mkdir -p $params.dir/3.Pear/
    pear -f $reads1 -r $reads2 -j $params.threads \\
    -o $params.dir/3.Pear/${samplename} > $params.dir/3.Pear/${samplename}.log
    """
}
