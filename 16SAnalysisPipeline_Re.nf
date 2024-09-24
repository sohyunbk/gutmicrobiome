params.dir = "/scratch/sb14489/10.Metagenome/"
params.reads = "/scratch/sb14489/10.Metagenome/1.RawData/Leaf0.5_Re1_{1,2}.fastq.gz"
params.threads = 4

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    FastQC(read_pairs_ch) 
    trimmed_reads_ch = Trimmomatic(read_pairs_ch) 
    MergeReads(trimmed_reads_ch)

    }

process FastQC {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(reads)

    output:
    stdout

    script:
    """
    fastqc --threads $params.threads ${reads[0]} -o $params.dir/1.RawData/
    fastqc --threads $params.threads ${reads[1]} -o $params.dir/1.RawData/
    """
}
process Trimmomatic {
    tag "$pair_id"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_1_trimmed_paired.fq.gz"), path("${pair_id}_2_trimmed_paired.fq.gz")

    publishDir "$params.dir/2.Trimmomatic/", mode: 'copy'  // Move the outputs to the desired directory

    script:
    """
    mkdir -p $params.dir/2.Trimmomatic/
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        $params.dir/1.RawData/${reads[0]} $params.dir/1.RawData/${reads[1]} \\
        ${pair_id}_1_trimmed_paired.fq.gz ${pair_id}_1_trimmed_unpaired.fq.gz \\
        ${pair_id}_2_trimmed_paired.fq.gz ${pair_id}_2_trimmed_unpaired.fq.gz \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    """
}

process MergeReads {
    input:
    tuple val(pair_id), path(trimmed1), path(trimmed2)

    output:
    stdout

    script:
    def samplename = params.reads1.replace("_1.fastq.gz", "")
    """
    mkdir -p $params.dir/3.Pear/
    pear -f $trimmed1 -r $trimmed2 -j $params.threads \\
        -o $params.dir/3.Pear/${pair_id} > $params.dir/3.Pear/${pair_id}.log
    """
}
