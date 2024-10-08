params.dir = "/scratch/sb14489/10.Metagenome/"
params.reads = "/scratch/sb14489/10.Metagenome/1.RawData/*_{1,2}.fastq.gz"
params.threads = 4
params.qiime = "ON"
workflow {
     // Upto Merging Reads
    if [ $params.qiime != "ON" ]; then
        read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        FastQC(read_pairs_ch)
        trimmed_reads_ch = Trimmomatic(read_pairs_ch)
        MergeReads(trimmed_reads_ch)

    // Qiime2
    else
        Writing_fastqManifest


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

    publishDir "$params.dir/2.Trimmomatic/", mode: 'copy'  

    script:
    """
    mkdir -p $params.dir/2.Trimmomatic/
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \\
        -threads $params.threads \\
        ${reads[0]} ${reads[1]} \\
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
    """
    mkdir -p $params.dir/3.Pear/
    pear -f $trimmed1 -r $trimmed2 -j $params.threads \\
        -o $params.dir/3.Pear/${pair_id} > $params.dir/3.Pear/${pair_id}.log
    """
}

process Writing_fastqManifest{

    """
    #!/home/sb14489/miniconda3/envs/Spatial/bin/python


    """
}