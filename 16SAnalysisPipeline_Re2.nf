params.dir = "/scratch/sb14489/10.Metagenome/"

params.threads = 4

params.reads = "${params.dir}/1.RawData/*_{1,2}.fastq.gz"
params.mergedFiles = "${params.dir}/3.Pear/*.assembled.fastq"
workflow {
     // Upto Merging Reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    FastQC(read_pairs_ch)
    trimmed_reads_ch = Trimmomatic(read_pairs_ch)
    MergeReads(trimmed_reads_ch)

    // Qiime2
    all_files_ch = Channel.fromPath(params.mergedFiles, checkIfExists: true)
    Writing_fastqManifest(all_files_ch)


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
    input:
    path assembled_files

    output:
    path "manifest_33.txt" 
    publishDir "$params.dir/4.Importing/", mode: 'copy'  
    
    script:
    """
    echo $assembled_files
    mkdir -p $params.dir/4.Importing/
    echo "# single-end PHRED 33 fastq manifest file for forward reads" > "manifest_33.txt"
    echo "sample-id,absolute-filepath,direction" >> "manifest_33.txt"
    
    for sFile in ${assembled_files}; do
        # Extract the base name and strip the '.assembled.fastq' extension
        FileName=\$(basename "\$sFile")
        FileName="\${FileName%.assembled.fastq}"

        # Write to the manifest file
        echo "\$FileName,\$sFile,forward" >> "manifest_33.txt"
    done
    """
}
