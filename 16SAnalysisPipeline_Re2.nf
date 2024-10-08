params.dir = "/scratch/sb14489/10.Metagenome/"
params.metadata = "/scratch/sb14489/10.Metagenome/1.RawData/Metadata.txt"
params.threads = 10
params.commName = "ChickenGut"
params.reads = "${params.dir}/1.RawData/*_{1,2}.fastq.gz"
params.mergedFiles = "${params.dir}/3.Pear/*.assembled.fastq"
params.preprocessing=true
params.featuretable=false
params.qiimeDiversity=false
params.figures=false

workflow {
    if (params.preprocessing) {
        read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        FastQC(read_pairs_ch)
        trimmed_reads_ch = Trimmomatic(read_pairs_ch)
        MergeReads(trimmed_reads_ch)
    }
    if (params.featuretable){
        all_files_ch = Channel.fromPath(params.mergedFiles, checkIfExists: true).collect()
        manifestfile = Writing_fastqManifest(all_files_ch)
        demultiplexed_QZA = Making_MultiflexedQZAFile(manifestfile)
        outputQZAs = Denoising_ChimeraRemov_ASV_dada2_QZAFile(demultiplexed_QZA)
        tableQZA = outputQZAs[0]
        repseqQZA = outputQZAs[1]
        filteredQZA = MakeQZVs_FilteringMissingSamples(tableQZA,repseqQZA)
    }

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

process Writing_fastqManifest {
    input:
    path assembled_files

    output:
    path "${params.commName}_manifest_33.txt"
    publishDir "$params.dir/4.Importing/", mode: 'copy'

    script:
    """
    for sFile in ${assembled_files}; do
        echo "\${sFile}"
    done
    
    mkdir -p ${params.dir}/4.Importing/
    
    echo "# single-end PHRED 33 fastq manifest file for forward reads" > "${params.commName}_manifest_33.txt"
    echo "sample-id,absolute-filepath,direction" >> "${params.commName}_manifest_33.txt"
    
    # Loop through each file in assembled_files and write to the manifest
    for sFile in ${assembled_files}; do
        # Extract the base name and strip the '.assembled.fastq' extension
        FileName=\$(basename "\${sFile}" | sed 's/\\.assembled\\.fastq//')
        # Write to the manifest file
        echo "\${FileName},\$(realpath \${sFile}),forward" >> "${params.commName}_manifest_33.txt"
    done
    """
}


process Making_MultiflexedQZAFile{
    input:
    path "${params.commName}_manifest_33.txt"

    output:
    path "${params.commName}_demultiplexed.qza"
    publishDir "$params.dir/4.Importing/", mode: 'copy'

    script:
    """
    qiime tools import --type 'SampleData[SequencesWithQuality]' \\
    --input-path "${params.commName}_manifest_33.txt" \\
     --output-path "${params.commName}_demultiplexed.qza" \\
     --input-format SingleEndFastqManifestPhred33
    """

}

process Denoising_ChimeraRemov_ASV_dada2_QZAFile{
    input:
    path "${params.commName}_demultiplexed.qza"

    output:
    path "${params.commName}_table.qza"
    path "${params.commName}_rep_seqs.qza"
    path "${params.commName}_stats-dada2.qza"
    publishDir "$params.dir/5.AfterQC/", mode: 'copy'

    script:
    """
    mkdir -p ${params.dir}/5.AfterQC/
    qiime dada2 denoise-single --p-n-threads ${params.threads} \\
     --i-demultiplexed-seqs "${params.commName}_demultiplexed.qza"  \\
     --p-trunc-len 0  --p-trim-left 0 --o-representative-sequences \\
      ${params.commName}_rep_seqs.qza --o-table \\
       ${params.commName}_table.qza  \\
     --o-denoising-stats ${params.commName}_stats-dada2.qza
    """

}

process MakeQZVs_FilteringMissingSamples{
    input:
    path "${params.commName}_table.qza"
    path "${params.commName}_rep_seqs.qza"

    output:
    path "${params.commName}_filtered-table.qza"
    publishDir "$params.dir/5.AfterQC/", mode: 'copy'

    script:
    """
    qiime feature-table summarize --i-table  "${params.commName}_table.qza" \\
     --o-visualization "${params.commName}_table.qzv" \\
      --m-sample-metadata-file $params.metadata 

    qiime feature-table tabulate-seqs --i-data "${params.commName}_rep_seqs.qza" \\ 
    --o-visualization "${params.commName}_rep_seqs.qzv"

    ##ignore missing samples
    qiime feature-table filter-samples --i-table ./4.QualityControl/table.qza  \\
    --m-metadata-file ./0.MetaInfo/sample_metadata.tsv \\
    --o-filtered-table ./4.QualityControl/filtered-table.qza
    """
}

process DifferentialAbundanceTest_EdgeR{
    input:
    stdin

    output:
    stdout

    """
    python "${params.ScriptDir}/DifferentialTest_EdgeR.py --input_table --output_name "
    """
    }


process StatckedBarPlot{
    input:
    stdin

    output:
    stdout

    """
    python "${params.ScriptDir}/StatckedBarPlot_ggplot.py --input_table --output_name "
    """
    }
