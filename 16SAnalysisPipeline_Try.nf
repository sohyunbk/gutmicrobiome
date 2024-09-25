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
params.

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
    path "${params.commName}_featureTable"

    output:
    stdout

    """
    Rscript "${params.ScriptDir}/DifferentialTest_EdgeR.R \\
     --SampleName  "${params.commName}"   
     --inputPath "${params.commName}_featureTable" \\
     --outputPath "${params.dir}/9.AbundanceTest" \\
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


process DifferentialAbundanceTest_EdgeR{
    input:
    path "${params.commName}_featureTable.tsv"
    path 

    output:
    stdout

    """
    #!/home/sb14489/miniconda3/envs/Metagenomics_qiime/bin/R
    library(edgeR)
    library(ggplot2)
    library(limma)
    library(DESeq2)
    library(edgeR)
    library(ggplot2)
    library(gridExtra)

    option_list = list(
        make_option(c("--SampleName"), type="character",
                    help="Type SampleName", metavar="character"),
        make_option(c("--inputPath"), type="character",
                    help="Type InputPath"),
        make_option(c("--outputPath"), type="character",
                    help=Type OutputPath"),
        make_option(c("--MetaData"), type="character",
                    help="MetaData is needed for target object")       

        );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);
    Datalist <- list()
    InputDir <- 
    Taxa_All <- c("Phylum","Class","Order","Family","Genus","Minimum_Species")
    file_names <- paste0(opt$inputPath, "/FeatureTable_L", 2:7, ".tsv")

    Input <- Datalist[[nNumb]]
    Input_subset <- subset(Input,subset=c(Input$SpecificTaxa!="Unassigned(Kingdom)"))
    Data <- subset(Input_subset,subset=c(Input2$SpecificTaxa!="Bacteria(Kingdom)"))

    CountData_Raw <- Data[,3:ncol(Data)]
    head(CountData_Raw)
    head(Data)

    Taxa <- Data[,1]
    Taxonomy <- Data[,1]
    unique(Taxa)
    Taxa_all <- Data[,2]
    Taxonomy_all <- Data[,2]

    CountData_Raw <- CountData_Raw[,as.character(NameOrder)]
    head(CountData_Raw)




###################################################
###################################################
######### 1) Leaf 
###################################################
###################################################
head(CountData_Raw)


CountData <-CountData_Raw[,c(1:16)]
head(CountData)

########MetaData############


Combination<- c(rep("Negative control",4),
rep("Positive control",4),
rep("Leaf_0.3",4),
rep("Leaf_0.5",4))

Combination <- factor(Combination,
levels=c("Negative control","Positive control","Leaf_0.3","Leaf_0.5"))

Tissue <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Leaf",8))

Tissue <- factor(Tissue,
levels=c("Negative control","Positive control","Leaf"))


MetaData <- data.frame(colnames(CountData),Combination,Tissue)
str(MetaData)


targets <- MetaData ##Metadata
design <- model.matrix(~0+Tissue, data=targets)
head(CountData)

##################################################################
######## DAO : Fitting Treatment
##################################################################


y <- DGEList(counts=CountData, gene=Taxa)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
plotBCV(y)
fit <- glmFit(y, design)
TMM <- cpm(y, normalized.lib.sizes=TRUE,log=T)
TMMcolName <- colnames(TMM)
TMMcolName <- paste("TMM_",TMMcolName,sep="")
str(TMMcolName)
colnames(TMM) <- as.factor(TMMcolName)



##################################################################
######## DAO1 : Negative Control vs leaf
##################################################################
head(design)


lrt <- glmLRT(fit, contrast=c(1,0,-1))
Result1 <- topTags(lrt, n=dim(CountData)[1], sort.by="none")$table
colnames(Result1)<-c("Taxa","logFC","logCPM","LR","PValue","FDR")

head(Result1)
nrow(CountData)[1]
print(sum(Result1$PValue < 0.05))

Result_table1 <- cbind(Taxa_all,Result1,CountData,TMM)
head(Result_table1)

OutputName <- paste(Taxa_All[nNumb],"Negative_Control_vs_Leaf.csv",sep="_") 
write.csv(Result_table1, OutputName, quote=F, row.names=F)


### Drawing box-plot for top 10 DEGs
print(length(which(Result_table1$PValue < 0.05)))

Index <- order(Result_table1$PValue)[1:12]
OrderedTable <- Result_table1[Index,]
OrderedTMM <- TMM[Index,]

Index <-which(OrderedTable$PValue < 0.05)

plot_data <- data.frame(Abundance = OrderedTMM[Index[1],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[1]])

for(i in 2:length(Index)){
  temp <- data.frame(Abundance = OrderedTMM[Index[i],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[i]])

  
  plot_data <- rbind(plot_data, temp)
}


str(plot_data)
dim(plot_data)

plot_data_subset <- subset(plot_data,subset=(Condition == "Negative control" | 
Condition == "Leaf_0.3" | Condition == "Leaf_0.5"))
plot_data_subset$Taxonomy <- gsub("_","\n_",plot_data_subset$Taxonomy)
ggplot(data=plot_data_subset, aes(x=Condition, y=Abundance)) +
  geom_boxplot(aes(fill = Condition)) +
  facet_wrap(~ Taxonomy, nrow=3, ncol=4, scale = "free") +
  guides(fill=FALSE) +
  ylab("log2 TMM normalized values") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

OutputName <- paste(Taxa_All[nNumb],"Negative_Control_vs_Leaf.tiff",sep="_")
ggsave(OutputName, width=8, height=8)




##################################################################
######## DAO2 : Positive Control vs leaf
##################################################################
head(design)


lrt <- glmLRT(fit, contrast=c(0,1,-1))
Result1 <- topTags(lrt, n=dim(CountData)[1], sort.by="none")$table
colnames(Result1)<-c("Taxa","logFC","logCPM","LR","PValue","FDR")

head(Result1)
nrow(CountData)[1]
print(sum(Result1$PValue < 0.05))

Result_table1 <- cbind(Taxa_all,Result1,CountData,TMM)
head(Result_table1)

OutputName <- paste(Taxa_All[nNumb],"Positive_Control_vs_Leaf.csv",sep="_") 
write.csv(Result_table1, OutputName, quote=F, row.names=F)


### Drawing box-plot for top 10 DEGs
print(length(which(Result_table1$PValue < 0.05)))

Index <- order(Result_table1$PValue)[1:12]
OrderedTable <- Result_table1[Index,]
OrderedTMM <- TMM[Index,]

Index <-which(OrderedTable$PValue < 0.05)

plot_data <- data.frame(Abundance = OrderedTMM[Index[1],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[1]])

for(i in 2:length(Index)){
  temp <- data.frame(Abundance = OrderedTMM[Index[i],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[i]])

  
  plot_data <- rbind(plot_data, temp)
}


str(plot_data)
dim(plot_data)

plot_data_subset <- subset(plot_data,subset=(Condition == "Positive control" | 
Condition == "Leaf_0.3" | Condition == "Leaf_0.5"))
plot_data_subset$Taxonomy <- gsub("_","\n_",plot_data_subset$Taxonomy)
ggplot(data=plot_data_subset, aes(x=Condition, y=Abundance)) +
  geom_boxplot(aes(fill = Condition)) +
  facet_wrap(~ Taxonomy, nrow=3, ncol=4, scale = "free") +
  guides(fill=FALSE) +
  ylab("log2 TMM normalized values") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

OutputName <- paste(Taxa_All[nNumb],"Positive_Control_vs_Leaf.tiff",sep="_")
ggsave(OutputName, width=8, height=8)


###################################################
###################################################
######### 2) Root 
###################################################
###################################################
head(CountData_Raw)
dim(CountData_Raw)


CountData <-CountData_Raw[,c(1:8,17:24)]
head(CountData)

########MetaData############


Combination<- c(rep("Negative control",4),
rep("Positive control",4),
rep("Root_0.3",4),
rep("Root_0.5",4))

Combination <- factor(Combination,
levels=c("Negative control","Positive control","Root_0.3","Root_0.5"))

Tissue <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Root",8))

Tissue <- factor(Tissue,
levels=c("Negative control","Positive control","Root"))


MetaData <- data.frame(colnames(CountData),Combination,Tissue)
str(MetaData)


targets <- MetaData ##Metadata
design <- model.matrix(~0+Tissue, data=targets)
head(CountData)



##################################################################
######## DAO : Fitting Treatment
##################################################################


y <- DGEList(counts=CountData, gene=Taxa)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
plotBCV(y)
fit <- glmFit(y, design)
TMM <- cpm(y, normalized.lib.sizes=TRUE,log=T)
TMMcolName <- colnames(TMM)
TMMcolName <- paste("TMM_",TMMcolName,sep="")
str(TMMcolName)
colnames(TMM) <- as.factor(TMMcolName)



##################################################################
######## DAO1 : Negative Control vs leaf
##################################################################
head(design)


lrt <- glmLRT(fit, contrast=c(1,0,-1))
Result1 <- topTags(lrt, n=dim(CountData)[1], sort.by="none")$table
colnames(Result1)<-c("Taxa","logFC","logCPM","LR","PValue","FDR")

head(Result1)
nrow(CountData)[1]
print(sum(Result1$PValue < 0.05))

Result_table1 <- cbind(Taxa_all,Result1,CountData,TMM)
head(Result_table1)

OutputName <- paste(Taxa_All[nNumb],"Negative_Control_vs_Root.csv",sep="_") 
write.csv(Result_table1, OutputName, quote=F, row.names=F)


### Drawing box-plot for top 10 DEGs
print(length(which(Result_table1$PValue < 0.05)))

Index <- order(Result_table1$PValue)[1:12]
OrderedTable <- Result_table1[Index,]
OrderedTMM <- TMM[Index,]

Index <-which(OrderedTable$PValue < 0.05)

plot_data <- data.frame(Abundance = OrderedTMM[Index[1],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[1]])

for(i in 2:length(Index)){
  temp <- data.frame(Abundance = OrderedTMM[Index[i],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[i]])

  
  plot_data <- rbind(plot_data, temp)
}


str(plot_data)
dim(plot_data)

plot_data_subset <- subset(plot_data,subset=(Condition == "Negative control" | 
Condition == "Root_0.3" | Condition == "Root_0.5"))
plot_data_subset$Taxonomy <- gsub("_","\n_",plot_data_subset$Taxonomy)
ggplot(data=plot_data_subset, aes(x=Condition, y=Abundance)) +
  geom_boxplot(aes(fill = Condition)) +
  facet_wrap(~ Taxonomy, nrow=3, ncol=4, scale = "free") +
  guides(fill=FALSE) +
  ylab("log2 TMM normalized values") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

OutputName <- paste(Taxa_All[nNumb],"Negative_Control_vs_Root.tiff",sep="_")
ggsave(OutputName, width=8, height=8)




##################################################################
######## DAO2 : Positive Control vs leaf
##################################################################
head(design)


lrt <- glmLRT(fit, contrast=c(0,1,-1))
Result1 <- topTags(lrt, n=dim(CountData)[1], sort.by="none")$table
colnames(Result1)<-c("Taxa","logFC","logCPM","LR","PValue","FDR")

head(Result1)
nrow(CountData)[1]
print(sum(Result1$PValue < 0.05))

Result_table1 <- cbind(Taxa_all,Result1,CountData,TMM)
head(Result_table1)

OutputName <- paste(Taxa_All[nNumb],"Positive_Control_vs_Root.csv",sep="_") 
write.csv(Result_table1, OutputName, quote=F, row.names=F)


### Drawing box-plot for top 10 DEGs
print(length(which(Result_table1$PValue < 0.05)))

Index <- order(Result_table1$PValue)[1:12]
OrderedTable <- Result_table1[Index,]
OrderedTMM <- TMM[Index,]

Index <-which(OrderedTable$PValue < 0.05)

plot_data <- data.frame(Abundance = OrderedTMM[Index[1],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[1]])

for(i in 2:length(Index)){
  temp <- data.frame(Abundance = OrderedTMM[Index[i],],
                        Condition = Combination,
                        Taxonomy = OrderedTable$Taxa[Index[i]])

  
  plot_data <- rbind(plot_data, temp)
}
    """
    }