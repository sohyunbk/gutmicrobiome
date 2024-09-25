 #!/home/sb14489/miniconda3/envs/Metagenomics_qiime/bin/R
    library(edgeR)
    library(ggplot2)
    library(limma)
    library(gridExtra)

    option_list = list(
        make_option(c("--SampleName"), type="character",
                    help="Type SampleName", metavar="character"),
        make_option(c("--inputPath"), type="character",
                    help="Type InputPath"),
        make_option(c("--outputPath"), type="character",
                    help="Type OutputPath"),
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
    Taxa <- Data[,1]
    Taxonomy <- Data[,1]
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