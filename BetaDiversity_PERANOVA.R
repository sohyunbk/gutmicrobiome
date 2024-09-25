####################
##	18/01/ 12	##
##	By Sohyun   ##
####################
setwd("D:/BIOPOP_sohyun/Project/4.Metagenome_Chicken/2.Analysis/7.BetaDiversity_R")

DistanceTable  <- read.table("unifrac_weight_distance-matrix.tsv",header=T)
head(DistanceTable)

####################################################

Condition <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Leaf_0.3",4),
rep("Leaf_0.5",4),
rep("Root_0.3",4),
rep("Root_0.5",4))


Treatment <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Leaf",8),
rep("Root",8))
Treatment <- factor(Treatment,levels=c("Negative control","Positive control","Leaf","Root"))

Ratio <- c(rep("Negative control",4),
rep("Positive control",4),
rep("0.3",4),
rep("0.5",4),
rep("0.3",4),
rep("0.5",4)
)


Allium_hookeri <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Allium_hookeri",16))

Group <- c(rep("Negative control",4),
rep("Positive control",4),
rep("Leaf_0.3",4),
rep("Leaf_0.5",4),
rep("Root_0.3",4),
rep("Root_0.5",4))

##################################################



Distance <- as.matrix(DistanceTable[,2:ncol(DistanceTable)])
head(Distance)
Header <- colnames(Distance)
dim(Distance)
rownames(Distance) <- Header

dim(Distance)
str(Distance)
as.dist(Distance)
str(dist(Distance))

colnames(Distance) <- Treatment
rownames(Distance) <- Treatment


hClustering <-  hclust(as.dist(Distance), method = 'complete')
plot(hClustering, hang = -2)

Cut_2_group_hc <- cutree(hClustering, 4)

my_color = rainbow(4)
my_color = c("#FF7F00","#B983FF")

library(ape)

tiff(paste("Clustering_WeightedUnifrac_Tissue",".tiff",sep=""), width = 3200, height = 3200, units = "px", res = 400)
par(mfrow = c(1,2))
plot(as.phylo(hClustering), type = "cladogram", tip.color = my_color[Cut_2_group_hc])
plot(as.phylo(hClustering), type = "radial", tip.color = my_color[Cut_2_group_hc])
dev.off()

################################################################


Distance <- as.matrix(DistanceTable[,2:ncol(DistanceTable)])
head(Distance)
Header <- colnames(Distance)
dim(Distance)
rownames(Distance) <- Header


library(ggplot2)

Color <- c( "#377EB8","#8DD3C7",   "#E41A1C" , "#FF7F00","#7FC97F"  ,"#984EA3","#F781BF", "#F781BF")

PCA_plot_Class <- cmdscale(Distance,eig=TRUE, k=2)
plotData <- data.frame(PCA_plot_Class$points, Combination = Treatment)

ggplot(plotData, aes(x=X1, y=X2, color=Treatment)) +
	geom_point(size=2) +
	scale_color_manual(breaks = levels(as.factor(Treatment)), values = Color) +
	stat_ellipse(aes(x=X1,y=X2, fill=Treatment), geom="polygon", level=0.9, alpha=0.01) +
	theme_bw()+
	theme(
	#axis.text.x = element_blank(),
	axis.title.x=element_text(size= 20),
	axis.title.y=element_text(size= 20),
	#legend.position='none',
	plot.title = element_text(size = 30, face = "bold"),
	axis.text.y = element_text(size= 20),
	axis.text.x = element_text(size= 20))+
   xlab("PC1") +
   ylab("PC2")
DataName <- "Combination"
ggsave(paste("PCoA_",DataName,".tiff",sep=""), width=10, height=8)



################################################################
## perm manova
##################################################################
library(RVAideMemoire)
require(vegan)
data(iris)



Condition

PairWise <- pairwise.perm.manova(as.dist(Distance),Condition,nperm=10000,p.method="none")
PairWise


PairWise2 <- pairwise.perm.manova(as.dist(Distance),Treatment,nperm=10000,p.method="none")
PairWise2


MetaData <- data.frame(Treatment,Ratio,Allium_hookeri,Group)


#adonis(as.dist(Distance) ~ Treatment * Ratio , data=MetaData, permutations=1000)
adonis(as.dist(Distance) ~ Treatment  , data=MetaData, permutations=1000)
adonis(as.dist(Distance) ~ Ratio  , data=MetaData, permutations=1000)
adonis(as.dist(Distance) ~ Group  , data=MetaData, permutations=1000)
