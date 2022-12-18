
BiocManager::install("DESeq2")
BiocManager::install("WGCNA")
library(devtools)
install_github("ctlab/fgsea")

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(WGCNA)
library(wesanderson)
library(Rmisc)
library(factoextra)
library(scales)
library(fgsea)

#Original file


matrixFile <- read.delim("your_matrix_file.txt",header = F,row.names = 1)
sampleFile <- read.csv("your_sample_file.csv")

#setting up coldata

coldata = sampleFile
cts = matrixFile
colnames(cts) = coldata$SampleID
row.names(coldata) = coldata$SampleID

cts = cts[,colnames(cts) %in% row.names(coldata)]

#Deseq2 model: design = ~  ~  Batch  + group)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~  accession + longshort + time + longshort:time)


#This is to filter lowly-expressed genes
keep <- rowSums(counts(dds)) >= ncol(cts)
dds <- dds[keep,]


#Deseq analysis
dds_out <- DESeq(dds)

dds_out_normalized <- data.frame(counts(dds_out,normalized=TRUE))
vsd <- vst(dds_out, blind=FALSE)
head(assay(vsd),3)

cts_normalized <- assay(vsd)
cts_normalized <- as.data.frame(cts_normalized)


#Outlier_filtering
PCA_All_DEGs <- prcomp(t(cts_normalized), scale = TRUE)
fviz_eig(PCA_All_DEGs)
eig.val_PCA_All_DEGs <- get_eigenvalue(PCA_All_DEGs)
res_idv_PCA_All_DEGs <- get_pca_ind(PCA_All_DEGs)
head(res_idv_PCA_All_DEGs$coord[,1:4])
PCA_All_DEGs_figure_data <- data.frame(filter(coldata, row.names(coldata) %in% colnames(cts_normalized)),res_idv_PCA_All_DEGs$coord[,1:4])

percentVar_PCA_All_DEGs <- round(eig.val_PCA_All_DEGs$variance.percent,1)

library(ggplot2)

#This is just for visualization
#PC1 and PC2
PCA_All_DEGs_figure_1 <- ggplot(PCA_All_DEGs_figure_data, aes(x=Dim.1, y=Dim.2, col = your_interest, shape = your_interest2)) +
  geom_point( size = 3) +
  xlab(paste0("PC1: ",percentVar_PCA_All_DEGs[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_PCA_All_DEGs[2],"% variance")) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 9)) + 
  theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
  labs(col = "your_col", shape = "your_shape")  +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold"))

PCA_All_DEGs_figure_1

#PC2 and PC3

PCA_All_DEGs_figure_2 <- ggplot(PCA_All_DEGs_figure_data, aes(x=Dim.2, y=Dim.3, col = your_interest, shape = your_interest2)) +
  geom_point( size = 3) +
  xlab(paste0("PC2: ",percentVar_PCA_All_DEGs[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar_PCA_All_DEGs[3],"% variance")) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 9)) + 
  theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
  labs(col = "your_col", shape = "your_shape")  +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold"))


PCA_All_DEGs_figure_2

#Further filtering
#This is all based on your PC1/2/3 Visualization

further_exclude_PCA <- row.names(filter(PCA_All_DEGs_figure_data, Dim.2 < -200))

coldata <- filter(coldata, !row.names(coldata) %in% further_exclude_PCA)
cts <- cts[,colnames(cts) %in% row.names(coldata)]

#Test whether rownames of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))


#New Deseq2 with filtered coldata and cts


dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ your design)
keep <- rowSums(counts(dds)) >= ncol(cts)
dds <- dds[keep,]
dds_out <- DESeq(dds)

dds_out_normalized <- data.frame(counts(dds_out,normalized=TRUE))

#cts_normalized_filtering

vsd <- vst(dds_out, blind=FALSE)
head(assay(vsd),3)

cts_normalized <- assay(vsd)
cts_normalized <- as.data.frame(cts_normalized)


#PCA for general characterization

PCA_All_DEGs <- prcomp(t(cts_normalized), scale = TRUE)
fviz_eig(PCA_All_DEGs)
eig.val_PCA_All_DEGs <- get_eigenvalue(PCA_All_DEGs)
res_idv_PCA_All_DEGs <- get_pca_ind(PCA_All_DEGs)
head(res_idv_PCA_All_DEGs$coord[,1:4])
PCA_All_DEGs_figure_data <- data.frame(filter(coldata,row.names(coldata) %in% colnames(cts_normalized)),res_idv_PCA_All_DEGs$coord[,1:4])

percentVar_PCA_All_DEGs <- round(eig.val_PCA_All_DEGs$variance.percent,1)


#New PCA
#This is just for visualization
#PC1 and PC2
PCA_All_DEGs_figure_1 <- ggplot(PCA_All_DEGs_figure_data, aes(x=Dim.1, y=Dim.2, col = your_interest, shape = your_interest2)) +
  geom_point( size = 3) +
  xlab(paste0("PC1: ",percentVar_PCA_All_DEGs[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_PCA_All_DEGs[2],"% variance")) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 9)) + 
  theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
  labs(col = "your_col", shape = "your_shape")  +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold"))

PCA_All_DEGs_figure_1

#PC2 and PC3

PCA_All_DEGs_figure_2 <- ggplot(PCA_All_DEGs_figure_data, aes(x=Dim.2, y=Dim.3, col = your_interest, shape = your_interest2)) +
  geom_point( size = 3) +
  xlab(paste0("PC2: ",percentVar_PCA_All_DEGs[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar_PCA_All_DEGs[3],"% variance")) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 9)) + 
  theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
  labs(col = "your_col", shape = "your_shape")  +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold"))


PCA_All_DEGs_figure_2






#WGCNA

library(WGCNA)
library(dplyr)


datExpr0 = as.data.frame(t(cts_normalized))

sampleTree = hclust(dist(datExpr0), method = "average")

library(ggplot2)
library(ggdendro)

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#WGCNA setting

options(stringsAsFactors = FALSE)

filtered_coldata <- coldata

#trait data

allTraits = filtered_coldata
filtered_samples <- row.names(datExpr)
datTraits <- filter(allTraits, row.names(allTraits) %in% filtered_samples)
row.names(datTraits) = filtered_samples

#Automatic network construction and module detection

#Soft threshholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power

cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# One-step network construction and module detection
net = blockwiseModules(datExpr, power = 12, networkType = "signed",
                       TOMType = "signed", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "supervised_norm_count_TOM",
                       verbose = 3)

#Caution: WGCNA might conflict with other packages and resulted in "Error" by the end of blockwiseModules.
#to solve this, do "cor <- WGCNA::cor" before running blockwiseModules, and do "cor<-stats::cor" after the run is done.

table(net$colors)

mergedColors = labels2colors(net$colors)
table(mergedColors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "unsupervise-02-networkConstruction-auto_n50=.RData")

#Check eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

modNames = substring(names(MEs), 3)

#Module membership
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMpval = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
geneModuleMembership$gene_name = row.names(geneModuleMembership)

MMpval$gene_name = row.names(MMpval)


names(geneModuleMembership) = modNames;
geneModuleMembership$final_module = moduleColors

names(MMpval) = modNames
MMpval$final_module = moduleColors

gmm_long = gather(geneModuleMembership,key = 'module',value = 'MM', modNames[1]:modNames[length(modNames)])
MMpval_long = gather(MMpval,key = 'module',value = 'MMpval', modNames[1]:modNames[length(modNames)])

colnames(gmm_long)[1] = 'gene_name'
colnames(MMpval_long)[1] = 'gene_name'

gmm_MMpval_long = data.frame(gmm_long,
                             MMpval = MMpval_long$MMpval)

##

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datTraits_1 <- datTraits

#Here the number in c() should be the column number of your variable of interest

moduleTraitCor = cor(MEs, datTraits_1[,c(4,5)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


#Module_trait relationships plot

tiff('Module_trait relationships.tiff', units="in", width = 4, height=8, res=600, compression = 'lzw')
par(mar = c(4,6, 4, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits_1)[c(4,5)],
               ySymbols = names(MEs),
               yLabels = names(MEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(100),
               setStdMargins = FALSE,
               cex.text = 0.5,
               textMatrix = textMatrix,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               cex.lab.x = 0.7,
               cex.lab.y = 0.5,
               xLabelsAngle = 30,
               xLabelsPosition = "bottom",
               xLabelsAdj = 0.9)

dev.off()


# Plot 3: summary of number of gene per module

library(dplyr)
library(tidyverse)

table(mergedColors)
N_gene_per_module <- as.data.frame(table(mergedColors))
names(N_gene_per_module) = c("module_name", "Gene_number")
N_gene_per_module <- N_gene_per_module[order(-N_gene_per_module$Gene_number),]
module_name <- N_gene_per_module$module_name
cols <- c("grey"="grey","turquoise"="turquoise","blue"="blue","brown"="brown","yellow"="yellow","green"="green","red"="red","black"="black","pink"="pink","magenta"="magenta","purple"="purple","greenyellow"="greenyellow","tan"="tan","salmon"="salmon","cyan"="cyan","midnightblue"="midnightblue","lightcyan"="lightcyan","grey60"="grey60","lightgreen"="lightgreen","lightyellow"="lightyellow","royalblue"="royalblue","red"="red","green"="green","turquoise"="turquoise","grey"="grey") 
N_gene_per_module$name_and_number <- paste("ME",N_gene_per_module$module_name," (",N_gene_per_module$Gene_number,")", sep ="")


library(ggplot2)
library(ggforce)


# Plot 4: Model eigengene trend

library(ggplot2)
library(Rmisc)
library(tidyr)

sample_name <- row.names(datTraits)

MEs_visualization <- data.frame(sample_name,MEs,datTraits)

#A and B should be the first and last column name with with ME, e.g. MEturquoise

MEs_visualization_long <- gather(MEs_visualization,Module,ME, A:B, factor_key=TRUE)


#You can use the df of MEs_visualization_long  to visualize based on your interest
#Here is just an example. THIS WILL NOT WORK WITH YOUR DATASET


MEs_visualization_long$field_chilling_factor <- factor(MEs_visualization_long$field_chilling)
MEs_visualization_long$Module <- factor(MEs_visualization_long$Module,levels = module_levels,labels  = N_gene_per_module$name_and_number)



surpervised_MEs_visualization_1 <- ggplot (MEs_visualization_long,aes(x=GC_collection, y=ME, col = field_chilling_factor,group = field_chilling_factor)) +
  geom_jitter(size = 1.5) +
  geom_smooth(aes(x=GC_collection, y=ME,fill =  field_chilling_factor), method = 'loess',alpha = 0.2) +
  facet_wrap( ~ Module, scales = "free_y",ncol = 5) +
  xlab("Time in forcing (h)") +
  ylab("Module eigengene") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(color = "black",size = 10),
        axis.text.x = element_text(angle = 45, size = 10)) + 
  theme(strip.text=element_text(color="Black", size = 10, face = "bold"), strip.background = element_blank()) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  theme(axis.text.y = element_blank(),axis.ticks.y =  element_blank()) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
  labs(col = "Chilling unit \n (NC") +
  scale_color_manual(values = field_chilling_cols) +
  scale_fill_manual(values = field_chilling_cols) +
  guides(fill = FALSE) +
  #theme(legend.direction = "horizontal") + 
  guides(color = guide_legend(override.aes = list(fill = field_chilling_cols)))


surpervised_MEs_visualization_1


library(ggpubr)

ggsave("MEs_visualization_figure_1.png", plot = surpervised_MEs_visualization_1, width = 28, height = 28, unit = "cm", dpi = 700)

