#Normalization,WGCNA and figures

#Packages

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(WGCNA)
library(wesanderson)
library(Rmisc)

#Input files

matrixFile <- "Hanna_ABA_gene_count.txt"
sampleFile <- "Hanna_sample_metadata.txt"

originalmatrixFile <- read.delim(matrixFile,header=FALSE,sep = "\t")
cleanedmatrixFile <- originalmatrixFile[!duplicated(originalmatrixFile$V1), ]

coldata <- read.delim(sampleFile, sep="\t", row.names=1, header=TRUE)
cts <- data.frame(cleanedmatrixFile[,-c(1,36)], row.names=cleanedmatrixFile[,1])
colnames(cts) <- rownames(coldata)

#Test whether rownames of coldata match colnames of cts
all(rownames(coldata) == colnames(cts))


#Set up factor type:

coldata$Time <-as.numeric(coldata$Time)
coldata$Treatment <- factor(coldata$Treatment, levels = c("ABA","Control"))

#Deseq2 model: design = ~ Treatment + Time)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Treatment + Time)

#Low count genefiltering. Criterium: Averagely, at least one count per sample

keep <- rowSums(counts(dds)) >= ncol(cts)
dds <- dds[keep,]

#DESeq2
dds_out <- DESeq(dds)
resultsNames(dds_out)

ddsout_normalized <- data.frame(counts(dds_out,normalized=TRUE))

summary(results(dds_out, name = "Time"))
summary(results(dds_out, contrast=c("Treatment","ABA","Control")))

res_Treatment <- data.frame(results(dds_out, contrast=c("Treatment","ABA","Control")))

#Using padj <0.05 as DEG criterium

res_Treatment_padj0.05 <- filter(res_Treatment, padj < 0.05)
res_Treatment_padj0.05_LFC1 <- filter(res_Treatment_padj0.05, log2FoldChange >1 | log2FoldChange < -1) 

#3010 genes differently expressed in the overall comparison of "Treatment_JA_vs_Control"

#Unsupervised (without any filtering) WGCNA
#Use vst (variance stabilization transformation) count as input

vsd <- vst(dds, blind=FALSE)
head(assay(vsd),3)

cts_vst <- assay(vsd)
cts_vst <- as.data.frame(cts_vst)

#WGCNA

library(WGCNA)
library(dplyr)
options(stringsAsFactors = FALSE)


datExpr0 = as.data.frame(t(cts_vst))
datExpr0 <- as.matrix(datExpr0)

datExpr0_check <- datExpr0

row.names(datExpr0_check) <- paste0(coldata$Treatment,'_',coldata$Time)

#Cluster samples (hclust) to remove outlier samples
sampleTree = hclust(dist(datExpr0_check), method = "average")

tiff('Sample dendrogram.tiff', units="in", width=18, height=10, res=1000, compression = 'lzw')

par(mar = c(2,8,5,3))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.axis = 2, cex = 1.5,font.axis = 2,cex.main = 2.5,font.lab =2, cex.lab = 3,ylab = "Height", hang = 0.01, lwd = 2.5, lend = 'round')

#Don't think there is any outliers

dev.off()
datExpr <- datExpr0

#trait data

datTraits = coldata

datTraits <- datTraits %>% 
  mutate(Treatment_num = if_else(Treatment == "ABA", 1, 0))


# Re-clustering samples 
collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits[,c(2,4)],signed = FALSE, colors = blueWhiteRed(100))

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = c("Time post-treatment","Treatment"),
                    sub="", xlab="", cex.lab = 1,main = "",
                    cex.axis = 0.7, cex.main = 1, cex.dendroLabels = 1, cex.colorLabels = 1)

save(datExpr, datTraits, file = "unsupervised WGCNA.RData")

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
                       TOMType = "signed", minModuleSize = 30,
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

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleTraitCor = cor(MEs, datTraits[,c(2,4)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


# Module_trait relationships

tiff('Module_trait relationships.tiff', units="in", width=4, height=6, res=1000, compression = 'lzw')
par(mar = c(4,6, 4, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Time post-treatment","ABA"),
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
               xLabelsAdj = 0.7)

dev.off()

#Visualization of MEs on timepoint

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


# Model eigengene visualization

library(ggplot2)
library(Rmisc)
library(tidyr)

sample_name <- row.names(datTraits)

MEs_visualization <- data.frame(sample_name,MEs,datTraits)

MEs_visualization$Treatment <-factor(MEs_visualization$Treatment)
MEs_visualization_long <- gather(MEs_visualization,Module,ME,MEbrown:MEgrey, factor_key=TRUE)

#Summarying

MEs_visualization_long_summary <- summarySE(MEs_visualization_long, measurevar="ME", groupvars=c("Module","Time","Treatment"))

MEs_visualization_long_summary_new <- MEs_visualization_long_summary

pd <- position_dodge(0.2)

module_levels = paste0("ME",N_gene_per_module$module_name)

MEs_visualization_long_summary_new$Module <- factor(MEs_visualization_long_summary_new$Module,levels = module_levels,labels  = N_gene_per_module$name_and_number)


MEs_visualization_long_summary_new$Treatment <- factor(MEs_visualization_long_summary_new$Treatment, levels = c("Control","ABA"))
MEs_visualization_long_summary_new$Time <- factor(MEs_visualization_long_summary_new$Time, levels = c("0","6","12","24","48","72"))




MEs_visualization_1 <- ggplot (MEs_visualization_long_summary_new,aes(x=Time, y=ME, group = Treatment, color = Treatment)) +
  geom_line(position=pd) +
  geom_point(position=pd, size = 1.7) +
  facet_wrap( ~ Module,scales = "free_y",ncol = 4) + 
  xlab("Time post-treatment (h)") +
  ylab("Module eigengene") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 10)) + 
  theme(strip.text=element_text(color="Black", size = 10, face = "bold"), strip.background = element_blank()) +
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  labs(color = "Treatment") +
  scale_color_manual(values = c("red","blue")) +
  theme(axis.text.y = element_blank(),axis.ticks.y =  element_blank()) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) 

MEs_visualization_1 


library(grid)
library(gridExtra)
library(lemon)

gtable_show_names(MEs_visualization_1)


MEs_visualization_1 <- reposition_legend(MEs_visualization_1 + facet_wrap(~Module, ncol=4,scales = "free_y"),
                                                     'center', panel=c('panel-3-4'))

library(ggpubr)

ggsave("MEs_visualization_figure.png", plot = MEs_visualization_1, width = 18, height = 18, unit = "cm", dpi = 1500)



#Characterization of transcriptome (all gene but exclude low count genes and 'grey' genes)

WGCNA_input_gene <- data.frame(t(datExpr))

WGCNA_input_gene <- data.frame(WGCNA_input_gene, module = moduleColors)                          
WGCNA_input_gene_filtered <- filter(WGCNA_input_gene, !module == "grey")

#16815 genes left after low count filtering and noise filtering

library(factoextra)

PCA_All_genes <- prcomp(t(WGCNA_input_gene_filtered[,1:34]), scale = TRUE)
fviz_eig(PCA_All_genes)
eig.val_PCA_All_genes <- get_eigenvalue(PCA_All_genes)
res_idv_PCA_All_genes <- get_pca_ind(PCA_All_genes)
head(res_idv_PCA_All_genes$coord[,1:4])

PCA_All_genes_figure_data <- data.frame(coldata,res_idv_PCA_All_genes$coord[,1:4])
PCA_All_genes_figure_data$Time <- factor(PCA_All_genes_figure_data$Time,
                                         levels = c("0","6","12","24","48","72"),
                                         labels = c("0 h", "6 h","12 h","24 h","48 h","72 h"))

percentVar_PCA_All_genes <- round(eig.val_PCA_All_genes$variance.percent,1)

library(ggplot2)

library(RColorBrewer)
display.brewer.all()

brewer.pal(n = 9, name = "YlOrRd")

PCA_All_genes_figure <- ggplot(PCA_All_genes_figure_data, aes(x=Dim.1, y=Dim.2, shape = Treatment)) +
  scale_fill_manual(values = c("#FFFFCC","#FED976","#FEB24C", "#FD8D3C", "#FC4E2A", "#BD0026"), aesthetics = "fill") +
  scale_shape_manual(values = c(21,24),labels = c("Control","ABA")) +
  geom_point(aes(fill = Time), color= "Black", size = 4) +
  xlab(paste0("PC1: ",percentVar_PCA_All_genes[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_PCA_All_genes[2],"% variance")) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(color = "black",size = 16, face = "bold" )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text.y  = element_text(color = "black",size = 12),
        axis.text.x  = element_text(color = "black",size = 12)) + 
  theme(plot.title= element_text(color="Black", size=18, face="bold.italic", hjust=0.5)) + 
  theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
  labs(fill = "Time post-treatment")  +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
  guides(fill = guide_legend(override.aes=list(shape=21),order = 2),shape = guide_legend(order = 1)) 


  PCA_All_genes_figure


  PCA_All_genes_figure_2 <- ggplot(PCA_All_genes_figure_data, aes(x=Dim.2, y=Dim.3, shape = Treatment)) +
    scale_fill_manual(values = c("#FFFFCC","#FED976","#FEB24C", "#FD8D3C", "#FC4E2A", "#BD0026"), aesthetics = "fill") +
    scale_shape_manual(values = c(21,24),labels = c("Control","ABA")) +
    geom_point(aes(fill = Time), color= "Black", size = 4) +
    xlab(paste0("PC2: ",percentVar_PCA_All_genes[2],"% variance")) +
    ylab(paste0("PC3: ",percentVar_PCA_All_genes[3],"% variance")) +
    theme_bw() + 
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(color = "black",size = 16, face = "bold" )) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "vertical", legend.box = "vertical") + 
    theme(axis.text.y  = element_text(color = "black",size = 12),
          axis.text.x  = element_text(color = "black",size = 12)) + 
    theme(plot.title= element_text(color="Black", size=18, face="bold.italic", hjust=0.5)) + 
    theme(legend.text = element_text(color="Black", size = 14, face = "bold")) +
    theme(legend.title = element_text(color="Black", size = 14, face = "bold")) +
    labs(fill = "Time post-treatment")  +
    theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
    theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
    guides(fill = guide_legend(override.aes=list(shape=21),order = 2),shape = guide_legend(order = 1)) 
  
  
  PCA_All_genes_figure_2

PCA_figure <- ggarrange(PCA_All_genes_figure, PCA_All_genes_figure_2, ncol = 2, 
                        font.label = list(size = 20, color ="black", face = "bold"),
                        labels = c("A","B"),common.legend = TRUE, legend = "right")

PCA_figure

ggsave("PCA_figure.png", plot = PCA_figure, width = 24, height = 10, unit = "cm", dpi = 1500)


write.table(as.data.frame(PCA_All_genes_figure_data),file = "PCA_All_genes_figure_data_deacclimation.txt",
            row.names = TRUE, col.names = TRUE,sep="\t", quote = FALSE)


##Gene select based on consistent, late and vanished

summary(results(dds_out, name = "Time"))
summary(results(dds_out, contrast=c("Treatment","ABA","Control")))

res_Treatment <- data.frame(results(dds_out, contrast=c("Treatment","ABA","Control")))
res_Treatment <- data.frame(res_Treatment,module = moduleColors)

#Using padj <0.05 as DEG criterium

res_Treatment_padj0.05 <- filter(res_Treatment, padj < 0.05)
res_Treatment_padj0.05_LFC1 <- filter(res_Treatment_padj0.05, log2FoldChange >1 | log2FoldChange < -1) 
res_Treatment_padj0.05_LFC1 <- filter(res_Treatment_padj0.05_LFC1, module %in% c("turquoise","yellow","green"))

#p-val output

list <- read.delim("list_genes_vitis_correspondencesV3_1.txt", header = TRUE)

probes = row.names(res_Treatment_padj0.05_LFC1)
probes2annot = match(probes, list$Final.v3.name..tentative.creation.of.new.genes)

sum(is.na(probes2annot))

res_Treatment_padj0.05_LFC1_V1_added = data.frame(res_Treatment_padj0.05_LFC1,
                                                    V1_genename = list$v1.name[probes2annot])

res_Treatment_padj0.05_LFC1_V1_added_pos = filter(res_Treatment_padj0.05_LFC1_V1_added, log2FoldChange > 0)
res_Treatment_padj0.05_LFC1_V1_added_neg = filter(res_Treatment_padj0.05_LFC1_V1_added, log2FoldChange < 0)


write.table(as.data.frame(res_Treatment_padj0.05_LFC1_V1_added),file = "res_Treatment_padj0.05_LFC1_V1_added_deacclimation.txt",
            row.names = TRUE, col.names = TRUE,sep="\t", quote = FALSE)

#GSEA input preparation

GSEA_input <- data.frame(Gene = res_Treatment_padj0.05_LFC1_V1_added$V1_genename, 
                         padj = res_Treatment_padj0.05_LFC1_V1_added$padj)

GSEA_input <- filter(GSEA_input, !Gene =="")

write.table(as.data.frame(GSEA_input),file = "GSEA_input_all.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)


GSEA_input_pos <- data.frame(Gene = res_Treatment_padj0.05_LFC1_V1_added_pos$V1_genename, 
                             padj = res_Treatment_padj0.05_LFC1_V1_added_pos$padj)

GSEA_input_pos  <- filter(GSEA_input_pos, !Gene =="")

write.table(as.data.frame(GSEA_input_pos),file = "GSEA_input_pos.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)

GSEA_input_neg <- data.frame(Gene = res_Treatment_padj0.05_LFC1_V1_added_neg$V1_genename, 
                             padj = res_Treatment_padj0.05_LFC1_V1_added_neg$padj)

GSEA_input_neg <- filter(GSEA_input_neg, !Gene =="")

write.table(as.data.frame(GSEA_input_neg),file = "GSEA_input_neg.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)


#Late response gene

coldata <- coldata %>% 
  mutate(Within_day = if_else(Time %in% c(6,12,24), "Yes", "No"))

coldata$Trt_time <- paste0(coldata$Treatment, coldata$Time)
coldata$Trt_within_day <- paste0(coldata$Treatment, coldata$Within_day)

dds_1 <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Trt_within_day)

#Low count genefiltering. Criterium: Averagely, at least one count per sample

keep <- rowSums(counts(dds_1)) >= ncol(cts)
dds_1 <- dds_1[keep,]

#DESeq2
dds_1_out <- DESeq(dds_1)
resultsNames(dds_1_out)

res_late <- data.frame(results(dds_1_out, contrast=c("Trt_within_day","ABANo","ControlNo")))

res_late <- data.frame(res_late,module = moduleColors)

#Using padj <0.05 as DEG criterium

res_late_padj0.05 <- filter(res_late, padj < 0.05)
res_late_padj0.05_LFC1 <- filter(res_late_padj0.05, log2FoldChange >1 | log2FoldChange < -1) 
res_late_padj0.05_LFC1 <- filter(res_late_padj0.05_LFC1, module %in% c("black","purple","midnightblue","lightcyan"))

#p-val output

list <- read.delim("list_genes_vitis_correspondencesV3_1.txt", header = TRUE)

probes = row.names(res_late_padj0.05_LFC1)
probes2annot = match(probes, list$Final.v3.name..tentative.creation.of.new.genes)

sum(is.na(probes2annot))

res_late_padj0.05_LFC1_V1_added = data.frame(res_late_padj0.05_LFC1,
                                                  V1_genename = list$v1.name[probes2annot])

res_late_padj0.05_LFC1_V1_added_pos = filter(res_late_padj0.05_LFC1_V1_added, log2FoldChange > 0)
res_late_padj0.05_LFC1_V1_added_neg = filter(res_late_padj0.05_LFC1_V1_added, log2FoldChange < 0)

#GSEA input preparation

GSEA_input_late <- data.frame(Gene = res_late_padj0.05_LFC1_V1_added$V1_genename, 
                         padj = res_late_padj0.05_LFC1_V1_added$padj)

GSEA_input_late <- filter(GSEA_input_late, !Gene =="")

write.table(as.data.frame(GSEA_input_late),file = "GSEA_input_late_all.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)


GSEA_input_late_pos <- data.frame(Gene = res_late_padj0.05_LFC1_V1_added_pos$V1_genename, 
                             padj = res_late_padj0.05_LFC1_V1_added_pos$padj)

GSEA_input_late_pos  <- filter(GSEA_input_late_pos, !Gene =="")

write.table(as.data.frame(GSEA_input_late_pos),file = "GSEA_input_late_pos.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)

GSEA_input_late_neg <- data.frame(Gene = res_late_padj0.05_LFC1_V1_added_neg$V1_genename, 
                             padj = res_late_padj0.05_LFC1_V1_added_neg$padj)

GSEA_input_late_neg <- filter(GSEA_input_late_neg, !Gene =="")

write.table(as.data.frame(GSEA_input_late_neg),file = "GSEA_input_late_neg.txt",
            row.names = FALSE, col.names = FALSE,sep="\t", quote = FALSE)


#GSEA summary

GSEA_summary <- read.delim(file = "Consistent_GSEA_input.txt", header = T)

GSEA_summary$VitisNet_Category <- factor(GSEA_summary$VitisNet_Category, 
                                         levels = c("Metabolism","Genetic information processing networks",
                                                    "Environmental information processing","Cellular processes",
                                                    "Transport","Transcription factors"),labels =c("Metabolism","Genetic information \nprocessing networks",
                                                                                                   "Environmental information \nprocessing","Cellular processes",
                                                                                                   "Transport","Transcription factors"))


GSEA_summary$expression_cat <-factor(GSEA_summary$expression_cat,
                                     levels = c("All gene","Positive genes","Negative genes"),
                                     labels = c("All \n(2264)","Upregulation \n(1748)","Downregulation \n(516)"))


GSEA_summary_figure <-ggplot(GSEA_summary,aes(x=expression_cat, y=NAME,size =NOM.p.val,col = VitisNet_Category, label = SIZE)) +
  geom_point(aes(y = reorder(NAME,-NOM.p.val))) +
  ylab("Enriched patywah in GSEA") +
  xlab("Regulation category with ABA") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 10)) + 
  theme(axis.title = element_text(color = "black",size = 18, face = "bold" )) +
  labs(size = "NOM p-value",col = "VitisNet pathway category") +
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 10),
        axis.text.y = element_text(size = 7,angle = 20)) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
  scale_size_continuous(trans = 'reverse',range = c(0.5,5)) +
  guides(color = guide_legend(override.aes=list(size=4),order = 1)) +
  geom_label(hjust = 0 , nudge_x = 0.1, size =4, show.legend = FALSE) +
  theme(legend.key.size = unit(2, 'lines')) +
  scale_color_manual(values = c("orchid","cyan3","chocolate1","cornflowerblue","olivedrab3")) +
  scale_y_discrete(labels = c("Nac" = "NAC","Bhlh" = "BHLH",
                              "Ap2 Erebp" = "AP2 EREBP",
                              "Aba Biosynthesis" = "ABA Biosynthesis",
                              "R Proteins From Plant-Pathogen Interaction" = "R Proteins from Plant-Pathogen Interaction",
                              "Porphyrin And Chlorophyll Metabolism" = "Porphyrin and Chlorophyll Metabolism",
                              "Starch And Sucrose Metabolism" = "Starch and Sucrose Metabolism",
                              "Arginine And Proline Metabolism" ="Arginine and Proline Metabolism",
                              "Alanine And Aspartate Metabolism"="Alanine and Aspartate Metabolism",
                              "Abc Transporters"="ABC Transporters",
                              "Bhlh"="bHLH",
                              "Abi3Vp1" = "ABI3VP1",
                              "Hsf"="HSF",
                              "Urea Cycle And Metabolism Of Amino Groups"="Urea Cycle and Metabolism of Amino Groups",
                              "Other ZF-C3HC4" = "Other ZF-C3HC4",
                              "Channels And Pores" = "Channels and Pores",
                              "Regulation Of Actin Cytoskeleton"="Regulation of Actin Cytoskeleton",
                              "Phenylalanine Tyrosine And Tryptophan Biosynthesis" = "Phenylalanine Tyrosine and Tryptophan Biosynthesis",
                              "Protein Processing In Endoplasmic Reticulum" ="Protein Processing in Endoplasmic Reticulum",
                              "Snare Interactions In Vesicular Transport" ="Snare Interactions in Vesicular Transport",
                              "Mrna Surveillance Pathway" ="mRNA Surveillance Pathway",
                              "Phenylalanine Tyrosine And Tryptophan Biosynthesis"="Phenylalanine Tyrosine and Tryptophan Biosynthesis",
                              "Urea Cycle And Metabolism Of Amino Groups" ="Urea Cycle and Metabolism of Amino Groups",
                              "Rna Transport" ="RNA Transport",
                              "Hsf"="HSF",
                              "Primary Active Transporter Cat A5 To A8"="Primary Active Transporter Cat A5 to A8",
                              "Other Zf-C3Hc4" = "Other Zf-C3HC4",
                              "Platz" = "PLATZ",
                              "Myb" = "MYB",
                              "Valine Leucine And Isoleucine Biosynthesis" = "Valine Leucine and Isoleucine Biosynthesis",
                              "Rna Degradation" = "RNA Degradation",
                              "Biosynthesis Of Unsaturated Fatty Acids" ="Biosynthesis of Unsaturated Fatty Acids"))
                            

GSEA_summary_figure

ggsave("GSEA.png", plot = GSEA_summary_figure, width = 22, height = 16, unit = "cm", dpi = 1500)


#Gene expression check

gene_matrix_expression <- as.data.frame(ddsout_normalized[,row.names(datExpr)])
gene_matrix_expression <- data.frame(gene_name = row.names(gene_matrix_expression), gene_matrix_expression)
gene_matrix_expression_long = gather(gene_matrix_expression,Sample,Expression,ABA_T12_R1 : ControL_T72_R3, factor_key=TRUE)

gene_matrix_expression_long = na.omit(gene_matrix_expression_long)
datTraits_1 = data.frame(t(datTraits))
datTraits_1 = data.frame(t(datTraits_1))

probes = gene_matrix_expression_long$Sample
probes2annot = match(probes, row.names(datTraits_1))

gene_matrix_expression_long = data.frame(gene_matrix_expression_long,
                                         Treatment = datTraits_1$Treatment[probes2annot],
                                         Time = datTraits_1$Time[probes2annot])


probes = gene_matrix_expression_long$gene_name
probes2annot = match(probes, list$Final.v3.name..tentative.creation.of.new.genes)
gene_matrix_expression_long = data.frame(gene_matrix_expression_long,
                                         V1_name = list$v1.name[probes2annot])

library(Rmisc)
gene_matrix_expression_long_summary <- summarySE(gene_matrix_expression_long, 
                                                 measurevar="Expression", groupvars=c("V1_name","Treatment","Time","gene_name"))


gene_expression_long_check <- filter(gene_matrix_expression_long, V1_name == "VIT_02s0033g00770")
gene_matrix_expression_long_summary_demo<- gene_matrix_expression_long_summary
gene_matrix_expression_long_summary_demo$Time <- as.numeric(gene_matrix_expression_long_summary_demo$Time)
gene_matrix_expression_long_summary_demo <- filter(gene_matrix_expression_long_summary_demo, !Time >0)
gene_matrix_expression_long_summary_demo$Time <- factor(gene_matrix_expression_long_summary_demo$Time)

gene_matrix_expression_long_summary_demo <- gene_matrix_expression_long_summary_demo %>% 
  mutate(Treatment = if_else(Treatment == 'Control', 'ABA', 'ABA'))

gene_matrix_expression_long_summary_new <- rbind(gene_matrix_expression_long_summary,gene_matrix_expression_long_summary_demo)
gene_matrix_expression_long_summary_new$Time <- as.numeric(gene_matrix_expression_long_summary_new$Time)
gene_matrix_expression_long_summary_new$Time <- factor(gene_matrix_expression_long_summary_new$Time)


##RACK1 gene expression check figures

select_gene = c("VIT_17s0000g02750")

select = filter(gene_matrix_expression_long_summary_new, V1_name == select_gene)
select$Treatment <- factor(select$Treatment,
                             levels = c("Control","ABA"), labels = c("Deacclimation\n Control","Deacclimation\n ABA"))
select$Time <- factor(select$Time)
                           
                             
gene_expression_fig <- ggplot (select,aes(x=Time, y=Expression, color = Treatment, group =Treatment)) +
  geom_line(position=pd, size = 2.2) +
  geom_point(position=pd, size = 6.4) +
  xlab("Time post-treatment (h)") +
  ylab("Normalized expression") +
  theme_bw() + 
  labs(subtitle = paste0(select_gene,'\n A. thaliana RACK1 homolog')) +
  theme(plot.title = element_text(size = 23, face = "bold", hjust = -0.1),
        plot.subtitle = element_text(hjust = 0.5,size = 18, face = "bold")) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 14,face = "bold")) + 
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  labs(color = "Treatment") +
  scale_color_manual(values = c("#D95F02","#1B9E77")) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(colour = "black", size = 2))

gene_expression_fig

ggsave("Grape_RACK1_deacc_new_new.png", plot =gene_expression_fig , width = 17.5, height = 12, unit = "cm", dpi = 1500)


##theme(axis.title = element_blank()) +



##Single gene expression

#Shared_ABA (size 10x8 cm, color: Control:"#D95F02" ABA = "#1B9E77")

#VIT_02s0087g00930
#VIT_06s0004g05460
#VIT_07s0031g01770
#VIT_11s0016g03180
#VIT_13s0019g02200
#VIT_18s0001g10450
#VIT_19s0093g00550

#Auxin (size 10x8 cm)
#VIT_02s0025g04560
#VIT_02s0033g00660
#VIT_02s0033g00670
#VIT_02s0033g00770
#VIT_02s0033g00870
#VIT_03s0091g00310
#VIT_03s0180g00280
#VIT_05s0020g03280
#VIT_05s0094g01020
#VIT_06s0004g07230
#VIT_08s0040g01520
#VIT_18s0001g02570
#VIT_19s0014g04690


select_gene = c("VIT_03s0091g00310")

select = filter(gene_matrix_expression_long_summary_new, V1_name == select_gene)
select$Treatment <- factor(select$Treatment,
                           levels = c("Control","ABA"))
select$Time <- factor(select$Time)


gene_expression_fig <- ggplot (select,aes(x=Time, y=Expression, color = Treatment, group =Treatment)) +
  geom_line(position=pd, size = 2.2) +
  geom_point(position=pd, size = 6.4) +
  xlab("Time post-treatment (h)") +
  ylab("Normalized expression") +
  theme_bw() + 
  labs(title = select_gene) +
  theme(plot.title = element_text(size = 23, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "vertical", legend.box = "vertical") + 
  theme(axis.text = element_text(color = "black",size = 14,face = "bold")) + 
  theme(axis.title = element_text(color = "black",size = 14, face = "bold" )) +
  labs(color = "Treatment") +
  scale_color_manual(values = c("#D95F02","#1B9E77")) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2),
        axis.ticks = element_line(colour = "black", size = 2)) +
  theme(legend.position= "none",axis.title = element_blank()) 
  
gene_expression_fig

ggsave("VIT_03s0091g00310_deacc_new_new.png", plot =gene_expression_fig , width = 12.5, height = 10, unit = "cm", dpi = 1500)

#GEO submission file preparation

cts_GEO_colnames = paste0("Deacclimation_",row.names(datExpr0_check),"_",coldata$Rep)

Deacclimation_cts_GEO = cts
colnames(Deacclimation_cts_GEO) = cts_GEO_colnames

write.table(as.data.frame(Deacclimation_cts_GEO),file = "Deacclimation_cts_GEO.txt",
            row.names = T, col.names = T,sep="\t", quote = FALSE)

write.table(as.data.frame(cts_GEO_colnames),file = "Deacclimation_GEO_Sample_names.txt",
            row.names = F, col.names = F,sep="\t", quote = FALSE)
