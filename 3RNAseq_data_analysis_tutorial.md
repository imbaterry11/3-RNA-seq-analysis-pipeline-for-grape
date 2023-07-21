RNA-seq data analysis for grapevine tutorial
================
Hongrui Wang
2023-07-21

- <a href="#rna-seq-data-analysis-for-grapevine"
  id="toc-rna-seq-data-analysis-for-grapevine">RNA-seq data analysis for
  grapevine</a>
  - <a href="#input-files" id="toc-input-files">Input files</a>
  - <a href="#generate-normalized-gene-count"
    id="toc-generate-normalized-gene-count">Generate normalized gene
    count</a>
  - <a href="#pca-to-detect-outliers" id="toc-pca-to-detect-outliers">PCA to
    detect outliers</a>
  - <a href="#wgcna-for-co-expression-network-analysis"
    id="toc-wgcna-for-co-expression-network-analysis">WGCNA for
    co-expression network analysis</a>
  - <a href="#gsea-for-pathway-enrichment-analysis"
    id="toc-gsea-for-pathway-enrichment-analysis">GSEA for pathway
    enrichment analysis</a>
  - <a href="#next-steps" id="toc-next-steps">Next steps</a>

# RNA-seq data analysis for grapevine

This is an RNA-seq analysis tutorial aiming to provide the entire
analysis pipeline and the biological reasoning of all the steps
included. The pipeline includes PCA-based outlier filtering, PCA-based
general transcriptome characterization, DeSeq2-based normalization,
WGCNA-based co-expression network analysis, DeSeq2-based differential
expression analysis and GSEA-based pathway enrichment analysis.

The tutorial is based on the method in Wang, H. et al. (2022)
‘Transcriptomic analysis of grapevine in response to ABA application
reveals its diverse regulations during cold acclimation and
deacclimation’, Fruit Research, 2(1), pp. 1–12. Available at:
[here](https://doi.org/10.48130/FruRes-2022-0001)

The demo data used in this tutorial is available at
[here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184114)

## Input files

sample_metadata.txt is a metadata file providing the treatment of each
sample (library) in the experiment. Here, we have two treatments:
control and ABA (a plant hormone that tends to paly a role during winter
dormancy). In this experiment, we are testing if the application of ABA
on woolly grapevine buds would delay budbreak and the transcriptome
mechanism associated with the delay. In the metadata file, ‘Treatment’
indicates the treatment during the experiment, and ‘Time’ indicates
number of hours that the samples has been in treatment

``` r
#cold data is the metadata file
coldata <- fread("sample_metadata.txt")
```

gene_count.txt is a gene count matrix with row.names as gene names and
each column as a library.

``` r
#cts data is the gene count matrix file
cts_original <- fread("gene_count.txt")
#remove duplications
cts_original = cts_original[!duplicated(cts_original[,1]),]
cts <- cts_original[,-c(1,36)]
row.names(cts) = cts_original$V1
colnames(cts) <- coldata$Sample_name
```

## Generate normalized gene count

## PCA to detect outliers

``` r
#Here we use vsd count for PCA
PCA_All_genes <- prcomp(t(cts_normalized), scale = TRUE)
fviz_eig(PCA_All_genes)
```

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/PCA%20to%20detect%20outlier-1.png)<!-- -->

``` r
eig.val_PCA_All_genes <- get_eigenvalue(PCA_All_genes)
res_idv_PCA_All_genes <- get_pca_ind(PCA_All_genes)
```

PCA1, 2 and 3 (to detect outliers)

``` r
PCA_All_genes_figure_data <- data.frame(coldata,res_idv_PCA_All_genes$coord[,1:4])
PCA_All_genes_figure_data$Time <- factor(PCA_All_genes_figure_data$Time,
                                         levels = c("0","6","12","24","48","72"),
                                         labels = c("0 h", "6 h","12 h","24 h","48 h","72 h"))

percentVar_PCA_All_genes <- round(eig.val_PCA_All_genes$variance.percent,1)

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
```

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/PCA%20to%20detect%20outlier%20plotting-1.png)<!-- -->

``` r
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
```

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/PCA%20to%20detect%20outlier%20plotting-2.png)<!-- -->
None of these samples seem to be an outlier, so we move to the next step
to do WGCNA for co-expression network analysis.

## WGCNA for co-expression network analysis

Now we have constructed the gene co-expression network and processed
module eigengenes (MEs). We have to do some analysis and visualize the
result to make decisions for the next step.

``` r
#Correlation analysis of MEs and your variables of interest. In the design of this experiment, we should care about if there are genes responding to 1) Treatment and 2) Time under treatment (columns 1 and 2)
#Here the number in c() should be the column number of your variable of interest

#WGCNA can only handle correlation analysis using numeric variables. Thus, we have to transform our treatment to boolean type values (0 or 1) to facilitate the analysis.

datTraits_analysis = datTraits
datTraits_analysis$Treatment = ifelse(datTraits_analysis$Treatment == 'ABA', 1, 0)
#So ABA is 1 and control is 0. A positive correlation of ME with Treatment means the genes in the module tend to be up-regulated by ABA compared to control.

moduleTraitCor = cor(MEs, datTraits_analysis[,c(1,2)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


#Module_trait relationships plot

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits_analysis)[c(1,2)],
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
```

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/WGCNA%20visualization%20correlation%20analysis-1.png)<!-- -->

``` r
##There might be some issues with showing the plot in Rmd
```

Based on the correlation analysis, we can easily tell that MEbrown and
MEyellow are the modules that show most significant positive correlation
with treatment. This result indicates that the genes in these two
modules are generally upregulated by ABA. In contrast, MEred and MEgreen
are the top two module that show negative correlation with treatment,
indicating that the genes in these modules are downregulated by ABA.

Since our firsts goal is the identify ABA’s impact on transcriptome, we
are ignoring the ‘Time’ variable. During the analysis of a project in
real world, we have to consider to effect of each variable.

``` r
#Summary of number of gene per module
N_gene_per_module <- as.data.frame(table(mergedColors))
names(N_gene_per_module) = c("module_name", "Gene_number")
N_gene_per_module <- N_gene_per_module[order(-N_gene_per_module$Gene_number),]
module_name <- N_gene_per_module$module_name
N_gene_per_module$name_and_number <- paste("ME",N_gene_per_module$module_name," (",N_gene_per_module$Gene_number,")", sep ="")

#Visualization of Model eigengenes (MEs)
sample_name <- row.names(datTraits)
MEs_visualization <- data.frame(sample_name,MEs,datTraits)
MEs_visualization$Treatment <-factor(MEs_visualization$Treatment)
MEs_visualization_long <- gather(MEs_visualization,Module,ME,colnames(MEs)[1]:MEgrey, factor_key=TRUE)

#Summarizing
module_levels = paste0("ME",N_gene_per_module$module_name)
MEs_visualization_long$Module = factor(MEs_visualization_long$Module,levels = module_levels,labels  = N_gene_per_module$name_and_number)

MEs_visualization_long$Treatment <- factor(MEs_visualization_long$Treatment, levels = c("Control","ABA"))
MEs_visualization_long$Time <- factor(MEs_visualization_long$Time, levels = c("0","6","12","24","48","72"))

pd <- position_dodge(0.2)

MEs_visualization_1 <- ggplot (MEs_visualization_long,aes(x=Time, y=ME, group = Treatment, color = Treatment)) +
  geom_jitter(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(x=Time, y=ME,fill =  Treatment), method = 'loess',alpha = 0.2) +
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
  scale_fill_manual(values = c("red","blue")) +
  theme(axis.text.y = element_blank(),axis.ticks.y =  element_blank()) +
  theme(legend.text = element_text(color="Black", size = 12, face = "bold")) +
  theme(legend.title = element_text(color="Black", size = 12, face = "bold")) 

MEs_visualization_1 
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/WGCNA%20visualization%20MEs%20ploting-1.png)<!-- -->

``` r
#Note that module grey contains all the noise genes that can'bt be categorized to any co-expression modules.
```

The goal of the above steps (correlation analysis and the visualization
of MEs) is to help us form the hypothesis for differential expression
analysis in a complex factorial design experiment (especially when it’s
a time serial experiment).

Here, we can see that the MEs of module brown and module yellow are
always higher in ABA treatment, whereas the MEs of module green and red
are always lower in ABA treatment. Interestingly, the fitted lines using
‘loess’ function of MEbrown and MEyellow seem exactly opposite to the
lines of MEred and MEgreen. This result might indicate that the genes in
these modules are functional in similar pathways but are negative or
positive regulator. To achieve a phenotypical change, it sometimes
requires the coordination of negative regulators and positive regulators
(e.g. upregulation of positive regulators along with downregulation of
negative regulator). Based on this observation, one hypothesis we can
form here is that, modules brown, yellow, green and red contains the
genes whose expressions are consistently altered by ABA treatment during
the experiment. To test tge hypothesis, we should first make sure if all
the genes in these four modules are consistently responsive to ABA
during the experiment. To do that, we should do a constrasting in
DeSeq2.

## GSEA for pathway enrichment analysis

``` r
#To test of hypothesis, the DeSeq2 model is simply just design = ~ Treatment
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Treatment)
```

    ##   it appears that the last variable in the design formula, 'Treatment',
    ##   has a factor level, 'Control', which is not the reference level. we recommend
    ##   to use factor(...,levels=...) or relevel() to set this as the reference level
    ##   before proceeding. for more information, please see the 'Note on factor levels'
    ##   in vignette('DESeq2').

``` r
#Low count genefiltering. Criterium: Averagely, at least one count per sample
keep <- rowSums(counts(dds)) >= ncol(cts)
dds <- dds[keep,]

#DeSeq2
dds_out <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 1 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
resultsNames(dds_out)
```

    ## [1] "Intercept"                "Treatment_Control_vs_ABA"

``` r
ddsout_normalized <- data.frame(counts(dds_out,normalized=TRUE))
summary(results(dds_out, contrast=c("Treatment","ABA","Control")))
```

    ## 
    ## out of 18906 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 5052, 27%
    ## LFC < 0 (down)     : 4749, 25%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#res_Treatment is a df containing all the statistics of the contrast
res_Treatment <- data.frame(results(dds_out, contrast=c("Treatment","ABA","Control")))
res_Treatment$V3_gene_name = row.names(cts)[keep]

#Assign the module of each gene
res_Treatment$module = moduleColors

#Filter the whole gene list based on padj LFC and module
res_Treatment_target = filter(res_Treatment, padj < 0.05, abs(log2FoldChange) > 1, module %in% c('brown','yellow','red','green'))
```

After filtering, 2381 genes were left. These genes are our target genes
to identify if any pathways were significantly impacted by ABA
treatment. To do that, we do some pathway enrichment analysis. Here, we
use GSEA (gene set enrichment analysis).

In each biological pathway, there are positive regulators and negative
regulators. To test if the pathway is enriched or if a pathway is
impacted, we have combined the impact of both types of regulators. This
means that, the upregulated genes and the downregulated genes should be
analyzed separately and together to reveal all the impacted pathways.

``` r
#We use GSEA for the pathway enrichment analysis

#Since the current gene name do not have proper functional annotation, we transform it to another gene name system that is properly annotated. 

list <- read.delim("list_genes_vitis_correspondencesV3_1.txt", header = TRUE)
Gene_function <- read.delim(file =  'Vitis_Net_gene_function.txt', header = T)

res_Treatment_target$V1_gene_name = list$v1.name[match(res_Treatment_target$V3_gene_name, list$Final.v3.name..tentative.creation.of.new.genes)]
res_Treatment_target$encoding_protein = Gene_function$Functional.annotation[match(res_Treatment_target$V1_gene_name, Gene_function$Unique.ID)]
res_Treatment_target$pathway = Gene_function$Network[match(res_Treatment_target$V1_gene_name, Gene_function$Unique.ID)]

#Some filtering, e.g., if there is no pathway assigned to the gene, we are not using the gene to do the pathway enrichment analysis

res_Treatment_target = filter(res_Treatment_target, !is.na(pathway), !pathway == '')

#Dataframe preparation for GSEA (We rank only based on absolute LFC here, but you can also rank based on pval or padj)

#GSEA using all target DEGs
res_Treatment_target_gsea = data.frame(gene = res_Treatment_target$V1_gene_name, LFC = abs (res_Treatment_target$log2FoldChange))
res_Treatment_target_gsea = res_Treatment_target_gsea[order(res_Treatment_target_gsea$LFC, decreasing = T),]
res_Treatment_target_gsea$rank = 1:nrow(res_Treatment_target_gsea)
res_Treatment_target_gsea = res_Treatment_target_gsea[,-2]
  
ranks = deframe(res_Treatment_target_gsea)

vitis_pathway <- gmtPathways("GSEA_geneset_from_VitisNet.gmt")

fgsearesult_all = fgsea(pathways = vitis_pathway,
                              stats = ranks,
                              minSize = 2,
                              maxSize = 500)

#Which pathway were enriched?
head(fgsearesult_all[order(padj), ])
```

    ##                         pathway         pval       padj   log2err         ES
    ## 1:              vv23010Ribosome 0.0002363368 0.03521418 0.5188481  0.6070175
    ## 2:        vv10195Photosynthesis 0.0473790323 0.78438620 0.2065879  0.4743766
    ## 3:   vv10620Pyruvate_metabolism 0.0128077756 0.78438620 0.3807304  0.6260227
    ## 4: vv30010Gibberellin_signaling 0.0469111325 0.78438620 0.3217759 -0.6234888
    ## 5:             vv40006Cell_wall 0.0410410410 0.78438620 0.2220560  0.4302652
    ## 6:     vv50111Tethering_factors 0.0387931034 0.78438620 0.2765006  0.8045608
    ##          NES size
    ## 1:  1.826690   23
    ## 2:  1.396700   20
    ## 3:  1.614642   10
    ## 4: -1.527009    5
    ## 5:  1.358568   33
    ## 6:  1.464980    3
    ##                                                                                                        leadingEdge
    ## 1: VIT_16s0115g00310,VIT_05s0020g01140,VIT_06s0061g00770,VIT_06s0004g06140,VIT_02s0025g02890,VIT_18s0001g12850,...
    ## 2: VIT_07s0005g04400,VIT_04s0023g00410,VIT_19s0015g01760,VIT_04s0008g01730,VIT_04s0008g04310,VIT_05s0020g03490,...
    ## 3:                       VIT_06s0009g01930,VIT_04s0008g03560,VIT_00s1995g00010,VIT_14s0060g00420,VIT_05s0049g01130
    ## 4:                       VIT_18s0001g14270,VIT_01s0026g00620,VIT_17s0000g06210,VIT_07s0104g00930,VIT_18s0072g01110
    ## 5: VIT_18s0001g10040,VIT_19s0015g00730,VIT_00s1349g00010,VIT_15s0046g01630,VIT_12s0035g01900,VIT_04s0023g02980,...
    ## 6:                                                           VIT_08s0058g00880,VIT_06s0004g04010,VIT_08s0040g02730

``` r
#It seems that Ribosome pathway was significantly enriched
plotEnrichment(vitis_pathway[["vv23010Ribosome"]],
               ranks) + labs(title="Ribosomal proteins")
```

![](3RNAseq_data_analysis_tutorial_files/figure-gfm/decision%20and%20GSEA-1.png)<!-- -->

``` r
#The next step is to do gsea with only upregulated genes and only downregulated genes.
```

## Next steps

The next step of the pipeline will be determined by the result from
pathway enrichment analysis. Typically, the top pathways will be
carefully scouted along with gene expression visualization to determine
what is happening in the pathway. This step needs some help from
external tools (e.g. Cytoscape).
