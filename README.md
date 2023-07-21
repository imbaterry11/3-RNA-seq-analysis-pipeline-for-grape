# 3-RNA-seq-analysis-pipeline-for-grape
This repo includes the whole pipeline for 3'RNA-seq analysis of grape tissue.

## Files description
__Part1_alignment.sh__ includes all the code for 3'RNA-seq data trimming, alignment, QC and result output. These codes are for the operation under Linux. <br>
<br>
__3RNAseq_data_analysis_tutorial.md__ includes all the steps from outlier detection to pathway enrichment analysis. <br>
<br>
__gene_count.txt__ is a gene count files from Linux-based RNA-seq alignment. It is the output from __Part1_alignment.sh__ and one of the input for __3RNAseq_data_analysis_tutorial.md__. <br>
<br>
__sample_metadata.txt__ is the sample description file needed for differential expression analysis. It is another input for __3RNAseq_data_analysis_tutorial.md__. <br>
<br>
__All the other files__ are for the analysis in __3RNAseq_data_analysis_tutorial.md__. <br>


## Example RNA-seq data
The data used in this tutorial is from [this paper](https://www.maxapress.com/article/doi/10.48130/FruRes-2022-0001). The analysis pipeline is also a modification of the pipeline described in the paper. Please cite the paper if you develop your pipeline based on this tutorial <br>
The raw RNA-seq data used in this tutotial is available [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184114)
