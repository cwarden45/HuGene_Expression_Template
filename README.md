Template for normalization, QC, differential expression, and preparation of files for downstream analysis

### Order to Run Scripts ###

0)

1) 

### Dependencies (some optional) ###

*Normalization*

Affy Power Tools (apt-probeset-summarize): http://www.affymetrix.com/estore/partners_programs/programs/developer/tools/powertools.affx

*Array Annotation Files*:

HuGene Arrary Files: http://www.affymetrix.com/catalog/131453/AFFY/Human+Gene+ST+Arrays

HuGene-2.0 (.pgf, .clf, .bgp, Transcript Cluster .mps files): [Library Files](http://www.affymetrix.com/Auth/analysis/downloads/lf/wt/HuGene-2_0-st/HuGene-2_0-st_rev1.zip)

HuGene-2.0 (gene annotations): [probeset](http://www.affymetrix.com/Auth/analysis/downloads/na36/wtgene/HuGene-2_0-st-v1.na36.hg19.probeset.csv.zip) and [transcript cluster](http://www.affymetrix.com/Auth/analysis/downloads/na36/wtgene/HuGene-2_0-st-v1.na36.hg19.transcript.csv.zip) annotations are under `Current NetAffx Annotation Files`

*Differential Expression*

limma: https://bioconductor.org/packages/release/bioc/html/limma.html

qvalue: https://bioconductor.org/packages/release/bioc/html/qvalue.html

*Visualization*

gplots: https://cran.r-project.org/web/packages/gplots/index.html

RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

*Gene Set Enrichment*

IPA: http://www.ingenuity.com/

PANTHER: http://www.pantherdb.org/

DAVID: https://david.ncifcrf.gov/summary.jsp


### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value.|
|CEL_Input|List of CEL files for normalization.  Can be created using `create_cel_list.py`|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|summary_method|Probeset summarization level.  Can be *transcript_cluster* or *gene_symbol|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|sample_description_file|Name of Sample Description File|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|RMA_expression_cutoff|Rounding Value for RMA expression (minimum reliable expression level)|
|minimum_fraction_expressed|Minimum fraction of samples with expression above *rpkm_expression_cutoff*. Filter for differential expression anaylsis.|
|fold_change_cutoff|Minimum fold-change difference to consider a gene differentially expressed|
|cor_cutoff|If using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally expressed|
|fdr_cutoff|Maximum FDR to consider a gene differentially expressed|
|sec_fold_change_cutoff|If comparing two gene lists, fold-change threshold for list you want to filter out|
|sec_cor_cutoff|If comparing two gene lists and using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|sec_pvalue_cutoff|If comparing two gene lists, p-value threshold for list you want to filter out|
|sec_fdr_cutoff|If comparing two gene lists, FDR threshold for list you want to filter out|
|interaction| Method for comparing an interaction of two variables.  Can be *model*, *filter-overlap*, or *no*|
|secondary_trt| If comparing two gene lists, this is treatment group for the list that you want to filter out; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value (also converts second variable from factor to numeric, even if interaction is set to *no*)|
