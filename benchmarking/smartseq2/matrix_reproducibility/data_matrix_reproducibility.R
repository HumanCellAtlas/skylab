#' ---
#' title: HTML report Pipeline Reproducibility Test
#' author: Jishu Xu
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' 
#' # Load  functions and input data
#' ## Input parameters
#' 
#' * --matrix1 data matrix from pipeline, such as base pipeline.
#' * --matrix2 data matrix from pipeline, such as updated pipeline.
#' * --gtf_file gene annotation file in gtf format.
#' * --metadata_file meta information about cells, such as cell type. 
#' * --output_prefix output prefix or output name.
#' * --output_prefix output prefix or output name.
#' 
#' # Inputs
library(knitr)
source('/usr/local/scripts/analysis_functions.R')
# color palette
palette(c("#00AFBB", "#E7B800"))
option_list <- list(
  make_option(
    "--matrix1",
    type = "character",
    default = NULL,
    help = "data matrix file name",
    metavar = "character"
  ),
  make_option(
    "--matrix2",
    type = "character",
    default = NULL,
    help = "updated data matrix file name",
    metavar = "character"
  ),
  make_option(
    "--gtf_file",
    type = "character",
    default = NULL,
    help = "gtf annotation file",
    metavar = "character"
  ),
  make_option(
    "--metadata_file",
    type = "character",
    default = NULL,
    help = "metadata file",
    metavar = "character"
  ),
  make_option(
    "--output_prefix",
    type = "character",
    default = NULL,
    help = "output prefix",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
matrixfile1 <- opt$matrix1
matrixfile2 <- opt$matrix2
gtf_file <-  opt$gtf_file
meta_file <-  opt$metadata_file
grps <- opt$groups
output_name <- opt$output_prefix
#'
#' # Parse Inputs
#' ## metadata
#' metadata records information related to this dataset,
#' such as 'cell lineage','celltype','is population'
# parse metadata
meta <- read.delim(meta_file, sep = ',', header = T)
knitr::kable(head(meta))
#' ## Load Data Matrix
#' Data matrix files are *.csv file and can be count, TPM or FPKM matrix. 
#' For benchmarking purpose, TPM would be optimal choice.
# load matrix file 
mat1 <- read.csv(matrixfile1)
mat2 <- read.csv(matrixfile2)
rownames(mat1)<-mat1[,1]
rownames(mat2)<-mat2[,1]
#' # Data Transformation
#' First, take log2 transformation then filter out cells with too low overall expression level,
#' such as over all sum expression <10000
# log2 transformation
mat1.log2 <- takelog2(mat1[,-c(1:2)])
mat2.log2 <- takelog2(mat2[,-c(1:2)])
# filter cells/samples if total expression less than 10000
mat1.d <- FilterCellsbyExp(mat1.log2, 10000)
mat2.d <- FilterCellsbyExp(mat2.log2, 10000)
#' Then we need to match data matrix with metadata by cell ID(SRAID). 
# parse meta data to match filtered data matrix
mlist1 <- match(colnames(mat1.d), meta$sra)
mlist2 <- match(colnames(mat2.d), meta$sra)
grp1 <- meta[mlist1, 'cell_type']
grp2 <- meta[mlist2, 'cell_type']
#'
#' Finally, we have log2-transformed and filtered data matrix, which is matched with metadata.
#' # Run Analysis
#' ## Group t-test: SC vs Bulk
#' Here we carry out t-test between overall mean of single cells vs bulk samples
#'
#' Run ttest on base pipeline results  
test.res1 <- apply(mat1.d, 1, function(x) {
  RunTest(x = x, y = grp1, 'ttest')
})
# print out header of table
knitr::kable(head(t(test.res1)))
#' Run ttest on updated pipeline results
test.res2 <- apply(mat2.d, 1, function(x) {
  RunTest(x = x, y = grp2, 'ttest')
})
# print out header of table
knitr::kable(head(t(test.res2)))
#'
#' # Visualize t-test results
#' Now we have two sets of t-test results from two pipeline results. 
#' Next we will visualize them by using scatterplot.
#' ## Comparison results
#' We test the reproducibility of two pipelines by checking whether 
#' two pipelines can produce consisitent signatures gene. For example, 
#' t-test of gene A return a significant p-value from both pipeline, 
#' then we call this gene A is consistent signature gene, otherwise it 
#' is inconsistent signature. 
#' ## Fold Changes
#' To visualize this comparison, we first plot fold change of two pipeline into scatterplot, 
#' x-axis represents foldchange from base pipeline 
#' y-axis represents foldchange form updated pipeline. 
# extract foldchange from data matrix
res1 <- data.frame(t(test.res1))
res2 <- data.frame(t(test.res2))
fc1 <- res1$a1 - res1$a2
names(fc1) <- rownames(res1)
fc2 <- res2$a1 - res2$a2
names(fc2) <- rownames(res2)
#' ## Color Schema Of Signatures
#' We define the following color schema:
#' * Consistent signatures in both pipeline
#' * Inconsistent signatures in base but not in updated
#' * Inconsistent signatures not in base but in updated
#' * Consistent signatures not in base or updated
#' 
#' Next we convert pvalue to color schema.
# convert p value
res1$pvalue <-
  factor(ifelse(res1$pvalue > 0.01 | is.na(res1$pvalue), "Base.insig", "Base.sig"))
res2$pvalue <-
  factor(ifelse(res2$pvalue > 0.01 | is.na(res2$pvalue), "Updated.insig", "Updated.sig"))
gcols <- paste(res1$pvalue, res2$pvalue, sep = ',')
# create data frame
df <- merge(fc1, fc2, by = 'row.names')
colnames(df) <- c('ensID', 'Base', 'Updated')
df['groups'] <- gcols
#' Then use 'ggpubr' package to visualize the comparison.
sp <-
  ggscatter(
    df,
    x = 'Base',
    y = 'Updated',
    color = "groups",
    size = 3,
    alpha = 0.6,
    palette = "jco",
    ggtheme = theme_minimal(),
    legend.title = "t-test SC vs Bulk",
    xlab = "Fold Changes: Base Bulk -SC",
    ylab = "Fold Changes: Updated Bulk - SC"
  ) + border()
sp
#' ## Inconsistent Signature Genes
#' It would be more interesting to find out what type of genes are inconsistent between pipelines in term of their biotype.
#' First we group these inconsistent signature genes.
# parse out two inconsistent signatures categories
df.base <- subset(df, df$groups %in% c("Base.sig,Updated.insig"))
df.updated <- subset(df, df$groups %in% c("Base.insig,Updated.sig"))
#' Then query gene annotation based on ensID of inconsistent signatures.
# parse out gene annotation
genes <- ParseGene(gtf_file)
annt1 <- subset(genes, genes$gene_id %in% df.base$ensID)
annt2 <- subset(genes, genes$gene_id %in% df.updated$ensID)
#' Generate total counts table by gene annotation/biotype.
tb1 <- table(annt1$gene_type)
tb2 <- table(annt2$gene_type)
#' Now we have two vectors of counts, one is for 'Base.sig,Updated.insig' and another is 'Base.insig,Updated.sig'.
#' For visualization purpose, we would give negative values to one of the count vectors.
#'
# assign negetive value to one of the vectors
tb <- rbind(tb1, -1 * tb2)
rownames(tb) <-
  c('Base.sig,Updated.insig', 'Base.insig,Updated.sig')
#' Then we can visualize inconsisitent signatures by barplot. 
mtb <- melt(tb)
pb <- ggbarplot(
  mtb,
  x = 'Var2',
  y = 'value',
  color = "Var1",
  fill = "Var1",
  # change fill color by FC
  x.text.angle = 90,
  # Rotate vertically x axis texts
  xlab = "# Of Inconsistent Signature Genes",
  ylab = "",
  rotate = TRUE,
  ggtheme = theme_minimal(),
  palette = "jco",
  legend.title = "Inconsistent Signature Genes"
) + border()
pb
#' Finally, save t-test results into '.csv' files.
write.csv(
  res1,
  file = paste(
    output_name,
    "_base_ttest_res.csv",
    sep = ''),
    quote = F,
    row.names = F
)
write.csv(
  res2,
  file = paste(
    output_name,
    "_updated_ttest_res.csv",
    sep = ''),
  quote = F,
  row.names = F
)
#' ## Converting R Script to HTML/PDF  
#' If you did everyhing right, above this is the easy part.  Simply render the script as desired with the `render` 
#' function from `rmarkdown`.  
#' rmarkdown::render('/usr/local/scripts/data_matrix_reproducibility2.R')
