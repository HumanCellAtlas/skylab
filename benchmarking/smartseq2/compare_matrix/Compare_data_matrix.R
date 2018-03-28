# This script is to compare data matrix
# from two pipelines or pre and after pipeline updated
# In this task, we will calculate and test
# the correlation of two matrix, cell vs cell
# We will also do clustering analysis then compare two clustering results
source('analysis_functions.R')
set.seed(42)

# inputs
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
    "--metadata_file",
    type = "character",
    default = NULL,
    help = "sample metadata file",
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
metadata_file <- opt$metadata_file
output_name <- opt$output_prefix
# color palette
palette(c("#00AFBB", "#E7B800"))
# load data matrix
mat1 <- read.csv(matrixfile1)
mat2 <- read.csv(matrixfile2)
rownames(mat1)<-mat1[,1]
rownames(mat2)<-mat2[,1]
mat1.d<-mat1[,-c(1:2)]
mat2.d<-mat2[,-c(1:2)]
#match row
mlist <- match(rownames(mat1.d), rownames(mat2.d))
# match sample column
nlist <- match(colnames(mat1.d), colnames(mat2.d))
mat2.d <- mat2.d[mlist, nlist]
# meta
meta <- read.csv(metadata_file,header=T)
mlist <- match(colnames(mat1.d), meta$sra)
grps <- meta[mlist, 'cell_type']

#log2 transformation
mat1.log2<-takelog2(mat1.d)
mat2.log2<-takelog2(mat2.d)

# correlation
# return two matrix of corraltion matrix
# take diag of correlation matrix
# plot into histogram
mcor<-cor(mat1.log2,mat2.log2)
cxy<-diag(mcor)
phist <- gghistogram(
  cxy,
  fill = "lightgray",
  add = "mean",
  rug = TRUE,
  ggtheme = theme_minimal(),
  xlab = "Correlation: Expression of Base vs Updated pipeline"
)
phist <- addTheme(phist)
# SNN cluster
mat1.snn <- RunSNNCluster(mat1)
mat2.snn <- RunSNNCluster(mat2)
# Calculate adjusted rand index between two cluster results
# to figure out how similar the results are
# ideal, if we habe pre-labeled data, we can calculate RAND index
# between clustering result and labels
rand <- adjustedRandIndex(mat1.snn$membership, mat2.snn$membership)

# visualize SNN clustering
snn.p1 <- ggscatter(
  mat1.snn,
  x = "PC1",
  y = "PC2",
  color = "membership",
  palette = "jco",
  size = 8,
  alpha = 0.6,
  ggtheme = theme_minimal(),
  title = 'Base Pipeline'
) + border()
snn.p2 <- ggscatter(
  mat2.snn,
  x = "PC1",
  y = "PC2",
  color = "membership",
  palette = "jco",
  size = 8,
  alpha = 0.6,
  ggtheme = theme_minimal(),
  title = 'Updated Pipeline'
) + border()
snn.p1 <- addTheme(snn.p1)
snn.p2 <- addTheme(snn.p2)
pcl <- ggarrange(snn.p1,
                snn.p2,
                ncol = 2,
                nrow = 1)

pdf(paste(output_name, '_cluster_data_matrix_comparison.pdf', sep = ''),
    30,
    60)
print(phist)
print(pcl)
dev.off()


grps <-
  merge(mat1.snn,
        mat2.snn,
        by = 'row.names',
        suffixes = c(".base", ".updated"))
write.csv(
  grps,
  file = paste(output_name, '_SNN_cluster_data_matrix_comparison.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)