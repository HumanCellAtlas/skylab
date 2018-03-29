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
print(opt)
matrixfile1 <- opt$matrix1
matrixfile2 <- opt$matrix2
metadata_file <- opt$metadata_file
output_name <- opt$output_prefix
# color palette
palette(c("#00AFBB", "#E7B800"))
# load data matrix
mat1 <- read.csv(matrixfile1)
mat2 <- read.csv(matrixfile2)
rownames(mat1) <- mat1[, 1]
rownames(mat2) <- mat2[, 1]
mat1.d <- mat1[, -c(1:2)]
mat2.d <- mat2[, -c(1:2)]
#match row
mlist <- match(rownames(mat1.d), rownames(mat2.d))
# match sample column
nlist <- match(colnames(mat1.d), colnames(mat2.d))
mat2.d <- mat2.d[mlist, nlist]
# meta
meta <- read.csv(metadata_file, header = T)
mlist <- match(colnames(mat1.d), meta$sra)
meta <- meta[mlist, ]

#log2 transformation
mat1.log2 <- takelog2(mat1.d)
mat2.log2 <- takelog2(mat2.d)

# correlation
# return two matrix of corraltion matrix
# take diag of correlation matrix
# plot into histogram
mcor <- cor(mat1.log2, mat2.log2)
cxy <- diag(mcor)
phist <- gghistogram(
  cxy,
  fill = "lightgray",
  add = "mean",
  rug = TRUE,
  ggtheme = theme_minimal(),
  xlab = "Correlation: Expression of Base vs Updated pipeline"
) + border()
phist <- addTheme(phist)
# SNN cluster
mat1.snn <- RunSNNCluster(mat1.log2)
mat2.snn <- RunSNNCluster(mat2.log2)
# Calculate adjusted rand index between clustering results to metadata
labels <- paste(meta$lineage, meta$cell, sep = '-')
r1 <- adjustedRandIndex(mat1.snn$membership,
                        labels)
r2 <- adjustedRandIndex(mat2.snn$membership,
                        labels)
t1 <- table(mat1.snn$membership, labels)
t2 <- table(mat2.snn$membership, labels)

dtable.1 <-
  ggtexttable(t1,
              theme = ttheme('mBlue', base_size = 18),
              rows = rownames(t1))
dtable.2 <-
  ggtexttable(t2,
              theme = ttheme('mBlue', base_size = 18),
              rows = rownames(t2))

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
p1 <- ggarrange(
  snn.p1,
  dtable.1,
  ncol = 1,
  nrow = 2,
  heights = c(2, 1),
  labels = c(
    paste('A: Clustering of base pipeline:Rand Index:', r1),
    "B: clustering index vs sample labels"
  )
)
p2 <- ggarrange(
  snn.p2,
  dtable.2,
  ncol = 1,
  nrow = 2,
  heights = c(2, 1),
  labels = c(
    paste('A: Clustering of updated pipeline:Rand Index:', r2),
    "B: clustering index vs sample labels"
  )
)
pdf(paste(output_name, '_cluster_data_matrix_comparison.pdf', sep = ''),
    20,
    15)
print(phist)
print(p1)
print(p2)
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