# R script to do data matrix comparison
# loading library
library(ggplot2)
library(MASS)
library(plyr)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(ggpmisc)
library(cowplot)
library(corrplot)
library(ggrepel)
library(optparse)
library(rsvd)
library(scran)
library(igraph)
library(mclust)
library('rtracklayer')
set.seed(42)

addTheme <- function(p) {
  p <- p + theme(legend.position = "top")
  p <-
    p + theme(plot.title = element_text(
      hjust = 0.5,
      size = 20,
      face = 'bold'
    ))
  p <-
    p + theme(
      axis.title.x = element_text(
        color = 'black',
        size = 18,
        face = 'bold'
      ),
      axis.title.y = element_text(
        size = 18,
        color = 'black',
        face = 'bold'
      )
    )
  p <-
    p + theme(
      axis.text.y = element_text(
        color = 'black',
        size = 18,
        face = 'bold'
      ),
      axis.text.x = element_text(
        color = 'black',
        size = 18,
        face = 'bold'
      )
    )
  p <-
    p + theme(legend.text = element_text(size = 15, face = 'bold'),
              legend.title = element_blank())
  return(p)
}
# parse GTF file
ParseGene <- function(gtf_file) {
  gtf_gencode <-
    readGFF(
      gtf_file,
      version = 2L,
      tags = c("gene_name", "gene_id", "transcript_id", "gene_type")
    )
  genes <- subset(gtf_gencode, gtf_gencode$type == "gene")
  return(genes)
}
# cover to foldchnages
# matrix1 and matrix2 are two data matrix
foldchanges <- function(mat1, mat2) {
  # log2 transformation
  # first 2 columns are gene ID and length
  logmd1 <- log(mat1[, -c(1:2)] + 1, base = 2)
  # match gene ID
  mlist <- match(mat1[, 1], mat2[, 1])
  # match sample column
  nlist <- match(colnames(mat1), colnames(mat2))
  mat2 <- mat2[mlist, nlist]
  # log2 transformation
  logmd2 <- log(mat2[mlist, -c(1:2)] + 1, base = 2)
  fcmat <- logmd1 - logmd2
  out <- cbind(mat1[, 1:2], fcmat)
  return(out)
}
# log2 transformation
takelog2 <- function(mat) {
  dd <- log(mat[,-c(1:2)] + 1, base = 2)
  out <- cbind(mat[, c(1:2)], dd)
  return(out)
}
# run correlation test
# return stats
RunCorrTest <- function(x, y) {
  pval <- c()
  cval <- c()
  for (i in colnames(x)) {
    z <- cor.test(x[, i], y[, i])
    pval <- c(pval, z$p.value)
    cval <- c(cval, z$estimate)
  }
  return(data.frame(
    'pvalue' = pval,
    'cor' = cval,
    'sample' = colnames(x)
  ))
}
# correlatin between two data matrix
CorrDataMatrix <- function(mat1, mat2, islog2) {
  # match gene ID
  if (islog2 == 1) {
    mlist <- match(mat1[, 1], mat2[, 1])
    # match sample column
    nlist <- match(colnames(mat1), colnames(mat2))
    mat2 <- mat2[mlist, nlist]
    # log2 transformation
    mat1.log <- takelog2(mat1)
    mat2.log <- takelog2(mat2)
  } else{
    mat1.log <- mat1
    mat2.log <- mat2
  }
  #mcor <-  cor(mat1.log[,-c(1:2)],mat2.log[,-c(1:2)])
  stat.cor <- RunCorrTest(mat1.log[,-c(1:2)], mat2.log[,-c(1:2)])
  return(stat.cor)
}
# summary fold changes by gene's biotype
# foldchanges: median or mean of FC cross dataset
# genes: gene annotation from gft file
# threshold: a cut-off to select significant foldchange
# values.
summaryFoldChanges <- function(foldchanges, genes, threshold) {
  upIDs <- names(which(foldchanges > threshold))
  downIDs <- names(which(foldchanges < -threshold))
  # parse up- down- genes
  upGene <- subset(genes, genes$gene_id %in% upIDs)
  downGene <- subset(genes, genes$gene_id %in% downIDs)
  upGene <- cbind('FC' = rep('UP', nrow(upGene)), upGene)
  downGene <- cbind('FC' = rep('DN', nrow(downGene)), downGene)
  # put together
  fcGene <- rbind(upGene, downGene)
  # summary up- down- FC
  fc_tb <- table(fcGene$FC, fcGene$gene_type)
  fc_tb['DN', ] <- (-1) * fc_tb['DN', ]
  return(fc_tb)
}
# Run SNN-Cliq to generate SNN graphic and
# then cluster cells based on graph to generate cluster/group
# then visualize the SNN clustering result in PCA
RunSNNCluster <- function(datamatrix) {
  mat.log <- takelog2(datamatrix)
  # First run randome PCA then extract top 2 PCs for visualization
  mat.pca <- rpca(mat.log[,-c(1:2)], scale = T)
  PCs <- mat.pca$rotation[, 1:2]
  # Run SNN graph
  grps <- buildSNNGraph(mat.log[,-c(1:2)])
  # extract memebership
  clusters <- cluster_fast_greedy(grps)
  fc <-  membership(clusters)
  out <-
    data.frame('PC1' = PCs[, 1],
               'PC2' = PCs[, 2],
               'membership' = as.factor(fc))
  return(out)
}

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
    "--gtf_file",
    type = "character",
    default = NULL,
    help = "gtf annotation file",
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
output_name <- opt$output_prefix
# color palette
palette(c("#00AFBB", "#E7B800"))
# load data matrix
mat1 <- read.csv(matrixfile1)
mat2 <- read.csv(matrixfile2)
# calculate folder changes
fc <- foldchanges(mat1, mat2)
# calucate mean or median FC cross samples
fc.dd <- fc[,-c(1:2)]
fc.m <- apply(fc.dd, 1, median)
# boxplot
# Box plot (bp)
bxp <- ggviolin(
  fc.m,
  ylab = 'log FoldChanges',
  xlab = '',
  ggtheme = theme_minimal(),
  add = "jitter"
)
bxp <- addTheme(bxp)
# correlation
# return two matrix of corraltion matrix
# take diag of correlation matrix
# plot into histogram
mcor <- CorrDataMatrix(mat1, mat2, 0)
phist <- gghistogram(
  mcor,
  x = "cor",
  fill = "lightgray",
  add = "mean",
  rug = TRUE,
  xlab = "Correlation"
)
phist <- addTheme(phist)
# parse GTF
genes <- ParseGene(gtf_file)
# summary foldchanges by biotypes
names(fc.m) <- fc[, 1]
fc_tb <- summaryFoldChanges(fc.m, genes, 0.5)
# SNN cluster
mat1.snn <- RunSNNCluster(mat1)
mat2.snn <- RunSNNCluster(mat2)
# Calculate adjusted rand index between two cluster results
# to figure out how similar the results are
# ideal, if we habe pre-labeled data, we can calculate RAND index
# between clustering result and labels
rand <- adjustedRandIndex(mat1.snn$membership,mat2.snn$membership)
# visualization
mtb <- melt(fc_tb)
pb <- ggbarplot(
  mtb,
  x = 'Var2',
  y = 'value',
  color = "Var1",
  fill = "Var1",
  # change fill color by FC
  x.text.angle = 90,
  # Rotate vertically x axis texts
  xlab = "",
  ylab = "",
  rotate = TRUE,
  ggtheme = theme_minimal(),
  legend.title = "Summary of FC genes"
)
pb <- addTheme(pb)
# combine plots
text <- paste(
  "A: Violinplot of log-foldchanges.",
  "B: Histogram of correlation between data matrix",
  "C: Barplot to summary log-foldchanges by gene biotyps",
  sep = ' '
)
# visualize SNN clustering
snn.p1 <- ggscatter(
  mat1.snn,
  x = "PC1",
  y = "PC2",
  color = "membership",
  palette = "jco",
  size = 8,
  alpha = 0.6,
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
  title = 'Updated Pipeline'
) + border()
snn.p1 <- addTheme(snn.p1)
snn.p2 <- addTheme(snn.p2)

p1 <- ggarrange(bxp,
                phist,
                labels = c("A", "B"),
                ncol = 1,
                nrow = 2)
p2 <- ggarrange(pb, labels = c("C"))
p3 <- ggarrange(p1,
                p2,
                ncol = 2,
                nrow = 1)
p4 <- ggparagraph(text,
                  face = "bold",
                  size = 20,
                  color = "black")
p5 <- ggarrange(p3,
                p4,
                ncol = 1,
                nrow = 2,
                heights = c(1, 0.1))
text <- paste(
  "A: SNN-Cliq clustering using counts matrix generated by Pipline1",
  "B: SNN-Cliq clustering using counts matrix generated by Pipline2",
  "Adjusted Rand index is",
  rand,
  sep = ' '
)
p6 <- ggarrange(snn.p1,
                snn.p2,
                labels = c("A","B"),
                ncol = 2,
                nrow = 1)
p7 <- ggparagraph(text,
                  face = "bold",
                  size = 20,
                  color = "black")
p8 <- ggarrange(p6,
                p7,
                ncol = 1,
                nrow = 2,
                heights = c(1, 0.1))
pdf(paste(output_name, '_data_matrix_comparison.pdf', sep = ''),
    30,
    30)
print(p5)
dev.off()
pdf(paste(output_name, '_cluster_data_matrix_comparison.pdf', sep = ''),
    30,
    60)
print(p8)
dev.off()
write.csv(
  mcor,
  file = paste(output_name, '_correlation_data_matrix_comparison.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
write.csv(
  fc_tb,
  file = paste(output_name, '_data_matrix_comparison.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
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