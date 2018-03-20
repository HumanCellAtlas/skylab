# R script to do metrics test
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
library('rtracklayer')
set.seed(42)
# This R script is to evaluate the changes/difference in two data matrix
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
gtf_file <- gtf_file
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

# parse GTF
genes <- ParseGene(gtf_file)
# summary foldchanges by biotypes
names(fc.m) <- fc[,1]
fc_tb <- summaryFoldChanges(fc.m,genes,0.5)
# plot
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
text <- paste("A: Violin of log-foldchanges.", sep = ' ')
text.a <-
  ggparagraph(
    text = text,
    face = "italic",
    size = 20,
    color = "black"
  )
text <-
  paste("B: Barplot to summary log-foldchanges by gene biotyps", sep = ' ')
text.b <-
  ggparagraph(
    text = text,
    face = "italic",
    size = 20,
    color = "black"
  )
p1 <- ggarrange(bxp, pb, labels = c("A", "B"))
p2 <- ggarrange(text.a, text.b, ncol = 2, nrow = 1)
p3 <- ggarrange(p1,
                p2,
                ncol = 1,
                nrow = 2,
                heights = c(1, 0.1))
ggsave(
  p3,
  file = paste(output_name, '_data_matrix_comparison.pdf', sep = ''),
  width = 30,
  height = 30,
  limitsize = FALSE
)
write.csv(
  fc_tb,
  file = paste(output_name, '_data_matrix_comparison.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)