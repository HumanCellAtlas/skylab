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
set.seed(42)

# if p<0.05, return 1 value.
convert2SigLevels <- function(P) {
  x <- floor(-1 * log(P, base = 10))
  x[x > 1] <- 1
  return(x)
}
# cover to foldchnages
# mat1 and mat2 are two data matrix
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
  dd <- log(mat[, -c(1:2)] + 1, base = 2)
  out <- cbind(mat[, c(1:2)], dd)
  return(out)
}
# correlation between qc metrics and data matrics.
# first run PCA on data matrix and top top N PCs
# Then calculate correlation and correlation test between
# QC metrics and N PCs
# Visualize this correlation and test results in corrplot
# Y-axis rank by total number of significant correlation test
# QCmetrics with High rank(on the top) indicates these metrics have
# high impact on data matrix
# input cmd: wither rpca or pca. rpca is randomized PCA
CorrQCvsPCs <- function(qc_mets, cnts, npcs, output_name, cmd) {
  if (cmd == "rpca") {
    cnt.pca <- rpca(cnts[, -c(1:2)], scale = T)
  } else{
    cnt.pca <- prcomp(cnts[, -c(1:2)], scale = T)
  }
  # take top npcs from rotation matrix
  pcs <- cnt.pca$rotation[, 1:npcs]
  nc1 <- ncol(qc_mets)
  # merge qc metrics with pcs by rowname which is sample ID
  dt <- merge(qc_mets, pcs, by = 0)
  # run correlation test
  res <- cor.mtest(dt[, -1])
  pmat <- res$p # extract p values
  # pmat is square matrix
  # each row and column represent qc metrics + top pcs
  colnames(pmat) <- colnames(dt[, -1])
  rownames(pmat) <- colnames(dt[, -1])
  # calculate the correlation
  M <- cor(dt[, -1])
  # nc2 should be # pcs + # qc metrics
  nc2 <- ncol(M)
  # take subset of M, only use correlation matrix between qc metrics vs PCs
  M.sub <- M[(nc1 + 1):nc2, 1:nc1]
  p.sub <- pmat[(nc1 + 1):nc2, 1:nc1]
  # rank column by # of significant correlation test between qc metrics and PCs
  rank.p <- order(colSums(p.sub < 0.05, na.rm = T), decreasing = T)
  png(
    paste(output_name, '_corrplot_top_', npcs, '_vs_qc_mets.png', sep = ''),
    width = 5000,
    height = 5000
  )
  corrplot(
    t(M.sub[, rank.p]),
    p.mat = t(p.sub[, rank.p]),
    is.corr = TRUE,
    insig = "label_sig",
    sig.level = c(.001, .01, .05),
    pch.cex = 1.5,
    pch.col = "white",
    na.label = " ",
    tl.cex = 5,
    cl.cex = 5,
    tl.col = "black",
    method = "color",
    col = colorRampPalette(c("blue", "white", "red"))(200)
  )
  dev.off()
  # use the full correlation matrix
  # rank column/row by # of significant correlation test between qc metrics vs PCs
  rank.p <- order(rowSums(pmat < 0.05, na.rm = T), decreasing = T)
  png(
    paste(
      output_name,
      '_corrplot_top_',
      npcs,
      '_vs_qc_mets_all.png',
      sep = ''
    ),
    width = 8000,
    height = 8000
  )
  corrplot(
    t(M[rank.p, rank.p]),
    p.mat = t(pmat[rank.p, rank.p]),
    is.corr = TRUE,
    insig = "label_sig",
    sig.level = c(.001, .01, .05),
    pch.cex = 1.5,
    pch.col = "white",
    na.label = " ",
    tl.cex = 5,
    cl.cex = 5,
    tl.cor = 'black',
    method = "color",
    col = colorRampPalette(c("blue", "white", "red"))(200)
  )
  dev.off()
  # only return the subset of correlation test results, which is qc metrics vs PCs only
  return(p.sub)
}
# add extra theme
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
# xy-plot: # significant correlation between QC and data
# x-axis is the # significant test from one pipeline
# y-axis is the # significant test form another pipeline
# inputs are correlation tesxt matrix, such as output of CorrQCvsPCs()
plotSignCorr <- function(base_mat, updated_mat, output_name) {
  # sum significant test
  sum.sig1 <- colSums(base_mat < 0.05, na.rm = T)
  sum.sig2 <- colSums(updated_mat < 0.05, na.rm = T)
  # merge by qc metrics names
  df <- merge(sum.sig1, sum.sig2, by = 0)
  colnames(df) <- c('QC', 'Base', 'Updated')
  # scatter plot and 1x1 line
  pl <- ggplot(df) +
    geom_point(aes(Base, Updated), color = 'red') +
    geom_text_repel(aes(Base, Updated, label = QC),
                    size = 6.0,
                    alpha = 0.85) +
    theme_classic(base_size = 24) +
    geom_abline(slope = 1,
                intercept = 0,
                color = "blue")
  pl <-
    pl + ggtitle(paste(
      "# of significant correlation between QC metrics and Expression PCs"
    ))
  pl <- addTheme(pl)
  ggsave(
    pl,
    file = paste(output_name, '_sum_nsig_p_xyplot.png', sep = ''),
    width = 25,
    height = 25,
    limitsize = FALSE
  )
}
# python style input args
option_list <- list(
  make_option(
    "--bdatafile",
    type = "character",
    default = NULL,
    help = "data matrix file name",
    metavar = "character"
  ),
  make_option(
    "--udatafile",
    type = "character",
    default = NULL,
    help = "updated data matrix file name",
    metavar = "character"
  ),
  make_option(
    "--bmetrics",
    type = "character",
    default = NULL,
    help = " base metrics file name",
    metavar = "character"
  ),
  make_option(
    "--umetrics",
    type = "character",
    default = NULL,
    help = " updated metrics file name",
    metavar = "character"
  ),
  make_option(
    "--npcs",
    type = "integer",
    default = 10,
    help = " number of PCs to collect",
    metavar = "number"
  ),
  make_option(
    "--out",
    type = "character",
    default = "out",
    help = "output file name [default= %default]",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
## load files and variables
met1 <- read.csv(opt$bmetrics, row.names = 1)
met2 <- read.csv(opt$umetrics, row.names = 1)
cnt1 <- read.csv(opt$bdatafile)
cnt2 <- read.csv(opt$udatafile)
npcs <- opt$npcs
output_name <- opt$out
## before run PCA, do log2 transformation
logcnt1<-takelog2(cnt1)
logcnt2<-takelog2(cnt2)
# Do correlation analysis between data matrix and QC metrics.
# Use *rpca* to run PCA analysis on data matrix
# Run correlation analysis by using two pipeline results separately
pmat1 <-
  CorrQCvsPCs(t(met1),
              logcnt1,
              npcs,
              paste(output_name, '_rpca_base', sep = ''),
              cmd = 'rpca')
pmat2 <-
  CorrQCvsPCs(t(met2),
              logcnt2,
              npcs,
              paste(output_name, '_rpca_updated', sep = ''),
              cmd = 'rpca')
write.csv(
  pmat1,
  file = paste(output_name, '_base_rpca_cor_test.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
write.csv(
  pmat2,
  file = paste(output_name, '_updated_rpca_cor_test.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
# Then compare two pipelines correlation test results
# in scatter plot, should look for data points move away from 1x1 line
plotSignCorr(pmat1, pmat2, paste(output_name, '_rpcs', sep = ''))

# second analysis
# delta vs delta kind of analysis
# delta values in data matrix vs delta value in qc metrics
# delta vs delta analysis can tell us the impact of changes in QC metrics to changes in data matrix
# do fold changes
## folder changes of data file(matrix) and do rPCA
logfc <- foldchanges(cnt1, cnt2)
## difference in metrics
delta <- met1 - met2
## correlation to randomized PCs
pmat3 <-
  CorrQCvsPCs(t(delta),
              logfc,
              npcs,
              paste(output_name, '_fc_rpca_delta', sep = ''),
              cmd = 'rpca')
write.csv(
  pmat3,
  file = paste(output_name, '_fc_rpca_delta_cor_test.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
