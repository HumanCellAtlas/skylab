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
library(rtracklayer)
# R functions used in analysis
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
    p + theme(
      legend.text = element_text(size = 15, face = 'bold'),
      legend.title = element_text(size = 15, face = 'bold')
    )
  return(p)
}
plotlm <- function(df) {
  # lm test
  x <- df$Base
  y <- df$Updated
  fit <- lm(y ~ x)
  # summary
  sfit <- summary(fit)
  # r2
  r2  <- round(sfit$adj.r.squared, 3)
  # intercept and slope
  beta <- round(sfit$coefficients[2], 3)
  a <- round(sfit$coefficients[1], 3)
  # f-stats and p values
  f <- sfit$fstatistic
  fpval <- pf(f['value'], f['numdf'], f['dendf'], lower.tail = F)
  fpval.s <- convert2star(fpval)
  # start to plot
  coefs <-
    data.frame(
      'ic' = c(0, a),
      's' = c(1, beta),
      'tl' = c('1x1', 'fitted')
    )
  my.formula <- y ~ x
  p <-
    ggplot(data = df, aes(x = Base, y = Updated)) + geom_point(shape = 1) +
    theme_bw()
  p <- p + stat_poly_eq(
    formula = my.formula,
    eq.with.lhs = "italic(hat(y))~`=`~",
    aes(
      label = paste(..eq.label.., ..rr.label.., sep = "~~~"),
      size = 8
    ),
    parse = TRUE,
    size = 8
  )
  p <-
    p + stat_fit_glance(
      method = 'lm',
      method.args = list(formula = my.formula),
      geom = 'text',
      aes(
        label = paste("P-value:", fpval , sep = ""),
        size = 12
      ),
      label.x.npc = 'left',
      label.y.npc = 0.85,
      size = 8
    )
  p <- p + theme(legend.position = "top")
  p <- addTheme(p)
  p <- p + xlab(paste('Base')) + ylab(paste('Updated'))
  p <-
    p + geom_abline(
      data = coefs,
      mapping = aes(
        slope = s,
        intercept = ic,
        linetype = factor(tl),
        color = factor(tl)
      )
    )
  return(list(
    'p' = p,
    'beta' = beta,
    'a' = a,
    'r2' = r2,
    'fpval' = fpval
  ))
}
# plot violin plot
# run mean comparison test
plotBox <- function(df) {
  mz <- melt(df)
  p <-
    ggviolin(
      mz,
      x = 'variable',
      y = 'value',
      color = "variable",
      fill = "variable",
      add = "Violin plot",
      palette = "jco",
      add.params = list(fill = "white")
    )
  p <-
    p + stat_compare_means(comparisons = list(c('Base', 'Updated')), label = "p.signif") +
    stat_compare_means(size = 10) + # Add pairwise comparisons p-value
    grids(linetype = "dashed") +
    theme_bw()
  p$layers[[2]]$aes_params$size = 1
  p$layers[[2]]$aes_params$textsize <- 10
  p <- addTheme(p)
  stat <- wilcox.test(df[, 'Base'], df[, 'Updated'])
  return(list('p' = p, 'pval' = stat$p.value))
}
# plot histogram/density plot
plothist <- function(df) {
  # melt dataframe
  mz <- melt(df)
  # group mean values
  mu <- apply(df, 2, mean)
  density.p <- ggdensity(
    mz,
    x = "value",
    add = "mean",
    rug = TRUE,
    color = "variable",
    fill = "variable"
  )
  density.p <- density.p + grids(linetype = "dashed") + theme_bw()
  density.p <- addTheme(density.p)
  return(list('p' = density.p))
}
# ks test
plotKS <- function(df) {
  mz <- melt(df)
  x <- df$Base
  y <- df$Updated
  ks <- ks.test(x, y)
  # run cdf
  cdf1 <- ecdf(x)
  cdf2 <- ecdf(y)
  # find max distance D-stats
  minMax <- seq(min(x, y), max(x, y), length.out = length(x))
  x0 <-
    minMax[which(abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))))]
  y0 <- cdf1(x0)
  y1 <- cdf2(x0)
  # plot ecdf
  ks.p <-
    ggplot(mz, aes(x = value, group = variable, color = variable)) +
    stat_ecdf(size = 1) +
    theme_bw() +
    theme(legend.position = "top") +
    ylab("ECDF")
  # do KS test
  dtable <-
    data.frame(
      'test' = 'K-S',
      'D-stats' = round(ks$statistic, 4),
      'P-value' = convert2star(ks$p.value)
    )
  dtable.p <-
    ggtexttable(dtable,
                row = NULL,
                theme = ttheme('mBlue', base_size = 18))
  # plot segement to represent the D statistics
  ks.p <-
    ks.p + geom_segment(aes(
      x = x0[1],
      y = y0[1],
      xend = x0[1],
      yend = y1[1]
    ),
    linetype = "dashed",
    color = "red")
  ks.p <-
    ks.p + geom_point(aes(x = x0[1] , y = y0[1]), color = "red", size = 4)
  ks.p <-
    ks.p + geom_point(aes(x = x0[1] , y = y1[1]), color = "red", size = 4)
  ks.p <- ks.p + ggtitle(paste("K-S Test"))
  ks.p <- addTheme(ks.p)
  ks.p <-
    ks.p + theme(legend.title = element_blank()) + annotation_custom(ggplotGrob(dtable.p))
  return(list(
    'p' = ks.p,
    'D' = ks$statistic,
    'pval' = ks$p.value
  ))
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
# calculate correlation matrix btw matrix data
CorrDataMatrix <- function(mat1, mat2, isnotlog2) {
  # match gene ID
  if (isnotlog2) {
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
  stat.cor <- RunCorrTest(mat1.log[,-c(1:2)], mat2.log[,-c(1:2)])
  return(stat.cor)
}
# remove cells if total expr < threshold
# threshold value default 10000
FilterCellsbyExp <- function(matrixdata, threshold) {
  d <- matrixdata
  tot.exp <- apply(d, 2, sum)
  rmlist <- which(tot.exp < threshold)
  sub.d <- d[,-c(rmlist)]
  return(sub.d)
}
# calculate group mean
AggregateMeanbyGroup <- function(x, y) {
  cnt <- aggregate(x, list(y), function(x) {
    sum(x > 0)
  })
  m <- aggregate(x, list(y), mean)
  out <- c(
    'n1' = cnt[1, 2],
    'n2' = cnt[2, 2],
    'a1' = m[1, 2],
    'a2' = m[2, 2]
  )
  return(out)
}
# run statistical test between data matrix with conditions
# condition can be sc vs bulk, celltype or phenotype.
RunTest <- function(x, y, testname) {
  cnt <- aggregate(x, list(y), function(x) {
    sum(x > 0)
  })
  if (testname == "ttest") {
    res <-
      t.test(x ~ y, alternative = "two.sided", var.equal = FALSE)
  } else if (testname == "wilcox") {
    res <- wilcox.test(x ~ y)
  }
  m <- aggregate(x, list(y), mean)
  pval <- res$p.value
  out <- c(
    'pvalue' = pval,
    "a1" = m[1, 2],
    "a2" = m[2, 2],
    "n1" = cnt[1, 2],
    "n2" = cnt[2, 2]
  )
  return(out)
}

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
  fc_tb['DN',] <- (-1) * fc_tb['DN',]
  return(fc_tb)
}
