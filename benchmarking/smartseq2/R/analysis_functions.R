liblist <-
  c(
    "ggplot2",
    "MASS",
    "plyr",
    "reshape2",
    "ggpubr",
    "gridExtra",
    "ggpmisc",
    "cowplot",
    "corrplot",
    "ggrepel",
    "optparse",
    "rsvd",
    "scran",
    "igraph",
    "rtracklayer",
    "knitr",
    "mclust",
    "Rtsne",
    "factoextra"
  )
for (libname in liblist) {
  suppressWarnings(suppressMessages(library(libname, character.only = TRUE)))
}

options(verbose = FALSE)
options(warn = 0)
# color palette
palette(c("#00AFBB", "#E7B800"))
knitr::opts_chunk$set(message = FALSE)
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
  if (sd(x) > 0 && sd(y) > 0) {
    beta <- round(sfit$coefficients[2, 1], 3)
    a <- round(sfit$coefficients[1, 1], 3)
    # f-stats and p values
    f <- sfit$fstatistic
    fpval <- pf(f['value'], f['numdf'], f['dendf'], lower.tail = F)
    fpval.s <- convert2star(fpval)
  }else if (sd(x) == 0 && sd(y) == 0){
    a <- 'NA'
    beta <- 'NA'
    fpval <- 'NA'
    fpval.s <- 'NA'
  }else{
    a <- round(sfit$coefficients[1], 3)
    beta <- 'NA'
    fpval <- sfit$coefficients[4]
    fpval.s <- convert2star(fpval)
  }
  
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
    theme_bw() +
    labs(caption = "linear regression test of metric between two pipeline")
  p <- p + stat_poly_eq(
    formula = my.formula,
    eq.with.lhs = "italic(hat(y))~`=`~",
    aes(
      label = paste(..eq.label.., ..rr.label.., sep = "~~~"),
      size = 4
    ),
    parse = TRUE,
    size = 4
  )
  p <-
    p + stat_fit_glance(
      method = 'lm',
      method.args = list(formula = my.formula),
      geom = 'text',
      aes(
        label = paste("P-value:", fpval , sep = ""),
        size = 8
      ),
      label.x.npc = 'left',
      label.y.npc = 0.85,
      size = 4
    )
  p <- p + theme(legend.position = "top")
  
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
      add.params = list(fill = "white")
    )
  p <-
    p + stat_compare_means(comparisons = list(c('Base', 'Updated')), label = "p.signif") +
    stat_compare_means() + # Add pairwise comparisons p-value
    grids(linetype = "dashed") +
    theme_bw() + theme(legend.position = "top") +
    labs(caption = "violin plot of metric, labeled with Wilcoxon test results")
  #p$layers[[2]]$aes_params$size = 1
  #p$layers[[2]]$aes_params$textsize <- 4
  
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
  density.p <- density.p +
    grids(linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(caption = "Density plot of metric of two pipelines")
  
  return(list('p' = density.p))
}
# ks test
plotKS <- function(df) {
  mz <- melt(df)
  x <- df$Base
  y <- df$Updated
  ks <- suppressWarnings(suppressMessages(ks.test(x, y)))
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
    ylab("ECDF") + labs(caption = paste(
      "D-stats: ",
      round(ks$statistic, 2),
      "significant level: ",
      convert2star(ks$p.value),
      sep = ""
    ))
  # do KS test
  
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
  
  ks.p <-
    ks.p + theme(legend.title = element_blank())
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
foldchanges <- function(mat1.log2, mat2.log2) {
  # log2 transformation
  # first 2 columns are gene ID and length
  mat.fc <- mat1.log2 - mat2.log2
  return(mat.fc)
}
# log2 transformation
takelog2 <- function(mat) {
  dd <- log(mat + 1, base = 2)
  return(dd)
}
RunCorrTest <- function(x, y) {
  pval <- c()
  cval <- c()
  for (i in colnames(x)) {
    z <- cor.test(x[, i], y[, i], method = 'spearman')
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
    # log2 transformation
    mat1.log <- takelog2(mat1)
    mat2.log <- takelog2(mat2)
  } else{
    mat1.log <- mat1
    mat2.log <- mat2
  }
  stat.cor <- RunCorrTest(mat1.log, mat2.log)
  return(stat.cor)
}
# remove cells if total expr < threshold
# threshold value default 10000
FilterCellsbyExp <- function(matrixdata, threshold) {
  d <- matrixdata
  tot.exp <- apply(d, 2, sum)
  rmlist <- which(tot.exp < threshold)
  sub.d <- d[, -c(rmlist)]
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
  fc_tb['DN', ] <- (-1) * fc_tb['DN', ]
  return(fc_tb)
}
# Run tsne to reduce dimensions and then
RunClustering <- function(mat.log) {
  mat.pcs <- rpca(mat.log, scale = T)
  x <- mat.pcs$rotation[, c(1:50)]
  tsne <-
    Rtsne(
      t(mat.log),
      check_duplicates = FALSE,
      dims = 2,
      perplexity = 30,
      max_iter = 500,
      verbose = FALSE
    )
  grps <- buildSNNGraph(t(x),
                        rand.seed = 1000)
  clusters <- cluster_fast_greedy(grps)
  rownames(tsne$Y) <- colnames(mat.log)
  names(clusters$membership) <- colnames(mat.log)
  return(list(
    'tsne' = tsne,
    'clusters' = clusters,
    'graph' = grps
  ))
}
# Run SNN-Cliq to generate SNN graphic and
# then cluster cells based on graph to generate cluster/group
# then visualize the SNN clustering result in PCA
RunSNNCluster <- function(mat.log) {
  # First run randome PCA then extract top 2 PCs for visualization
  mat.pca <- rpca(mat.log, scale = T)
  PCs <- mat.pca$rotation[, 1:2]
  # Run SNN graph
  grps <- buildSNNGraph(mat.log)
  # extract memebership
  clusters <- cluster_fast_greedy(grps)
  fc <-  membership(clusters)
  out <-
    data.frame('PC1' = PCs[, 1],
               'PC2' = PCs[, 2],
               'membership' = as.factor(fc))
  return(out)
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
##convert p-value to star char
convert2star <- function(x) {
  if (x > 0.05) {
    p <- 'ns'
  } else if (x > 0.01) {
    p <- '*'
  } else if (x > 0.001) {
    p <- '**'
  } else if (x > 0.0001) {
    p <- '***'
  } else{
    p <- '****'
  }
  return(p)
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
    cnt.pca <- rpca(cnts, scale = T)
  } else{
    cnt.pca <- prcomp(cnts, scale = T)
  }
  # take top npcs from rotation matrix
  pcs <- cnt.pca$rotation[, 1:npcs]
  nc1 <- ncol(qc_mets)
  # merge qc metrics with pcs by rowname which is sample ID
  dt <- merge(qc_mets, pcs, by = 0)
  # run correlation test
  res <- cor.mtest(dt[,-1])
  pmat <- res$p # extract p values
  # pmat is square matrix
  # each row and column represent qc metrics + top pcs
  colnames(pmat) <- colnames(dt[,-1])
  rownames(pmat) <- colnames(dt[,-1])
  # calculate the correlation
  M <- cor(dt[,-1])
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
SummaryPerColumn <- function(cnts, threshold) {
  cnt.dd <- cnts[,-c(1:2)]
  detected <- apply(cnt.dd, 2, function(x) {
    sum(x > threshold)
  })
  ratio <- detected / nrow(cnts)
  return(ratio)
}
# calculate MT contents
# input gtf_file(gencode annotation, version v2 of gff file)
# input cnt, the data matrix, can be either TPM or counts
# return the ratio of
# total reads/TPM in MT genes vs total reads/TPM per sample
ParseMTGene <- function(gtf_file, cnt) {
  gtf_gencode <-
    readGFF(
      gtf_file,
      version = 2L,
      tags = c("gene_name", "gene_id", "transcript_id", "gene_type")
    )
  genes <- subset(gtf_gencode, gtf_gencode$type == "gene")
  mt.genes <-
    subset(genes, genes$gene_type %in% c('Mt_tRNA', 'Mt_rRNA'))
  # select MT gene IDs
  x <- subset(cnt[,-c(1:2)], cnt$gene_id %in% mt.genes$gene_id)
  # MT gene read counts
  cnt.mt <- apply(x, 2, sum)
  cnt.tot <- apply(cnt[,-c(1:2)], 2, sum)
  mt.ratio <- cnt.mt / cnt.tot
  return(mt.ratio)
}
# Combine Picard metrics with MT and
# detectable gene ratio into single file
CombineMetrics <- function(cnt, met, gtf_file, nthreshold) {
  ## blacklist of metrics
  met.core <- subset(met,!(met$metrics %in% BLACKLIST))
  rownames(met.core) <- make.names(met.core$metrics, unique = TRUE)
  met.core <- met.core[,-1]
  ## combine QC,summary of quantification
  cnt.ratio <- round(SummaryPerColumn(cnt, nthreshold), 5)
  mt.ratio <- round(ParseMTGene(gtf_file, cnt), 5)
  ## combine wiht QC
  mlist <- match(colnames(met.core), names(cnt.ratio))
  x <- data.frame(cnt.ratio[mlist])
  colnames(x) <- 'detected_ratio'
  mlist <- match(colnames(met.core), names(mt.ratio))
  y <- data.frame(mt.ratio[mlist])
  colnames(y) <- 'MT_ratio'
  met.core <- rbind(t(x), t(y), met.core)
  return(met.core)
}

RunVarianceAnalysis <- function(expn, df) {
  metKeys <- colnames(df)
  fo <- as.formula(paste('expn', "~", paste(metKeys, collapse = "+")))
  fit <- lm(fo, data = df)
  m <- anova(fit)
  # total sum of sq
  m.tss <- sum(m[, 2])
  # % of total sum of sq of each predictor
  m.ss <- m[, 2] / m.tss
  names(m.ss) <- c(colnames(df), 'residual')
  return(m.ss)
}

SelectGeneByVars <- function(x, k) {
  vlist <- apply(x, 1, var)
  olist <- order(vlist, decreasing = T)
  x.sub <- x[olist[1:k], ]
  return(x.sub)
}
findProportion <- function(x, v) {
  m <- apply(x, 1, sum)
  p <- m / sum(m)
  pv <- aggregate(p,
                  by = list(v),
                  FUN = sum,
                  na.rm = TRUE)
  return(pv)
}
CumulateVar <- function(expn.pca) {
  # Eigenvalues
  eig <- (expn.pca$sdev) ^ 2
  # Variances in percentage
  variance <- eig * 100 / sum(eig)
  # Cumulative variances
  cumvar <- cumsum(variance)
  expn.var <- data.frame(eig = eig,
                         variance = variance,
                         cumvariance = cumvar)
  
  return(expn.var)
}
RunPCAMetrix <- function(mets) {
  mets.x <- apply(mets, 1, function(x) {
    (x - mean(x)) / sd(x)
  })
  mets.x <- data.frame(t(na.omit(t(mets.x))))
  mets.pca <- prcomp(t(mets.x), retx = T)
  mets.var <- CumulateVar(mets.pca)
  return(list('PCs' = mets.pca, 'CumuVar' = mets.var))
}
digitalizeMat <- function(data, n, m) {
  d <- matrix(0, nrow(data), ncol(data))
  for (i in 1:ncol(data)) {
    x <- data[, i]
    d[, i] <-
      cut(
        x,
        breaks = c(-Inf, n, m, Inf),
        labels = c(-1, 0, 1),
        right = FALSE
      )
  }
  return(d)
}
SummaryQuantificationByType <- function(data, genes, nthreshold) {
  biotypes <- genes[match(rownames(data), genes$gene_id), 'gene_type']
  headers <- colnames(data)
  summ <- data.frame()
  for (i in 1:ncol(data)) {
    x <- data[, i]
    y <-
      data.frame(aggregate(x == nthreshold, by = list(biotypes), FUN = sum))
    colnames(y) <- c('biotype', headers[i])
    if (i == 1) {
      summ <- y
    } else{
      summ <- merge(summ, y, by = 'biotype')
    }
  }
  return(summ)
}
