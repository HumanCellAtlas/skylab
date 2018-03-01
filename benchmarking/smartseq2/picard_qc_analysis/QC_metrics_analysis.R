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
library(optparse)

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
## standard theme
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
  fpval <- pf(f[1], f[2], f[3], lower.tail = F)
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
    ggplot(data = df, aes(x = Base, y = Updated)) + geom_point(shape = 1) + theme_bw(base_size = 20)
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
        label = paste("P-value:", fpval.s , sep = ""),
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
  outs <-
    list(
      'p' = p,
      'beta' = beta,
      'a' = a,
      'r2' = r2,
      'fpval' = fpval
    )
}
plotBox <- function(df) {
  mz <- melt(df)
  p <-
    ggviolin(
      mz,
      x = 'variable',
      y = 'value',
      color = "variable",
      palette = c("#00AFBB", "#E7B800"),
      fill = "variable",
      add = "boxplot",
      add.params = list(fill = "white")
    )
  p <-
    p + stat_compare_means(comparisons = list(c('Base', 'Updated')), label = "p.signif") +
    stat_compare_means() + # Add pairwise comparisons p-value
    grids(linetype = "dashed")
  p <- addTheme(p)
  stat <- wilcox.test(x, y)
  
  return(list('p' = p, 'pval' = stat$p.value))
}
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
    fill = "variable",
    palette = c("#00AFBB", "#E7B800")
  )
  density.p <- density.p + grids(linetype = "dashed")
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
    theme_bw(base_size = 20) +
    theme(legend.position = "top") +
    ylab("ECDF") +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))
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
# metrics files
# metfile1 is the production pipeline results
# metfile2 is the updated pipeline results
## params,python style
option_list <- list(
  make_option(
    "--bmetrics",
    type = "character",
    default = NULL,
    help = " First/Base metrics file name",
    metavar = "character"
  ),
  make_option(
    "--umetrics",
    type = "character",
    default = NULL,
    help = " Second/Updated metrics file name",
    metavar = "character"
  ),
  make_option(
    "--out",
    type = "character",
    default = "out",
    help = "output file name [default= %default]",
    metavar = "character"
  ),
  make_option(
    "--metKeys",
    type = "character",
    default = "out",
    help = "A list of metrics to analyze",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
## load params
output_name <- opt$out
metKeys <- strsplit(opt$metKeys, split = ',')[[1]]
print(metKeys)
met1 <- read.csv(opt$bmetrics, row.names = 1)
met2 <- read.csv(opt$umetrics, row.names = 1)
# match column names between two metrics tables
colnames1 <- colnames(met1)
colnames2 <- colnames(met2)
mlist <- match(colnames1, colnames2)
nalist <- is.na(mlist)
if (sum(nalist) > 0) {
  print("input files have different column elements")
  exit()
}
# re-order
met2 <- met2[, mlist]
# select subset of metrics
met1.core <- met1[metKeys, ]
met2.core <- met2[metKeys, ]
out <- c()
pouts <- list()
for (ii in 1:length(metKeys)) {
  x <- as.numeric(met1.core[metKeys[ii], ])
  y <- as.numeric(met2.core[metKeys[ii], ])
  z <- data.frame('Base' = x, 'Updated' = y)
  # linear regression model
  p1 <- plotlm(z)
  # plot histogram
  p2 <- plothist(z)
  # plot box
  p3 <- plotBox(z)
  # plot KS
  p4 <- plotKS(z)
  # arranage plots in one figure
  gp <-
    ggarrange(
      p2$p,
      p3$p,
      p1$p,
      p4$p,
      labels = c("A", "B", "C", "D"),
      common.legend = TRUE, legend = "bottom"
    )
  gp <-
    gp + ggtitle(paste(metKeys[ii])) + theme(plot.title = element_text(
      hjust = 0.5,
      size = 20,
      face = 'bold'
    ))
  text <-
    paste(
      "A: Density plot of",
      metKeys[ii],
      " from two pipeline.",
      "B: Boxplot of ",
      metKeys[ii],
      ", labeled with Wilcoxon test result.",
      "C: Linear regression test of ",
      metKeys[ii],
      "between two pipelines.",
      "D: K-S test, D statistics is labeled",
      sep = " "
    )
  text.p <-
    ggparagraph(
      text = text,
      face = "bold",
      size = 20,
      color = "black"
    )
  gp<-ggarrange(gp,text.p,
            ncol = 1, nrow = 2,
            heights = c(1, 0.1),
            common.legend = TRUE)
  pouts[[ii]] <- gp
  out <-
    rbind(out,
          c(
            metKeys[ii],
            p1$beta,
            p1$a,
            p1$r2,
            p1$fpval,
            p3$pval,
            p4$D,
            p4$pval
          ))
}

# save multiple page into one pdf
pdf(paste(output_name, '_group_plots_all.pdf', sep = ''), 25, 25)
for (ii in 1:nrow(met1.core)) {
  print(pouts[[ii]])
}
dev.off()
colnames(out) <-
  c('metrics',
    'beta',
    'a',
    'r2',
    'pvalue',
    'wilcoxon',
    'ks-D-stats',
    'ks-Pvalue')
write.table(
  out,
  file = paste(output_name, '_tests_stats.csv', sep = ''),
  quote = F,
  row.names = F,
  col.names = T,
  sep = ','
)
