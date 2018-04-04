source('/usr/local/scripts/analysis_functions.R')
## inputs, python style
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
palette(c("#00AFBB", "#E7B800"))
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
      labels = c("A", "B", "C", "D")
    ) # arrange 4 plots
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
    ) # text block for each test
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
          )) # test, statistics outputs
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
