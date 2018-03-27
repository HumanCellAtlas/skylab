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
library(rtracklayer)
source('analysis_functions.R')

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
# load genes annotation
genes <- ParseGene(gtf_file)
meta_file <-  opt$metadata_file
grps <- opt$groups
output_name <- opt$output_prefix
# load meta/pheno data
meta <- read.delim(meta_file, sep = ',', header = T)
# color palette
palette(c("#00AFBB", "#E7B800"))
# load data matrix
mat1 <- read.csv(matrixfile1)
mat2 <- read.csv(matrixfile2)
mat1.log2 <- takelog2(mat1)
mat2.log2 <- takelog2(mat2)
# filter cells/samples if total expression less than 10000
mat1.log2.d <- mat1.log2[,-c(1:2)]
rownames(mat1.log2.d) <- mat1.log2[, 1]
mat2.log2.d <- mat2.log2[,-c(1:2)]
rownames(mat2.log2.d) <- mat2.log2[, 1]
mat1.d <- FilterCellsbyExp(mat1.log2.d, 10000)
mat2.d <- FilterCellsbyExp(mat2.log2.d, 10000)
# parse meta data to match filtered data matrix
mlist1 <- match(colnames(mat1.d), meta$sra)
mlist2 <- match(colnames(mat2.d), meta$sra)
grp1 <- meta[mlist1, 'cell_type']
grp2 <- meta[mlist2, 'cell_type']

# run t-test
test.res1 <- apply(mat1.d, 1, function(x) {
  RunTest(x = x, y = grp1, 'ttest')
})
test.res2 <- apply(mat2.d, 1, function(x) {
  RunTest(x = x, y = grp2, 'ttest')
})

# results
res1 <- data.frame(t(test.res1))
res2 <- data.frame(t(test.res2))
fc1 <- res1$a1 - res1$a2
names(fc1) <- rownames(res1)
fc2 <- res2$a1 - res2$a2
names(fc2) <- rownames(res2)
# refactor pvalues
res1$pvalue <-
  factor(ifelse(res1$pvalue > 0.01 | is.na(res1$pvalue), "Base.insig", "Base.sig"))
res2$pvalue <-
  factor(ifelse(res2$pvalue > 0.01 | is.na(res2$pvalue), "Updated.insig", "Updated.sig"))
gcols <- paste(res1$pvalue, res2$pvalue, sep = ',')
# join data
x <- merge(fc1, fc2, by = 'row.names')
colnames(x) <- c('ensID', 'Base', 'Updated')
x['groups'] <- gcols
sp <-
  ggscatter(
    x,
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
sp <- addTheme(sp)

# summary sig gene in biotype
x.base <- subset(x, x$groups %in% c("Base.sig,Updated.insig"))
x.updated <- subset(x, x$groups %in% c("Base.insig,Updated.sig"))
annt1 <- subset(genes, genes$gene_id %in% x.base$ensID)
annt2 <- subset(genes, genes$gene_id %in% x.updated$ensID)
tb1 <- table(annt1$gene_type)
tb2 <- table(annt2$gene_type)
tb <- rbind(tb1, -1 * tb2)
rownames(tb) <-
  c('Base.sig,Updated.insig', 'Base.insig,Updated.sig')
# visualize by barplot
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
  legend.title = ""
) + border()
pb <- addTheme(pb)
# arrange plots
gp <- ggarrange(
  sp,
  pb,
  labels = c("A: T-test Bulk -SCs", "B: Signature Genes Summary in Biotypes"),
  common.legend = FALSE,
  ncol = 2,
  nrow = 1
)
pdf(paste(output_name, '_matrix_reproducibility.pdf', sep = ''),
    30,
    30)
print(gp)
dev.off()
write.csv(
  res1,
  file = paste(
    output_name,
    "_base_ttest_res.csv",
    sep = ''),
    quote = F,
    row.names = F,
    col.names = T
)
write.csv(
  res2,
  file = paste(
    output_name,
    "_updated_ttest_res.csv",
    sep = ''),
    quote = F,
    row.names = F,
    col.names = T
)
