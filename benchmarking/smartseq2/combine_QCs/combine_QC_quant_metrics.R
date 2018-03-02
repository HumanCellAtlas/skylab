# R script to do metrics test
# loading library
library('rtracklayer')
library(optparse)
# blacklist includes are non-numerics QC metrics
BLACKLIST <-
  c(
    'BAD_CYCLES',
    'CATEGORY',
    'LIBRARY',
    'SAMPLE',
    'READ_GROUP',
    'READ_PAIRS',
    'MEAN_READ_LENGTH',
    'PAIR_ORIENTATION'
  )
# calcaute detectable genes ratio
# input cnts is data matrix
# threshold is minimum counts/tpm
SummaryPerColumn <- function(cnts, threshold) {
  cnt.dd <- cnts[, -c(1:2)]
  detected <- apply(cnt.dd, 2, function(x) {
    sum(x > threshold)
  })
  ratio <- detected / nrow(cnts)
  return(ratio)
}
# calculate MT contents
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
  x <- subset(cnt[, -c(1:2)], cnt$gene_id %in% mt.genes$gene_id)
  # MT gene read counts
  cnt.mt <- apply(x, 2, sum)
  cnt.tot <- apply(cnt[, -c(1:2)], 2, sum)
  mt.ratio <- cnt.mt / cnt.tot
  return(mt.ratio)
}
# Combine Picard metrics and MT, detectable gene ratio into single file
CombineMetrics <- function(cnt, met, gtf_file, nthreshold) {
  ## blacklist of metrics
  met.core <- subset(met, !(met$metrics %in% BLACKLIST))
  print(dim(met.core))
  rownames(met.core) <- make.names(met.core[, 1], unique = TRUE)
  met.core <- met.core[, -1]
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
# python style input args
option_list <- list(
  make_option(
    "--datafile",
    type = "character",
    default = NULL,
    help = "dataset file name",
    metavar = "character"
  ),
  make_option(
    "--gtf",
    type = "character",
    default = NULL,
    help = "gtf annotation file",
    metavar = "character"
  ),
  make_option(
    "--metrics",
    type = "character",
    default = NULL,
    help = " metricsfile name",
    metavar = "character"
  ),
  make_option(
    "--nthreshold",
    type = "integer",
    default = 10,
    help = " threshold cut off in datafile",
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
# threshold cut off and output
nthreshold <- opt$nthreshold
output_name <- opt$out
# load gtf
gtf_file <- opt$gtf
# load quantification and QC metrics
metfile <- opt$metrics
cntfile <- opt$datafile
met <- read.csv(metfile)
cnt <- read.csv(cntfile)
# rename first column by metrics name
colnames(met)[1] <- 'metrics'
met.core <- CombineMetrics(cnt, met, gtf_file, nthreshold)
write.csv(
  met.core,
  file = paste(output_name, '_metrics_combined.csv', sep = ''),
  sep = ',',
  quote = F,
  row.names = T,
  col.names = T
)
