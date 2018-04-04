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
# input threshold is minimum counts/tpm
# return the ratio of detected genes
# over total number genes
source('/usr/local/scripts/analysis_functions.R')
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
# parse input parameters and load data
gtf_file <- opt$gtf
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
