source('/usr/local/scripts/analysis_functions.R')
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
mat1 <- read.csv(opt$bdatafile)
mat2 <- read.csv(opt$udatafile)
npcs <- opt$npcs
output_name <- opt$out
## before run PCA, do log2 transformation
rownames(mat1)<-mat1[,1]
rownames(mat2)<-mat2[,1]
mat1.d <- mat1[, -c(1:2)]
mat2.d <- mat2[, -c(1:2)]
mlist <- match(rownames(mat1.d), rownames(mat2.d))
# match sample column
nlist <- match(colnames(mat1.d), colnames(mat2.d))
mat2.d <- mat2.d[mlist, nlist]
mat1.log2 <- takelog2(mat1.d)
mat2.log2 <- takelog2(mat2.d)

# Do correlation analysis between data matrix and QC metrics.
# Use *rpca* to run PCA analysis on data matrix
# Run correlation analysis by using two pipeline results separately
pmat1 <-
  CorrQCvsPCs(t(met1),
              mat1.log2,
              npcs,
              paste(output_name, '_rpca_base', sep = ''),
              cmd = 'rpca')
pmat2 <-
  CorrQCvsPCs(t(met2),
              mat2.log2,
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
## fold changes of data file(matrix) and do rPCA
logfc <- foldchanges(mat1.log2, mat2.log2)
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
