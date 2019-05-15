#!/usr/bin/env Rscript

library(Matrix)
library(data.table)
library(numDeriv)
library(optparse)

## Define helper functions

#' Prints a message to stderr and exits R with error code 1
#' @param msg message to standard error
errorExit <- function(msg) {
  cat(msg,file=stderr());
  quit(save='no',status=1);
}

## Parse the input arguments
option_list <- list(
  make_option(c('-i','--input-rds'),
              type='character',
              default=NULL, ## required
              dest='input_rds',
              help='input RDS file containing the data matrix in dgCMatrix format gene x droplet orientation'),
  make_option(c('-n','--n_cells_expected'),
              type='integer',
              default=3000,
              dest='n_cells_expected',
              help='number of cells expected in experiment'),
  make_option(c('-m','--filter_mode'),
              type='character',
              default='cellranger',
              dest='filter_mode',
              help='mode to use to select cell.  Options are cellranger, inflection, both, either'
              ),
  make_option(c('-o','--output_csv'),
              type='character',
              default = NULL,
              dest='output_csv',
              help='output csv file'
  )
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

## Check the parsed arguments
if(is.null(opt$input_rds))  errorExit("Input RDS is not specified\n");
if(is.null(opt$output_csv)) errorExit("Output CSV is not specified\n");
if(!file.exists(opt$input_rds)) errorExit("Input RDS doesn't exist!\n");
if(file.exists(opt$output_csv)) errorExit("Output CSV file exists!\n");
if(!opt$filter_mod %in% c("cellranger","inflection","both","either")) errorExit("Filter mode not a valid option\n");
  
feature_barcode_mat <- readRDS(opt$input_rds)

sorted_umis_per_barcode <- sort(colSums(feature_barcode_mat),decreasing = TRUE)
sorted_umis_per_barcode <- as.data.table(sorted_umis_per_barcode,keep.rownames = TRUE)
sorted_umis_per_barcode[,idx:=1:.N]

cutoff_count_cr <- sorted_umis_per_barcode[1:opt$n_cells_expected,quantile(sorted_umis_per_barcode,.99)]/10
sorted_umis_per_barcode[,is_cell_cr:=sorted_umis_per_barcode>=cutoff_count_cr]

#get inflection point by minimizing derivative (largest negative value)
smooth_spline <- smooth.spline(sorted_umis_per_barcode[sorted_umis_per_barcode>1,idx],
                               sorted_umis_per_barcode[sorted_umis_per_barcode>1,sorted_umis_per_barcode],
                               spar=0.005)

func <- function(x) {
  return(predict(smooth_spline,x)$y)
}

deriv <- function(x) {
  return(grad(func,x))
}

cut_off_point_inflection <- optimize(interval=c(cutoff_count_cr/5,cutoff_count_cr*5),f=deriv)
cut_off_count_inflection <- sorted_umis_per_barcode[idx<=cut_off_point_inflection$minimum,min(sorted_umis_per_barcode)]
sorted_umis_per_barcode[,is_cell_inflection:=sorted_umis_per_barcode>=cut_off_count_inflection]

if (opt$filter_mode=='cellranger') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr]
} else if (opt$filter_mode=='inflection') {
  sorted_umis_per_barcode[,is_cell:=is_cell_inflection]
} else if (opt$filter_mode=='both') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr & is_cell_inflection]
} else if (opt$filter_mode=='either') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr | is_cell_inflection]
}

output_table <- sorted_umis_per_barcode[,.("barcode"="rn","is_cell_cr","is_cell_inflection","is_cell")]
fwrite(sorted_umis_per_barcode,opt$output_csv)

