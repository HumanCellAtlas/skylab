#!/usr/bin/env Rscript

# Load libraries
library(Matrix)
library(data.table)
library(numDeriv)
library(optparse)
library(ggplot2)

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

## Read the input matrix
feature_barcode_mat <- readRDS(opt$input_rds)

## Sum over genes and sort
sorted_umis_per_barcode <- sort(colSums(feature_barcode_mat),decreasing = TRUE)
sorted_umis_per_barcode <- as.data.table(sorted_umis_per_barcode,keep.rownames = TRUE)
sorted_umis_per_barcode[,idx:=1:.N]

## Replicate cellranger 2.0 functionality
## Get 1/10th 99th percentile of cells within the expected cutoff
cutoff_count_cr <- sorted_umis_per_barcode[1:opt$n_cells_expected,quantile(sorted_umis_per_barcode,.99)]/10
sorted_umis_per_barcode[,is_cell_cr:=sorted_umis_per_barcode>=cutoff_count_cr]

## Alternative way use above as starting point
## get inflection point by minimizing derivative (largest negative value)
smooth_spline <- smooth.spline(sorted_umis_per_barcode[sorted_umis_per_barcode>1,idx],
                               sorted_umis_per_barcode[sorted_umis_per_barcode>1,sorted_umis_per_barcode],
                               spar=0.005)

## Helper functions
smoother <- function(x) {
  return(predict(smooth_spline,x)$y)
}

deriv <- function(x) {
  return(grad(smoother,x))
}

## Get the inflection point
cut_off_point_inflection <- optimize(interval=c(cutoff_count_cr/5,cutoff_count_cr*5),f=deriv)
cut_off_count_inflection <- sorted_umis_per_barcode[idx <= cut_off_point_inflection$minimum,min(sorted_umis_per_barcode)]

## In place assignment in data.table
sorted_umis_per_barcode[,is_cell_inflection:=sorted_umis_per_barcode>=cut_off_count_inflection]

## This is the one used for populating the is_cell column
if (opt$filter_mode=='cellranger') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr]
} else if (opt$filter_mode=='inflection') {
  sorted_umis_per_barcode[,is_cell:=is_cell_inflection]
} else if (opt$filter_mode=='both') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr & is_cell_inflection]
} else if (opt$filter_mode=='either') {
  sorted_umis_per_barcode[,is_cell:=is_cell_cr | is_cell_inflection]
}

## Output the table
output_table <- sorted_umis_per_barcode[,.("barcode"="rn","is_cell_cr","is_cell_inflection","is_cell")]
fwrite(sorted_umis_per_barcode,opt$output_csv)

## Plots
p1 <- ggplot(sorted_umis_per_barcode[sorted_umis_per_barcode>1],aes(x=idx, y=sorted_umis_per_barcode)) +
    geom_line(aes(color=is_cell_inflection)) + scale_y_log10() + scale_x_log10() + xlab("barcode rank") + ylab("UMIs per Barcode") +
    geom_point(data=last(sorted_umis_per_barcode[is_cell_cr==TRUE,]),aes(x=idx, y=sorted_umis_per_barcode)) +
    geom_text(data=last(sorted_umis_per_barcode[is_cell_cr==TRUE,]),label="Last Real Cell (cellranger)",hjust=0, vjust=0)

xmin_zoomed <- min(last(sorted_umis_per_barcode[is_cell_cr==TRUE,idx])-100,last(sorted_umis_per_barcode[is_cell_inflection==TRUE,idx])-100)
xmax_zoomed <- max(last(sorted_umis_per_barcode[is_cell_cr==TRUE,idx])+100,last(sorted_umis_per_barcode[is_cell_inflection==TRUE,idx])+100)

ymin_zoomed <- min(last(sorted_umis_per_barcode[is_cell_cr==TRUE,sorted_umis_per_barcode])-100,last(sorted_umis_per_barcode[is_cell_inflection==TRUE,sorted_umis_per_barcode])-100)
ymax_zoomed <- max(last(sorted_umis_per_barcode[is_cell_cr==TRUE,sorted_umis_per_barcode])+100,last(sorted_umis_per_barcode[is_cell_inflection==TRUE,sorted_umis_per_barcode])+100)

p2 <- ggplot(sorted_umis_per_barcode[sorted_umis_per_barcode>1],aes(x=idx, y=sorted_umis_per_barcode)) +
    geom_line(aes(color=is_cell_inflection)) + xlim(xmin_zoomed,xmax_zoomed) + ylim(ymin_zoomed,ymax_zoomed) + xlab("barcode rank") + ylab("UMIs per Barcode") +
    geom_point(data=last(sorted_umis_per_barcode[is_cell_cr==TRUE,]),aes(x=idx, y=sorted_umis_per_barcode)) +
    geom_text(data=last(sorted_umis_per_barcode[is_cell_cr==TRUE,]),label="Last Real Cell (cellranger)",hjust=0, vjust=0)

ggsave("umis_per_barcode.png",p1)
ggsave("umis_per_barcode_zoomed.png",p2)

