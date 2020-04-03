#!/usr/bin/env R

## Load Matrix Library
library(Matrix)

## Read the matrix
matrix <- readRDS('matrix.rds')

## Print info on matrix dimensions
matrix.dims <- dim(matrix)
print(paste0('Input matrix dimensions: ', matrix.dims[1] , ' x ', matrix.dims[2]))

## Put rows and columns in a defined order
matrix <- matrix[ sort(rownames(matrix)), sort(colnames(matrix)) ]

## We generate a new reference matrix here in case we want to replace
newReferenceMatrix <- drop0(matrix)
saveRDS(newReferenceMatrix,'newReferenceMatrix.rds')
rm(newReferenceMatrix)

## Generate simple diagnostic plots
png('reads_per_cell_histogram.png')
hist(rowSums(matrix),main="Reads per Cell")
dev.off()

png('reads_per_gene_histogram.png')
hist(colSums(matrix), main="Reads per Gene")
dev.off()

png('number_of_genes_per_cell.png')
hist(colSums(matrix > 1), main="Genes per Cell")
dev.off()

## Read in the reference matrix
referenceMatrix <- readRDS('referenceMatrix.rds')

## Compare matrices
matrix <- drop0(matrix)
str(matrix)

## Here we are checking the matrices for equality by looking at the
## element contents which will be identical after drop0() has been run for
## identical matrices
if(all(matrix@i == referenceMatrix@i) && all(matrix@p == referenceMatrix@p) &&
    all(matrix@Dimnames[[1]] == referenceMatrix@Dimnames[[1]]) &&
    all(matrix@Dimnames[[2]] == referenceMatrix@Dimnames[[2]]) &&
   all(matrix@x == referenceMatrix@x)) {
    print('PASS: Matrices are identical')
    quit(status=0)
} else {
    print('FAIL: Matrices differ')
    quit(status=1)
}
