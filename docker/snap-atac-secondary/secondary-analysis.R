library("optparse")

args = commandArgs(trailingOnly=TRUE)

path_input_snap <- args[1]
path_outputs <- args[2]
path_data <- args[3]
sample_name <- args[4]
num_cores <- 1

# option_list = list(
#   make_option(c("-i", "--input"), type="character", default=NULL,
#               help="snap filename", metavar="character"),
#   make_option(c("-o", "--path_output"), type="character", default=NULL,
#               help="path to output", metavar="character"),
#   make_option(c("-d", "--path_data"), type="character", default=NULL,
#               help="path to data", metavar="character"),
#   make_option(c("-s", "--sample"), type="character", default=NULL,
#               help="sample name", metavar="character")
# );

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

# fixme:
# do some validation

# fixme: use optparse
# path_outputs <- opt$path_output
# path_input_snap <- opt$input
# path_data <- opt$path_data
# sample_name <- opt$sample
# num_cores <- 1

print(path_input_snap)
print(path_outputs)
print(path_data)
print(sample_name)

# stop("Test")

library(SnapATAC);

# -----------------------------------------------------
# Step 1. Barcode selection (SnapATAC)
print("Step 1");

x.sp = createSnap(
  file=path_input_snap,
  sample=sample_name,
  num.cores=num_cores
);

# fixme
# path_barcode_pdf <- file.path(path_outputs, "barcode.pdf")
path_barcode_pdf <- "barcode.pdf"

plt <- plotBarcode(
  obj=x.sp,
  pdf.file.name=path_barcode_pdf,
  pdf.width=7,
  pdf.height=7,
  col="grey",
  border="grey",
  breaks=50
);

# -----------------------------------------------------
# filter cells only using number of fragments and UMI with the following cutoffs

x.sp = filterCells(
  obj=x.sp,
  subset.names=c("fragment.num", "UMI"),
  low.thresholds=c(1000,1000),
  high.thresholds=c(Inf, Inf)
);

# -----------------------------------------------------
# Step 2. Bin size selection (SnapATAC)
print("Step 2");

# show what bin sizes exist in atac_v1_adult_brain_fresh_5k.snap file
showBinSizes(path_input_snap);

x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=1);

# -----------------------------------------------------
# Step 3. Fragments-in-promoter ratio.
print("Step 3");

library(GenomicRanges);
promoter.df = read.table(file.path(path_data, "promoter.bed"));
promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]));
ov = findOverlaps(x.sp@feature, promoter.gr);
idy = queryHits(ov);
promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat");

# fixme
# path_bin_size_png <- file.path(path_outputs, "bin-size.png")
path_bin_size_png <- "bin-size.png"
plot(log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10), promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ));
# dev.copy(png, filename=path_bin_size_png);
# dev.off();

idx = which(promoter_ratio > 0.2 & promoter_ratio < 0.8);
x.sp = x.sp[idx,];

# fixme: write to file
x.sp;

# fixme: write to file
summarySnap(x.sp);

# -----------------------------------------------------
# Step 4. Matrix binarization (SnapATAC)
print("Step 4");

# We next convert the cell-by-bin count matrix to a binary matrix. We found some items in the matrix have abnormally high coverage perhaps due to the alignment error. Therefore, we first remove top 0.1% items in the count matrix followed by converting the rest of the values into binary.
x.sp = makeBinary(x.sp, mat="bmat");

# -----------------------------------------------------
# Step 5. Bin filtration (SnapATAC)
# We next filter out any bins overlapping with the ENCODE blacklist and bins belonging to chrM or random chromsomes to prevent from any potential artifacts.
print("Step 5");

library(GenomicRanges);
black_list = read.table(file.path(path_data, "mm10.blacklist.bed.gz"));
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
);
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];
x.sp

# fixme
# path_bin_coverage_pdf <- file.path(path_outputs, "bin-coverage.pdf")
path_bin_coverage_pdf <- "bin-coverage.pdf"
plotBinCoverage(
  x.sp,
  pdf.file.name=path_bin_coverage_pdf,
  col="grey",
  border="grey",
  breaks=10,
  xlim=c(-6,6)
);

x.sp = filterBins(
  x.sp,
  low.threshold=-1.5,
  high.threshold=1.5,
  mat="bmat"
);
x.sp

# -----------------------------------------------------
# Step 6. Jaccard matrix (SnapATAC)
print("Step 6");

x.sp = runJaccard(
  obj = x.sp,
  tmp.folder=tempdir(),
  mat = "bmat",
  max.var=2000
);

# -----------------------------------------------------
# Step 7. Normalization (SnapATAC)
print("Step 7");

x.sp = runNormJaccard(
  obj = x.sp,
  tmp.folder=tempdir(),
  ncell.chunk=1000,
  method="normOVE",
  row.center=TRUE,
  row.scale=TRUE,
  low.threshold=-5,
  high.threshold=5,
  do.par=TRUE,
  num.cores=5,
  seed.use=10
);

# -----------------------------------------------------
# Step 8. Linear Dimentionality Reduction (SnapATAC)
print("Step 8");

x.sp = runDimReduct(
  x.sp,
  pc.num=50,
  input.mat="jmat",
  method="svd",
  center=TRUE,
  scale=FALSE,
  seed.use=10
);

# -----------------------------------------------------
# Step 9. Determine statistically significant principal components (SnapATAC)
print("Step 9");

# fixme
path_pca_elbow_pdf <- "pca-elbow.pdf"
plotDimReductElbow(
  obj=x.sp,
  point.size=1.5,
  point.shape=19,
  point.color="red",
  point.alpha=1,
  pdf.file.name=path_pca_elbow_pdf,
  pdf.height=7,
  pdf.width=7,
  labs.title="PCA Elbow plot",
  labs.subtitle=NULL
);

# fixme
path_pca_pw_pdf <- "pca-pw.pdf"
plotDimReductPW(
  obj=x.sp,
  pca.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=path_pca_pw_pdf,
  pdf.height=7,
  pdf.width=7
);

# -----------------------------------------------------
# Step 10. KNN Graph Construction (SnapATAC)
print("Step 10");

x.sp = runKNN(
  obj=x.sp,
  pca.dims=1:40,
  weight.by.sd=FALSE,
  k=15
);

# -----------------------------------------------------
# Step 11. Clustering (SnapATAC)
print("Step 11");

x.sp = runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);

# alternative (requires leiden, leidenalg)
# library(leiden);
# x.sp = runCluster(
#   obj=x.sp,
#   tmp.folder=tempdir(),
#   louvain.lib="leiden",
#   seed.use=10,
#   resolution=1
# );

# -----------------------------------------------------
# Step 12. Non-linear dimentionality reduction (SnapATAC)
print("Step 12");

x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  pca.dims=1:40,
  weight.by.sd=FALSE,
  method="Rtsne",
  fast_tsne_path=NULL,
  Y.init=NULL,
  seed.use=10,
  num.cores=5
);

# this requires umap
# https://cran.r-project.org/web/packages/umap/index.html
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  pca.dims=1:40,
  weight.by.sd=FALSE,
  method="umap",
  fast_tsne_path=NULL,
  Y.init=NULL,
  seed.use=10,
  num.cores=5
);

# -----------------------------------------------------
# Step 13. Visulization (SnapATAC)
print("Step 13");

# fixme:
path_tsne_pdf <- "tsne.pdf"
plotViz(
  obj=x.sp,
  method="tsne",
  point.size=0.5,
  point.shape=19,
  point.alpha=0.8,
  point.color="cluster",
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  pdf.file.name=path_tsne_pdf,
  pdf.width=7,
  pdf.height=7
);

# fixme:
path_umap_pdf <- "umap.pdf"
plotViz(
  obj=x.sp,
  method="umap",
  point.size=0.5,
  point.shape=19,
  point.alpha=0.8,
  point.color="cluster",
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=TRUE,
  pdf.file.name=path_umap_pdf,
  pdf.width=7,
  pdf.height=7
);

# feature.value = SnapATAC::rowSums(x.sp@bmat);
# feature.value = pmin(feature.value, quantile(feature.value, 0.99));
# feature.value = pmax(feature.value, 0);
# feature.value = (feature.value-min(feature.value))/(max(feature.value)-min(feature.value));

# PlotFeatureSingle(
#   obj=x.sp,
#   feature.value=feature.value,
#   method="tsne",
#   point.size=0.3,
#   point.shape=19,
#   point.color="red",
#   down.sample=10000,
#   pdf.file.name=NULL,
#   pdf.width=7,
#   pdf.height==7
# );

# PlotFeatureSingle(
#   obj=x.sp,
#   feature.value=feature.value,
#   method="umap",
#   point.size=0.2,
#   point.shape=19,
#   point.color="red",
#   down.sample=10000,
#   pdf.file.name=NULL,
#   pdf.width=7,
#   pdf.height==7
# );

# -----------------------------------------------------
# Step 14. Gene-body based annotation for expected cell types (SnapATAC)
print("Step 14");

genes = read.table(file.path(path_data, "gencode.vM16.gene.bed"));
genes.gr = GRanges(genes[,1],
                     IRanges(genes[,2], genes[,3]),
                     name=genes[,4]
);
marker.genes = c(
  "Snap25", "Gad2", "Apoe",
  "C1qb", "Pvalb", "Vip",
  "Sst", "Lamp5", "Slc17a7",
  "Mog", "Pdgfra", "Cspg4",
  "Cx3cr1","F3","Aqp4",
  "Rorb"
);
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
x.sp = createGmat(
  obj=x.sp,
  genes= genes.sel.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=10
);

# normalize the matrix by cell coverage
x.sp = scaleCountMatrix(
  obj=x.sp,
  cov=SnapATAC::rowSums(x.sp, mat="bmat"),
  mat="gmat",
  method = "RPM"
);

# plot enrichment for marker genes
# fixme:
path_marker_genes_pdf <- "marker-genes.pdf"
plotGene(
  obj=x.sp,
  gene.names=marker.genes,
  viz.method="tsne",
  point.size=0.3,
  point.color="red",
  point.shape=19,
  background.point=TRUE,
  background.point.color="grey",
  background.point.alpha=0.3,
  background.point.size=0.1,
  background.point.shape=19,
  low.value=0.0,
  high.value=0.95,
  down.sample=5000,
  seed.use=10,
  plot.nrow=4,
  plot.ncol=4,
  pdf.file.name=path_marker_genes_pdf,
  pdf.height=7,
  pdf.width=7
);

# -----------------------------------------------------
# Step 15. Heretical clustering of the clusters (SnapATAC)
print("Step 15");

# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})

# cluster using 1-cor as distance
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
plot(hc, hang=-1, xlab="");

# -----------------------------------------------------
# Step 16. Gene-body based annotation for excitatory neurons
print("Step 16");

idx = which(x.sp@cluster %in% c(3, 9, 27, 14, 8, 17, 22, 24, 15, 4, 12));
length(idx) # 2449 56% of total population
marker.genes = c(
  "Cux2", "Rorb", "Deptor",
  "Vat1l", "Sulf1", "Tle4",
  "Foxp2", "Tshz2", "Grik3"
);
x.exc.sp = x.sp[idx,];
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
x.exc.sp = createGmat(
  obj=x.exc.sp,
  genes=genes.sel.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=10
);

# normalize the matrix by cell coverage
x.exc.sp = scaleCountMatrix(
  obj=x.exc.sp,
  cov=SnapATAC::rowSums(x.exc.sp, mat="bmat"),
  mat="gmat",
  method = "RPM"
);

# fixme:
path_marker_genes_pdf <- "marker-genes.pdf"
plotGene(
  obj=x.exc.sp,
  gene.names=marker.genes,
  viz.method="tsne",
  point.size=0.2,
  point.color="red",
  point.shape=19,
  background.point=TRUE,
  background.point.color="grey",
  background.point.alpha=0.3,
  background.point.size=0.1,
  background.point.shape=19,
  low.value=0.0,
  high.value=1.0,
  down.sample=5000,
  seed.use=10,
  plot.nrow=3,
  plot.ncol=3,
  pdf.file.name=NULL,
  pdf.height=7,
  pdf.width=7
);

# -----------------------------------------------------
# Step 17. Change the cluster label to cell type
print("Step 17");

library(plyr);
current.cluster.ids <- 1:27;
new.cluster.ids  <- c(
  "Other", "Ogc.a", "Exc.a", "Exc.b", "Gaba.a",
  "Opc", "Gaba.b", "Exc.c", "Exc.d", "Ogc.b",
  "Asc.b", "Exc.e", "Gaba.c", "Exc.f", "Exc.g",
  "Gaba.d", "Exc.h", "Gaba.e", "Ogc.c", "Gaba.f",
  "Asc.c", "Exc.i", "Mgc.a", "Exc.j", "Mgc.b",
  "Asc.a", "Exc.k"
);
x.sp@cluster <- plyr::mapvalues(
  x = x.sp@cluster,
  from = current.cluster.ids,
  to = new.cluster.ids
);
plotViz(
  obj=x.sp,
  method="tsne",
  point.size=0.5,
  point.shape=19,
  point.alpha=0.8,
  point.color="cluster",
  text.add=TRUE,
  text.size=1.2,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  pdf.file.name=NULL,
  pdf.width=7,
  pdf.height=7
);

# -----------------------------------------------------
# Step 18. Identify cis-elements for each cluster seperately
print("Step 18")

# peaks_sst.df = runMACS(
#   obj=x.sp[which(x.sp@cluster=="Gaba.a"),],
#   output.prefix="atac_v1_adult_brain_fresh_5k.Sst",
#   path.to.snaptools="/usr/local/bin/snaptools",
#   path.to.macs="/usr/local/bin/macs2",
#   gsize="mm",
#   buffer.size=500,
#   num.cores=5,
#   macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
#   tmp.folder=tempdir()
# );

# nrow(peaks_sst.df);

peaks.gr = runMACSForAll(
  obj=x.sp,
  path.to.snaptools="/usr/local/bin/snaptools",
  path.to.macs="/usr/local/bin/macs2",
  output.prefix="macs2_out",
  num.cores=16,
  min.cells=100,
  gsize="mm",
  buffer.size=500,
  macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
  tmp.folder=tempdir()
);



print("DONE.");
