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
library(corrplot)
library(ggrepel)
library(optparse)
set.seed(42)
# take top N high variance genes
TopVarGenes<-function(cnts,topvar){
  cnt.dd<-cnts[,-c(1:2)]
  cnt.label<-cnts[,c(1:2)]
  sd.list<-order(apply(cnt.dd,1,sd),decreasing = T)
  kept<-cnts[sd.list[1:topvar],]
  return(kept)
}
# if p<0.05, return 1 value.
convert2SigLevels<-function(P){
  x<-floor(-1*log(P,base=10))
  x[x>1]<-1
  return(x)
}
# correlation between qc metrics and data matrics. 
# first run PCA on data matrix and top top N PCs 
# Then calculate correlation and correlation test between
# QC metrics and N PCs 
# Visualize this correlation and test results in corrplot
# Y-axis rank by total number of significant correlation test
# QCmetrics with High rank(on the top) indicates these metrics have
# high impact on data matrix
CorrQCvsPCs<-function(qc_mets, cnts,npcs,output_name){
  cnt.pca<-prcomp(cnts[,-c(1:2)],scale=T)
  pcs<-cnt.pca$rotation[,1:npcs]
  nc1<-ncol(qc_mets)
  dt<-merge(qc_mets,pcs,by=0)
  res<-cor.mtest(dt[,-1])
  pmat<-res$p
  colnames(pmat)<-colnames(dt[,-1])
  rownames(pmat)<-colnames(dt[,-1])
  M<-cor(dt[,-1])
  nc2<-ncol(M)
  M.sub<-M[(nc1+1):nc2,1:nc1]
  p.sub<-pmat[(nc1+1):nc2,1:nc1]
  rank.p<-order(colSums(p.sub<0.05,na.rm = T),decreasing = T)
  png(paste(output_name,'_corrplot_top_',npcs,'_vs_qc_mets.png',sep=''),width=5000,height=5000)
  corrplot(t(M.sub[,rank.p]), p.mat = t(p.sub[,rank.p]),is.corr=TRUE,insig = "label_sig",
           sig.level = c(.001, .01, .05), 
           pch.cex = 1.5, 
           pch.col = "white",
           na.label = " ",
           tl.cex =5,
           cl.cex=5,
           tl.col="black",
           method = "color")
  dev.off()
  ## all metrics and PCs
  rank.p<-order(rowSums(pmat<0.05,na.rm=T),decreasing = T)
  png(paste(output_name,'_corrplot_top_',npcs,'_vs_qc_mets_all.png',sep=''),width=8000,height=8000)
  corrplot(t(M[rank.p,]), p.mat = t(pmat[rank.p,]),is.corr=TRUE,insig = "label_sig",
           sig.level = c(.001, .01, .05), 
           pch.cex = 1.5, 
           pch.col = "white",
           na.label = " ",
           tl.cex =5,
           cl.cex=5,
           tl.cor='black',
           method = "color")
  dev.off()
  return(p.sub)
}
addTheme<-function(p){
  p<-p+theme(legend.position ="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=18,face='bold'),axis.title.y = element_text(size=18,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=18,face='bold'),axis.text.x = element_text(color='black',size=18,face='bold'))
  p<-p+theme(legend.text=element_text(size=15,face='bold'),legend.title = element_blank())
  return(p)
}

# params
# python style input args
option_list <- list(
  make_option("--bdatafile", type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option("--udatafile", type="character", default=NULL, 
              help="updated dataset file name", metavar="character"),
  make_option("--bmetrics", type="character", default=NULL, 
              help=" base metricsfile name", metavar="character"),
  make_option("--umetrics", type="character", default=NULL, 
              help=" updated metricsfile name", metavar="character"),
  make_option("--npcs", type="integer", default=10, 
              help=" number of PCs to collect", metavar="number"),
  make_option("--nvar", type="integer", default=1000, 
              help=" number of top variantion genes", metavar="number"),
  make_option("--out", type="character", default="out", 
              help="output file name [default= %default]", metavar="character")
)
opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)


## load files
met1<-read.csv(opt$bmetrics,row.names = 1)
met2<-read.csv(opt$umetrics,row.names = 1)
cnt1<-read.csv(opt$bdatafile)
cnt2<-read.csv(opt$udatafile)
nvar<-opt$nvar
npcs<-opt$npcs
output_name<-opt$out
## pre-PCA, take top 5000 variation genes
cnt1.kept<-TopVarGenes(cnt1,nvar)
cnt2.kept<-TopVarGenes(cnt2,nvar)

## run PCA analysis and then calculate correlation between QC and PCs
pmat1<-CorrQCvsPCs(t(met1), cnt1.kept,npcs,paste(output_name,'_base',sep=''))
write.csv(pmat1,file=paste(output_name,'_base_cor_test.csv',sep=''),
          sep=',',quote=F,row.names=T,col.names=T)
pmat2<-CorrQCvsPCs(t(met2), cnt2.kept,npcs,paste(output_name,'_updated',sep=''))
write.csv(pmat2,file=paste(output_name,'_updated_cor_test.csv',sep=''),
          sep=',',quote=F,row.names=T,col.names=T)
## plot two pipeline QC rank based on correlation to matrix
sum.sig1<-colSums(pmat1<0.05,na.rm = T)
sum.sig2<-colSums(pmat2<0.05,na.rm = T)
df<-merge(sum.sig1,sum.sig2,by=0)
colnames(df)<-c('QC','Base','Updated')
pl<-ggplot(df) +
  geom_point(aes(Base, Updated), color = 'red') +
  geom_text_repel(aes(Base, Updated, label = QC)) +
  theme_classic(base_size = 16)+
  geom_abline(slope = 1, intercept = 0,color="blue")
pl<-pl+ggtitle(paste("# of significant correlation between QC metrics amd Expression PCs"))
pl<-addTheme(pl)
ggsave(pl,file=paste(output_name,'_sum_nsig_p_xyplot.png',sep=''),
       width=20,height=20,limitsize = FALSE)
