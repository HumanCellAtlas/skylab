## ------------------------------------------------------------------------
library(ggplot2)
library(rjson)
library(dplyr)
library(MASS)
library('rtracklayer')
require(gplots) 
library(plyr)
library(reshape2)
library(VennDiagram)

## function


## ------------------------------------------------------------------------
load('~/metrics.RData')
ls()

## ------------------------------------------------------------------------
## alignment
col1<-paste('star',colnames(star.aln),sep='.')
col2<-paste('hisat2',colnames(hisat2.aln),sep='.')
colnames(star.aln)<-col1
colnames(hisat2.aln)<-col2
dd<-merge(star.aln,hisat2.aln,by.x='star.sraID',by.y='hisat2.sraID')
colnames(dd)[1]<-'Run_s'
## matching sraID
mlist<-match(dd$Run_s,metadata$Run_s)
dd<-data.frame(dd,'rlen'=factor(metadata$AvgSpotLen_l[mlist]))
for(m in c('TOTAL_READS','PF_READS','PF_READS_ALIGNED','PCT_PF_READS_ALIGNED','PF_ALIGNED_BASES','PF_MISMATCH_RATE','PF_HQ_ERROR_RATE','PF_READS_IMPROPER_PAIRS','PCT_PF_READS_IMPROPER_PAIRS')){
  i<-paste('star',m,sep='.')
  j<-paste('hisat2',m,sep='.')
  xmax<-max(dd[,i])*1.1
  ymax<-max(dd[,j])*1.1
  amax<-max(xmax,ymax)
  p1<-ggplot(dd,aes_string(x=i,y=j,color='rlen'))+geom_point()+xlim(0,amax)+ylim(0,amax)
  p1<-p1+geom_abline(intercept=0, slope=1,color='blue')
  p1
  xx<-data.frame('Run_s'=dd$Run_s,'hisat2'=dd[,j],'star'=dd[,i],'rlen'=dd$rlen)
  xx.m<-melt(xx)
  p<-ggplot(data=xx.m,aes(x=value,colour=paste(variable,paste(rlen,'bp',sep='')),linetype=paste(rlen,'bp',sep='')))
  if(length(grep('PCT' , m))>1){
    p<-p+geom_density(alpha=0.2)+xlim(c(0,1))
  }else{
    p<-p+geom_density(alpha=0.2)
  }

  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=12,face='bold'),legend.title = element_blank())
  p<-p+xlab('values')+ylab(paste('Density of ',m,sep=''))+ggtitle(paste('HISAT2 vs STAR',m,sep=' '))
  p
}

## ------------------------------------------------------------------------
## rna metrics
col1<-paste('star',colnames(star.rna),sep='.')
col2<-paste('hisat2',colnames(hisat2.rna),sep='.')
colnames(star.rna)<-col1
colnames(hisat2.rna)<-col2
dd<-merge(star.rna,hisat2.rna,by.x='star.sraID',by.y='hisat2.sraID')
dd<-data.frame(dd,'rlen'=factor(metadata$AvgSpotLen_l[mlist]))
colnames(dd)[1]<-'Run_s'

for(m in c('PCT_CODING_BASES','PCT_UTR_BASES','PCT_USABLE_BASES','PCT_INTERGENIC_BASES','PCT_INTRONIC_BASES','MEDIAN_5PRIME_BIAS','MEDIAN_3PRIME_BIAS','MEDIAN_5PRIME_TO_3PRIME_BIAS')){
  i<-paste('star',m,sep='.')
  j<-paste('hisat2',m,sep='.')
  xmax<-max(dd[,i])*1.1
  ymax<-max(dd[,j])*1.1
  amax<-max(xmax,ymax)
  p1<-ggplot(dd,aes_string(x=i,y=j,color='rlen'))+geom_point()+xlim(0,amax)+ylim(0,amax)
  p1<-p1+geom_abline(intercept=0, slope=1,color='blue')
  p
  ##hist
  xx<-data.frame('Run_s'=dd$Run_s,'hisat2'=dd[,j],'star'=dd[,i],'rlen'=dd$rlen)
  xx.m<-melt(xx)
  p<-ggplot(data=xx.m,aes(x=value,colour=paste(variable,paste(rlen,'bp',sep='')),linetype=paste(rlen,'bp',sep='')))
  if(length(grep('PCT' , m))>1){
    p<-p+geom_density(alpha=0.2)+xlim(c(0,1))
  }else{
    p<-p+geom_density(alpha=0.2)
  }
  
  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=12,face='bold'),legend.title = element_blank())
  p<-p+xlab('values')+ylab(paste('Density of ',m,sep=''))+ggtitle(paste('HISAT2 vs STAR',m,sep=' '))
  p
}

## ------------------------------------------------------------------------
## duplication
col1<-paste('star',colnames(star.dup),sep='.')
col2<-paste('hisat2',colnames(hisat2.dup),sep='.')
colnames(star.dup)<-col1
colnames(hisat2.dup)<-col2
dd<-merge(star.dup,hisat2.dup,by.x='star.sraID',by.y='hisat2.sraID')
dd<-data.frame(dd,'rlen'=factor(metadata$AvgSpotLen_l[mlist]))
colnames(dd)[1]<-'Run_s'

for(m in c('PERCENT_DUPLICATION','ESTIMATED_LIBRARY_SIZE','SECONDARY_OR_SUPPLEMENTARY_RDS','UNMAPPED_READS','UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES')){
  i<-paste('star',m,sep='.')
  j<-paste('hisat2',m,sep='.')
  xmax<-max(dd[,i])*1.1
  ymax<-max(dd[,j])*1.1
  amax<-max(xmax,ymax)
  p1<-ggplot(dd,aes_string(x=i,y=j,color='rlen'))+geom_point()+xlim(0,amax)+ylim(0,amax)
  p1<-p1+geom_abline(intercept=0, slope=1,color='blue')
  p1
  xx<-data.frame('Run_s'=dd$Run_s,'hisat2'=dd[,j],'star'=dd[,i],'rlen'=dd$rlen)
  xx.m<-melt(xx)
  p<-ggplot(data=xx.m,aes(x=value,colour=paste(variable,paste(rlen,'bp',sep='')),linetype=paste(rlen,'bp',sep='')))
  if(length(grep('PCT' , m))>1){
    p<-p+geom_density(alpha=0.2)+xlim(c(0,1))
  }else{
    p<-p+geom_density(alpha=0.2)
  }

  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=12,face='bold'),legend.title = element_blank())
  p<-p+xlab('values')+ylab(paste('Density of ',m,sep=''))+ggtitle(paste('HISAT2 vs STAR',m,sep=' '))
  p
}

## ------------------------------------------------------------------------
## insertion size
col1<-paste('star',colnames(star.int),sep='.')
col2<-paste('hisat2',colnames(hisat2.int),sep='.')
colnames(star.int)<-col1
colnames(hisat2.int)<-col2
dd<-merge(star.int,hisat2.int,by.x='star.sraID',by.y='hisat2.sraID')
colnames(dd)[1] <-'Run_s'
mlist<-match(dd$Run_s,metadata$Run_s)
dd<-data.frame(dd,'rlen'=factor(metadata$AvgSpotLen_l[mlist]))
colnames(dd)[1]<-'Run_s'
for(m in c('MEAN_INSERT_SIZE','MEDIAN_INSERT_SIZE')){
  i<-paste('star',m,sep='.')
  j<-paste('hisat2',m,sep='.')
  xmax<-max(dd[,i])*1.1
  ymax<-max(dd[,j])*1.1
  amax<-max(xmax,ymax)
  p1<-ggplot(dd,aes_string(x=i,y=j,color='rlen'))+geom_point()+xlim(0,500)+ylim(0,500)
  p1<-p1+geom_abline(intercept=0, slope=1,color='blue')
  p1
  xx<-data.frame('Run_s'=dd$Run_s,'hisat2'=dd[,j],'star'=dd[,i],'rlen'=dd$rlen)
  xx.m<-melt(xx)
  p<-ggplot(data=xx.m,aes(x=value,colour=paste(variable,paste(rlen,'bp',sep='')),linetype=paste(rlen,'bp',sep='')))
  if(length(grep('PCT' , m))>1){
    p<-p+geom_density(alpha=0.2)+xlim(c(0,1))
  }else{
    p<-p+geom_density(alpha=0.2)
  }
  
  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=12,face='bold'),legend.title = element_blank())
  p<-p+xlab('values')+ylab(paste('Density of ',m,sep=''))+ggtitle(paste('HISAT2 vs STAR',m,sep=' '))
  p
}

