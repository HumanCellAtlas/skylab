## ------------------------------------------------------------------------
library('ggplot2')
library('rtracklayer')
require(gplots) 
library(plyr)
library(reshape2)


options(stringsAsFactors = FALSE)
gtf_gencode <- readGFF("gencode.v27.primary_assembly.annotation.gtf", version=2L, tags = c("gene_name","gene_id", "transcript_id","gene_type"))
genes<-subset(gtf_gencode,gtf_gencode$type == "gene")
#load cnt and tpm data
load('cnts.RData')
load('tpm.RData')
load('metrics.RData')

## plot function

sumBiotype<-function(x,labelname,ratio){
  x.miRNA<-x$miRNA+x$misc_RNA+x$snoRNA+x$snRNA
  x.pseudo<-x$pseudogene+x$transcribed_processed_pseudogene+x$translated_processed_pseudogene+x$processed_pseudogene+x$transcribed_unprocessed_pseudogene+x$unprocessed_pseudogene
  x.protein_coding<-x$protein_coding
  x.antisense<-x$antisense_RNA
  x.lncRNA<-x$lincRNA
  if(ratio == 'ratio'){
    x.others<- 1- (x.pseudo+x.protein_coding+x.antisense+x.lncRNA+x.miRNA)
  }else{
    x.tot<-apply(x,1,sum)
    x.others<- x.tot- (x.pseudo+x.protein_coding+x.antisense+x.lncRNA+x.miRNA)
  }
  d1<-data.frame('Discordant Genes'=rep(labelname,nrow(x)),'misc_RNA'=x.miRNA,'protein_coding'=x.protein_coding,'pseudo'=x.pseudo,'linRNA'=x.lncRNA,'antisensse_RNA'=x.antisense,'others'=x.others)
  m1<-apply(d1[,-1],2,mean,na.rm=T)
  s1<-apply(d1[,-1],2,sd,na.rm=T)
  mdd<-data.frame('discordance'=rep(labelname,6),'mean'=m1,'sd'=s1,'biotype'=names(m1))
  return(mdd)
}
plotBarplot<-function(mdd,labels){
  
  p<- ggplot(mdd, aes(x=biotype, y=mean, fill=discordance)) + geom_bar(stat="identity", color="black", position=position_dodge()) 
  p<-p+geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) 
  p<-p+xlab(labels[['xlab']])+ylab(labels[['ylab']])+ggtitle(labels[['title']])
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(legend.position="top")
  p<-p+theme(axis.title.x = element_text(color='black',size=18,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=18,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=15,face='bold'),legend.title =element_text(size=15,face='bold') )+scale_fill_discrete(name = labels[['legend']])
  return(p)
}


## ------------------------------------------------------------------------
metadata$AvgSpotLen_l<-paste(metadata$AvgSpotLen_l,'bp',sep='')
nalist<-which(metadata$patient_id == "")
metadata$patient_id_s[nalist]<-metadata$cell_line_s[nalist]
biotypes<-names(table(genes$gene_type))

## ------------------------------------------------------------------------
## compare star.mult.tpm vs hisat2.tpm
db<-list('star'=star.tpm,'hisat2'=hisat2.mult.tpm)
## overlapping samples 
sraIDs<-intersect(colnames(star.mult.tpm[,-c(1:2)]),colnames(hisat2.mult.tpm[,-c(1:2)]))
stats.lab<-c()
stats<-c()
for(n in c(5,10,50,100)){
  print(n)
  star<-db[['star']]
  biolist<-list()
  hisat2<-db[['hisat2']]
  l1<-match(as.matrix(sraIDs),colnames(star))
  l2<-match(as.matrix(sraIDs),colnames(hisat2))
  star.dd<-star[,l1]
  hisat2.dd<-hisat2[,l2]
  star.lab<-star[,c(1)]
  hisat2.lab<-hisat2[,c(1)]
  for(ss in sraIDs){
    ##print(ss)
    x1<-star.dd[,ss]
    x2<-hisat2.dd[,ss]
    g1<-subset(genes$gene_type,genes$gene_id %in% star.lab[x1 >=n & x2 <n] ) ## high in star but low in hisat2
    g2<-subset(genes$gene_type,genes$gene_id %in% hisat2.lab[x1 <n & x2 >=n]) ## high in hisat2 but low in star
    t1<-table(factor(g1,levels=biotypes)) ## high in star but low in hisat2
    t2<-table(factor(g2,levels=biotypes)) ## high in hisat2 but low in star
    if(length(stats)==0){
        stats<-rbind(t1,t2)
        stats.lab<-rbind(c(ss,n,'hisat2'),c(ss,n,'star'))
    }else{
      stats<-rbind(stats,t1,t2)
      stats.lab<-rbind(stats.lab,rbind(c(ss,n,'hisat2'),c(ss,n,'star')))
    }
    
    }
}
colnames(stats.lab)<-c('sra','cutoff','aligner')
stats.lab<-data.frame(stats.lab)
outdata<-cbind(stats.lab,as.data.frame(stats))
## ------------------------------------------------------------------------
##plot outdata
star.drop<-subset(outdata,outdata$aligner == 'star')
hisat2.drop<-subset(outdata,outdata$aligner == "hisat2")
star.tot<-apply(star.drop[,-c(1:3)],1,sum)
hisat2.tot<-apply(hisat2.drop[,-c(1:3)],1,sum)
rmlist<-which(star.tot==0)
if(length(rmlist)>0){
  star.rr<-star.drop[-c(rmlist),-c(1:3)]/star.tot[-c(rmlist)]
  hisat2.rr<-hisat2.drop[-c(rmlist),-c(1:3)]/hisat2.tot[-c(rmlist)]
  star.dis<-star.drop[-c(rmlist),-c(1:3)]
  hisat2.dis<-hisat2.drop[-c(rmlist),-c(1:3)]
}else{
  star.rr<-star.drop[,-c(1:3)]/star.tot
  hisat2.rr<-hisat2.drop[,-c(1:3)]/hisat2.tot
  star.dis<-star.drop[,-c(1:3)]
  hisat2.dis<-hisat2.drop[,-c(1:3)]
}
## different threshold 
for(n  in c(5,10,50,100))
{
  x<-subset(star.rr,star.drop$cutoff ==n)
  y<-subset(hisat2.rr,hisat2.drop$cutoff ==n)
  m1<-sumBiotype(x,'I','ratio')
  m2<-sumBiotype(y,'II','ratio')
  mdd<-rbind(m1,m2)
  labels<-list('xlab'= 'Biotype',
               'ylab'= "Avg % of discordant genes in each biotype groups",
               'title'= paste('TPM <',n,sep=''),
               'legend'='Discordance Group')
  plotBarplot(mdd,labels)
  ## total number
  x<-subset(star.dis,star.drop$cutoff ==n)
  y<-subset(hisat2.dis,hisat2.drop$cutoff ==n)
  m1<-sumBiotype(x,'I','total')
  m2<-sumBiotype(y,'II','total')
  mdd<-rbind(m1,m2)
  labels<-list('xlab'= 'Biotype',
               'ylab'= "Avg # of discordant genes in each biotype groups",
               'title'= paste('TPM <',n,sep=''),
               'legend'='Discordance Group')
  plotBarplot(mdd,labels)
}


## ------------------------------------------------------------------------
## extreme case such as pipeline A have very low expression but pipeline B have very high expression
db<-list('star'=star.mult.tpm,'hisat2'=hisat2.mult.tpm)
stats.lab<-c()
stats<-c()
n1<-5
n2<-20
star<-db[['star']]
biolist<-list()
hisat2<-db[['hisat2']]
l1<-match(as.matrix(sraIDs),colnames(star))
l2<-match(as.matrix(sraIDs),colnames(hisat2))
star.dd<-star[,l1]
hisat2.dd<-hisat2[,l2]
star.lab<-star[,c(1)]
hisat2.lab<-hisat2[,c(1)]
## pairwise and sample to sample 
for(ss in sraIDs){
  ##print(ss)
  x1<-star.dd[,ss]
  x2<-hisat2.dd[,ss]
  g1<-subset(genes$gene_type,genes$gene_id %in% star.lab[x2 >=n2 & x1 <n1] ) ## high in hisat2 but low in star
  g2<-subset(genes$gene_type,genes$gene_id %in% hisat2.lab[x2 <n1 & x1 >=n2]) ## high in star but low in hisat2
  t1<-table(factor(g1,levels=biotypes)) ## high in hsiat2 but low in star
  t2<-table(factor(g2,levels=biotypes)) ## high in hisat2 but low in hisat2
  if(length(stats)==0){
    stats<-rbind(t1,t2)
    stats.lab<-rbind(c(ss,'star'),c(ss,'hisat2'))
  }else{
    stats<-rbind(stats,t1,t2)
    stats.lab<-rbind(stats.lab,rbind(c(ss,'star'),c(ss,'hisat2')))
  }
}
colnames(stats.lab)<-c('sra','aligner')
stats.lab<-data.frame(stats.lab)
outdata<-cbind(stats.lab,as.data.frame(stats))

## ------------------------------------------------------------------------
##
star.drop<-subset(outdata,outdata$aligner == 'star')
hisat2.drop<-subset(outdata,outdata$aligner == "hisat2")
star.tot<-apply(star.drop[,-c(1:2)],1,sum)
hisat2.tot<-apply(hisat2.drop[,-c(1:2)],1,sum)
rmlist<-which(star.tot==0)
if(length(rmlist)>0){
  star.rr<-star.drop[-c(rmlist),-c(1:2)]/star.tot[-c(rmlist)]
  hisat2.rr<-hisat2.drop[-c(rmlist),-c(1:2)]/hisat2.tot[-c(rmlist)]
  star.dis<-star.drop[-c(rmlist),-c(1:3)]
  hisat2.dis<-hisat2.drop[-c(rmlist),-c(1:2)]
}else{
  star.rr<-star.drop[,-c(1:2)]/star.tot
  hisat2.rr<-hisat2.drop[,-c(1:2)]/hisat2.tot
  star.dis<-star.drop[,-c(1:2)]
  hisat2.dis<-hisat2.drop[,-c(1:2)]
}
##barplot of ratio, with error bar
x<-star.rr
y<-hisat2.rr
## group biotypes into protein_coding, pseudo gene, miscRNA, lincRNA groups.
m1<-sumBiotype(x,'I','ratio')
m2<-sumBiotype(y,'II','ratio')
mdd<-rbind(m1,m2)
labels<-list('xlab'= 'Biotype',
             'ylab'= "Avg % of discordant genes in each biotype groups",
             'title'= paste('TPM <5',sep=''),
             'legend'='Discordance Group')
plotBarplot(mdd,labels)
## total number
x<-subset(star.dis,star.drop$cutoff ==n)
y<-subset(hisat2.dis,hisat2.drop$cutoff ==n)
m1<-sumBiotype(x,'I','total')
m2<-sumBiotype(y,'II','total')
mdd<-rbind(m1,m2)
labels<-list('xlab'= 'Biotype',
             'ylab'= "Avg # of discordant genes in each biotype groups",
             'title'= paste('TPM <5 ',sep=''),
             'legend'='Discordance Group')
plotBarplot(mdd,labels)




