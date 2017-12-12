
## ------------------------------------------------------------------------
library('ggplot2')
require(gplots) 
library(plyr)
library(reshape2)

##-------------------------------------------------------------------------
load('cnts.RData')
load('tpm.RData')


## ------------------------------------------------------------------------

tot<-data.frame()
for(bname in names(cnts)){
  genes<-cnts[[bname]]
  dd<-cnt[,-c(1:2)]
  gt0cnt<-apply(dd,2,function(x){sum(x>0)})
  gt5cnt<-apply(dd,2,function(x){sum(x>5)})
  gt10cnt<-apply(dd,2,function(x){sum(x>10)})
  gt50cnt<-apply(dd,2,function(x){sum(x>50)})
  gt100cnt<-apply(dd,2,function(x){sum(x>100)})
  stats<-data.frame('pipeline'=rep(bname,ncol(dd)),gt0cnt,gt5cnt,gt10cnt,gt50cnt,gt100cnt)
  tot<-rbind(tot,stats)
}

##plot histogram
## ------------------------------------------------------------------------
for(cc in c(0,5,10,50,100)){
  cname<-paste('gt',cc,'cnt',sep='')
  p<-ggplot(data=tot,aes_string(x=cname,colour='pipeline'))+geom_density(alpha=0.5)
  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=10,face='bold'),legend.title = element_blank())
  p<-p+xlab('# of Detected Genes')+ylab('Density')+ggtitle(paste('# of genes with counts >',cc,sep=' ' ))
  p
}

## ------------------------------------------------------------------------
## load tpm

for(bname in names(cnts)){
  genes<-cnts[[bname]]
  tpm<-genes[,-c(1:2)]
  gt0tpm<-apply(tpm,2,function(x){sum(x>0)})
  gt5tpm<-apply(tpm,2,function(x){sum(x>5)})
  gt10tpm<-apply(tpm,2,function(x){sum(x>10)})
  gt50tpm<-apply(tpm,2,function(x){sum(x>50)})
  gt100tpm<-apply(tpm,2,function(x){sum(x>100)})
  cnttpm<-data.frame('pipeline'=rep(bname,ncol(tpm)),gt0tpm,gt5tpm,gt10tpm,gt50tpm,gt100tpm)
  tpmstats<-rbind(tpmstats,cnttpm)
}

## ------------------------------------------------------------------------
for(cc in c(0,5,10,50,100)){
  cname<-paste('gt',cc,'tpm',sep='')
  p<-ggplot(data=tpmstats,aes_string(x=cname,colour='pipeline'))+geom_density(alpha=0.5)
  p<-p+theme(legend.position="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=15,face='bold'),axis.title.y = element_text(size=15,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=15,face='bold'),axis.text.x = element_text(color='black',size=15,face='bold'))
  p<-p+theme(legend.text=element_text(size=10,face='bold'),legend.title = element_blank())
  p<-p+xlab('# of Detected Genes')+ylab('Density')+ggtitle(paste('# of genes with TPM >',cc,sep=' ' ))
  p
  
}

