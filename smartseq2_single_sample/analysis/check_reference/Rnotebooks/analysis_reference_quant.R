## ------------------------------------------------------------------------
library('ggplot2')
library('rtracklayer')
require(gplots) 
library(plyr)
library(reshape2)
library(plotly)
library(VennDiagram)
system.time(gtf_gencode_comp <- readGFF("~/Documents/HCA/reference/refs/gencode.v27.chr_patch_hapl_scaff.annotation.gtf", version=2L, tags = c("gene_name","gene_id", "transcript_id","gene_type")))
system.time(gtf_gencode_basic <- readGFF("~/Documents/HCA/reference/refs/gencode.v27.chr_patch_hapl_scaff.basic.annotation.gtf", version=2L, tags = c("gene_name","gene_id", "transcript_id","gene_type")))
system.time(gtf_ensembl <- readGFF("~/Documents/HCA/reference/refs/Homo_sapiens.GRCh38.90.gtf", version=2L,tags = c("gene_name","gene_id", "transcript_id","gene_biotype")))
system.time(gtf_refseq <- readGFF("~/Documents/HCA/reference/refs/ncbi-genomes-2017-10-05/GCF_000001405.37_GRCh38.p11_genomic.gff.gz", tags = c("ID", "Name","gbkey","gene","gene_biotype","Parent")))

## ------------------------------------------------------------------------
tags<-c('GeneID','Chr','Start','End','Strand','Length','Counts')
head(gtf_ensembl)

## ------------------------------------------------------------------------
g4<-unique(gtf_gencode_basic$gene_name)
g3<-unique(gtf_gencode_comp$gene_name)
g2<-unique(gtf_ensembl$gene_name)
g1<-unique(na.omit(gtf_refseq$gene))
## two packages can plot venn diagram
##venn.plot<-venn.diagram(x = list('RefSeq'=g1,'Ensembl'=g2,'Gencode'=g3,'GencodeBasic'=g4), imagetype='png',filename = "~/Documents/HCA/reference/refs/gtf_gene_id_venn.png",col = "transparent", fill = c("cornflowerblue","green","yellow","darkorchid1"),label.col = c("orange", "white", "darkorchid4", "white", "white",  "white",    "white", "white", "darkblue", "white", "white", "white", "white",  "darkgreen", "white"),alpha = 0.50,cex = 1.5, fontfamily = "serif", fontface = "bold",cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 270,margin = 0.2)
venn(list('RefSeq'=g1,'Ensembl'=g2,'Gencode'=g3,'GencodeBasic'=g4),show.plot = T)

## ------------------------------------------------------------------------
g1.gene<-subset(gtf_refseq,gtf_refseq$type == "gene")
g1.sub<-subset(g1.gene, !(g1.gene$gene %in% g3))
unq.gene.refseq<-data.frame(table(g1.sub$gene_biotype))
piepercent<- round(100*unq.gene.refseq$Freq/sum(unq.gene.refseq$Freq), 1)
## only show >5%
x<-subset(unq.gene.refseq, piepercent >5)
y<-subset(unq.gene.refseq,piepercent <5)
x<-rbind(x,data.frame('Var1'='others','Freq'=sum(y$Freq)))
piepercent<- round(100*x$Freq/sum(x$Freq), 1)
pie(x$Freq,piepercent,main='Unq gene in refSeq',col=rainbow(nrow(x)))
legend("topright",legend=x$Var1 ,cex = 0.8,fill = rainbow(nrow(x)))

## ------------------------------------------------------------------------
ensembl.counts<-read.delim('/Users/jishuxu/Documents/HCA/reference/counts/SRR1294900_25_GRCh38_Ensembl.gene.unq.counts.txt',sep='\t',header=T,skip=1)
colnames(ensembl.counts)<-tags
mlist<-match(ensembl.counts$GeneID,gtf_ensembl$gene_id)
geneName<-gtf_ensembl$gene_name[mlist]
btype<-gtf_ensembl$gene_biotype[mlist]
c1<-data.frame("geneName"=geneName,'Ensembl.Counts'=log(ensembl.counts$Counts+1,base=2),'Ensembl.type'=gtf_ensembl$gene_biotype[mlist])
head(c1)

## ------------------------------------------------------------------------
gene.detects<-list()
t1<-sort(table(c1[c1$Ensembl.Counts>0,3]),decreasing=T)
t2<-sort(table(c1[c1$Ensembl.Counts>5,3]),decreasing=T)
gene.detects[['ensembl']]<-cbind('detected.genes'=t1[1:10],'good.genes'=t2[1:10])

## ------------------------------------------------------------------------
refseq.counts<-read.delim('/Users/jishuxu/Documents/HCA/reference/counts/SRR1294900_25_GRCh38_RefSeq.gene.unq.counts.txt',sep='\t',header=T,skip=1)
colnames(refseq.counts)<-tags
mlist<-match(refseq.counts$GeneID,gtf_refseq$ID)
geneName<-gtf_refseq$Name[mlist]
c2<-data.frame("geneName"=geneName,"RefSeq.Counts"=log(refseq.counts$Counts+1,base=2),'RefSeq.type'=gtf_refseq$gene_biotype[mlist])
head(c2)

## ------------------------------------------------------------------------
t1<-sort(table(c2[c2$RefSeq.Counts>0,3]),decreasing=T)
t2<-sort(table(c2[c2$RefSeq.Counts>5,3]),decreasing=T)
gene.detects[['refseq']]<-cbind('detected.genes'=t1[1:10],'good.genes'=t2[1:10])

## ------------------------------------------------------------------------
gencode.counts.comp<-read.delim('/Users/jishuxu/Documents/HCA/reference/counts/SRR1294900_25_GRCh38_GencodeV27.gene.unq.counts.txt',sep='\t',header=T,skip=1)
colnames(gencode.counts.comp)<-tags
mlist<-match(gencode.counts.comp$GeneID,gtf_gencode_comp$gene_id)
geneName<-gtf_gencode_comp$gene_name[mlist]
c3<-data.frame("geneName"=geneName,"GencodeComp.Counts"=log(gencode.counts.comp$Counts+1,base=2),'Gencode.type'=gtf_gencode_comp$gene_type[mlist])
dim(c3)

## ------------------------------------------------------------------------
gencode.counts.basic<-read.delim('/Users/jishuxu/Documents/HCA/reference/counts/SRR1294900_25_GRCh38_GencodeV27_basic.gene.unq.counts.txt',sep='\t',header=T,skip=1)
colnames(gencode.counts.basic)<-tags
mlist<-match(gencode.counts.basic$GeneID,gtf_gencode_basic$gene_id)
geneName<-gtf_gencode_basic$gene_name[mlist]
c4<-data.frame("geneName"=geneName,"GencodeBasic.Counts"=log(gencode.counts.basic$Counts+1,base=2),'Gencode.type'=gtf_gencode_basic$gene_type[mlist])
dim(c4)

## ------------------------------------------------------------------------
t1<-sort(table(c3[c3$GencodeComp.Counts>0,3]),decreasing=T)
t2<-sort(table(c3[c3$GencodeComp.Counts >5,3]),decreasing=T)
gene.detects[['GencodeComp']]<-cbind('detected.genes'=t1[1:10],'good.genes'=t2[1:10])

## ------------------------------------------------------------------------
t1<-sort(table(c4[c4$GencodeBasic.Counts>0,3]),decreasing=T)
t2<-sort(table(c4[c4$GencodeBasic.Counts >5,3]),decreasing=T)
gene.detects[["GencodeBasic"]]<-cbind('detected.genes'=t1[1:10],'good.genes'=t2[1:10])

## ------------------------------------------------------------------------
##top 3 biotype
library(reshape2)
bios<-c('protein_coding','processed_pseudogene','pseudogene','lincRNA','lncRNA')
g.ensembl<-subset(gene.detects$ensembl,rownames(gene.detects$ensembl) %in% bios)
g.refseq<-subset(gene.detects$refseq,rownames(gene.detects$refseq) %in% bios)
g.gencode<-subset(gene.detects$GencodeComp,rownames(gene.detects$GencodeComp) %in% bios)
g.gencodeBasic<-subset(gene.detects$GencodeBasic,rownames(gene.detects$GencodeBasic) %in% bios)
detects<-data.frame('biotype'=rep(c('protein_coding','pseudogene','lncRNA'),4),rbind(g.ensembl,g.refseq,g.gencode,g.gencodeBasic),'bundle'=rep(c('Ensembl','RefSeq','Gencode','GencodeV27'),c(3,3,3,3)))
m.detects<-melt(detects)
##detected >0X
ggplot(data=m.detects,aes(x=biotype,y=value,fill=bundle))+geom_bar(stat="identity",position=position_dodge())+facet_grid(variable~.,scales="free_y")+geom_text(aes(biotype, value, label = value),position = position_dodge(width = 1)) 


## ------------------------------------------------------------------------
ensembl.refseq<-merge(c1,c2,by=c('geneName'),all=F,all.x=F,all.y=F)
dim(ensembl.refseq)

## ------------------------------------------------------------------------
subdata<-subset(ensembl.refseq,Ensembl.Counts>0 & RefSeq.Counts >0)
fit<-lm(RefSeq.Counts~Ensembl.Counts, data=subdata)
summary(fit)

## ------------------------------------------------------------------------
pred<-predict(fit,data=ensembl.refseq,level=0)
conf_interval_3 <- predict(fit, newdata=data.frame(Ensembl.Counts=100), interval="prediction",level = 0.95)
conf_interval_3

## ------------------------------------------------------------------------
conf_interval <- predict(fit, newdata=ensembl.refseq, interval="prediction",level = 0.95)

## ------------------------------------------------------------------------
##this function just label model in plot
equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(coef(x)[2], digits = 2),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));                 
}
##label with color
type<-rep('normal',nrow(ensembl.refseq))
type[ensembl.refseq$Ensembl.Counts==0 & ensembl.refseq$RefSeq.Counts>0]<-'dropout.ensembl'
type[ensembl.refseq$Ensembl.Counts>0 & ensembl.refseq$RefSeq.Counts==0]<-'dropout.refseq'
type[(ensembl.refseq$RefSeq.Counts > conf_interval[,3]|ensembl.refseq$RefSeq.Counts<conf_interval[,2])&(ensembl.refseq$Ensembl.Counts>0 & ensembl.refseq$RefSeq.Counts>0)]<-'overdispersion' 
ensembl.refseq$type<-type
table(ensembl.refseq$type)
## label prefiction model in plot
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(fit)[1], digits = 2), b = format(coef(fit)[2], digits = 2), r2 = format(summary(fit)$r.squared, digits = 3)))
## scatter plot
p1 <- ggplot(ensembl.refseq,aes(x=Ensembl.Counts,y=RefSeq.Counts),color="grey")
p1<-p1+geom_ribbon(aes(ymin=conf_interval[,2],ymax=conf_interval[,3]),alpha=0.2,fill="red")
p1<-p1+geom_point(data=ensembl.refseq,aes(x=Ensembl.Counts,y=RefSeq.Counts,color=type))
p1<-p1+annotate("text", x = 4.1, y = 12, label = equation(fit), parse = TRUE)
p1<-p1+xlab("log2(Ensembl.Counts+1)")+ylab("log2(RefSeq+1)")
p1

## ------------------------------------------------------------------------
dropout.ensembl<-subset(ensembl.refseq,Ensembl.Counts==0 & RefSeq.Counts>0)
dim(dropout.ensembl)
## barplot
x<-data.frame(sort(table(dropout.ensembl$Ensembl.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout genes in Ensembl')

## ------------------------------------------------------------------------
library(GeneOverlap)
library(MSigDB)
## create function to do overlapping test
gsOverlap<-function(gs_cureated,gs){
  ##msigdb_c2<-MSigDB[["C2_CURATED"]]
  pvals<-c()
  for(i in names(gs_curated)){
    go.obj <- newGeneOverlap(gs_curated[[i]], gs, genome.size=60000)
    ##fisher exact
    go.obj <- testGeneOverlap(go.obj)
    overlap<-length(getIntersection(go.obj))
    PopSize<-60000
    f.pval<-getPval(go.obj)
    list1<-length(gs)
    list2<-length(gs_curated[[i]])
    ##hyperGT
    if(overlap>1){
     h.pval<-phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
  }else{
    h.pval<-1
    }
    pvals<-rbind(pvals,c(i,overlap,list1,list2,f.pval,h.pval))
  }
  colnames(pvals)<-c('gsName','overlap','dropout.gs.length','msigdb.gs.length','fisher','hyper')
  pvals<-data.frame(pvals)
  pvals$fisher.adj<-p.adjust(pvals$fisher)
  pvals$hyper.adj<-p.adjust(pvals$hyper)
  return(pvals)
}

## do test : ensembl
gs_curated<-MSigDB[["C2_CURATED"]]
gs<-dropout.ensembl$geneName
res<-gsOverlap(gs_curated,gs)
head(res)
sum(res$fisher.adj<0.05)
sum(res$hyper.adj<0.05)
## hypherGtest

## ------------------------------------------------------------------------
dropout.refseq<-subset(ensembl.refseq,Ensembl.Counts>0 & RefSeq.Counts==0)
dim(dropout.refseq)
##barplot(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10],main="Dropout genes: RefSeq")
x<-data.frame(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout genes in RefSeq')

## ------------------------------------------------------------------------
gs<-dropout.refseq$geneName
res<-gsOverlap(gs_curated,gs)
head(res)
sum(res$fisher.adj<0.05)
sum(res$hyper.adj<0.05)

## ------------------------------------------------------------------------
gencode.refseq<-merge(c2,c4,by=c('geneName'),all=F, all.x=F,all.y=F)
dim(gencode.refseq)

## ------------------------------------------------------------------------
subdata<-subset(gencode.refseq,GencodeBasic.Counts>0 & RefSeq.Counts >0)
fit<-lm(RefSeq.Counts~GencodeBasic.Counts, data=subdata)
summary(fit)

## ------------------------------------------------------------------------
pred<-predict(fit,data=gencode.refseq)
conf_interval_3 <- predict(fit, newdata=data.frame(GencodeBasic.Counts=100), interval="prediction",level = 0.95)
conf_interval_3

## ------------------------------------------------------------------------
conf_interval <- predict(fit, newdata=gencode.refseq, interval="prediction",level = 0.95)

## ------------------------------------------------------------------------
##label with color
type<-rep('normal',nrow(gencode.refseq))
type[gencode.refseq$GencodeBasic.Counts==0 & gencode.refseq$RefSeq.Counts>0]<-'dropout.gencode'
type[gencode.refseq$GencodeBasic.Counts>0 & gencode.refseq$RefSeq.Counts==0]<-'dropout.refseq'
type[(gencode.refseq$RefSeq.Counts > conf_interval[,3] | gencode.refseq$RefSeq.Counts<conf_interval[,2])&(gencode.refseq$GencodeBasic.Counts>0 & gencode.refseq$RefSeq.Counts>0)]<-'overdispersion'
gencode.refseq$type<-type
table(gencode.refseq$type)
## label prefiction model in plot
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(fit)[1], digits = 2), b = format(coef(fit)[2], digits = 2), r2 = format(summary(fit)$r.squared, digits = 3)))
## scatter plot
p1 <- ggplot(gencode.refseq,aes(x=GencodeBasic.Counts,y=RefSeq.Counts),color="grey")
p1<-p1+geom_ribbon(aes(ymin=conf_interval[,2],ymax=conf_interval[,3]),alpha=0.2,fill="red")
p1<-p1+geom_point(data=gencode.refseq,aes(x=GencodeBasic.Counts,y=RefSeq.Counts,color=type))
p1<-p1+annotate("text", x = 4.1, y = 12, label = equation(fit), parse = TRUE)
p1<-p1+xlab("log2(GencodeBasic.Counts+1)")+ylab("log2(RefSeq.Counts+1)")
p1

## ------------------------------------------------------------------------
dropout.gencode<-subset(gencode.refseq,GencodeBasic.Counts==0 & RefSeq.Counts>0)
dim(dropout.gencode)
#barplot(sort(table(dropout.gencode$Gencode.type),decreasing=T)[1:10],main="Dropout genes: GencodeBasic")
x<-data.frame(sort(table(dropout.gencode$Gencode.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout gene in Gencode')

## ------------------------------------------------------------------------
gs<-dropout.gencode$geneName
gs_curated<-MSigDB[["C2_CURATED"]]
res<-gsOverlap(gs_curated,gs)
head(res)
sum(res$fisher.adj<0.05)
sum(res$hyper.adj<0.05)

## ------------------------------------------------------------------------
dropout.refseq<-subset(gencode.refseq,GencodeBasic.Counts>0 & RefSeq.Counts==0)
dim(dropout.refseq)
##barplot(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10],main="Dropout genes: RefSeq")
x<-data.frame(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout genes in refseq')

## ------------------------------------------------------------------------
gs<-dropout.refseq$geneName
gs_curated<-MSigDB[["C2_CURATED"]]
res<-gsOverlap(gs_curated,gs)
head(res)
sum(res$fisher.adj<0.05)
sum(res$hyper.adj<0.05)

## ------------------------------------------------------------------------
gencodes<-cbind(c3,c4)
dim(gencodes)
subdata<-subset(gencodes,GencodeBasic.Counts>0 & GencodeComp.Counts >0)
fit<-lm(GencodeComp.Counts ~ GencodeBasic.Counts, data=subdata)
summary(fit)

## ------------------------------------------------------------------------
pred<-predict(fit,data=ensembl.refseq,level=0)
conf_interval_3 <- predict(fit, newdata=data.frame(GencodeBasic.Counts=100), interval="prediction",level = 0.95)
conf_interval_3
conf_interval <- predict(fit, newdata=gencodes, interval="prediction",level = 0.95)

## ------------------------------------------------------------------------
##label with color
type<-rep('normal',nrow(gencodes))
type[gencodes$GencodeBasic.Counts==0 & gencodes$GencodeComp.Counts>0]<-'dropout.basic'
type[gencodes$GencodeBasic.Counts>0 & gencodes$GencodeComp.Counts==0]<-'dropout.comp'
type[(gencodes$GencodeComp.Counts > conf_interval[,3]|gencodes$GencodeComp.Counts<conf_interval[,2])&(gencodes$GencodeComp.Counts>0 & gencodes$GencodeBasic.Counts>0)]<-'overdispersion' 
gencodes$type<-type
table(gencodes$type)
## label prefiction model in plot
eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, list(a = format(coef(fit)[1], digits = 2), b = format(coef(fit)[2], digits = 2), r2 = format(summary(fit)$r.squared, digits = 3)))
## scatter plot
p1 <- ggplot(gencodes,aes(x=GencodeBasic.Counts,y=GencodeComp.Counts),color="grey")
p1<-p1+geom_ribbon(aes(ymin=conf_interval[,2],ymax=conf_interval[,3]),alpha=0.2,fill="red")
p1<-p1+geom_point(data=gencodes,aes(x=GencodeBasic.Counts,y=GencodeComp.Counts,color=type))
p1<-p1+annotate("text", x = 4.1, y = 12, label = equation(fit), parse = TRUE)
p1<-p1+xlab("log2(GencodeBasic.Counts+1)")+ylab("log2(GencodeComp.Counts+1)")
p1

## ------------------------------------------------------------------------
dropout.basic<-subset(gencodes,GencodeComp.Counts>0 & GencodeBasic.Counts==0)
dim(dropout.basic)
##barplot(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10],main="Dropout genes: RefSeq")
x<-data.frame(sort(table(dropout.basic$Gencode.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout genes in GencodeBasic')

## ------------------------------------------------------------------------
dropout.comp<-subset(gencodes,GencodeBasic.Counts>0 & GencodeComp.Counts==0)
dim(dropout.comp)
##barplot(sort(table(dropout.refseq$RefSeq.type),decreasing=T)[1:10],main="Dropout genes: RefSeq")
x<-data.frame(sort(table(dropout.comp$Gencode.type),decreasing=T)[1:10])
colnames(x)<-c('Biotype','Counts')
ggplot(data.frame(x),aes(x=Biotype,y=Counts))+geom_bar(aes(fill = Biotype), position = "dodge", stat="identity")+xlab('Biotype')+ylab("# of dropout genes")+geom_text(aes(label=Counts), vjust=1.6, color="black", size=3.5) +theme(axis.text.x=element_blank())+ggtitle('dropout genes in GencodeComp')

## ------------------------------------------------------------------------

g1<-c1[c1$Ensembl.Counts>0,1]
g2<-c2[c2$RefSeq.Counts>0,1]
g3<-c3[c3$GencodeComp.Counts>0,1]
g4<-c4[c4$GencodeBasic.Counts>0,1]
##venn.plot<-venn.diagram(x = list('Ensembl'=g1,'RefSeq'=g2,'Gencode'=g3,'GencodeBasic'=g4), imagetype='png',filename = "Venn.png",
##col = "transparent", fill = c("cornflowerblue","green","yellow","darkorchid1"),
##label.col = c("orange", "white", "darkorchid4", "white", "white", 
##"white",    "white", "white", "darkblue", "white", "white", "white", "white", 
##"darkgreen", "white"),alpha = 0.50,cex = 1.5, fontfamily = "serif", fontface = "bold",
##cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,
##cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 270,
##margin = 0.2)
venn(list('Ensembl'=g1,'RefSeq'=g2,'Gencode'=g3,'GencodeBasic'=g4),show.plot = T)

## ------------------------------------------------------------------------
g1<-c1[c1$Ensembl.Counts>5,1]
g2<-c2[c2$RefSeq.Counts>5,1]
g3<-c3[c3$GencodeComp.Counts>5,1]
g4<-c4[c4$GencodeBasic.Counts>5,1]
##venn.plot<-venn.diagram(x = list('Ensembl'=g1,'RefSeq'=g2,'Gencode'=g3,'GencodeBasic'=g4), imagetype='png',filename = "Venn_good_genes.png",
##col = "transparent", fill = c("cornflowerblue","green","yellow","darkorchid1"),
##label.col = c("orange", "white", "darkorchid4", "white", "white", 
##"white",    "white", "white", "darkblue", "white", "white", "white", "white", 
##"darkgreen", "white"),alpha = 0.50,cex = 1.5, fontfamily = "serif", fontface = "bold",
##cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,
##cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 270,
##margin = 0.2)
venn(list('Ensembl'=g1,'RefSeq'=g2,'Gencode'=g3,'GencodeBasic'=g4),show.plot = T)

