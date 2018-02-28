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

##convert p-value to star char
convert2star<-function(x){
  if(x>0.05){
    p='ns'
  }else if(x> 0.01){
    p='*'
  }else if(x>0.001){
    p="**"
  }else if(x>0.0001){
    p='***'
  }else{
    p='****'
  }
  return(p)
}
## standard theme
addTheme<-function(p){
  p<-p+theme(legend.position ="top")
  p<-p+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  p<-p+theme(axis.title.x = element_text(color='black',size=18,face='bold'),axis.title.y = element_text(size=18,color='black',face='bold'))
  p<-p+theme(axis.text.y = element_text(color='black',size=18,face='bold'),axis.text.x = element_text(color='black',size=18,face='bold'))
  p<-p+theme(legend.text=element_text(size=15,face='bold'),legend.title = element_blank())
  return(p)
}

## params
option_list <- list(
  
  make_option("--bmetrics", type="character", default=NULL, 
              help=" base metricsfile name", metavar="character"),
  make_option("--umetrics", type="character", default=NULL, 
              help=" updated metricsfile name", metavar="character"),
  make_option("--metKeys", type="character",default=NULL, 
              help=" a list of metrics name ", metavar="character"),
  make_option("--out", type="character", default="out", 
              help="output file name [default= %default]", metavar="character")
)
opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)
## load files
metKeys<-opt$metKeys
output_name<-opt$out
# checking data format and header
met1<-read.table(opt$bmetrics,header=T,sep=',',stringsAsFactors=F)
met2<-read.table(opt$umetrics,header=T,sep=',',stringsAsFactors=F)
colnames(met1)[1]<-'metrics'
colnames(met2)[1]<-'metrics'
colnames1<-colnames(met1)
colnames2<-colnames(met2)
mlist<-match(colnames1,colnames2)
nalist<-is.na(mlist)
if(sum(nalist)>0){
  print("input files have discrepency columns")
  exit()
}
# re-order 
met2<-met2[,mlist]
# select subset of metrics 
mlist1<-match(metKeys,met1$metrics)
met1.core<-met1[mlist1,]
mlist2<-match(metKeys,met2$metrics)
met2.core<-met2[mlist2,]
out<-c()
pouts<-list()
for(ii in 1:length(metKeys)){
  x<-as.numeric(met1.core[ii,-1])
  y<-as.numeric(met2.core[ii,-1])
  z<-data.frame('Base'=x,'Updated'=y)
  # linear regression model
  fit<-lm(y~x)
  sfit<- summary(fit)
  r2<-round(sfit$adj.r.squared,3)
  beta<-round(sfit$coefficients[2],3)
  a<-round(sfit$coefficients[1],3)
  f<-sfit$fstatistic
  fpval<-pf(f[1],f[2],f[3],lower.tail=F)
  fpval.s<-convert2star(fpval)
  coefs<-data.frame('ic'=c(0,a),'s'=c(1,beta),'tl'=c('1x1','fitted'))
  my.formula <- y ~ x
  p <- ggplot(data = z, aes(x = Base, y = Updated)) +geom_point(shape=1)+theme_bw(base_size = 20)
  p<-p+  stat_poly_eq(formula = my.formula,
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"),size=8), 
                 parse = TRUE,size = 8) 
  p<-p+stat_fit_glance(method = 'lm', method.args = list(formula = my.formula),geom = 'text', aes(label = paste("P-value:",fpval.s , sep = ""),size=12),label.x.npc = 'left',
                       label.y.npc = 0.85, size = 8)
  p<-p+theme(legend.position ="top")
  p<-addTheme(p)
  p<-p+xlab(paste('Base'))+ylab(paste('Updated'))
  p<-p+geom_abline(data=coefs,mapping=aes(slope=s, intercept=ic, linetype=factor(tl),color=factor(tl)))

  # hist 
  z<-data.frame('Base'=x,'Updated'=y)
  mz<-melt(z)
  mu <- ddply(mz, "variable", summarise, grp.mean=mean(value))
  stable <- desc_statby(mz,measure.var = 'value',grps='variable')
  stable <- stable[, c("variable", "length", "mean", "sd")]
  stable$mean<-round(stable$mean,3)
  stable$sd<-round(stable$sd,3)
  stable.p<-ggtexttable(stable,row=NULL, theme=ttheme('mOrange',base_size = 18))
  density.p<-ggplot(data=mz,aes(x=value,color=variable,fill=variable))
  density.p<-density.p+theme_bw(base_size = 20)
  density.p<-density.p+geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),linetype="dashed")
  density.p<-density.p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  density.p<-density.p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  density.p<-density.p+theme(legend.position="top")
  density.p<-density.p+geom_histogram(aes(y=..density..), position="identity", alpha=0.5)
  density.p<-density.p+geom_density(alpha=0.6)
  density.p<-addTheme(density.p)
  density.p<-density.p+annotation_custom(ggplotGrob(stable.p))
  
  # ks test
  ks<-ks.test(x,y)
  cdf1<-ecdf(x)
  cdf2<-ecdf(y)
  minMax <- seq(min(x, y), max(x, y), length.out=length(x)) 
  x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
  y0 <- cdf1(x0) 
  y1 <- cdf2(x0)
  ks.p<-ggplot(mz, aes(x = value, group = variable, color = variable))+
    stat_ecdf(size=1) +
    theme_bw(base_size = 20) +
    theme(legend.position ="top") +
    ylab("ECDF")
  
  dtable<-data.frame('test'='K-S','D-stats'=round(ks$statistic,4),'P-value'=convert2star(ks$p.value))
  dtable.p<-ggtexttable(dtable,row=NULL, theme=ttheme('mBlue',base_size = 18))
  ks.p<-ks.p+geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),linetype = "dashed", color = "red") 
  ks.p<-ks.p+geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=4) 
  ks.p<-ks.p+geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=4)
  ks.p<-ks.p+ggtitle(paste("K-S Test"))
  ks.p<-addTheme(ks.p)
  ks.p<-ks.p+theme(legend.title=element_blank())+annotation_custom(ggplotGrob(dtable.p))
  # arrange layout 
  gt<-arrangeGrob(density.p, p, ks.p, ncol = 2, nrow = 2, layout_matrix = rbind(c(1,1), c(2,3)))
  gp <- as_ggplot(gt) + draw_plot_label(label = c("A", "B", "C"), size = 20,x = c(0, 0, 0.5), y = c(1, 0.5, 0.5)) # Add labels
  gp<-gp+ggtitle(paste(metKeys[ii]))+theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))
  ggsave(gp,file=paste(output_name,'/group_plots_',metKeys[ii],'.png',sep=''),width=20,height=20)
  pouts[[ii]]<-gp
  out<-rbind(out,c(metKeys[ii],beta,a,r2,fpval,ks$statistic,ks$p.value))
}
# save multiple page into one pdf
pdf(paste(output_name,'/group_plots_all.pdf',sep=''),25,25)
for(ii in 1:nrow(met1.core)){print(pouts[[ii]])}
dev.off()
colnames(out)<-c('metrics','beta','a','r2','pvalue','ks-D-stats','ks-Pvalue')
write.table(out,file=paste(output_name,'/tests_stats.csv',sep=''),quote=F,row.names=F,col.names=T,sep=',')
