library(ggplot2)
library(ggpubr)
library(stringr)
library(here)
set_here("..")

H_score<-function(numVec){
  cmp.numVec<-numVec[!is.na(numVec)]
  if(length(cmp.numVec)<3 | sd(cmp.numVec)==0){
    return(NA)
  }
  vec<- (cmp.numVec - mean(cmp.numVec))/sd(cmp.numVec)
  p.val<-ks.test(vec,"pnorm")$p.value
  return(2*pnorm(log(1-p.val)))
}

cancer.vec<-c("LUAD","OV")
Hscore.ls<-list()
g.ls<-list()
# pro.names<-character()
for(cIdx in 1:length(cancer.vec)){
  tmp.ls<-list()
  cancer.type<-cancer.vec[cIdx]
  Y.df<-readRDS(here("TCGA_Data","Proteomics_Data_Y",paste0("Y_",cancer.type,".rds")))
  Y.mean<-apply(Y.df,2,mean)
  Y.mean.mat<-matrix(rep(Y.mean,nrow(Y.df)),byrow=T, nrow=nrow(Y.df))
  Y.df.cen<-Y.df - Y.mean.mat
  Y.df<-Y.df.cen
  rm(Y.df.cen)
  ProHS<-rep(NA,ncol(Y.df))
  for(gIdx in 1:length(ProHS)){
    ProHS[gIdx]<-H_score(Y.df[,gIdx])
    names(ProHS)[gIdx]<-colnames(Y.df)[gIdx]
  }
  Hscore.ls[[cIdx]]<-ProHS
  if(cancer.type=="LUAD"){
    protein1<-paste0("M_","Akt")
    protein2<-paste0("M_","PTEN")
  }
  if(cancer.type=="OV"){
    protein1<-paste0("M_","E-Cadherin")
    protein2<-paste0("M_","Rb")
  }
  pro.names<-c(protein1,protein2)
  Idx1<-which(names(ProHS)==protein1)
  Idx2<-which(names(ProHS)==protein2)
  pr.df<-Y.df[,c(Idx1,Idx2)]
  colnames(pr.df)<-c(protein1,protein2)
  pr.plt.df<-reshape2::melt(pr.df)
  g1<-ggplot(pr.plt.df, aes(sample = value,color=variable)) + geom_qq() + geom_qq_line()+ theme(text = element_text(size = 20))
  g2<-ggplot(pr.plt.df[pr.plt.df$variable==protein1,], aes(x=value)) + geom_density() + 
    stat_function(fun = dnorm,args = list(mean = mean(pr.df[,protein1]), sd = sd(pr.df[,protein1])),col="#1b98e0") +
    theme(text = element_text(size = 20))
  g3<-ggplot(pr.plt.df[pr.plt.df$variable==protein2,], aes(x=value)) + geom_density() + 
    stat_function(fun = dnorm,args = list(mean = mean(pr.df[,protein2]), sd = sd(pr.df[,protein2])),col="#1b98e0") +
    theme(text = element_text(size = 20))
  tmp.ls[[1]]<-g1
  tmp.ls[[2]]<-g2
  tmp.ls[[3]]<-g3
  g.ls[[cIdx]]<-tmp.ls
  names(g.ls)[cIdx]<-cancer.type
}

# sort(Hscore.ls[[2]])
Hscore.df<-data.frame(LUAD=Hscore.ls[[1]],OV=Hscore.ls[[2]])
Hscore.plt<-reshape2::melt(Hscore.df)
g.HS<-ggplot(Hscore.plt,aes(x=variable, y=value)) + geom_violin(scale = "width") + geom_boxplot(width=0.1) + 
  labs(x="Cancer",y="H score") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g.den<-ggarrange(g.ls$LUAD[[2]],g.ls$LUAD[[3]],g.ls$OV[[2]],g.ls$OV[[3]],nrow=1,ncol=4)
g.QQ<-ggarrange(g.ls$LUAD[[1]],g.ls$OV[[1]],nrow=1,ncol=2)
g.ind<-ggarrange(g.den,g.QQ,nrow=2,ncol=1)
g.final<-ggarrange(g.ind,g.HS,nrow=1,ncol=2,widths = c(1, 0.5))
# g.LUAD.den<-ggarrange(g.ls$LUAD[[2]],g.ls$LUAD[[3]],nrow=1,ncol=2,labels = pro.names[1:2],font.label=list(size=25))
# g.LUAD<-ggarrange(g.ls$LUAD[[1]],g.LUAD.den,nrow=2,ncol=1)
# g.OV.den<-ggarrange(g.ls$OV[[2]],g.ls$OV[[3]],nrow=1,ncol=2,labels = pro.names[3:4],font.label=list(size=25))
# g.OV<-ggarrange(g.ls$OV[[1]],g.OV.den,nrow=2,ncol=1)
# g.final<-ggarrange(g.LUAD,g.OV,g.HS,nrow=1,ncol=3, labels=c("(A) LUAD","(B) OV","(C)"),font.label=list(size=30))
ggsave(plot=g.final,file=here("Plots","Hscore.png"),width=32,height=18)