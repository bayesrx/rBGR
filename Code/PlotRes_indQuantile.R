library(ggplot2)
library(ggpubr)
library(stringr)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(here)


f.path<-here("..","MCMC_Res")
immune.vec.tmp<-c("T cells","Monocytes","Neutrophils")
immune.vec<-c("Intercept",immune.vec.tmp)
cancer.vec<-c("LUAD","OV")
q<-length(immune.vec)

match.ls<-Y.ls<-X.ls<-ind.res.ls<-pop.res.ls<-pop.samePathway<-list()
for(cIdx in 1:length(cancer.vec)){
  print(cIdx)
  cancer.type<-cancer.vec[cIdx]
  pop.coef<-data.frame(readRDS(here(f.path,cancer.type,"sym_intSS_PopRes.RDS"))$coef)
  
  Y.df<-as.matrix(readRDS(here("..","TCGA_Data","Proteomics_Data_Y",paste0("Y_",cancer.type,".rds"))))
  Y.mean<-apply(Y.df,2,mean)
  Y.mean.mat<-matrix(rep(Y.mean,nrow(Y.df)),byrow=T, nrow=nrow(Y.df))
  Y.df.cen<-Y.df - Y.mean.mat
  Y.df<-Y.df.cen
  rm(Y.df.cen)
  Y.ls[[cIdx]]<-Y.df
  
  X.df<-as.matrix(readRDS(here("..","TCGA_Data","mRNA_Immune_X",paste0("X_",cancer.type,".rds"))))
  X.df<-log(X.df)
  X.sd<-apply(X.df,2,sd)
  X.sd.mat<-matrix(rep(X.sd,nrow(X.df)),nrow=nrow(X.df),byrow=T)
  X.mean<-apply(X.df,2,mean)
  X.mean.mat<-matrix(rep(X.mean,nrow(X.df)),nrow=nrow(X.df),byrow=T)
  X.df.reScl<-(X.df-X.mean.mat)/X.sd.mat
  X.df<-X.df.reScl[,c("T cells","Monocytes","Neutrophils")]
  rm(X.df.reScl)
  X.int<-cbind(1,X.df)
  rm(X.df)
  X.ls[[cIdx]]<-X.int
  
  match.df<-readRDS(file=here("..","MCMC_Res","match_df.rds"))
  match.df$Id<-match(match.df$protein,str_sub(colnames(Y.df),3))
  match.df$Pathway<-ifelse(match.df$multi.Path==0,match.df$PRECISE.cat,"Core reactive")
  match.df<-match.df[,c("Id","protein","Pathway")]
  match.ls[[cIdx]]<-match.df
  
  pop.coef$k.pathway<-pop.coef$j.pathway<-pop.coef$k.protein<-pop.coef$j.protein<-NA
  for(rIdx in 1:nrow(pop.coef)){
    pop.coef$j.protein[rIdx]<-match.df$protein[match.df$Id==pop.coef$j[rIdx]]
    pop.coef$k.protein[rIdx]<-match.df$protein[match.df$Id==pop.coef$k[rIdx]]
    pop.coef$j.pathway[rIdx]<-match.df$Pathway[match.df$Id==pop.coef$j[rIdx]]
    pop.coef$k.pathway[rIdx]<-match.df$Pathway[match.df$Id==pop.coef$k[rIdx]]
  }
  pop.res.ls[[cIdx]]<-pop.coef
  res<-pop.coef[pop.coef$j.pathway==pop.coef$k.pathway,]
  pop.samePathway[[cIdx]]<-res
}

quan.ls.noZ<-list()
for(cIdx in 1:length(cancer.vec)){
  cancer.type<-cancer.vec[cIdx]
  print(cIdx)
  Y.df<-Y.ls[[cIdx]]
  p<-ncol(Y.df)
  X.df<-X.ls[[cIdx]]
  q<-ncol(X.df)
  # immu.comp.edge<-data.frame(matrix(NA,nr=p*(p-1)*q*5/2,nc=15))
  # colnames(immu.comp.edge)<-c("cancer","j","k","h","X.value",
  #                             "PrNeg_jk","PrNeg_kj","PrZero_jk","PrZero_kj","PrPos_jk","PrPos_kj",
  #                             "j.protein","k.protein","j.pathway","k.pathway")
  immu.comp.edge<-data.frame(matrix(NA,nr=p*(p-1)*(q-1)*7/2,nc=13))
  colnames(immu.comp.edge)<-c("cancer","j","k","h","X.value",
                              "PrNeg_jk","PrNeg_kj","PrZero_jk","PrZero_kj","PrPos_jk","PrPos_kj",
                              "Sign","ePP")
  immu.comp.edge$cancer<-cancer.type

  ind.res<-readRDS(here(f.path,cancer.type,"/IteWise_intSS_SymParam_noZ.RDS"))
  alpha.mcmc<-ind.res$alpha
  t.mcmc<-ind.res$t
  iteNum<-dim(t.mcmc)[1]
  rIdx<-0
  for(h in 2:q){
    X.quan<-quantile(X.df[,h],c(0,0.05,0.25,0.5,0.75,0.95,1))
    X.df.pred<-cbind(1,matrix(0,nrow=length(X.quan),ncol=q-1))
    X.df.pred[,h]<-X.quan
    n.pred<-nrow(X.df.pred)
    ePos<-eNeg<-eZero<-array(0, dim=c(p,n.pred,p))
    
    print("Calculate nonZero Posterior Probability")
    for(j in 1:p){
      # print(j)
      negTmp<-zeroTmp<-posTmp<-array(0,dim=c(n.pred,p-1))
      for(l in 1:iteNum){
        alpha.k<-alpha.mcmc[l,j,-j,]
        theta.j.X<-X.df.pred %*% t(alpha.k)
        beta.j.X<-theta.j.X * as.numeric(abs(theta.j.X)>t.mcmc[l,j])
        posTmp<-posTmp+as.numeric(beta.j.X>0)
        negTmp<-negTmp+as.numeric(beta.j.X<0)
        zeroTmp<-zeroTmp+as.numeric(beta.j.X==0)
      }
      posTmp<-posTmp/iteNum
      negTmp<-negTmp/iteNum
      zeroTmp<-zeroTmp/iteNum
      ePos[j,,-j]<-posTmp
      eNeg[j,,-j]<-negTmp
      eZero[j,,-j]<-zeroTmp
    }
    
    
    for(i in 1:n.pred){
      net.i.pos<-ePos[,i,]
      net.i.neg<-eNeg[,i,]
      net.i.zero<-eZero[,i,]
      for(j in 1:(p-1)){
        for(k in (j+1):p){
          rIdx<-rIdx+1
          # j.protein<-match.df$protein[match.df$Id==j]
          # k.protein<-match.df$protein[match.df$Id==k]
          # j.pathway<-match.df$Pathway[match.df$Id==j]
          # k.pathway<-match.df$Pathway[match.df$Id==k]
          immu.comp.edge[rIdx,2:11]<-c(j,k,h,X.df.pred[i,h],
                                       net.i.neg[j,k], net.i.neg[k,j],
                                       net.i.zero[j,k], net.i.zero[k,j],
                                       net.i.pos[j,k], net.i.pos[k,j])
          nonZeroProb<-1-immu.comp.edge[rIdx,c("PrZero_jk","PrZero_kj")]
          selIdx<-which.max(nonZeroProb) # Use Max Rules
          immu.comp.edge$ePP[rIdx]<-as.numeric(nonZeroProb[selIdx])
          if(immu.comp.edge$ePP[rIdx]<0.5){
            immu.comp.edge$Sign[rIdx]=0
          }else{
            if(selIdx==1){
              immu.comp.edge$Sign[rIdx]<-2*(which.max(immu.comp.edge[rIdx,c("PrNeg_jk","PrPos_jk")])-1)-1
            }else{
              immu.comp.edge$Sign[rIdx]<-2*(which.max(immu.comp.edge[rIdx,c("PrNeg_kj","PrPos_kj")])-1)-1
            }
          }
        }
      }
    }
  }
  quan.ls.noZ[[cIdx]]<-immu.comp.edge
}
saveRDS(quan.ls.noZ,file=here(f.path,"quanPred_intSS_res_noZ.RDS"))
quan.ls.noZ<-readRDS(file=here(f.path,"quanPred_intSS_res_noZ.RDS"))

cIdx<-2
cancer.type<-cancer.vec[cIdx]
PrBeta<-quan.ls.noZ[[cIdx]]
PrBeta$Sign<-ifelse(PrBeta$Sign==1,"Positive",ifelse(PrBeta$Sign==-1,"Negative","Zero"))
PrBeta$Sign.fac<-factor(PrBeta$Sign,levels = c("Negative","Zero","Positive"))
g.ls<-list()
for(h in 2:q){
  PrBeta.h<-PrBeta[PrBeta$h==h,]
  quan.val.vec<-unique(PrBeta.h$X.value)
  g.quan.h<-list()
  for(qIdx in 1:length(quan.val.vec)){
    quan.val<-quan.val.vec[qIdx]
    PrBeta.h.quan<-PrBeta.h[PrBeta.h$X.value==quan.val,]
    
    # ggplot(PrBeta.h.quan,aes(x=j,y=k,fill=ePP))+geom_tile()
    
    match.df<-match.ls[[cIdx]]
    match.df<-match.df[order(match.df$protein,match.df$Pathway),]
    plt.df<-PrBeta.h.quan[PrBeta.h.quan$ePP>0.5,c("j","k","ePP","Sign.fac")]
    ig<-graph_from_data_frame(plt.df, directed = FALSE, vertices = match.df)
    # ig.sign<-graph_from_data_frame(plt.sign.df, directed = FALSE, vertices = match.df)
    if(nrow(plt.df)==0){
      g<-ggraph(ig,layout = "circle") +
        geom_node_point()+
        geom_node_text(aes(label = protein),size=6, repel=TRUE)+
      # guides(alpha="none",color="none") +
        theme_graph() + theme(legend.title=element_text(size=30),legend.text = element_text(size = 30))
    }else{
      g<-ggraph(ig,layout = "circle") +
        geom_edge_link(aes(color=factor(Sign.fac,levels = c("Negative","Zero","Positive")))) +
        geom_node_point() +
        scale_edge_colour_manual(name="Sign",values=c("Negative"="#F8766D","Zero"="black","Positive"="#00BA38"),drop=T)+
        geom_node_text(aes(label = protein),size=6, repel=TRUE)+
        theme_graph()
    }
    if(h==2){
      g<-g+theme(legend.title=element_text(size=30),legend.text = element_text(size = 30),
                 plot.background = element_rect(fill = brewer.pal(n = 9, name = 'Reds')[1]))
    }
    if(h==3){
      g<-g+theme(plot.background = element_rect(fill = brewer.pal(n = 9, name = 'YlGn')[1]),
                 legend.position = "none")
    }
    if(h==4){
      g<-g+theme(plot.background = element_rect(fill = brewer.pal(n = 9, name = 'Greens')[1]),
                 legend.position = "none")
    }
    g.quan.h[[qIdx]]<-g
  }
  g.h<-ggarrange(g.quan.h[[1]],g.quan.h[[3]],g.quan.h[[4]],g.quan.h[[5]],g.quan.h[[7]],nrow=1,common.legend = T)
  # ggsave(plot=g.h,file=paste0("./Plots/RPPA/hotcold_cmp/LUAD_intSS_quan_",immune.vec[h],".png"),width=32,heigh=9)
  g.ls[[(h-1)]]<-g.quan.h
}

g.T<-ggarrange(plotlist = g.ls[[1]][2:6],common.legend = T,ncol=5,nrow=1)
g.Mono<-ggarrange(plotlist = g.ls[[2]][2:6],common.legend = F,ncol=5,nrow=1)
g.Neu<-ggarrange(plotlist = g.ls[[3]][2:6],common.legend = F,ncol=5,nrow=1)
gLab<-paste0(rep("(",3),LETTERS[1:3],rep(") ",3),immune.vec[2:4])
g.network<-ggarrange(g.T,g.Mono,g.Neu,labels=gLab,nrow=3,ncol=1,
                     hjust=nchar(gLab)*0.0001,font.label = list(size=30),common.legend = F)

n <- 1800
df.triangle <- data.frame(xstart = 1:n,xend=1:n, ystart=0,yend=seq(0,1,length.out=n))
g.triangle<-ggplot(df.triangle) +
  geom_segment(aes(x=xstart,xend=xend,y=ystart,yend=yend,color=xstart)) +
  scale_color_gradient(limits=c(1, n), low = brewer.pal(n=9,name="Blues")[2], high = brewer.pal(n=9,name="Blues")[9])+
  theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 

g.final<-ggarrange(g.network,g.triangle,nrow=2,ncol=1,heights=c(3,0.5),common.legend = F) +bgcolor("white")

# g.final<-ggarrange(g.T,g.Mono,g.Neu,g.triangle,nrow=4,ncol=1,heights=c(1,1,1,0.5),
#                    labels=c("(A) T Cells","(B) Monocytes","(C) Neutrophils","Immune Components Abundance"),
#                    font.label = list(size=30),common.legend = T) + 
#   bgcolor("white")
ggsave(plot=g.final,file=here("..","Plots",paste0(cancer.type,"_quan_intSS_immuAll.png")),width=32,heigh=27)
