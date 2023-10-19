library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(gridExtra)
library(here)

immune.vec<-c("Intercept","T cells","Monocytes","Neutrophils")
cancer.vec<-c("LUAD","OV")
Z.plt.ls<-list()
g.net<-list()
deg.ls<-pop.coef.ls<-list()
for(cIdx in 1:length(cancer.vec)){
  print(cIdx)
  cancer.type<-cancer.vec[cIdx]
  f.path<-here("..","MCMC_Res",cancer.type)
  pop.res<-readRDS(here(f.path,"sym_intSS_PopRes.RDS"))
  
  match.df<-readRDS(file=here("..","MCMC_Res","match_df.rds"))
  Y.df<-as.matrix(readRDS(here("..","TCGA_Data","Proteomics_Data_Y",paste0("Y_",cancer.type,".rds"))))
  match.df$Id<-match(match.df$protein,str_sub(colnames(Y.df),3))
  match.df$Pathway<-ifelse(match.df$multi.Path==0,match.df$PRECISE.cat,"Core reactive")
  match.df<-match.df[,c("Id","protein","Pathway")]
  p<-ncol(Y.df)
  q<-length(immune.vec)
  Z.all<-pop.res$Z
  Z.plt<-matrix(NA,nr=0,nc=4)
  rIdx<-0
  for(j in 1:(p-1)){
    for(k in (j+1):p){
      for(h in 1:q){
        rIdx<-rIdx+1
        tmp<-c(j,k,h,Z.all[j,k,h])
        Z.plt<-rbind(Z.plt,tmp)
      }
    }
  }
  colnames(Z.plt)<-c("j","k","h","PIP")
  rownames(Z.plt)<-NULL
  Z.plt<-data.frame(Z.plt)
  Z.plt$Cancer<-paste0(cancer.type)
  Z.plt.ls[[cIdx]]<-Z.plt
  pop.coef<-data.frame(pop.res$coef)
  pop.coef$Immune.type<-ifelse(pop.coef$h==1,immune.vec[1],
                               ifelse(pop.coef$h==2,immune.vec[2],
                                      ifelse(pop.coef$h==3,immune.vec[3],immune.vec[4])))
  pop.coef.ls[[cIdx]]<-pop.coef
  g.h<-list()
  deg.df<-match.df
  deg.df$Cancer<-cancer.type
  tmp<-matrix(0,nrow=nrow(deg.df),ncol=length(immune.vec))
  colnames(tmp)<-immune.vec
  deg.df<-cbind(deg.df,tmp)
  match.df<-match.df[order(match.df$Pathway,match.df$protein),]
  for(h in 1:q){
    tmp.df<-pop.coef[pop.coef$h==h,]
    deg<-table(c(tmp.df$j,tmp.df$k))
    match.df$Degree<-1
    for(dIdx in 1:length(deg)){
      match.df$Degree[match.df$Id==names(deg)[dIdx]]<-(deg[dIdx]*50)
      deg.df[deg.df$Id==names(deg)[dIdx],colnames(deg.df)==immune.vec[h]]<-deg[dIdx]
    }
    ig<-graph_from_data_frame(tmp.df, directed = FALSE, vertices = match.df)
    g<-ggraph(ig,layout = "circle") +
      geom_edge_link() +
      # geom_node_point(aes(size=Degree,color=Pathway)) +
      geom_node_point(aes(size=Degree),color=brewer.pal(n = 3, name = 'Greys')[2]) +
      scale_size_continuous(range = c(0.5, 1.5*max(deg)))+
      geom_node_text(aes(label = protein),size=10, repel=TRUE)+
      guides(size="none") + labs(caption=immune.vec[h])+theme_graph() + 
      theme(legend.title=element_text(size=30),legend.text = element_text(size = 30),
            plot.caption = element_text(hjust=0.5,size=30,face='bold'))
    if(cIdx==1){
      g<-g+theme(plot.background = element_rect(fill = brewer.pal(n = 9, name = 'Reds')[1]))
    }
    if(cIdx==2){
      g<-g+theme(plot.background = element_rect(fill = brewer.pal(n = 9, name = 'Blues')[1]))
    }
    g.h[[h]]<-g
  }
  g.net[[cIdx]]<-g.h
  deg.ls[[cIdx]]<-deg.df
}
PIP.plt<-rbind(Z.plt.ls[[1]],Z.plt.ls[[2]])
PIP.plt$Immune.type<-ifelse(PIP.plt$h==1,immune.vec[1],
                            ifelse(PIP.plt$h==2,immune.vec[2],
                                   ifelse(PIP.plt$h==3,immune.vec[3],immune.vec[4])))
PIP.plt$Immune.type.fac<-factor(PIP.plt$Immune.type,levels = immune.vec)
g.PIP<-ggplot(PIP.plt[PIP.plt$h!=1,],aes(x=Immune.type.fac,y=PIP,fill=Cancer))+geom_boxplot()+labs(x="Immune Component")+theme(text=element_text(size=30))
PIP.plt[,c("PIP","Cancer","Immune.type")] %>% group_by(Cancer, Immune.type) %>%
  summarize_all(list(mean = ~ mean(., na.rm = TRUE), sd= ~ sd(.,na.rm=T), median=~ median(., na.rm=T))) 

g.LUAD<-ggarrange(plotlist = g.net[[1]][2:4],nrow=1,ncol=3)
g.OV<-ggarrange(plotlist = g.net[[2]][2:4],nrow=1,ncol=3)
g.network<-ggarrange(g.LUAD,g.OV,nrow=2,ncol=1,labels=c("(A) LUAD", "(B) OV"),font.label = list(size=30))
g.final<-ggarrange(g.network,g.PIP,widths = c(3,1),ncol=2,labels=c("","(C) PIP"),hjust=c(0,0.1),font.label = list(size=30))
ggsave(plot=g.final,file=here("..","Plots","intSSPop_net_all.png"),width=48,height=18)
