library(ggplot2)
setwd("~/Research/GR/rBGR/")
source("./code/expTail/rBGR_func.R")

args <- commandArgs(TRUE)
print(args)
j<-as.numeric(args[1])
nIdx<-as.numeric(args[2])
b<-as.integer(args[3])
n<-as.integer(args[4])
p<-as.integer(args[5])
q<-as.integer(args[6])

nonNormLvVec<-c(0,0.5,0.8)
nonNormLv<-nonNormLvVec[nIdx]
print(paste0("nonNormLv: ", nonNormLv, "; b: ", b))
df.ls<-readRDS(paste0("./data_sim/expTail/n",n,"p",p,"q",q,"/nonNormLv",nonNormLv,".rds"))[[b]]
Y.df<-df.ls$X.nonNorm
X.df<-df.ls$U
tail<-rep("Exp",ncol(Y.df))
D_dist<-"expTail"
iteNum=10000
# iteNum=3000
# iteNum=300
print(j)
system.time({
print(iteNum)
fit.mcmc<-rBGR_mcmc_Int(y=Y.df[,j], x=Y.df[,-j], u=X.df,
                        iteNum=iteNum, burnIn=(iteNum * 0.95), 
                        seed_=j *5679 + b*12 + nIdx, priorB = 1, sigProp = 0.2, tail = tail)
})

saveRDS(fit.mcmc,paste0("./MCMCRes/",D_dist,"/n",n,"p",p,"q",q,"/rBGR/nonNormLv",nonNormLv,"_b",b,"_j",j,".rds"))


