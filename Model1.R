# setwd("C:/Users/MyPC/Desktop/Model1")
##example 1: categorical covariates and binary response variable
stimu1 <- function(ind,J,R,P,N,ymethod,p,j,d0){
  #p: number of group 
  #j: number of predictors in each group
  if (ymethod=="balanced"){
    pr=rep(1/R,R)
    set.seed(ind)
    y=round(runif(N, min = 0.5, max = 2.5))
  }else if (ymethod=="unbalanced"){
    pr=2*(1+(R-c(1:R))/(R-1))/(3*R)
    set.seed(ind)
    y=rbinom(N,1,pr[2])+1
  }
  
  # group indices
  if (p*j<P) {group <- c(rep(1:p,each = j),c((p+1):(P-p*j+p))) } 
  if (p*j==P) {group <- c(rep(1:p,each = j))}
  
  #Bulid latent variable
  z=matrix(0,nrow = N,ncol = P)
  for (p in 1:P) {
    if (p<=d0) {
      set.seed(p);
      y1=which(y==1)
      z[y1,p]=rnorm(length(y1),-0.5,1)
      y2=which(y==2)
      z[y2,p]=rnorm(length(y2),0.5,1)
    }
    if (p>d0) {
      set.seed(p);
      z[,p]=rnorm(N,0,1)}
  }
  
  #Bulid x variable
  x=matrix(0,nrow = N,ncol = P)
  for (p in 1:p) {
    #odd number
    if (p%%2!=0) {
      for (i in 1:N) {
        x[i,p]=ifelse(z[i,p]>qnorm(0.5,0,1),2,1)
      }
    }
    #even number
    if (p%%2==0) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(0.2,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(0.2,0,1))&&(z[i,p]<=qnorm(0.4,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(0.4,0,1))&&(z[i,p]<=qnorm(0.6,0,1))) {x[i,p]=3}
        if ((z[i,p]>qnorm(0.6,0,1))&&(z[i,p]<=qnorm(0.8,0,1))) {x[i,p]=4}
        if (z[i,p]>qnorm(0.8,0,1)) {x[i,p]=5}
      }
    }
  }
  
  #da=data.frame(y,x);
  #colnames(da)=c("y",paste("x",group,sep = ""))
  re=list()
  re$y=y
  re$x=x
  re$group=group
  return(re)
}
realda=stimu1(ind=12,J=4,R=2,P=1500,N=200,ymethod = "balanced",p=3,j=3,d0=9)
realda1=stimu1(ind=12,J=4,R=2,P=1500,N=200,ymethod = "balanced",p=200,j=3,d0=9)
realda2=stimu1(ind=12,J=4,R=2,P=1500,N=200,ymethod = "balanced",p=500,j=3,d0=20)
head(realda2)
y=realda$y
x=realda$x
group=realda2$group
group

source("GIF-SIS.R")
source("GIGR-SIS.R")
source("GRPSIS.R")
source("IG-SIS.R")

GIGRSIS(cls = realda2$y,atr = realda2$x,method = "inforate",group)#GIGRSIS
IGSIS(cls = realda$y,atr = realda$x,method = "inforate")#IGSIS
GIFSIS(cls = realda$y,atr = realda$x,method = "inforate",group)#GIFSIS
GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gSIS")#gSIS
GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gHOLP")#gHOLP



testone <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  
  pc11=matrix(0,nrow = rep,ncol = length(activex))
  pc12=matrix(0,nrow = rep,ncol = length(activex))
  pc13=matrix(0,nrow = rep,ncol = length(activex))
  pc21=matrix(0,nrow = rep,ncol = d0u)
  pc22=matrix(0,nrow = rep,ncol = d0u)
  pc23=matrix(0,nrow = rep,ncol = d0u)
  pc31=matrix(0,nrow = rep,ncol = length(activex))
  pc32=matrix(0,nrow = rep,ncol = length(activex))
  pc33=matrix(0,nrow = rep,ncol = length(activex))
  pc41=matrix(0,nrow = rep,ncol = length(activex))
  pc42=matrix(0,nrow = rep,ncol = length(activex))
  pc43=matrix(0,nrow = rep,ncol = length(activex))
  pc51=matrix(0,nrow = rep,ncol = length(activex))
  pc52=matrix(0,nrow = rep,ncol = length(activex))
  pc53=matrix(0,nrow = rep,ncol = length(activex))
  ec1=c()
  ec2=c()
  ec3=c()
  ec4=c()
  ec5=c()
  mms1=c()
  mms2=c()
  mms3=c()
  mms4=c()
  mms5=c()
  
  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re1<-GIGRSIS(cls=realda$y,atr = realda$x,method = SISmethodu,group = realda$group)
    re2<-IGSIS(cls=realda$y,atr = realda$x,method = SISmethodu)
    re3<-GIFSIS(cls=realda$y,atr = realda$x,method = SISmethodu,group = realda$group)
    re4<-GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gSIS")
    re5<-GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gHOLP")
    for (j in 1:length(activex)) {
      pc11[i,j]=ifelse(j %in% (re1[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc12[i,j]=ifelse(j %in% (re1[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc13[i,j]=ifelse(j %in% (re1[1:d3]),1,0)
    }
    ec1[i]=ifelse(sum(pc13[i,])==length(activex),1,0)
    
    mm1=c()
    for (w in 1:length(activex)) {
      mm1[w]=which(re1==activex[w])
    }
    mms1[i]=max(mm1)
    
    for (j in 1:d0u) {
      pc21[i,j]=ifelse(j %in% (re2[1:d1]),1,0)
    }
    for (j in 1:d0u) {
      pc22[i,j]=ifelse(j %in% (re2[1:d2]),1,0)
    }
    for (j in 1:d0u) {
      pc23[i,j]=ifelse(j %in% (re2[1:d3]),1,0)
    }
    ec2[i]=ifelse(sum(pc23[i,])==d0u,1,0)
    
    mm2=c()
    for (w in 1:d0u) {
      mm2[w]=which(re2==w)
    }
    mms2[i]=max(mm2)
    
    for (j in 1:length(activex)) {
      pc31[i,j]=ifelse(j %in% (re3[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc32[i,j]=ifelse(j %in% (re3[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc33[i,j]=ifelse(j %in% (re3[1:d3]),1,0)
    }
    ec3[i]=ifelse(sum(pc33[i,])==length(activex),1,0)
    
    mm3=c()
    for (w in 1:length(activex)) {
      mm3[w]=which(re3==activex[w])
    }
    mms3[i]=max(mm3)
    
    for (j in 1:length(activex)) {
      pc41[i,j]=ifelse(j %in% (re4[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc42[i,j]=ifelse(j %in% (re4[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc43[i,j]=ifelse(j %in% (re4[1:d3]),1,0)
    }
    ec4[i]=ifelse(sum(pc43[i,])==length(activex),1,0)
    
    mm4=c()
    for (w in 1:length(activex)) {
      mm4[w]=which(re4==activex[w])
    }
    mms4[i]=max(mm4)
    
    for (j in 1:length(activex)) {
      pc51[i,j]=ifelse(j %in% (re5[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc52[i,j]=ifelse(j %in% (re5[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc53[i,j]=ifelse(j %in% (re5[1:d3]),1,0)
    }
    ec5[i]=ifelse(sum(pc53[i,])==length(activex),1,0)
    
    mm5=c()
    for (w in 1:length(activex)) {
      mm5[w]=which(re5==activex[w])
    }
    mms5[i]=max(mm5)
    
  }
  repc11=mean(apply(pc11, 1, sum)/length(activex))
  repc12=mean(apply(pc12, 1, sum)/length(activex))
  repc13=mean(apply(pc13, 1, sum)/length(activex))
  repc21=mean(apply(pc21, 1, sum)/d0u)
  repc22=mean(apply(pc22, 1, sum)/d0u)
  repc23=mean(apply(pc23, 1, sum)/d0u)
  repc31=mean(apply(pc31, 1, sum)/length(activex))
  repc32=mean(apply(pc32, 1, sum)/length(activex))
  repc33=mean(apply(pc33, 1, sum)/length(activex))
  repc41=mean(apply(pc41, 1, sum)/length(activex))
  repc42=mean(apply(pc42, 1, sum)/length(activex))
  repc43=mean(apply(pc43, 1, sum)/length(activex))
  repc51=mean(apply(pc51, 1, sum)/length(activex))
  repc52=mean(apply(pc52, 1, sum)/length(activex))
  repc53=mean(apply(pc53, 1, sum)/length(activex))
  reec1=mean(ec1)
  reec2=mean(ec2)
  reec3=mean(ec3)
  reec4=mean(ec4)
  reec5=mean(ec5)
  remms1=quantile(mms1,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms2=quantile(mms2,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms3=quantile(mms3,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms4=quantile(mms4,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms5=quantile(mms5,probs = c(0.05,0.25,0.5,0.75,0.95))

  re=list()
  re$pc=rbind(c(repc11,repc12,repc13),c(repc21,repc22,repc23),c(repc31,repc32,repc33),c(repc41,repc42,repc43),c(repc51,repc52,repc53))
  re$ec=rbind(reec1,reec2,reec3,reec4,reec5)
  re$mms=rbind(remms1,remms2,remms3,remms4,remms5)
  return(re)
}

starttime=Sys.time()
result1=testone(rep = 2,Ju=4,Ru=2,Pu=1500,Nu=200,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result1,file = "s1-200-1500-b-g.csv")

starttime=Sys.time()
result2=testone(rep = 100,Ju=4,Ru=2,Pu=1500,Nu=100,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result2,file = "s1-100-1500-b-g.csv")
#
starttime=Sys.time()
result3=testone(rep = 100,Ju=4,Ru=2,Pu=1500,Nu=120,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3,file = "s1-120-1500-b-g.csv")
#
starttime=Sys.time()
result4=testone(rep = 100,Ju=4,Ru=2,Pu=1500,Nu=80,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result4,file = "s1-80-1500-ub-g.csv")
#
starttime=Sys.time()
result5=testone(rep = 100,Ju=4,Ru=2,Pu=1500,Nu=100,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result5,file = "s1-100-1500-ub-g.csv")
#
starttime=Sys.time()
result6=testone(rep = 100,Ju=4,Ru=2,Pu=1500,Nu=120,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result6,file = "s1-120-1500-ub-g.csv")

