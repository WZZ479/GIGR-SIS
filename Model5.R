# setwd("C:/Users/MyPC/Desktop/Model5")
library(microbenchmark)#Calculate running time
##example 5: categorical covariates and binary response variable
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

source("GIF-SIS.R")
source("GIGR-SIS.R")
source("GRPSIS.R")
source("IG-SIS.R")

testone1 <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  pc11=matrix(0,nrow = rep,ncol = length(activex))
  pc12=matrix(0,nrow = rep,ncol = length(activex))
  pc13=matrix(0,nrow = rep,ncol = length(activex))
  ec1=c()
  mms1=c()
  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re1<-GIGRSIS(cls=realda$y,atr = realda$x,method = SISmethodu,group = realda$group)
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
  }
  repc11=mean(apply(pc11, 1, sum)/length(activex))
  repc12=mean(apply(pc12, 1, sum)/length(activex))
  repc13=mean(apply(pc13, 1, sum)/length(activex))
  reec1=mean(ec1)
  remms1=quantile(mms1,probs = c(0.05,0.25,0.5,0.75,0.95))
  re=list()
  re$pc=rbind(c(repc11,repc12,repc13))
  re$ec=rbind(reec1)
  re$mms=rbind(remms1)
  return(re)
}
testone2 <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  pc21=matrix(0,nrow = rep,ncol = d0u)
  pc22=matrix(0,nrow = rep,ncol = d0u)
  pc23=matrix(0,nrow = rep,ncol = d0u)
  ec2=c()
  mms2=c()
  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re2<-IGSIS(cls=realda$y,atr = realda$x,method = SISmethodu)
    
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
    # 
    
  }
  repc21=mean(apply(pc21, 1, sum)/d0u)
  repc22=mean(apply(pc22, 1, sum)/d0u)
  repc23=mean(apply(pc23, 1, sum)/d0u)
  reec2=mean(ec2)
  remms2=quantile(mms2,probs = c(0.05,0.25,0.5,0.75,0.95))
  re=list()
  re$pc=rbind(c(repc21,repc22,repc23))
  re$ec=rbind(reec2)
  re$mms=rbind(remms2)
  return(re)
}
testone3 <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  pc31=matrix(0,nrow = rep,ncol = length(activex))
  pc32=matrix(0,nrow = rep,ncol = length(activex))
  pc33=matrix(0,nrow = rep,ncol = length(activex))
  ec3=c()
  mms3=c()

  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re3<-GIFSIS(cls=realda$y,atr = realda$x,method = SISmethodu,group = realda$group)

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
  }
  repc31=mean(apply(pc31, 1, sum)/length(activex))
  repc32=mean(apply(pc32, 1, sum)/length(activex))
  repc33=mean(apply(pc33, 1, sum)/length(activex))
  reec3=mean(ec3)
  remms3=quantile(mms3,probs = c(0.05,0.25,0.5,0.75,0.95))
  
  re=list()
  re$pc=rbind(c(repc31,repc32,repc33))
  re$ec=rbind(reec3)
  re$mms=rbind(remms3)
  return(re)
}
testone4 <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  pc41=matrix(0,nrow = rep,ncol = length(activex))
  pc42=matrix(0,nrow = rep,ncol = length(activex))
  pc43=matrix(0,nrow = rep,ncol = length(activex))
  ec4=c()
  mms4=c()
  
  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re4<-GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gSIS")
    # 
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
  }
  repc41=mean(apply(pc41, 1, sum)/length(activex))
  repc42=mean(apply(pc42, 1, sum)/length(activex))
  repc43=mean(apply(pc43, 1, sum)/length(activex))
  reec4=mean(ec4)
  remms4=quantile(mms4,probs = c(0.05,0.25,0.5,0.75,0.95))
  
  re=list()
  re$pc=rbind(c(repc41,repc42,repc43))
  re$ec=rbind(reec4)
  re$mms=rbind(remms4)
  return(re)
}
testone5 <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  
  if (pu*ju<Pu) {group <- c(rep(1:pu,each = ju),c((pu+1):(Pu-pu*ju+pu))) } 
  if (pu*ju==Pu) {group <- c(rep(1:pu,each = ju))}
  activex=unique(group[1:d0u])
  pc51=matrix(0,nrow = rep,ncol = length(activex))
  pc52=matrix(0,nrow = rep,ncol = length(activex))
  pc53=matrix(0,nrow = rep,ncol = length(activex))
  ec5=c()
  mms5=c()
  
  for (i in 1:rep) {
    cat(i,",")
    realda=stimu1(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re5<-GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gHOLP")
    # 
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
  repc51=mean(apply(pc51, 1, sum)/length(activex))
  repc52=mean(apply(pc52, 1, sum)/length(activex))
  repc53=mean(apply(pc53, 1, sum)/length(activex))
  reec5=mean(ec5)
  remms5=quantile(mms5,probs = c(0.05,0.25,0.5,0.75,0.95))
  
  re=list()
  re$pc=rbind(c(repc51,repc52,repc53))
  re$ec=rbind(reec5)
  re$mms=rbind(remms5)
  return(re)
}
#n=150,P=1500
result1<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=1500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=1500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=1500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=1500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=1500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq     mean   median       uq      max    neval
# GIGR-SIS 3.289554 3.330623 3.368177 3.353802 3.375225 3.537076   100
# IG-SIS   2.881828 2.919752 2.941701 2.929888 2.952077 3.117824   100
# GIG-SIS  3.309124 3.365101 3.395856 3.382953 3.404948 3.570489   100
# gSIS     1.470090 1.507313 1.533108 1.518816 1.535622 1.719818   100
# gHOLP    1.557598 1.587965 1.614787 1.601231 1.621020 1.792080   100
#n=150,P=2500
result2<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=2500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=2500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=2500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=2500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=2500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq     mean   median       uq      max    neval
# GIGR-SIS 5.494299 5.570217 5.659486 5.672852 5.696671 5.868789   100
# IG-SIS   4.770752 4.852069 4.910827 4.901412 4.964495 5.117287   100
# GIG-SIS  5.513120 5.627651 5.687498 5.690967 5.731353 5.921356   100
# gSIS     2.427813 2.502660 2.543500 2.534560 2.558392 2.717746   100
# gHOLP    2.579073 2.642461 2.674566 2.669960 2.694595 2.858579   100
#n=150,P=3500
result3<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=3500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=3500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=3500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=3500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=3500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq     mean   median       uq      max    neval
# GIGR-SIS 7.699115 7.873410 7.947772 7.973713 7.996539 8.253035   100
# IG-SIS   6.717779 6.813813 6.928664 6.952105 6.982997 7.174963   100
# GIG-SIS  7.759647 7.941331 8.023771 8.028897 8.063906 8.401624   100
# gSIS     3.436946 3.514005 3.566716 3.550886 3.596695 3.878102   100
# gHOLP    3.608111 3.716761 3.763079 3.753511 3.785692 3.994942   100
#n=150,P=4500
result4<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=4500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=4500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=4500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=4500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=4500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq       mean      median       uq        max     neval
# GIGR-SIS 9.952425  10.166974 10.296161 10.320424  10.379040  10.550945   100
# IG-SIS   8.667074  8.879759  8.924028  8.927287   9.007771   9.178570    100
# GIG-SIS  10.059064 10.286687 10.383764 10.399807  10.467705  10.652790   100
# gSIS     4.483718  4.543444  4.593066  4.573370   4.627062   4.841683    100
# gHOLP    4.679809  4.786540  4.829767  4.815617   4.854455   5.044709    100
#n=150,P=5500
result5<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=5500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=5500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=5500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=5500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=5500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq        mean      median       uq      max     neval
# GIGR-SIS  12.466047 12.686935 12.741751 12.725734 12.840434 12.970272   100
# IG-SIS    10.719647 10.847841 10.936060 10.930933 11.031558 11.160138   100
# GIG-SIS   12.557232 12.704236 12.802820 12.809088 12.878315 13.068262   100
# gSIS      5.508760  5.583007  5.633349  5.606822  5.665670  5.808685    100
# gHOLP     5.814848  5.892932  5.941496  5.925079  5.984452  6.121633    100
#n=150,P=6500
result6<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=6500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=6500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=6500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=6500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=6500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq        mean     median       uq      max     neval
# GIGR-SIS 14.864796 14.994788 15.102586 15.074988 15.190366 15.394388   100
# IG-SIS   12.640783 12.774058 12.863136 12.860302 12.935702 13.122587   100
# GIG-SIS  14.972075 15.098336 15.203301 15.186373 15.278473 15.602892   100
# gSIS     6.477899  6.583475  6.643575  6.625081  6.687297  6.848894    100
# gHOLP    6.816021  6.950912  7.013216  7.006791  7.071962  7.249152    100
#n=150,P=7500
result7<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=7500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=7500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=7500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=7500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=7500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq        mean      median       uq      max       neval
# GIGR-SIS  17.048138 17.392129 17.581368 17.636154 17.765857   17.903361   100
# IG-SIS   14.390048 14.644654 14.841824  14.914827  15.008665  15.144522   100
# GIG-SIS  17.156993 17.516188 17.693176  17.760998  17.871590  18.041020   100
# gSIS     7.388287  7.584858  7.667323   7.662402   7.758945   7.914610    100
# gHOLP    7.789030  8.030974  8.093217   8.093088   8.162005   8.345568    100
#n=150,P=8500
result8<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=8500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=8500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=8500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=8500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=8500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq        mean     median       uq      max     neval
# GIGR-SIS 19.303486 19.451828 19.540657 19.533410 19.611892 19.956965   100
# IG-SIS   16.316742 16.459391 16.550965 16.543809 16.626750 16.954287   100
# GIG-SIS  19.417523 19.575395 19.676734 19.660238 19.755665 20.128257   100
# gSIS     8.329623  8.476473  8.550979  8.556243  8.642435  8.818589    100
# gHOLP    8.792883  8.948956  9.032511  9.037577  9.100933  9.364427    100
#n=150,P=9500
result9<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=9500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=9500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=9500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=9500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=9500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq       mean      median       uq      max     neval
# GIGR-SIS 21.578171 21.739334 21.852565 21.837005 21.934978 22.303793   100
# IG-SIS   18.210729 18.432923 18.514709 18.512405 18.613995 19.028119   100
# GIG-SIS  21.698115 21.875557 21.973475 21.972164 22.068806 22.559843   100
# gSIS     9.325462  9.509575  9.592189  9.587706  9.676351  9.861513    100
# gHOLP    9.879014 10.046571  10.140044 10.115457 10.221710 10.561903   100
#n=150,P=10500
result10<-microbenchmark::microbenchmark(testone1(rep = 1,Ju=4,Ru=2,Pu=10500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone2(rep = 1,Ju=4,Ru=2,Pu=10500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone3(rep = 1,Ju=4,Ru=2,Pu=10500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone4(rep = 1,Ju=4,Ru=2,Pu=10500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),testone5(rep = 1,Ju=4,Ru=2,Pu=10500,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15),times = 100)
#             min       lq     mean   median       uq      max    neval
# GIGR-SIS 23.88069 24.14105 24.26169 24.24170 24.36945 24.64091   100
# IG-SIS   20.24643 20.46833 20.58509 20.57952 20.69359 20.90776   100
# GIG-SIS  24.05017 24.28460 24.41086 24.40391 24.53103 24.78011   100
# gSIS     10.39366 10.56367 10.64973 10.64932 10.75677 11.03488   100
# gHOLP    10.95037 11.17037 11.27016 11.26927 11.35170 11.64168   100
# 