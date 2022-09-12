# setwd("C:/Users/MyPC/Desktop/Model2")
##example 2: categorical covariates and two response variable
stimu2 <- function(ind,J,R,P,N,ymethod,p,j,d0){
  #p: number of group 
  #j: number of predictors in each group
  if (ymethod=="balanced"){
    pr=rep(1/R,R)
    set.seed(ind)
    y=round(runif(N, min = 0.5, max = R+0.5))
  }else if (ymethod=="unbalanced"){
    pr2=2*(1+(R-c(1:R))/(R-1))/(3*R)
    set.seed(ind)
    y=NULL;
    for (i in 1:N) {
      if (i<=round(pr2[1]*N)) {y[i]=1}
      for (j in 2:(R-1)) {
        x0=j-1
        x1=j
        if ((i>round(sum(pr2[1:x0]*N)))&&(i<=round(sum(pr2[1:x1]*N)))) {y[i]=x1}
      }
      
      if (i>round(sum(pr2[1:(R-1)]*N))) {y[i]=R}
    }
  }
  
  # group indices
  if (p*j<P) {group <- c(rep(1:p,each = j),c((p+1):(P-p*j+p))) } 
  if (p*j==P) {group <- c(rep(1:p,each = j))}
  
  mu=NULL;
  for (r in 1:R) {
    mu[r]=1.5*(-0.9)^r
  }
  
  #Bulid latent variable
  z=matrix(0,nrow = N,ncol = P)
  for (p in 1:P) {
    if (p<=d0) {
      for (r in 1:R) {
        set.seed(p*r+2000)
        #z[which(y==r),p]=rnorm(length(which(y==r)),0,1)+mu[r];
        z[which(y==r),p]=rt(length(which(y==r)),4)+mu[r];
      }
    }
    if (p>d0) {
      set.seed(p)
      #z[,p]=rnorm(N,0,1)
      z[,p]=rt(N,4)
    }
  }
  
  #Bulid x variable
  x=matrix(0,nrow = N,ncol = P)
  for (p in 1:P) {
    if (p<=400) {x[,p]=ifelse(z[,p]>qnorm(0.5,0,1),2,1)}
    if (p<=800&&p>400) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(0.25,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(0.25,0,1))&&(z[i,p]<=qnorm(0.5,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(0.5,0,1))&&(z[i,p]<=qnorm(0.75,0,1))) {x[i,p]=3}
        if (z[i,p]>qnorm(0.75,0,1)) {x[i,p]=4}
      }
    }
    if (p<=1200&&p>800) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(1/6,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(1/6,0,1))&&(z[i,p]<=qnorm(2/6,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(2/6,0,1))&&(z[i,p]<=qnorm(3/6,0,1))) {x[i,p]=3}
        if ((z[i,p]>qnorm(3/6,0,1))&&(z[i,p]<=qnorm(4/6,0,1))) {x[i,p]=4}
        if ((z[i,p]>qnorm(4/6,0,1))&&(z[i,p]<=qnorm(5/6,0,1))) {x[i,p]=5}
        if (z[i,p]>qnorm(5/6,0,1)) {x[i,p]=6}
      }
    }
    if (p<=1600&&p>1200) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(1/8,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(1/8,0,1))&&(z[i,p]<=qnorm(2/8,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(2/8,0,1))&&(z[i,p]<=qnorm(3/8,0,1))) {x[i,p]=3}
        if ((z[i,p]>qnorm(3/8,0,1))&&(z[i,p]<=qnorm(4/8,0,1))) {x[i,p]=4}
        if ((z[i,p]>qnorm(4/8,0,1))&&(z[i,p]<=qnorm(5/8,0,1))) {x[i,p]=5}
        if ((z[i,p]>qnorm(5/8,0,1))&&(z[i,p]<=qnorm(6/8,0,1))) {x[i,p]=6}
        if ((z[i,p]>qnorm(6/8,0,1))&&(z[i,p]<=qnorm(7/8,0,1))) {x[i,p]=7}
        if (z[i,p]>qnorm(7/8,0,1)) {x[i,p]=8}
      }
    }
    if (p>1600) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(1/10,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(1/10,0,1))&&(z[i,p]<=qnorm(2/10,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(2/10,0,1))&&(z[i,p]<=qnorm(3/10,0,1))) {x[i,p]=3}
        if ((z[i,p]>qnorm(3/10,0,1))&&(z[i,p]<=qnorm(4/10,0,1))) {x[i,p]=4}
        if ((z[i,p]>qnorm(4/10,0,1))&&(z[i,p]<=qnorm(5/10,0,1))) {x[i,p]=5}
        if ((z[i,p]>qnorm(5/10,0,1))&&(z[i,p]<=qnorm(6/10,0,1))) {x[i,p]=6}
        if ((z[i,p]>qnorm(6/10,0,1))&&(z[i,p]<=qnorm(7/10,0,1))) {x[i,p]=7}
        if ((z[i,p]>qnorm(7/10,0,1))&&(z[i,p]<=qnorm(6/10,0,1))) {x[i,p]=8}
        if ((z[i,p]>qnorm(8/10,0,1))&&(z[i,p]<=qnorm(7/10,0,1))) {x[i,p]=9}
        if (z[i,p]>qnorm(9/10,0,1)) {x[i,p]=10}
      }
    }
  }
  re=list()
  re$y=y
  re$x=x
  re$group=group
  return(re)
}
realda2=stimu2(ind=1,R=4,P=2000,N=200,ymethod = "balanced",p=4,j=3,d0=12)
y=realda2$y
x=realda2$x
group=realda2$group
source("GIF-SIS.R")
source("GIGR-SIS.R")
source("GRPSIS.R")
source("IG-SIS.R")

GIGRSIS2(R=4,cls=y,atr = x,method = "inforate",group = group)#GIGRSIS
IGSIS(cls = realda$y,atr = realda$x,method = "inforate")#IGSIS
GIFSIS2(R=4,cls = y,atr = x,method = "inforate",group = group)#GIFSIS
GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gSIS")#gSIS
GRPSIS(X = realda$x,y=realda$y,group = realda$group,"gHOLP")#gHOLP


testtwo <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
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
    realda1=stimu2(ind = i,Ju,Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re1<-GIGRSIS2(R=Ru,cls=realda1$y,atr = realda1$x,method = SISmethodu,group = realda1$group)
    re2<-IGSIS(cls=realda$y,atr = realda$x,method = SISmethodu)
    re3<-GIFSIS2(R=Ru,cls=realda1$y,atr = realda1$x,method = SISmethodu,group = realda$group)
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
result21=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=100,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result21,file = "s2-100-2000-b-g-33s.csv")
# 
starttime=Sys.time()
result22=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=150,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result22,file = "s2-150-2000-b-g-33s.csv")
# 
starttime=Sys.time()
result23=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=200,ymethodu = "balanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result23,file = "s2-200-2000-b-g-33s.csv")

starttime=Sys.time()
result210=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=100,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result210,file = "s2-100-2000-ub-ug16.csv")
# # 
starttime=Sys.time()
result211=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=150,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result211,file = "s2-100-2000-ub-ug16.csv")

starttime=Sys.time()
result212=testtwo(rep = 100,Ju=4,Ru=10,Pu=2000,Nu=200,ymethodu = "unbalanced",SISmethodu = "inforate",pu=5,ju=3,d0u=15)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result212,file = "s2-200-2000-ub-ug.csv")