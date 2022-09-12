# setwd("C:/Users/MyPC/Desktop/Model3")
##example 3: categorical and continuous covariates--------------------------------
stimu3 <- function(ind,JK,R,P,N,ymethod,p,j,d0){
  #ind=12
  #N=300;#600
  #P=5000;
  #d0=20;
  #R=4;
  #p: number of group 
  #j: number of predictors in each group
  if (ymethod=="balanced"){
    pr=rep(1/R,R)
    set.seed(ind)
    y=round(runif(N, min = 0.5, max = R+0.5))
  }else if (ymethod=="unbalanced"){
    indf=function(x,xn){
      n=length(xn)
      ind=NULL;
      for (j in 1:n) {
        ind[j]=which(x==xn[j])
      }
      return(ind)
    }
    
    y=NULL;
    r=c(1:R)
    pr2=2*(1+(R-c(1:R))/(R-1))/(3*R)
    
    indc=list()
    y=x=1:N
    N=length(x)
    num=round(pr2*N)
    num[R]=N-sum(num)+num[R]
    for (i in 1:R) {
      set.seed(ind*i)
      xn=sample(x,num[i])
      indc[[i]]=xn
      y[xn]=i
      x=x[-indf(x,xn)]
    }
  }
  
  theta1=c(0.2,0.8,0.7,0.2,0.2,0.9,0.1,0.1,0.7,0.7,0.3,0.5);
  theta2=c(0.9,0.3,0.3,0.7,0.8,0.4,0.7,0.6,0.4,0.1,0.8,0.2);
  theta3=c(0.1,0.9,0.6,0.1,0.3,0.1,0.4,0.3,0.6,0.4,0.4,0.7);
  theta4=c(0.7,0.2,0.1,0.6,0.7,0.6,0.8,0.9,0.1,0.8,0.8,0.6);
  theta=rbind(theta1,theta2,theta3,theta4);
  colnames(theta)=c(1,2,3,4,5,6,7,8,9,10,11,12);
  
  # group indices # seperate two parts #p is *4
  if (p*j<P) { p01=0.25*p;
  p02=0.25*P-0.25*p*j+0.25*p;
  p03=p02+p01;
  p04=0.5*P-0.5*p*j+0.5*p;
  p05=p04+0.5*p;
  p06=P-p*j+p;
  group1 <- c(rep(1:p01,each = j),c((p01+1):p02));
  group2 <- c(rep((p02+1):p03,each = j),c((p03+1):p04));
  group3 <- c(rep((p04+1):p05,each = j),c((p05+1):p06));
  group <-c(group1,group2,group3);
  } 
  
  if (p*j==P) {group <- c(rep(1:p,each = j));}
  act1=c(1:(0.25*d0))
  act2=c((0.25*P+1):(0.25*P+0.25*d0))
  act3=c((0.5*P+1):(0.5*P+0.5*d0))
  act=c(act1,act2,act3)
  
  #Bulid latent variable
  z=matrix(0,nrow = N,ncol = P)
  for (p in 1:P) {
    if (p %in% act){
      for (r in 1:R) {
        set.seed(r*p)
        number=length(which(y==r))
        mu=((-1)^r)*theta[r,which(act==p)]
        z[which(y==r),p]=rnorm(number,mu,1)
      }
    }
    if ((p %in% act)==FALSE){
      set.seed(p)
      z[,p]=rnorm(N,0,1)
    }
  }
  
  #Bulid x variable
  x=matrix(0,nrow = N,ncol = P)
  for (p in 1:P) {
    if (p<=0.25*P) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(0.25,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(0.25,0,1))&&(z[i,p]<=qnorm(0.5,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(0.5,0,1))&&(z[i,p]<=qnorm(0.75,0,1))) {x[i,p]=3}
        if (z[i,p]>qnorm(0.75,0,1)) {x[i,p]=4}
      }
    }
    #
    if ((p>0.25*P)&&(p<=0.5*P)) {
      for (i in 1:N) {
        if (z[i,p]<=qnorm(1/10,0,1)) {x[i,p]=1}
        if ((z[i,p]>qnorm(1/10,0,1))&&(z[i,p]<=qnorm(2/10,0,1))) {x[i,p]=2}
        if ((z[i,p]>qnorm(2/10,0,1))&&(z[i,p]<=qnorm(3/10,0,1))) {x[i,p]=3}
        if ((z[i,p]>qnorm(3/10,0,1))&&(z[i,p]<=qnorm(4/10,0,1))) {x[i,p]=4}
        if ((z[i,p]>qnorm(4/10,0,1))&&(z[i,p]<=qnorm(5/10,0,1))) {x[i,p]=5}
        if ((z[i,p]>qnorm(5/10,0,1))&&(z[i,p]<=qnorm(6/10,0,1))) {x[i,p]=6}
        if ((z[i,p]>qnorm(6/10,0,1))&&(z[i,p]<=qnorm(7/10,0,1))) {x[i,p]=7}
        if ((z[i,p]>qnorm(7/10,0,1))&&(z[i,p]<=qnorm(8/10,0,1))) {x[i,p]=8}
        if ((z[i,p]>qnorm(8/10,0,1))&&(z[i,p]<=qnorm(9/10,0,1))) {x[i,p]=9}
        if (z[i,p]>qnorm(9/10,0,1)) {x[i,p]=10}
      }
    }
    if (p>0.5*P){
      
      if (JK==2){
        for (i in 1:N) {
          x[i,p]=ifelse(z[i,p]<=qnorm(0.5,0,1),1,2) 
        }
      }
      if (JK!=2){
        for (i in 1:N) {
          if (z[i,p]<=qnorm(1/JK,0,1)) {x[i,p]=1}
          for (j in 2:(JK-1)) {
            x0=j-1;
            x1=j;
            if ((z[i,p]>qnorm(x0/JK,0,1))&&(z[i,p]<=qnorm(x1/JK,0,1))) {x[i,p]=j}
          }
          if (z[i,p]>qnorm((JK-1)/JK,0,1)) {x[i,p]=JK}
        }
      }
    }
    
  }
  
  re=list()
  re$y=y
  re$x=x
  re$group=group
  return(re)
}
realda=stimu3(ind=14,JK=8,R=4,P=3000,N=100,ymethod = "balanced",p=4,j=3,d0=12)
y=realda$y
x=realda$x
group=realda$group
source("GIF-SIS.R")
source("GIGR-SIS.R")
source("GRPSIS.R")
source("IG-SIS.R")

oldtestthree <- function(rep,Ju,Ru,Pu,Nu,ymethodu,SISmethodu,pu,ju,d0u){
  d1=floor(Nu/log(Nu,2));
  d2=2*d1;
  d3=3*d1;
  
  # group indices # seperate two parts #p is *4
  p=pu;j=ju;P=Pu;
  if (p*j<P) { p01=0.25*p;
  p02=0.25*P-0.25*p*j+0.25*p;
  p03=p02+p01;
  p04=0.5*P-0.5*p*j+0.5*p;
  p05=p04+0.5*p;
  p06=P-p*j+p;
  group1 <- c(rep(1:p01,each = j),c((p01+1):p02));
  group2 <- c(rep((p02+1):p03,each = j),c((p03+1):p04));
  group3 <- c(rep((p04+1):p05,each = j),c((p05+1):p06));
  group <-c(group1,group2,group3);
  act=c(c(1:p01),c((p02+1):p03),c((p04+1):p05))
  oact=c(c(1:(0.25*d0u)),c((0.25*P+1):(0.25*P+0.25*d0u)),c((0.5*P+1):(0.5*P+0.5*d0u)))
  } 
  
  if (p*j==P) {group <- c(rep(1:p,each = j));
  p0=d0u/ju;
  act=c(c(1:(0.25*p0)),c((0.25*p+1):(0.25*p+0.25*p0)),c((0.5*p+1):(0.5*p+0.5*p0)))
  oact=c(c(1:(0.25*d0u)),c((0.25*P+1):(0.25*P+0.25*d0u)),c((0.5*P+1):(0.5*P+0.5*d0u)))
  }
  activex=act#active covariables range
  #4
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
  #8
  pc61=matrix(0,nrow = rep,ncol = length(activex))
  pc62=matrix(0,nrow = rep,ncol = length(activex))
  pc63=matrix(0,nrow = rep,ncol = length(activex))
  pc71=matrix(0,nrow = rep,ncol = d0u)
  pc72=matrix(0,nrow = rep,ncol = d0u)
  pc73=matrix(0,nrow = rep,ncol = d0u)
  pc81=matrix(0,nrow = rep,ncol = length(activex))
  pc82=matrix(0,nrow = rep,ncol = length(activex))
  pc83=matrix(0,nrow = rep,ncol = length(activex))
  pc91=matrix(0,nrow = rep,ncol = length(activex))
  pc92=matrix(0,nrow = rep,ncol = length(activex))
  pc93=matrix(0,nrow = rep,ncol = length(activex))
  pc101=matrix(0,nrow = rep,ncol = length(activex))
  pc102=matrix(0,nrow = rep,ncol = length(activex))
  pc103=matrix(0,nrow = rep,ncol = length(activex))
  #10
  pc111=matrix(0,nrow = rep,ncol = length(activex))
  pc112=matrix(0,nrow = rep,ncol = length(activex))
  pc113=matrix(0,nrow = rep,ncol = length(activex))
  pc121=matrix(0,nrow = rep,ncol = d0u)
  pc122=matrix(0,nrow = rep,ncol = d0u)
  pc123=matrix(0,nrow = rep,ncol = d0u)
  pc131=matrix(0,nrow = rep,ncol = length(activex))
  pc132=matrix(0,nrow = rep,ncol = length(activex))
  pc133=matrix(0,nrow = rep,ncol = length(activex))
  pc141=matrix(0,nrow = rep,ncol = length(activex))
  pc142=matrix(0,nrow = rep,ncol = length(activex))
  pc143=matrix(0,nrow = rep,ncol = length(activex))
  pc151=matrix(0,nrow = rep,ncol = length(activex))
  pc152=matrix(0,nrow = rep,ncol = length(activex))
  pc153=matrix(0,nrow = rep,ncol = length(activex))
  
  ec1=c()
  ec2=c()
  ec3=c()
  ec4=c()
  ec5=c()
  ec6=c()
  ec7=c()
  ec8=c()
  ec9=c()
  ec10=c()
  ec11=c()
  ec12=c()
  ec13=c()
  ec14=c()
  ec15=c()

 
  
  mms1=c()
  mms2=c()
  mms3=c()
  mms4=c()
  mms5=c()
  mms6=c()
  mms7=c()
  mms8=c()
  mms9=c()
  mms10=c()
  mms11=c()
  mms12=c()
  mms13=c()
  mms14=c()
  mms15=c()

  
  for (i in 1:rep) {
    cat(i,",")
    realda1=stimu3(ind = i,Ju[1],Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    realda2=stimu3(ind = i,Ju[2],Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    realda3=stimu3(ind = i,Ju[3],Ru,Pu,Nu,ymethodu,pu,ju,d0u)
    re1<-GIGRSIS3(cls = realda1$y,atr = realda1$x,method = "inforate",group=realda1$group)
    re2<-IGSIS(realda1$y,realda1$x,method = SISmethodu)
    re3<-GIFSIS3(realda1$y,realda1$x,method = SISmethodu,group = realda1$group)
    re4<-GRPSIS(X = realda1$x,y=realda1$y,group = realda1$group,"gSIS")
    re5<-GRPSIS(X = realda1$x,y=realda1$y,group = realda1$group,"gHOLP")
    re6<-GIGRSIS3(cls = realda2$y,atr = realda2$x,method = "inforate",group=realda2$group)
    re7<-IGSIS(realda2$y,realda2$x,method = SISmethodu)
    re8<-GIFSIS3(realda2$y,realda2$x,method = SISmethodu,group = realda2$group)
    re9<-GRPSIS(X = realda2$x,y=realda2$y,group = realda2$group,"gSIS")
    re10<-GRPSIS(X = realda2$x,y=realda2$y,group = realda2$group,"gHOLP")
    re11<-GIGRSIS3(cls = realda3$y,atr = realda3$x,method = "inforate",group=realda3$group)
    re12<-IGSIS(realda3$y,realda3$x,method = SISmethodu)
    re13<-GIFSIS3(realda3$y,realda3$x,method = SISmethodu,group = realda3$group)
    re14<-GRPSIS(X = realda3$x,y=realda3$y,group = realda3$group,"gSIS")
    re15<-GRPSIS(X = realda3$x,y=realda3$y,group = realda3$group,"gHOLP")
    
    ###4
    for (j in 1:length(activex)) {
      pc11[i,j]=ifelse(activex[j] %in% (re1[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc12[i,j]=ifelse(activex[j] %in% (re1[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc13[i,j]=ifelse(activex[j] %in% (re1[1:d3]),1,0)
    }
    ec1[i]=ifelse(sum(pc13[i,])==length(activex),1,0)
    
    mm1=c()
    for (w in 1:length(activex)) {
      mm1[w]=which(re1==activex[w])
    }
    mms1[i]=max(mm1)
    
    for (j in 1:d0u) {
      pc21[i,j]=ifelse(oact[j] %in% (re2[1:d1]),1,0)
    }
    for (j in 1:d0u) {
      pc22[i,j]=ifelse(oact[j] %in% (re2[1:d2]),1,0)
    }
    for (j in 1:d0u) {
      pc23[i,j]=ifelse(oact[j] %in% (re2[1:d3]),1,0)
    }
    ec2[i]=ifelse(sum(pc23[i,])==d0u,1,0)
    
    mm2=c()
    for (w in 1:d0u) {
      mm2[w]=which(re2==oact[w])
    }
    mms2[i]=max(mm2)
    
    for (j in 1:length(activex)) {
      pc31[i,j]=ifelse(activex[j] %in% (re3[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc32[i,j]=ifelse(activex[j] %in% (re3[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc33[i,j]=ifelse(activex[j] %in% (re3[1:d3]),1,0)
    }
    ec3[i]=ifelse(sum(pc33[i,])==length(activex),1,0)
    
    mm3=c()
    for (w in 1:length(activex)) {
      mm3[w]=which(re3==activex[w])
    }
    mms3[i]=max(mm3)
    
    for (j in 1:length(activex)) {
      pc41[i,j]=ifelse(activex[j] %in% (re4[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc42[i,j]=ifelse(activex[j] %in% (re4[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc43[i,j]=ifelse(activex[j] %in% (re4[1:d3]),1,0)
    }
    ec4[i]=ifelse(sum(pc43[i,])==length(activex),1,0)
    
    mm4=c()
    for (w in 1:length(activex)) {
      mm4[w]=which(re4==activex[w])
    }
    mms4[i]=max(mm4)
    
    for (j in 1:length(activex)) {
      pc51[i,j]=ifelse(activex[j] %in% (re5[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc52[i,j]=ifelse(activex[j] %in% (re5[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc53[i,j]=ifelse(activex[j] %in% (re5[1:d3]),1,0)
    }
    ec5[i]=ifelse(sum(pc53[i,])==length(activex),1,0)
    
    mm5=c()
    for (w in 1:length(activex)) {
      mm5[w]=which(re5==activex[w])
    }
    mms5[i]=max(mm5)
    
    ###8
    for (j in 1:length(activex)) {
      pc61[i,j]=ifelse(activex[j] %in% (re6[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc62[i,j]=ifelse(activex[j] %in% (re6[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc63[i,j]=ifelse(activex[j] %in% (re6[1:d3]),1,0)
    }
    ec6[i]=ifelse(sum(pc63[i,])==length(activex),1,0)
    
    mm6=c()
    for (w in 1:length(activex)) {
      mm6[w]=which(re6==activex[w])
    }
    mms6[i]=max(mm6)
    
    for (j in 1:d0u) {
      pc71[i,j]=ifelse(oact[j] %in% (re7[1:d1]),1,0)
    }
    for (j in 1:d0u) {
      pc72[i,j]=ifelse(oact[j] %in% (re7[1:d2]),1,0)
    }
    for (j in 1:d0u) {
      pc73[i,j]=ifelse(oact[j] %in% (re7[1:d3]),1,0)
    }
    ec7[i]=ifelse(sum(pc73[i,])==d0u,1,0)
    
    mm7=c()
    for (w in 1:d0u) {
      mm7[w]=which(re7==oact[w])
    }
    mms7[i]=max(mm7)
    
    for (j in 1:length(activex)) {
      pc81[i,j]=ifelse(activex[j] %in% (re8[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc82[i,j]=ifelse(activex[j] %in% (re8[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc83[i,j]=ifelse(activex[j] %in% (re8[1:d3]),1,0)
    }
    ec8[i]=ifelse(sum(pc83[i,])==length(activex),1,0)
    
    mm8=c()
    for (w in 1:length(activex)) {
      mm8[w]=which(re8==activex[w])
    }
    mms8[i]=max(mm8)
    
    for (j in 1:length(activex)) {
      pc91[i,j]=ifelse(activex[j] %in% (re9[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc92[i,j]=ifelse(activex[j] %in% (re9[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc93[i,j]=ifelse(activex[j] %in% (re9[1:d3]),1,0)
    }
    ec9[i]=ifelse(sum(pc93[i,])==length(activex),1,0)
    
    mm9=c()
    for (w in 1:length(activex)) {
      mm9[w]=which(re9==activex[w])
    }
    mms9[i]=max(mm9)
    
    for (j in 1:length(activex)) {
      pc101[i,j]=ifelse(activex[j] %in% (re10[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc102[i,j]=ifelse(activex[j] %in% (re10[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc103[i,j]=ifelse(activex[j] %in% (re10[1:d3]),1,0)
    }
    ec10[i]=ifelse(sum(pc103[i,])==length(activex),1,0)
    
    mm10=c()
    for (w in 1:length(activex)) {
      mm10[w]=which(re10==activex[w])
    }
    mms10[i]=max(mm10)
    
    
    ###10
    for (j in 1:length(activex)) {
      pc111[i,j]=ifelse(activex[j] %in% (re11[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc112[i,j]=ifelse(activex[j] %in% (re11[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc113[i,j]=ifelse(activex[j] %in% (re11[1:d3]),1,0)
    }
    ec11[i]=ifelse(sum(pc113[i,])==length(activex),1,0)
    
    mm11=c()
    for (w in 1:length(activex)) {
      mm11[w]=which(re11==activex[w])
    }
    mms11[i]=max(mm11)
    
    for (j in 1:d0u) {
      pc121[i,j]=ifelse(oact[j] %in% (re12[1:d1]),1,0)
    }
    for (j in 1:d0u) {
      pc122[i,j]=ifelse(oact[j] %in% (re12[1:d2]),1,0)
    }
    for (j in 1:d0u) {
      pc123[i,j]=ifelse(oact[j] %in% (re12[1:d3]),1,0)
    }
    ec12[i]=ifelse(sum(pc123[i,])==d0u,1,0)
    
    mm12=c()
    for (w in 1:d0u) {
      mm12[w]=which(re12==oact[w])
    }
    mms12[i]=max(mm12)
    
    for (j in 1:length(activex)) {
      pc131[i,j]=ifelse(activex[j] %in% (re13[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc132[i,j]=ifelse(activex[j] %in% (re13[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc133[i,j]=ifelse(activex[j] %in% (re13[1:d3]),1,0)
    }
    ec13[i]=ifelse(sum(pc133[i,])==length(activex),1,0)
    
    mm13=c()
    for (w in 1:length(activex)) {
      mm13[w]=which(re13==activex[w])
    }
    mms13[i]=max(mm13)
    
    for (j in 1:length(activex)) {
      pc141[i,j]=ifelse(activex[j] %in% (re14[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc142[i,j]=ifelse(activex[j] %in% (re14[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc143[i,j]=ifelse(activex[j] %in% (re14[1:d3]),1,0)
    }
    ec14[i]=ifelse(sum(pc143[i,])==length(activex),1,0)
    
    mm14=c()
    for (w in 1:length(activex)) {
      mm14[w]=which(re14==activex[w])
    }
    mms14[i]=max(mm14)
    
    for (j in 1:length(activex)) {
      pc151[i,j]=ifelse(activex[j] %in% (re15[1:d1]),1,0)
    }
    for (j in 1:length(activex)) {
      pc152[i,j]=ifelse(activex[j] %in% (re15[1:d2]),1,0)
    }
    for (j in 1:length(activex)) {
      pc153[i,j]=ifelse(activex[j] %in% (re15[1:d3]),1,0)
    }
    ec15[i]=ifelse(sum(pc153[i,])==length(activex),1,0)
    
    mm15=c()
    for (w in 1:length(activex)) {
      mm15[w]=which(re15==activex[w])
    }
    mms15[i]=max(mm15)
    
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
  
  repc61=mean(apply(pc61, 1, sum)/length(activex))
  repc62=mean(apply(pc62, 1, sum)/length(activex))
  repc63=mean(apply(pc63, 1, sum)/length(activex))
  repc71=mean(apply(pc71, 1, sum)/d0u)
  repc72=mean(apply(pc72, 1, sum)/d0u)
  repc73=mean(apply(pc73, 1, sum)/d0u)
  repc81=mean(apply(pc81, 1, sum)/length(activex))
  repc82=mean(apply(pc82, 1, sum)/length(activex))
  repc83=mean(apply(pc83, 1, sum)/length(activex))
  repc91=mean(apply(pc91, 1, sum)/length(activex))
  repc92=mean(apply(pc92, 1, sum)/length(activex))
  repc93=mean(apply(pc93, 1, sum)/length(activex))
  repc101=mean(apply(pc101, 1, sum)/length(activex))
  repc102=mean(apply(pc102, 1, sum)/length(activex))
  repc103=mean(apply(pc103, 1, sum)/length(activex))
  
  repc111=mean(apply(pc111, 1, sum)/length(activex))
  repc112=mean(apply(pc112, 1, sum)/length(activex))
  repc113=mean(apply(pc113, 1, sum)/length(activex))
  repc121=mean(apply(pc121, 1, sum)/d0u)
  repc122=mean(apply(pc122, 1, sum)/d0u)
  repc123=mean(apply(pc123, 1, sum)/d0u)
  repc131=mean(apply(pc131, 1, sum)/length(activex))
  repc132=mean(apply(pc132, 1, sum)/length(activex))
  repc133=mean(apply(pc133, 1, sum)/length(activex))
  repc141=mean(apply(pc141, 1, sum)/length(activex))
  repc142=mean(apply(pc142, 1, sum)/length(activex))
  repc143=mean(apply(pc143, 1, sum)/length(activex))
  repc151=mean(apply(pc151, 1, sum)/length(activex))
  repc152=mean(apply(pc152, 1, sum)/length(activex))
  repc153=mean(apply(pc153, 1, sum)/length(activex))
  
  reec1=mean(ec1)
  reec2=mean(ec2)
  reec3=mean(ec3)
  reec4=mean(ec4)
  reec5=mean(ec5)
  reec6=mean(ec6)
  reec7=mean(ec7)
  reec8=mean(ec8)
  reec9=mean(ec9)
  reec10=mean(ec10)
  reec11=mean(ec11)
  reec12=mean(ec12)
  reec13=mean(ec13)
  reec14=mean(ec14)
  reec15=mean(ec15)

  remms1=quantile(mms1,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms2=quantile(mms2,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms3=quantile(mms3,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms4=quantile(mms4,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms5=quantile(mms5,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms6=quantile(mms6,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms7=quantile(mms7,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms8=quantile(mms8,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms9=quantile(mms9,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms10=quantile(mms10,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms11=quantile(mms11,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms12=quantile(mms12,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms13=quantile(mms13,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms14=quantile(mms14,probs = c(0.05,0.25,0.5,0.75,0.95))
  remms15=quantile(mms15,probs = c(0.05,0.25,0.5,0.75,0.95))

  
  re=list()
  re$pc=rbind(c(repc11,repc12,repc13),c(repc21,repc22,repc23),
              c(repc31,repc32,repc33),c(repc41,repc42,repc43),
              c(repc51,repc52,repc53),c(repc61,repc62,repc63),c(repc71,repc72,repc73),c(repc81,repc82,repc83),c(repc91,repc92,repc93),c(repc101,repc102,repc103),c(repc111,repc112,repc113),c(repc121,repc122,repc123),c(repc131,repc132,repc133),c(repc141,repc142,repc143),c(repc151,repc152,repc153))
  re$ec=rbind(reec1,reec2,reec3,reec4,reec5,reec6,reec7,reec8,reec9,reec10,reec11,reec12,reec13,reec14,reec15)
  re$mms=rbind(remms1,remms2,remms3,remms4,remms5,remms6,remms7,remms8,remms9,remms10,remms11,remms12,remms13,remms14,remms15)
  return(re)
}
starttime=Sys.time()
result3_1=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=180,ymethodu = "balanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_1,file = "s3GIGSIS-180-3000-b-gn.csv")

starttime=Sys.time()
result3_2=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=200,ymethodu = "balanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_2,file = "s3GIGSIS-200-3000-b-gn.csv")

starttime=Sys.time()
result3_3=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=220,ymethodu = "balanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_3,file = "s3GIGSIS-220-3000-b-gn.csv")

starttime=Sys.time()
result3_4=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=180,ymethodu = "unbalanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_4,file = "s3GIGSIS-180-3000-ub-gn.csv")

starttime=Sys.time()
result3_5=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=200,ymethodu = "unbalanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_5,file = "s3GIGSIS-200-3000-ub-gn.csv")

starttime=Sys.time()
result3_6=oldtestthree(rep = 100,Ju=c(4,8,10),Ru=4,Pu=3000,Nu=220,ymethodu = "unbalanced",SISmethodu = "inforate",pu=4,ju=3,d0u=12)
endtime=Sys.time()
procestime=endtime-starttime
write.csv(result3_6,file = "s3GIGSIS-220-3000-ub-gn.csv")