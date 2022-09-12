#Group feature screening based on information gain ratio
GInfoGainmul <- function(cls,atr,method=c("info","inforate"),group){
  #cls;Y
  #atr:X
  #One dimensional empirical entropy function
  HDfuncs <- function(atrcls){#aatrcls is vector
    I=length(atrcls)
    tatrcls=table(atrcls)
    atrclspvec=as.vector(tatrcls)/I#
    logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
    HD=-as.vector(atrclspvec%*%logatrclspvec)
    return(HD)
  }
  #One dimensional empirical conditional entropy function
  Hatrs <- function(atrvec,clsvec){#Input as feature vector and category vector
    n=length(atrvec)
    tatr=table(atrvec)
    atrpvec=as.vector(tatr)/n
    txy=table(atrvec,clsvec)
    ptxy=prop.table(txy,1)
    logptxy=ifelse(ptxy==0,0,log(ptxy,2))
    loctab=ptxy*logptxy#Multiplication of corresponding position elements
    atr_clspvec=apply(loctab,1,sum)
    hatr=-as.vector(atrpvec%*%atr_clspvec)
    return(hatr)
  }
  
  #Joint entropy (multidimensional empirical entropy function)
  JEfunc <- function(atrcls,cls){
    if (ncol(atrcls)==1) {
      JE=HDfuncs(c(atrcls)$x)
    }
    if (ncol(atrcls)!=1) {#atrcls is matrix
      library(vcd)
      atrcls=data.frame(cls,atrcls)
      mytable2=xtabs(atrcls)
      re=ftable(prop.table(mytable2))
      atrclspvec=as.vector(re)
      logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
      JE=-as.vector(atrclspvec%*%logatrclspvec)
    }
    return(JE)
  }
  
  #Multidimensional empirical conditional entropy function
  Hatrm <- function(atrvec,clsvec){#Input as characteristic matrix and category matrix x,y
    if (ncol(atrvec)==1) {
      hatr=Hatrs(c(atrvec)$x,clsvec)
    }
    if (ncol(atrvec)!=1) {#atrvec is matrix, contains y,x
      library(vcd)
      all=data.frame(clsvec,atrvec)
      mytable=xtabs(all)
      atrpvec=as.vector(ftable(prop.table(mytable)))
      txy=table(all)
      ptxy=ftable(prop.table(txy,1))
      logptxy=ifelse(ptxy==0,0,log(ptxy,2))
      loctab=ptxy*logptxy
      atr_clspvec=matrix(0,nrow = nrow(loctab)/2,ncol = ncol(loctab))
      for (i in 1:(nrow(loctab)/2)) {
        for (j in 1:ncol(loctab)) {
          atr_clspvec[i,j]=loctab[i,j]+loctab[nrow(loctab)/2+i,j]
        }
      }
      atr_clspvec=as.vector(atr_clspvec)
      hatr=-as.vector(t(atrpvec)%*%atr_clspvec)
    }
    return(hatr)
  }
  HDcls=HDfuncs(cls)#H(Y)
  #HDatr=HDfuncs(atr)#H(X)
  #Ha(x)
  HDatr=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HDatr[i]=HDfuncs(atrcls=data.frame(xn))
  }
  
  #EN(Y|X)
  HatrVec=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HatrVec[i]=Hatrm(atrvec = data.frame(xn), clsvec = cls)
  }
  
  if (method=="info"){
    infogain=HDcls-HatrVec
  }else if(method=="inforate"){
    infogain=(HDcls-HatrVec)/HDatr
  }else stop("Please choose a useable method")
  names(infogain)=paste("x",c(1:max(group)),sep="")
  re=list(infogain=infogain,HDcls=HDcls,HatrVec=HatrVec,HDatr=HDatr)
  return(re)
}

GIGRSIS <- function(cls,atr,method,group){
  options(scipen = 200)
  #le=nrow(atr)
  res=GInfoGainmul(cls,atr,method,group)$infogain
  re=order(abs(res),decreasing = T)
  return(re)
}#GIGR-SIS

GInfoGainmul4 <- function(R,cls,atr,method=c("info","inforate"),group){
  #cls;Y
  #atr:X
  #One dimensional empirical entropy function
  HDfuncs <- function(atrcls){#aatrcls is vector
    I=length(atrcls)
    tatrcls=table(atrcls)
    atrclspvec=as.vector(tatrcls)/I#
    logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
    HD=-as.vector(atrclspvec%*%logatrclspvec)
    return(HD)
  }
  #One dimensional empirical conditional entropy function
  Hatrs <- function(atrvec,clsvec){#Input as feature vector and category vector
    n=length(atrvec)
    tatr=table(atrvec)
    atrpvec=as.vector(tatr)/n
    txy=table(atrvec,clsvec)
    ptxy=prop.table(txy,1)
    logptxy=ifelse(ptxy==0,0,log(ptxy,2))
    loctab=ptxy*logptxy#Multiplication of corresponding position elements
    atr_clspvec=apply(loctab,1,sum)
    hatr=-as.vector(atrpvec%*%atr_clspvec)
    return(hatr)
  }
  
  #Joint entropy (multidimensional empirical entropy function)
  JEfunc <- function(atrcls,cls){
    if (ncol(atrcls)==1) {
      JE=HDfuncs(c(atrcls)$x)
    }
    if (ncol(atrcls)!=1) {#atrcls is matrix
      library(vcd)
      atrcls=data.frame(cls,atrcls)
      mytable2=xtabs(atrcls)
      re=ftable(prop.table(mytable2))
      atrclspvec=as.vector(re)
      logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
      JE=-as.vector(atrclspvec%*%logatrclspvec)
    }
    return(JE)
  }
  
  #Multidimensional empirical conditional entropy function
  Hatrm <- function(atrvec,clsvec,R){#Input as characteristic matrix and category matrix x,y
    if (ncol(atrvec)==1) {
      hatr=Hatrs(c(atrvec)$x,clsvec)
    }
    if (ncol(atrvec)!=1) {#atrvec is matrix, contains y,x
      library(vcd)
      all=data.frame(clsvec,atrvec)
      mytable=xtabs(all)
      atrpvec=as.vector(ftable(prop.table(mytable)))
      
      txy=table(all)
      ptxy=ftable(prop.table(txy,1))
      logptxy=ifelse(ptxy==0,0,log(ptxy,2))
      
      loctab=ptxy*logptxy
      atr_clspvec=matrix(0,nrow = nrow(loctab)/R,ncol = ncol(loctab))
      for (i in 1:(nrow(loctab)/R)) {
        for (j in 1:ncol(loctab)) {
          for (w in 1:R) {
            atr_clspvec[i,j]=atr_clspvec[i,j]+loctab[4*(w-1)+i,j]
          }
        }
      }
      atr_clspvec=as.vector(atr_clspvec)
      hatr=-as.vector(t(atrpvec)%*%atr_clspvec)
    }
    
    return(hatr)
  }
  
  HDcls=HDfuncs(cls)#H(Y)
  #HDatr=HDfuncs(atr)#H(X)
  #Ha(x)
  HDatr=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HDatr[i]=HDfuncs(atrcls=data.frame(xn))
  }
  
  #EN(Y|X)
  HatrVec=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HatrVec[i]=Hatrm(atrvec = data.frame(xn), clsvec = cls,R)
  }
  
  if (method=="info"){
    infogain=HDcls-HatrVec
  }else if(method=="inforate"){
    infogain=(HDcls-HatrVec)/HDatr
  }else stop("Please choose a useable method")
  names(infogain)=paste("x",c(1:max(group)),sep="")
  re=list(infogain=infogain,HDcls=HDcls,HatrVec=HatrVec,HDatr=HDatr)
  return(re)
}

GIGRSIS2 <- function(R,cls,atr,method,group){
  options(scipen = 200)
  #le=nrow(atr)
  res=GInfoGainmul4(R,cls,atr,method,group)$infogain
  re=order(abs(res),decreasing = T)
  return(re)
}#GIGR-SIS

GInfoGainmul10 <- function(cls,atr,method=c("info","inforate"),group){
  #cls;Y
  #atr:X
  #One dimensional empirical entropy function
  HDfuncs <- function(atrcls){#aatrcls is vector
    I=length(atrcls)
    tatrcls=table(atrcls)
    atrclspvec=as.vector(tatrcls)/I#
    logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
    HD=-as.vector(atrclspvec%*%logatrclspvec)
    return(HD)
  }
  #One dimensional empirical conditional entropy function
  Hatrs <- function(atrvec,clsvec){#Input as feature vector and category vector
    n=length(atrvec)
    tatr=table(atrvec)
    atrpvec=as.vector(tatr)/n
    txy=table(atrvec,clsvec)
    ptxy=prop.table(txy,1)
    logptxy=ifelse(ptxy==0,0,log(ptxy,2))
    loctab=ptxy*logptxy#Multiplication of corresponding position elements
    atr_clspvec=apply(loctab,1,sum)
    hatr=-as.vector(atrpvec%*%atr_clspvec)
    return(hatr)
  }
  
  #Joint entropy (multidimensional empirical entropy function)
  JEfunc <- function(atrcls,cls){
    if (ncol(atrcls)==1) {
      JE=HDfuncs(c(atrcls)$x)
    }
    if (ncol(atrcls)!=1) {#atrcls is matrix
      library(vcd)
      atrcls=data.frame(cls,atrcls)
      mytable2=xtabs(atrcls)
      re=ftable(prop.table(mytable2))
      atrclspvec=as.vector(re)
      logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
      JE=-as.vector(atrclspvec%*%logatrclspvec)
    }
    return(JE)
  }
  
  #Multidimensional empirical conditional entropy function
  Hatrm <- function(atrvec,clsvec){#Input as characteristic matrix and category matrix x,y
    if (ncol(atrvec)==1) {
      hatr=Hatrs(c(atrvec)$x,clsvec)
    }
    if (ncol(atrvec)!=1) {
      library(vcd)
      all=data.frame(clsvec,atrvec)#atrvec is matrix, contains y,x
      mytable=xtabs(all)
      atrpvec=as.vector(ftable(prop.table(mytable)))
      
      txy=table(all)
      ptxy=ftable(prop.table(txy,1))
      logptxy=ifelse(ptxy==0,0,log(ptxy,2))
      
      loctab=ptxy*logptxy
      atr_clspvec=matrix(0,nrow = nrow(loctab)/4,ncol = ncol(loctab))
      for (i in 1:(nrow(loctab)/4)) {
        for (j in 1:ncol(loctab)) {
          for (w in 1:4) {
            atr_clspvec[i,j]=atr_clspvec[i,j]+loctab[(w-1)*nrow(loctab)/4+i,j]
          }
        }
      }
      atr_clspvec=as.vector(atr_clspvec)
      hatr=-as.vector(t(atrpvec)%*%atr_clspvec)
    }
    
    return(hatr)
  }
  HDcls=HDfuncs(cls)#H(Y)
  #HDatr=HDfuncs(atr)#H(X)
  #Ha(x)
  HDatr=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HDatr[i]=HDfuncs(atrcls=data.frame(xn))
  }
  
  #EN(Y|X)
  HatrVec=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HatrVec[i]=Hatrm(atrvec = data.frame(xn), clsvec = cls)
  }
  
  if (method=="info"){
    infogain=HDcls-HatrVec
  }else if(method=="inforate"){
    infogain=(HDcls-HatrVec)/HDatr
  }else stop("Please choose a useable method")
  names(infogain)=paste("x",c(1:max(group)),sep="")
  re=list(infogain=infogain,HDcls=HDcls,HatrVec=HatrVec,HDatr=HDatr)
  return(re)
}

GIGRSIS3 <- function(cls,atr,method,group){
  options(scipen = 200)
  #le=nrow(atr)
  res=GInfoGainmul10(cls,atr,method,group)$infogain
  re=order(abs(res),decreasing = T)
  return(re)
}#GIGR-SIS

#Group feature screening based on information gain ratio
GInfoGainmul3 <- function(cls,atr,method=c("info","inforate"),group){
  #cls;Y
  #atr:X
  #One dimensional empirical entropy function
  HDfuncs <- function(atrcls){#aatrcls is vector
    atrcls<-as.vector(atrcls)
    I=length(atrcls)
    tatrcls=table(atrcls)
    atrclspvec=as.vector(tatrcls)/I#
    logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
    HD=-as.vector(atrclspvec%*%logatrclspvec)
    return(HD)
  }
  #One dimensional empirical conditional entropy function
  Hatrs <- function(atrvec,clsvec){#Input as feature vector and category vector
    n=length(atrvec)
    tatr=table(atrvec)
    atrpvec=as.vector(tatr)/n
    atrcls<-matrix(atrvec,clsvec)
    txy=table(atrcls)
    ptxy=prop.table(txy,1)
    logptxy=ifelse(ptxy==0,0,log(ptxy,2))
    loctab=ptxy*logptxy#Multiplication of corresponding position elements
    atr_clspvec=apply(loctab,1,sum)
    hatr=-as.vector(atrpvec%*%atr_clspvec)
    return(hatr)
  }
  
  #Joint entropy (multidimensional empirical entropy function)
  JEfunc <- function(atrcls,cls){
    if (ncol(atrcls)==1) {
      JE=HDfuncs(c(atrcls)$x)
    }
    if (ncol(atrcls)!=1) {#atrcls is matrix
      library(vcd)
      atrcls=data.frame(cls,atrcls)
      mytable2=xtabs(atrcls)
      re=ftable(prop.table(mytable2))
      atrclspvec=as.vector(re)
      logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
      JE=-as.vector(atrclspvec%*%logatrclspvec)
    }
    return(JE)
  }
  
  #Multidimensional empirical conditional entropy function
  Hatrm <- function(atrvec,clsvec){#Input as characteristic matrix and category matrix x,y
    if (ncol(atrvec)==1) {
      hatr=Hatrs(atrvec,clsvec)
    }
    if (ncol(atrvec)!=1) {#atrvec is matrix, contains y,x
      library(vcd)
      all=data.frame(clsvec,atrvec)
      mytable=xtabs(all)
      atrpvec=as.vector(ftable(prop.table(mytable)))
      txy=table(all)
      ptxy=ftable(prop.table(txy,1))
      logptxy=ifelse(ptxy==0,0,log(ptxy,2))
      loctab=ptxy*logptxy
      atr_clspvec=matrix(0,nrow = nrow(loctab)/2,ncol = ncol(loctab))
      for (i in 1:(nrow(loctab)/2)) {
        for (j in 1:ncol(loctab)) {
          atr_clspvec[i,j]=loctab[i,j]+loctab[nrow(loctab)/2+i,j]
        }
      }
      atr_clspvec=as.vector(atr_clspvec)
      hatr=-as.vector(t(atrpvec)%*%atr_clspvec)
    }
    return(hatr)
  }
  HDcls=HDfuncs(cls)#H(Y)
  # HDatr=HDfuncs(atr)#H(X)
  #Ha(x)
  
  HDatr=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HDatr[i]=HDfuncs(atrcls=xn)
  }
  
  #EN(Y|X)
  HatrVec=c()
  for (i in 1:max(group)) {
    xn=atr[,which(group==i)]
    HatrVec[i]=Hatrm(atrvec = matrix(xn), clsvec = cls)
  }
  
  if (method=="info"){
    infogain=HDcls-HatrVec
  }else if(method=="inforate"){
    infogain=(HDcls-HatrVec)/HDatr
  }else stop("Please choose a useable method")
  names(infogain)=paste("x",c(1:max(group)),sep="")
  re=list(infogain=infogain,HDcls=HDcls,HatrVec=HatrVec,HDatr=HDatr)
  return(re)
}

GIGRSIS4 <- function(cls,atr,method,group){
  options(scipen = 200)
  #le=nrow(atr)
  res=GInfoGainmul3(cls,atr,method,group)$infogain
  re=order(abs(res),decreasing = T)
  return(re)
}#GIGR-SIS