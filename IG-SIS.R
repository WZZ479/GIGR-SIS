InfoGain <- function(cls,atr,method=c("info","inforate")){
  HDfunc <- function(atrcls){
    #aatrcls为向量
    I=length(atrcls)
    tatrcls=table(atrcls)
    atrclspvec=as.vector(tatrcls)/I
    logatrclspvec=ifelse(atrclspvec==0,0,log(atrclspvec,2))
    HD=-as.vector(atrclspvec%*%logatrclspvec)
    return(HD)
  }
  #计算经验条件熵
  Hatr <- function(atrvec,clsvec){
    #输入为特征向量及类别向量
    n=length(atrvec)
    tatr=table(atrvec)
    atrpvec=as.vector(tatr)/n
    txy=table(atrvec,clsvec)
    ptxy=prop.table(txy,1)
    logptxy=ifelse(ptxy==0,0,log(ptxy,2))
    loctab=ptxy*logptxy
    atr_clspvec=apply(loctab,1,sum)
    hatr=-as.vector(atrpvec%*%atr_clspvec)
    return(hatr)
  }
  
  HDcls=HDfunc(cls)
  HDatr=apply(atr, 2, HDfunc)
  HatrVec=apply(atr, 2, Hatr,clsvec=cls)
  if (method=="info"){
    infogain=HDcls-HatrVec
  }else if(method=="inforate"){
    infogain=(HDcls-HatrVec)/HDatr
  }else stop("Please choose a useable method")
  names(infogain)=paste("x",c(1:ncol(atr)),sep="")
  list(infogain=infogain,HDcls=HDcls,HatrVec=HatrVec,HDatr=HDatr)
}#ISIS
InfoGain(cls=realda$y,atr = realda$x,method = "inforate")
IGSIS <- function(cls,atr,method){
  options(scipen = 200)
  #le=nrow(atr)
  res=InfoGain(cls,atr,method)$infogain
  re=order(abs(res),decreasing = T)
  return(re)
}#IG-SIS