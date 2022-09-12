#setwd("C:/Users/MyPC/Desktop/Model4")
library(grpregOverlap)#Loading package
data("pathway.dat")#p53 gene data 
y<-pathway.dat$mutation
X<-pathway.dat$expression
group<-pathway.dat$pathways
length(group)#308
group1<-group
names(group1)<-c(1:308)
incidence.mat <- incidenceMatrix(X, group1)
group2<-incidence.mat@i
group3<-as.vector(group2[1:4301])
write.csv(group3,file = "group.csv")
write.csv(X,file = "pathwayX.csv")
write.csv(y,file = "pathwayy.csv")
y1<-read.table("pathwayy.csv",header = T,sep = ",")
x1<-read.table("pathwayX.csv",header = T,sep = ",")
g<-read.table("group.csv",header = T,sep = ",")
g<-g[,-1]
y2<-y1[,-1]
x2<-x1[,-1]
Norlisan <- function(data,JK){
  
  P=ncol(data)
  N=nrow(data)
  
  x=matrix(0,nrow = nrow(data),ncol = ncol(data))
  for (p in 1:P) {
    mu=mean(data[,p]);
    sigma=sd(data[,p]);
    
    for (i in 1:N) {
      if (data[i,p]<=qnorm(1/JK,mu,sigma)) {x[i,p]=1}
      for (j in 2:(JK-1)) {
        x0=j-1;
        x1=j;
        if ((data[i,p]>qnorm(x0/JK,mu,sigma))&&(data[i,p]<=qnorm(x1/JK,mu,sigma))) {x[i,p]=j}
      }
      if (data[i,p]>qnorm((JK-1)/JK,mu,sigma)) {x[i,p]=JK}
    }
  }
  
  return(x)
}
#4-Categories
NEWX1=Norlisan(data = x2,4)
#8-Categories
NEWX2=Norlisan(data = x2,8)
#10--Categories
NEWX3=Norlisan(data = x2,10)

source("GIGR-SIS.R")
source("GIF-SIS.R")
source("GRPSIS.R")
x3<-as.matrix(x2)

GIFSIS3(cls=y2,atr=x3,method = "inforate",group =t(g))[1:12]
#247  91 100 102 198 275  63  22 143 189 157 244
GIGRSIS4(cls=y2,atr=x3,method = "inforate",group =t(g))[1:12]
#247  91 100 102 198 275  63  22 143 189 157 244
#4-Categories
GIGRSIS4(cls=y2,atr=NEWX1,method = "inforate",group =t(g))[1:12]
#91 177 241  74  82 167 264 195 185  68 146  35
#8-Categories
GIGRSIS4(cls=y2,atr=NEWX2,method = "inforate",group =t(g))[1:12]
#91 177  74 241  82  68  63 167 185 264  79 195
#10--Categories
GIGRSIS(cls=y2,atr=NEWX3,method = "inforate",group =t(g))[1:12]
#91 177 241  74  82 167  63  68 195 141 185  79
GRPSIS(X=x3,y=y2,group = t(g),"gSIS")[1:12]
#153 152 151 150 149 148 147 146 145 144 143 142
GRPSIS(X=x3,y=y2,group = t(g),"gHOLP")[1:12]
#153 152 151 150 149 148 147 146 145 144 143 142




