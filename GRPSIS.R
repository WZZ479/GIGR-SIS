# if (!require("devtools"))
#   install.packages("devtools")
# devtools::install_github("debinqiu/grpss")
library(grpss)

GRPSIS<- function(X,y,group,criterion){
  options(scipen = 200)
  res=grp.criValues(X,y,group,criterion)
  res1<-as.data.frame(res)
  res2<-res1$grp.values
  re=order(res2,decreasing = T)
  return(re)
}