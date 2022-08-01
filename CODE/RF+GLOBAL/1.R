options("expressions"=500000)
library(caret)
library(DALEX)
library(ggplot2)
library(ISLR)
library(lattice)
library(sampling)
library(randomForest)
library(magrittr)
library(dplyr)
library(Matrix)
library(tree)
library(rpart)
source("/work1/u1210789/RFplot.R")


kk=proc.time()
data=as.data.frame(read.csv("/work1/u1210789/csv/chr1.csv",header=T,row.names=1,check.name=T))
proc.time()-kk

#(0,1,2)轉換
data[which(data>0.5,arr.ind = T)]=1


kk=proc.time()
p=ncol(data)
for (i in 1:p) {
  k=p-i+1
  if(sum(data[,k])==0|sum(data[,k])==nrow(data)){
    data=data[,-k]
  }
}
proc.time()-kk
k=vector(length=(ncol(data)-1))
for (i in 1:(ncol(data)-1)) {
  k[i]=sum(data[,i]==data[,ncol(data)])
}
which(k==nrow(data))
which(k==0)
if(not(identical(union(which(k==nrow(data)),which(k==0)),integer(0))))
{data=data[,-c(which(k==0),which(k==nrow(data)))]}

for (i in 1:ncol(data)) {
  data[,i]=as.factor(data[,i])
}

x_p=ncol(data)-1
set.seed(1)
GN.val=round(table(data$disease)*0.1,0)
val_get=strata(data,"disease",size=c(GN.val[[1]],GN.val[[2]]),method="srswor")
val_index=val_get$ID_unit
val_data=data[val_index,]
train_index=c(1:nrow(data))[-val_index]
train_data=data[train_index,]

model= model=randomForest(disease~., 
                     data=train_data,
                     ntree=500,
                     importance=T)
rf_test_yhat=predict(model,val_data)

ptm <- proc.time()
score_table=matrix(0,nrow=500,ncol=(ncol(data)-1))
names(score_table)=names(data)[1:(ncol(data)-1)]

for(kk in 1:500) {
  x_name=names(data)[1:(ncol(data)-1)]
  little_tree=as.data.frame(getTree(model,kk,labelVar=TRUE))
  scorelist=c()
  k=1
  i=1
  while (max(i)<=nrow(little_tree)&sum(i)!=0)
  {
    scorelist[i]=k
    i=c(as.matrix(little_tree[i,1:2]))
    k=k+1
  }
  #有分支位置
  scoreindex=which(little_tree$'split var'!="NA")
  #分支名(重複)
  names_feature=little_tree$'split var'[scoreindex]
  
  #分支名在data中位置
  feature_index=vector(length=length(names_feature))
  for(l in 1:length(feature_index))
  {
    feature_index[l]=which(x_name==names_feature[l])
  }
  
  #分支名(無重複)
  level=unique(names_feature)
  row_data=rep(0,length=(ncol(data)-1))
  for(j in 1:length(level))
  {score_table[kk,unique(feature_index[which(names_feature==level[j])])]=mean(scorelist[scoreindex][which(names_feature==level[j])])}
  
}


forTime <- proc.time() - ptm


zero=matrix(nrow=2,ncol=ncol(score_table))

for(n in 1:ncol(score_table)){
  zero[1,n]=sum(score_table[,n])/sum(score_table[,n]!=0)
  zero[2,n]=sum(score_table[,n]!=0)
}




train_temp=train_data[,1:x_p]
train_temp$y_hat=predict(model,train_data)



err_time=proc.time()
err.rate=c()
for (i in 1:min(100,sum(is.na(zero[1,])==0))) {
  aa=paste("y_hat~",paste(names(data)[order(zero[1,])[1:i]],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate[i]=err/length(tree_test_yhat)
  
}
err_plot=plot(err.rate,type="l")
err_time=proc.time()-err_time

zero=as.data.frame(zero)
names(zero)=names(data)[1:x_p]
zero=zero[,order(as.numeric(zero[1,]),-as.numeric(zero[2,]))]
  
zero_up=zero[,which(zero[2,]>=mean(as.numeric(zero[2,])))]


err_time_up=proc.time()
err.rate_up=c()
for (i in 1:min(100,sum(is.na(zero_up[1,])==0))) {
  aa=paste("y_hat~",paste(names(zero_up)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_up[i]=err/length(tree_test_yhat)
  
}
err_plot_up=plot(err.rate_up,type="l")
err_time_up=proc.time()-err_time_up



zero_upup=zero[,which(zero[2,]>=(mean(as.numeric(zero[2,]))+sd(as.numeric(zero[2,]))))]

err_time_upup=proc.time()
err.rate_upup=c()
for (i in 1:min(100,sum(is.na(zero_upup[1,])==0))) {
  aa=paste("y_hat~",paste(names(zero_upup)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_upup[i]=err/length(tree_test_yhat)
  
}
err_plot_upup=plot(err.rate_upup,type="l")
err_time_upup=proc.time()-err_time_upup




new=zero
new[1,]=(as.numeric(new[1,])-mean(as.numeric(new[1,]),na.rm=T))/sd(as.numeric(new[1,]),na.rm=T)
new[3,]=-new[1,]
new[4,]=(new[2,]-mean(as.numeric(new[2,])))/sd(as.numeric(new[2,]))
new[5,]=new[4,]+new[3,]
max(new[5,],na.rm=T)
which.max(new[5,])
new=new[,order(-as.numeric(new[5,]))]

err_time_new=proc.time()
err.rate_new=c()
for (i in 1:min(100,sum(is.na(new[5,])==0))) {
  aa=paste("y_hat~",paste(names(new)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_new[i]=err/length(tree_test_yhat)
  
}
err_plot_new=plot(err.rate_new,type="l")
err_time_new=proc.time()-err_time_new
setwd("/work1/u1210789/RData/")
save.image("chr1_0_12.RData")


rm(list=ls())
options("expressions"=500000)
library(caret)
library(DALEX)
library(ggplot2)
library(ISLR)
library(lattice)
library(sampling)
library(randomForest)
library(magrittr)
library(dplyr)
library(Matrix)
library(tree)
library(rpart)
source("/work1/u1210789/RFplot.R")


kk=proc.time()
data=as.data.frame(read.csv("/work1/u1210789/csv/chr1.csv",header=T,row.names=1,check.name=T))
proc.time()-kk

#(0,1,2)轉換
data[which(data<1.5,arr.ind = T)]=0
data[which(data>1.5,arr.ind = T)]=1
data[1:97,ncol(data)]=1

kk=proc.time()
p=ncol(data)
for (i in 1:p) {
  k=p-i+1
  if(sum(data[,k])==0|sum(data[,k])==nrow(data)){
    data=data[,-k]
  }
}
proc.time()-kk
k=vector(length=(ncol(data)-1))
for (i in 1:(ncol(data)-1)) {
  k[i]=sum(data[,i]==data[,ncol(data)])
}
which(k==nrow(data))
which(k==0)
if(not(identical(union(which(k==nrow(data)),which(k==0)),integer(0))))
{data=data[,-c(which(k==0),which(k==nrow(data)))]}

for (i in 1:ncol(data)) {
  data[,i]=as.factor(data[,i])
}

x_p=ncol(data)-1
set.seed(1)
GN.val=round(table(data$disease)*0.1,0)
val_get=strata(data,"disease",size=c(GN.val[[1]],GN.val[[2]]),method="srswor")
val_index=val_get$ID_unit
val_data=data[val_index,]
train_index=c(1:nrow(data))[-val_index]
train_data=data[train_index,]

model= model=randomForest(disease~., 
                     data=train_data,
                     ntree=500,
                     importance=T)
rf_test_yhat=predict(model,val_data)

ptm <- proc.time()
score_table=matrix(0,nrow=500,ncol=(ncol(data)-1))
names(score_table)=names(data)[1:(ncol(data)-1)]

for(kk in 1:500) {
  x_name=names(data)[1:(ncol(data)-1)]
  little_tree=as.data.frame(getTree(model,kk,labelVar=TRUE))
  scorelist=c()
  k=1
  i=1
  while (max(i)<=nrow(little_tree)&sum(i)!=0)
  {
    scorelist[i]=k
    i=c(as.matrix(little_tree[i,1:2]))
    k=k+1
  }
  #有分支位置
  scoreindex=which(little_tree$'split var'!="NA")
  #分支名(重複)
  names_feature=little_tree$'split var'[scoreindex]
  
  #分支名在data中位置
  feature_index=vector(length=length(names_feature))
  for(l in 1:length(feature_index))
  {
    feature_index[l]=which(x_name==names_feature[l])
  }
  
  #分支名(無重複)
  level=unique(names_feature)
  row_data=rep(0,length=(ncol(data)-1))
  for(j in 1:length(level))
  {score_table[kk,unique(feature_index[which(names_feature==level[j])])]=mean(scorelist[scoreindex][which(names_feature==level[j])])}
  
}


forTime <- proc.time() - ptm


zero=matrix(nrow=2,ncol=ncol(score_table))

for(n in 1:ncol(score_table)){
  zero[1,n]=sum(score_table[,n])/sum(score_table[,n]!=0)
  zero[2,n]=sum(score_table[,n]!=0)
}




train_temp=train_data[,1:x_p]
train_temp$y_hat=predict(model,train_data)



err_time=proc.time()
err.rate=c()
for (i in 1:min(100,sum(is.na(zero[1,])==0))) {
  aa=paste("y_hat~",paste(names(data)[order(zero[1,])[1:i]],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate[i]=err/length(tree_test_yhat)
  
}
err_plot=plot(err.rate,type="l")
err_time=proc.time()-err_time

zero=as.data.frame(zero)
names(zero)=names(data)[1:x_p]
zero=zero[,order(as.numeric(zero[1,]),-as.numeric(zero[2,]))]
  
zero_up=zero[,which(zero[2,]>=mean(as.numeric(zero[2,])))]


err_time_up=proc.time()
err.rate_up=c()
for (i in 1:min(100,sum(is.na(zero_up[1,])==0))) {
  aa=paste("y_hat~",paste(names(zero_up)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_up[i]=err/length(tree_test_yhat)
  
}
err_plot_up=plot(err.rate_up,type="l")
err_time_up=proc.time()-err_time_up



zero_upup=zero[,which(zero[2,]>=(mean(as.numeric(zero[2,]))+sd(as.numeric(zero[2,]))))]

err_time_upup=proc.time()
err.rate_upup=c()
for (i in 1:min(100,sum(is.na(zero_upup[1,])==0))) {
  aa=paste("y_hat~",paste(names(zero_upup)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_upup[i]=err/length(tree_test_yhat)
  
}
err_plot_upup=plot(err.rate_upup,type="l")
err_time_upup=proc.time()-err_time_upup




new=zero
new[1,]=(as.numeric(new[1,])-mean(as.numeric(new[1,]),na.rm=T))/sd(as.numeric(new[1,]),na.rm=T)
new[3,]=-new[1,]
new[4,]=(new[2,]-mean(as.numeric(new[2,])))/sd(as.numeric(new[2,]))
new[5,]=new[4,]+new[3,]
max(new[5,],na.rm=T)
which.max(new[5,])
new=new[,order(-as.numeric(new[5,]))]

err_time_new=proc.time()
err.rate_new=c()
for (i in 1:min(100,sum(is.na(new[5,])==0))) {
  aa=paste("y_hat~",paste(names(new)[1:i],collapse="+"),sep="")
  one_tree= rpart(aa,data =train_temp, method="class")
  tree_test_yhat =predict(one_tree, val_data,type="class")
  err=sum(tree_test_yhat!=rf_test_yhat)
  err.rate_new[i]=err/length(tree_test_yhat)
  
}
err_plot_new=plot(err.rate_new,type="l")
err_time_new=proc.time()-err_time_new
setwd("/work1/u1210789/RData/")
save.image("chr1_01_2.RData")

