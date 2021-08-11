library(data.table)
library(DALEX)
library(ranger)

dist=fread(file="outputs/distribution.csv")
source("src/data_cleaining_scripts/cleaning_master.R")

#multivariate analysis - explain scc variance as a function of structural changes, parametric variance, SCC Year, discount rate

struc=dat%>%
  select("Carbon Cycle":"Learning")%>%
  replace(is.na(.),0)

for(i in 1:dim(struc)[2]) struc[,i]=fct_collapse(struc[,i],No=c("-1.0","0"),Yes=c("1.0","Calibrated"))
colnames(struc)=paste0(colnames(struc),"_struc")

param=dat%>%
  select("TFP Growth":"Risk Aversion (EZ Utility)")%>%
  replace(is.na(.),0)
for(i in 1:dim(param)[2]) param[,i]=fct_recode(param[,i],No="0",Yes="1.0")
colnames(param)=paste0(colnames(param),"_param")

#get other covariates
dicemodel=numeric(length=nrow(dat));dicemodel[c(grep("DICE",dat$`Base IAM (if applicable)`),grep("DICE",dat$`IAM Calibrated To (if applicable)`))]=1
fundmodel=numeric(length=nrow(dat));fundmodel[c(grep("FUND",dat$`Base IAM (if applicable)`),grep("FUND",dat$`IAM Calibrated To (if applicable)`))]=1
pagemodel=numeric(length=nrow(dat));pagemodel[c(grep("PAGE",dat$`Base IAM (if applicable)`),grep("PAGE",dat$`IAM Calibrated To (if applicable)`))]=1

backstop=numeric(length=nrow(dat));backstop[which(dat$`Backstop Price?`=="1.0")]=1
failure=numeric(length=nrow(dat));failure[which(dat$`Other Market Failure?`=="1.0")]=1
sccyear_from2020=as.numeric(dat$`SCC Year`)-2020
marketonly=numeric(length=nrow(dat));marketonly[which(dat$`Market Only Damages`=="1.0")]=1
declining=numeric(length=nrow(dat));declining[which(dat$`Declining Discounting?` =="1.0")]=1
discountrate=round(dat$discountrate,2)

covars=cbind(struc,param,dicemodel,fundmodel,pagemodel,backstop,sccyear_from2020,declining,discountrate)
distreg=cbind(dist,covars[dist$row,])

distreg=distreg[-which(distreg$draw<quantile(distreg$draw,0.01)|distreg$draw>quantile(distreg$draw,0.99)),]
distreg=distreg[complete.cases(distreg),]

#simplest linear regression to start
mod=lm(I(log(draw))~.,data=distreg[which(distreg$draw>0),-2])

#downsample distribution for random forest
distrf=distreg[which(distreg$draw>0),]
samp=1e5
distrf=distrf[sample(1:nrow(distrf),size=samp,replace=FALSE),]
distrf$y=log(distrf$draw)
distrf=as.data.frame(distrf)

distrf=distrf[,-which(colnames(distrf)%in%c("dicemodel","fundmodel","pagemodel"))]

rfmod=ranger(num.trees=100,min.node.size=500,max.depth=10,y=distrf$y,x=distrf[,-c(1:2,31)],verbose=TRUE,importance="permutation")

rfmod_explained=explain(rfmod,data=distrf[,-c(1:2,34)],y=distrf$y)
rfmod_diag=model_diagnostics(rfmod_explained)
rfmod_mod=model_parts(rfmod_explained)

