library(data.table)
library(DALEX)
library(ranger)
library(fixest)
library(modelsummary)
library(forcats)
library(lfe)
library(patchwork)
library(MetBrewer)
library(tidyverse)

dist=fread(file="outputs/distribution_v2.csv")
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

#multivariate analysis - explain scc variance as a function of structural changes, parametric variance, SCC Year, discount rate


#define sample type
#1. random sample - rand
#2. draw each structural change with equal probability -struc

# source("src/analysis/find_distribution.R")
# #if not equal-weighting by paper, then need to resample from original dataset
# #this code fits a distribution to each row in the dataset
# set.seed(12345)
# all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
# all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')
# 
# #start by generating distributions for each row
# dists=list()
# for (ii in 1:nrow(dat)) {
#   print(ii)
#   all.as <- t(dat[ii, all.as.cols])
#   qs <- all.qs[!is.na(all.as)]
#   as <- all.as[!is.na(all.as)]
#   mu <- dat$`Central Value ($ per ton CO2)`[ii]
#   if (is.na(mu) && length(qs) == 0) {
#     next
#   }
# 
#   dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
# }
# source("src/data_cleaining_scripts/cleaning_master.R")
# dat=dat%>%
#   mutate(across("Carbon Cycle":"Learning",~replace_na(.x, 0)))%>%
#   mutate(across("Carbon Cycle":"Learning",~fct_collapse(.x,No=c("-1.0","0"),Yes=c("1.0","Calibrated"))))
# cols=c(which(colnames(dat)=="Carbon Cycle"):which(colnames(dat)=="Learning"))
# colnames(dat)[cols]=paste0(colnames(dat)[cols],"_struc")
# 
# dat$Earth_system_struc=factor(ifelse(dat$"Carbon Cycle_struc"=="Yes","Yes",ifelse(dat$"Climate Model_struc"=="Yes","Yes","No")))
# #dat$Tipping_points_struc=factor(ifelse(dat$"Tipping Points_struc"=="Yes","Yes",ifelse(dat$"Tipping Points2_struc"=="Yes","Yes","No")))
# 
# type="struc" #possible values - rand, struc, struc2
# samp=1e6 #down-sample full 10e6 distribution to fit random forests
# 
# if(type=="rand"){
#   distrf=distreg[which(distreg$draw>0),]
#   distrf$Earth_system_struc=ifelse(distrf$Carbon.Cycle_struc=="Yes","Yes",ifelse(distrf$Climate.Model_struc=="Yes","Yes","No"))
#   distrf$Tipping_points_struc=ifelse(distrf$Tipping.Points_struc=="Yes","Yes",ifelse(distrf$Tipping.Points2_struc=="Yes","Yes","No"))
#   distrf=distrf[sample(1:nrow(distrf),size=samp,replace=FALSE),]
# }
# if(type=="struc"){
#   struc=grep("_struc",colnames(dat))[-c(1:4)]
# 
#   weights_struc=matrix(nrow=nrow(dat),ncol=length(struc))
#   for(i in 1:length(struc)){
#     #equally-weight rows with and without structural change
#     struc_yes=which(dat[,struc[i]]=="Yes");struc_no=which(dat[,struc[i]]=="No")
# 
#     #assign equal total weighting to the set of obserations with and without the structural change, with uniform sampling within each group
#     weights_struc[struc_yes,i]=1/length(struc_yes);weights_struc[struc_no,i]=1/length(struc_no)
#   }
#   #sum probability weights across rows and normalize to get probability weight for each row
#   weights=rowSums(weights_struc)/sum(rowSums(weights_struc))
# 
#   distrf=matrix(nrow=samp,ncol=2)
#   struc_samp=sample(1:nrow(dat),size=samp,replace=TRUE,prob=weights)
#   for(i in 1:samp){
#     if(i%%10000==0) print(i)
# 
#     distrf[i,1]=sample(dists[[struc_samp[i]]],1)
#     distrf[i,2]=struc_samp[i]
#   }
#   colnames(distrf)=c("draw","row")
# 
#   #bind in covariates
#   covars=cbind(dat[,struc],param,dicemodel,fundmodel,pagemodel,backstop,sccyear_from2020,declining,discountrate)
#   distreg=cbind(distrf,covars[as.matrix(distrf)[,2],])
# 
#   distreg=distreg[-which(distreg$draw<quantile(distreg$draw,0.01)|distreg$draw>quantile(distreg$draw,0.99)),]
#   distreg=distreg[complete.cases(distreg),]
# 
#   #fix column names
#   colnames(distreg) <- gsub(" ", ".", colnames(distreg))
#   colnames(distreg) <- gsub("/", ".", colnames(distreg))
#   colnames(distreg) <- gsub("-", ".", colnames(distreg))
#   colnames(distreg) <- gsub("\\(" ,".", colnames(distreg))
#   colnames(distreg) <- gsub(")", ".", colnames(distreg))
# 
#   distrf=distreg
#   fwrite(distrf,file="outputs/distribution_structuralchangeweighted_withcovars.csv")
# 
# }

distrf=fread(file="outputs/distribution_structuralchangeweighted_withcovars.csv")

distrf$y=log(distrf$draw)
distrf=as.data.frame(distrf)

#drop some variables
distrf=distrf[,-which(colnames(distrf)%in%c("dicemodel","fundmodel","pagemodel"))]
#limit to pre-2100 - vast majority of observations
distrf=distrf[-which(dat$`SCC Year`[distrf$row]>2100),]
#remove 2.6% of distribution with values <=0 that can't be logged
distrf=distrf[-which(is.na(distrf$y)|is.infinite(distrf$y)),]
#add in synthetic SCC variables from data frame
distrf$log.synth.scc=dat$log.scc.synth[distrf$row];distrf$missing.scc.synth=dat$missing.scc.synth[distrf$row]

rfmod=ranger(y~.,data=distrf%>%select(-c(draw,row)),num.trees=500,min.node.size=200,max.depth=12,verbose=TRUE,importance="permutation")

rfmod_explained=DALEX::explain(rfmod,data=distrf%>%select(-c(draw,row,y)),y=distrf$y)
rfmod_diag=model_diagnostics(rfmod_explained)
rfmod_mod=model_parts(rfmod_explained);plot(rfmod_mod)
rfmod_prof=model_profile(rfmod_explained,N=500,variable="sccyear_from2020",groups="Persistent...Growth.Damages_struc")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="discountrate")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="discountrate")
rfmod_prof_acc=model_profile(rfmod_explained,N=500,variable="discountrate",type="accumulated")



predat=distrf%>%select(-c(draw,row,y))

rfmod_pred=predict_parts(rfmod_explained, new_observation = predat[1,])
plot(rfmod_pred)
