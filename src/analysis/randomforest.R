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
library(ggridges)
library(plyr)
library(xtable)

dist=fread(file="outputs/distribution_v2.csv")
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

#multivariate analysis - explain scc variance as a function of structural changes, parametric variance, SCC Year, discount rate


#define sample type
#1. random sample - rand
#2. draw each structural change with equal probability -struc

source("src/analysis/find_distribution.R")
#if not equal-weighting by paper, then need to resample from original dataset
#this code fits a distribution to each row in the dataset
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
#   #deal with odd case where quantiles are the same value
#   if(ii==274) {qs=numeric(0);as=numeric(0)}
# 
#   dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
# }
# dat=dat%>%
#   mutate(dplyr::across("Carbon Cycle":"Learning",~replace_na(.x, "0")))%>%
#   mutate(dplyr::across("Carbon Cycle":"Learning",~as.factor(.x)))%>%
#   mutate(dplyr::across("Carbon Cycle":"Learning",~fct_collapse(.x,No=c("-1.0","0","-1"),Yes=c("1.0","Calibrated","1"))))
# cols=c(which(colnames(dat)=="Carbon Cycle"):which(colnames(dat)=="Learning"))
# colnames(dat)[cols]=paste0(colnames(dat)[cols],"_struc")
# 
# dat$Earth_system_struc=factor(ifelse(dat$"Carbon Cycle_struc"=="Yes","Yes",ifelse(dat$"Climate Model_struc"=="Yes","Yes","No")))
# #dat$Tipping_points_struc=factor(ifelse(dat$"Tipping Points_struc"=="Yes","Yes",ifelse(dat$"Tipping Points2_struc"=="Yes","Yes","No")))
# 
# samp=1e6 #down-sample full 10e6 distribution to fit random forests
# 
# struc=grep("_struc",colnames(dat))[-c(1:2)]
# 
# weights_struc=matrix(nrow=nrow(dat),ncol=length(struc))
# for(i in 1:length(struc)){
#   #equally-weight rows with and without structural change
#   struc_yes=which(dat[,struc[i]]=="Yes");struc_no=which(dat[,struc[i]]=="No")
#   
#   #assign equal total weighting to the set of obserations with and without the structural change, with uniform sampling within each group
#   weights_struc[struc_yes,i]=1/length(struc_yes);weights_struc[struc_no,i]=1/length(struc_no)
# }
# #sum probability weights across rows and normalize to get probability weight for each row
# weights=rowSums(weights_struc)/sum(rowSums(weights_struc))
# 
# distrf=matrix(nrow=samp,ncol=2)
# struc_samp=sample(1:nrow(dat),size=samp,replace=TRUE,prob=weights)
# for(i in 1:samp){
#   if(i%%10000==0) print(i)
#   
#   distrf[i,1]=sample(dists[[struc_samp[i]]],1)
#   distrf[i,2]=struc_samp[i]
# }
# colnames(distrf)=c("draw","row")
# 
# #bind in covariates
# param=dat%>%
#   select("TFP Growth":"Risk Aversion (EZ Utility)")%>%
#   replace(is.na(.),0)
# colnames(param)=paste0(colnames(param),"_param")
# 
# backstop=numeric(length=nrow(dat));backstop[which(dat$`Backstop Price?`=="1.0")]=1
# failure=numeric(length=nrow(dat));failure[which(!is.na(dat$`Other Market Failure?`))]=1
# sccyear_from2020=as.numeric(dat$`SCC Year`)-2020
# marketonly=numeric(length=nrow(dat));marketonly[which(dat$`Market Only Damages`=="1.0")]=1
# declining=numeric(length=nrow(dat));declining[which(dat$`Declining Discounting?` =="1.0")]=1
# discountrate=round(dat$discountrate,2)
# 
# covars=cbind(dat[,struc],sccyear_from2020,param,backstop,declining,marketonly,failure, discountrate,log.scc.synth=dat$log.scc.synth, missing.scc.synth=dat$missing.scc.synth)
# covars=covars%>%
#   mutate(dplyr::across("TFP Growth_param":"failure",~as.factor(.x)))%>%
#   mutate(dplyr::across("TFP Growth_param":"failure",~fct_recode(.x,No="0",Yes="1")))
# distrf=cbind(distrf,covars[distrf[,2],])
# 
# distrf=distrf[-which(distrf$draw<quantile(distrf$draw,0.01)|distrf$draw>quantile(distrf$draw,0.99)),]
# distrf=distrf[complete.cases(distrf),]
# 
# 
# #fix column names
# colnames(distrf) <- gsub(" ", ".", colnames(distrf))
# colnames(distrf) <- gsub("/", ".", colnames(distrf))
# colnames(distrf) <- gsub("-", ".", colnames(distrf))
# colnames(distrf) <- gsub("\\(" ,".", colnames(distrf))
# colnames(distrf) <- gsub(")", ".", colnames(distrf))
# 
# fwrite(distrf,file="outputs/distribution_structuralchangeweighted_withcovars_v2.csv")

distrf=fread(file="outputs/distribution_structuralchangeweighted_withcovars_v2.csv")

distrf$y=log(distrf$draw)
distrf=as.data.frame(distrf)

#limit to pre-2100 - vast majority of observations
distrf=distrf%>%filter(sccyear_from2020<=80)
#remove 2.6% of distribution with values <=0 that can't be logged
distrf=distrf[-which(is.na(distrf$y)|is.infinite(distrf$y)),]

rfmod=ranger(y~.,data=distrf%>%select(-c(draw,row)),num.trees=500,min.node.size=200,max.depth=12,verbose=TRUE,importance="impurity_corrected")

rfmod_explained=DALEX::explain(rfmod,data=distrf%>%select(-c(draw,row,y)),y=distrf$y)
rfmod_diag=model_diagnostics(rfmod_explained)
save(rfmod, rfmod_explained,rfmod_diag,file="outputs/randomforestmodel.Rdat")

rfmod_mod=model_parts(rfmod_explained);
rfmod_mod$variable=fct_recode(rfmod_mod$variable,"SCC Year"="sccyear_from2020","Discount Rate"="discountrate","Log Damage-based SCC"="log.scc.synth","Declining DR"="declining","Earth System"="Earth_system_struc","Climate Tipping Points"="Tipping.Points_struc","Damages Tipping Points"="Tipping.Points2_struc","Growth Damages"="Persistent...Growth.Damages_struc","Epstein Zin"="Epstein.Zin_struc","Ambiguity"="Ambiguity.Model.Uncertainty_struc","Limited-Substitutability"="Limitedly.Substitutable.Goods_struc","Inequality Aversion"="Inequality.Aversion_struc","Learning"="Learning_struc","TFP Growth"="TFP.Growth_param","Pop Growth"="Population.Growth_param","Emissions Growth"="Emissions.Growth_param","Trans. Climate Resp"="Transient.Climate.Response_param","Carbon Cycle (param)"="Carbon.Cycle2_param","Eqm. Climate Sens."="Equilibrium.Climate.Sensitivity_param","Tipping Point Size"="Tipping.Point.Magnitude_param","Damage Function"="Damage.Function_param","Adaptation Rates"="Adaptation.Rates_param","Income Elasticity"="Income.Elasticity_param","Const. Disc. Rate"="Constant.Discount.Rate_param","EMUC"="EMUC2_param","PRTP"="PRTP2_param","Risk Aversion"="Risk.Aversion..EZ.Utility._param","Backstop Price"="backstop","Other Market Failure"="failure","Market Only Damages"="marketonly","Missing Damage SCC"="missing.scc.synth")
#customize feature importance plot
moddat=rfmod_mod%>%
  group_by(variable)%>%
  dplyr::summarize(mean=mean(dropout_loss),min=min(dropout_loss),max=max(dropout_loss))%>%
  filter(variable!="_baseline_"&variable!="_full_model_")%>%
  arrange(desc(mean))
moddat$variable=fct_reorder(moddat$variable,moddat$mean)
moddat$type=c("Other","Other","Damage Func","Struc","Struc","Damage Func","Damage Func","Struc","Param","Struc","Struc","Struc","Param","Param","Struc","Param","Param","Struc","Param","Param","Struc","Param","Param","Param","Other","Param","Param","Param","Other","Other","Other")
moddat$type=fct_relevel(moddat$type,c("Struc","Param","Damage Func","Other"))

a=ggplot(moddat,aes(y=variable,yend=variable,x=1,xend=mean,xmin=min,xmax=max,color=type))+geom_segment(size=5)
a=a+theme_bw()+geom_errorbar(col="black",width=0)+labs(title="Feature Importance, Random Forest Model",y="",x="RMSE Loss After Permutations",color="")
a=a+scale_color_manual(values=met.brewer(name="Lakota", n=4, type="discrete"),labels=c("Structural","Parametric","Damage Func.","Other"))
a=a+theme(legend.position = c(0.92, 0.2))

rfmod_prof=model_profile(rfmod_explained,N=500,variable="sccyear_from2020",groups="Persistent...Growth.Damages_struc")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="discountrate")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="log.scc.synth")

plot(rfmod_prof,geom="profiles")

#add predictions from random forest

samppred=1e6

predcols=colnames(distrf[-which(colnames(distrf)%in%c("draw","row","y"))])
sampdat=matrix(nrow=samppred,ncol=length(predcols))
colnames(sampdat)=predcols
sampdat=as.data.frame(sampdat)

#set obvious ones
sampdat[,grep("param",predcols)]="Yes"
sampdat$backstop="No";sampdat$declining="Yes";sampdat$marketonly="No"
sampdat$failure="No";sampdat$missing.scc.synth=FALSE
#for synthetic scc - draw from literature values
sampdat$log.scc.synth=sample(dat$log.scc.synth[-which(dat$log.scc.synth==0)],samppred,replace=TRUE)

#for discount rate, sample from Moritz's expert survey results
discountsurvey=read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
sampdat$discountrate=sample(discountsurvey$SDR[-which(is.na(discountsurvey$SDR))],samppred,replace=TRUE)

#for structural changes, use data from expert survey on assessed quality of structural changes
load("outputs/expert_survey_data_products/fig2surveydata.Rdat")

#define mapping from agree-disagree to probabilities
level_key = c("Strongly Disagree" = "0", "Disagree" = "0.25", "Neither Agree nor Disagree" = "0.5","Agree"="0.75","Strongly Agree"="1")

#transform factor levels to probabilities
strucprobs=fig2dat_qual%>%
  dplyr::mutate(across(-ID,~recode(.x,!!!level_key)))%>%
  dplyr::mutate(across(-ID,~as.numeric(.x)))

averageprobs=colMeans(strucprobs[,-1],na.rm=T)

struclookup=data.frame(name1=predcols[grep("struc",predcols)],name2=c("Tipping Points: Climate","Tipping Points: Damages","Persistent / Growth Damages","Epstein-Zin","Ambiguity/Model Uncertainty","Limited Substitutability","Inequality Aversion","Learning","Earth System"))

for(i in 1:nrow(struclookup)){
  probname=which(names(averageprobs)==struclookup$name2[i])
  col=which(colnames(sampdat)==struclookup$name1[i])
  sampdat[,col]=runif(samppred);sampdat[,col]=ifelse(sampdat[,col]<averageprobs[probname],"Yes","No")
}

#loop through a set of years to get scc distribution over time
years=c(2020,2030,2040,2050,2060,2070,2080,2090,2100)

predictionyears=matrix(nrow=samppred,ncol=length(years))

for(i in 1:length(years)){
  sampdat$sccyear_from2020=years[i]-2020
  
  predictions=predict(rfmod,sampdat)
  
  #add in residuals from random forest to get full distribution
  fulldist=predictions$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)
  predictionyears[,i]=fulldist
  print(years[i])
}

means=colMeans(exp(predictionyears))
quants=apply(predictionyears,MARGIN=2,FUN=function(x) quantile(x,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)))
rf_table=rbind(exp(quants),means)
colnames(rf_table)=years
rownames(rf_table)=c(rownames(rf_table)[1:9],"Mean")

tab<-xtable(s, caption= "summary statistics of air pollution data", 
            align=c("|c","|c","|c","|c","|c|"))

colnames(predictionyears)=years
  
temp=as.tibble(predictionyears)%>%
  pivot_longer(everything(),names_to="year",values_to="scc")%>%
  filter(year%in%c(2020,2050,2100))%>%
  group_by(year)%>%
  filter(quantile(scc,0.01)<scc&quantile(scc,0.99)>scc)%>%
  group_modify(~ ggplot2:::compute_density(.x$scc, NULL,bw=0.1))

breaks_ln=c(0,2.5,5,7.5)
a=ggplot(temp,aes(x=x,y=as.factor(year),group=year,height=density))+theme_ridges()+geom_density_ridges(stat="identity",lwd=0.75,fill="#f5dc62",scale=0.92)
a=a+labs(x="SCC ($ per ton CO2)",y="")+scale_x_continuous(limits=c(0,8),breaks=breaks_ln,labels=c(round_any(exp(breaks_ln[1:3]),10),round_any(exp(breaks_ln[4]),100)))
a
