## Prepare data for multivariate_analysis.R or lasso.R

library(data.table)

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

dist=fread(file="outputs/distribution_v2.csv")

struc=dat%>%
  select("Carbon Cycle":"Learning")%>%
  replace(is.na(.),0)

for(i in 1:dim(struc)[2]) struc[,i]=fct_collapse(struc[,i],No=c("-1.0","0","-1"),Yes=c("1.0","Calibrated","1"))
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

backstop=numeric(length=nrow(dat));backstop[which(dat$`Backstop Price?` %in% c("1.0", "1"))]=1
failure=numeric(length=nrow(dat));failure[which(!is.na(dat$`Other Market Failure?`))]=1
sccyear_from2020=as.numeric(dat$`SCC Year`)-2020
marketonly=numeric(length=nrow(dat));marketonly[which(dat$`Market Only Damages`%in% c("1.0", "1"))]=1
declining=numeric(length=nrow(dat));declining[which(dat$`Declining Discounting?` %in% c("1.0", "1"))]=1
discountrate=round(dat$discountrate,2)

covars=cbind(struc,param,dicemodel,fundmodel,pagemodel,backstop,sccyear_from2020,declining,discountrate,marketonly,failure, log.scc.synth=dat$log.scc.synth, missing.scc.synth=dat$missing.scc.synth)
distreg=cbind(dist,covars[dist$row,])

distreg=distreg[-which(distreg$draw<quantile(distreg$draw,0.01)|distreg$draw>quantile(distreg$draw,0.99)),]
distreg=distreg[complete.cases(distreg),]

#add info on paper
distreg$paper=as.factor(dat$ID_number[distreg$row])

#fix column names
colnames(distreg) <- gsub(" ", ".", colnames(distreg))
colnames(distreg) <- gsub("/", ".", colnames(distreg))
colnames(distreg) <- gsub("-", ".", colnames(distreg))
colnames(distreg) <- gsub("\\(" ,".", colnames(distreg))
colnames(distreg) <- gsub(")", ".", colnames(distreg))

#concatenate carbon cycle and climate model changes into a single "Earth System Change"
distreg$Earth_system_struc=factor(ifelse(distreg$"Carbon.Cycle_struc"=="Yes","Yes",ifelse(distreg$"Climate.Model_struc"=="Yes","Yes","No")))
