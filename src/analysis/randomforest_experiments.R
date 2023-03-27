library(data.table)
library(DALEX)
library(ranger)
library(tidyverse)

dist=fread(file="outputs/distribution_v2.csv")
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

load(file="outputs/randomforestmodel.Rdat")

#data used to fit random forest model
distrf=fread(file="outputs/distribution_structuralchangeweighted_withcovars_v2.csv")
distrf=as.data.frame(distrf)
distrf=distrf%>%filter(sccyear_from2020<=80)

#basic prediction matrix - from randomforest.R script
#includes scc year 2020, all parametric variation, discounting and structural changes from expert survey, random sample of damage functions
samppred=1e5
sampdat=fread(file="outputs/rf_experiments/basepredictiondata.csv")

#----A. No structural changes (classic DICE assumptions)-----
sampdat_DICE=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
#no structural or parametric changes
sampdat_DICE[,c(grep("param",predcols),grep("struc",predcols))]="No"
#discount rate matches DICE 2023 target discount rate of 4.6%
#see: https://yale.app.box.com/s/whlqcr7gtzdm4nxnrfhvap2hlzebuvvm/file/1150018404527 
sampdat_DICE$discountrate=4.6
#damage functions are sampled from DICE damage functions from post 2000
rel=which(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("DICE-2007","DICE-2013R","DICE-2016R2","DICE 2007","DICE 2010","DICE 2013","DICE 2013R","DICE 2016","DICE2007","DICE2010","DICE2013","DICE2016R"))
sampdat_DICE$log.scc.synth=sample(dat$log.scc.synth[rel],samppred,replace=TRUE)
sampdat_DICE$declining="No"
predictions_DICE=predict(rfmod,sampdat_DICE)
predictions_DICE=predictions_DICE$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)
fwrite(as.data.table(predictions_DICE),file="outputs/rf_experiments/A_DICE.csv")

#B. EPA assumptions
sampdat_EPA=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
#few structural or parametric changes
sampdat_EPA[,c(grep("param",predcols),grep("struc",predcols))]="No"
#structural changes to Earth System
sampdat_EPA$Earth_system_struc="Yes"
#parametric uncertainty in tfp growth, pop growth, earth system and damage functions
sampdat_EPA$TFP.Growth_param="Yes";sampdat_EPA$Population.Growth_param="Yes"
sampdat_EPA$Emissions.Growth_param="Yes"
sampdat_EPA$Transient.Climate.Response_param="Yes"; sampdat_EPA$Carbon.Cycle2_param="Yes"
sampdat_EPA$Equilibrium.Climate.Sensitivity_param='Yes'; sampdat_EPA$Damage.Function_param="Yes"
#central discount rate of 2%
sampdat_EPA$discountrate=2
#damage function from Howard and Sterner
rel=which(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("HowardSterner","HowardSterner (0.007438*T^2)"))
sampdat_EPA$log.scc.synth=sample(dat$log.scc.synth[rel],samppred,replace=TRUE)

predictions_EPA=predict(rfmod,sampdat_EPA)
predictions_EPA=predictions_EPA$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)

fwrite(as.data.table(predictions_EPA),file="outputs/rf_experiments/B_EPA.csv")

#C. Expert elicitation result - done in randomforest.R
#D. All structural changes
sampdat_allstruc=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
sampdat_allstruc[,grep("_struc",colnames(sampdat_allstruc))]="Yes"
predictions_allstruc=predict(rfmod,sampdat_allstruc)
predictions_allstruc=predictions_allstruc$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)

fwrite(as.data.table(predictions_allstruc),file="outputs/rf_experiments/C_allstruc.csv")

#E. Multiple discount rates (2, 3, 5, declining)
#F. Multiple publication years (2000, 2010, 2020)
#G. Multiple SCC pulse years (2020, 2050) - done in randomforest.R
#H. Multiple synthetic SCCs (DICE, FUND, PAGE, Howard & Sterner)


