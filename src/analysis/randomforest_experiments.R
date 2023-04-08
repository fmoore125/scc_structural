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
#predictions_DICE=predictions_DICE$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)
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
#predictions_EPA=predictions_EPA$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)

fwrite(as.data.table(predictions_EPA),file="outputs/rf_experiments/B_EPA.csv")

#C. Expert elicitation result - done in randomforest.R
#D. All structural changes
sampdat_allstruc=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
sampdat_allstruc[,grep("_struc",colnames(sampdat_allstruc))]="Yes"
predictions_allstruc=predict(rfmod,sampdat_allstruc)
#predictions_allstruc=predictions_allstruc$predictions+sample(rfmod_explained$residuals,samppred,replace=TRUE)

fwrite(as.data.table(predictions_allstruc),file="outputs/rf_experiments/C_allstruc.csv")

#E. Multiple constant discount rates (1, 1.5, 2.5, 3, 5)
discs=c(1,1.5,2.5,3,5)

for(i in 1:length(discs)){
  sampdat_discs=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
  sampdat_discs$discountrate=discs[i];sampdat_discs$declining="No"
  predictions_disc=predict(rfmod,sampdat_discs)
  fwrite(as.data.table(predictions_disc),file=paste0("outputs/rf_experiments/E_discountrates_",discs[i],".csv"))
  print(i)
}

#F. Multiple publication years (2000, 2010, 2020)
pubyears=c(2000,2010,2020)

for(i in 1:length(pubyears)){
  sampdat_pubs=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
  sampdat_pubs$PublicationYear=pubyears[i]
  predictions_pubs=predict(rfmod,sampdat_pubs)
  fwrite(as.data.table(predictions_pubs),file=paste0("outputs/rf_experiments/F_pubyears_",pubyears[i],".csv"))
  print(i)
}

#G. Multiple SCC pulse years (2020, 2050) - done in randomforest.R
#H. Multiple synthetic SCCs (DICE, FUND, PAGE, Howard & Sterner)

damages=c(11.72858,6.337801,14.02219,38.60881)
damagemods=c("DICE2016r2","FUND38","PAGE2009","HowardSterner")

#From James: DICE 2016r2: 11.72858
#FUND 3.8: 6.337801
#PAGE 2009: 14.02219
#Howard & Sterner: 38.60881

for(i in 1:length(damages)){
  sampdat_damages=sampdat[sample(1:nrow(sampdat),samppred,replace=FALSE),]
  sampdat_damages$log.scc.synth=log(damages[i])
  predictions_damages=predict(rfmod,sampdat_damages)
  fwrite(as.data.table(predictions_damages),file=paste0("outputs/rf_experiments/H_damages_",damagemods[i],".csv"))
  print(i)
}




