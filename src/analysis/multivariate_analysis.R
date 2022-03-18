library(data.table)
library(DALEX)
library(ranger)
library(fixest)
library(modelsummary)
library(forcats)

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


#simplest linear regression to start
mod=feols(fml=I(log(draw))~sccyear_from2020+I(sccyear_from2020^2)+discountrate+declining+Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+dicemodel+fundmodel+pagemodel+backstop,cluster="row",data=distreg[which(distreg$draw>0),])

mod_paperfe=feols(fml=I(log(draw))~sccyear_from2020+I(sccyear_from2020^2)+discountrate+declining+Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+dicemodel+fundmodel+pagemodel+backstop|paper,cluster="row",data=distreg[which(distreg$draw>0),])



varnames=c("sccyear_from2020" ="SCC Year","I(I(sccyear_from2020^2))" ="SCC Year^2","discountrate" ="Discount Rate","declining"="Declining DR","Carbon.Cycle_strucYes"="Carbon Cycle (Struc)","Climate.Model_strucYes"="Climate Model", "Tipping.Points_strucYes" ="Climate Tipping Points","Tipping.Points2_strucYes"="Damages Tipping Points", "Persistent...Growth.Damages_strucYes"="Growth Damages","Epstein.Zin_strucYes"="Epstein Zin","Ambiguity.Model.Uncertainty_strucYes"="Ambiguity","Limitedly.Substitutable.Goods_strucYes" ="Limited-Substitutability","Inequality.Aversion_strucYes"="Inequality AVersion","Learning_strucYes"="Learning","TFP.Growth_paramYes"="TFP Growth","Population.Growth_paramYes"="Pop Growth","Emissions.Growth_paramYes"="Emissions Growth", "Transient.Climate.Response_paramYes" ="Trans. Climate Resp.","Carbon.Cycle2_paramYes" ="Carbon Cycle (Param)","Equilibrium.Climate.Sensitivity_paramYes"="Eqm. Climate Sens.", "Tipping.Point.Magnitude_paramYes"="Tipping Point Size","Damage.Function_paramYes"="Damage Function","Adaptation.Rates_paramYes" ="Adaptation Rates","Income.Elasticity_paramYes"="Income Elasticity","Constant.Discount.Rate_paramYes" ="Const. Discount Rate","EMUC2_paramYes" ="EMUC", "PRTP2_paramYes"="PRTP","Risk.Aversion..EZ.Utility._paramYes" ="Risk Aversion", "dicemodel"="DICE","fundmodel"="FUND","pagemodel" ="PAGE","backstop"="Backstop Price")

modelsummary(models=list(mod,mod_paperfe),coef_rename = varnames,output="outputs\\sccolsanalysis.docx",stars=TRUE)

#impossible to show properly in table to make a plot of coefficient values instead
coefs=data.frame(coefs=c(mod$coefficients,mod_paperfe$coefficients),sd=c(mod$se,mod_paperfe$se),names=c(names(mod$coefficients),names(mod_paperfe$coefficients)),mod=c(rep("No Fixed-Effects",length(mod$coefficients)),rep("Paper Fixed-Effects",length(mod_paperfe$coefficients))))
coefs$type="Structural";coefs$type[grep("param",coefs$names)]="Parametric";coefs$type[c(grep("Intercept",coefs$names),grep("sccyear",coefs$names),grep("declining",coefs$names),grep("discountrate",coefs$names),grep("dicemodel",coefs$names),grep("fundmodel",coefs$names),grep("pagemodel",coefs$names),grep("backstop",coefs$names))]="Other"
coefs$names=fct_recode(coefs$names,'SCC Year'="sccyear_from2020", "SCC Year^2"="I(I(sccyear_from2020^2))","Discount Rate"="discountrate","Declining DR"="declining","Carbon Cycle (Struc)"="Carbon.Cycle_strucYes","Climate Model"="Climate.Model_strucYes",  "Climate Tipping Points"="Tipping.Points_strucYes","Damages Tipping Points"="Tipping.Points2_strucYes","Growth Damages"= "Persistent...Growth.Damages_strucYes","Epstein Zin"="Epstein.Zin_strucYes","Ambiguity"="Ambiguity.Model.Uncertainty_strucYes","Limited-Substitutability"="Limitedly.Substitutable.Goods_strucYes" ,"Inequality Aversion"="Inequality.Aversion_strucYes","Learning"="Learning_strucYes","TFP Growth"="TFP.Growth_paramYes","Pop Growth"="Population.Growth_paramYes","Emissions Growth"="Emissions.Growth_paramYes", "Trans. Climate Resp."="Transient.Climate.Response_paramYes" ,"Carbon Cycle (Param)"="Carbon.Cycle2_paramYes" ,"Eqm. Climate Sens."="Equilibrium.Climate.Sensitivity_paramYes", "Tipping Point Size"="Tipping.Point.Magnitude_paramYes","Damage Function"="Damage.Function_paramYes","Adaptation Rates"="Adaptation.Rates_paramYes" ,"Income Elasticity"="Income.Elasticity_paramYes","Const. Discount Rate"="Constant.Discount.Rate_paramYes" ,"EMUC"="EMUC2_paramYes", "PRTP"="PRTP2_paramYes","Risk Aversion"="Risk.Aversion..EZ.Utility._paramYes" , "DICE"="dicemodel","FUND"="fundmodel","PAGE"="pagemodel" ,"Backstop Price"="backstop")

coefs$names=ordered(coefs$names,levels=rev(c("Carbon Cycle (Struc)","Climate Model","Climate Tipping Points","Damages Tipping Points","Limited-Substitutability","Growth Damages","Epstein Zin","Inequality Aversion","Learning","Ambiguity","TFP Growth","Pop Growth","Emissions Growth","Trans. Climate Resp.","Carbon Cycle (Param)","Eqm. Climate Sens.","Tipping Point Size","Damage Function","Adaptation Rates","Income Elasticity","Const. Discount Rate","EMUC","PRTP","Risk Aversion","DICE","PAGE","FUND","Backstop Price","Discount Rate","Declining DR","SCC Year","SCC Year^2","(Intercept)")))

a=ggplot(coefs%>%filter(!names%in%c("(Intercept)","SCC Year","SCC Year^2")),aes(x=names,y=coefs,ymin=coefs-1.67*sd,ymax=coefs+1.67*sd,pch=mod,col=type))+geom_point(position=position_dodge(width = 0.9))+geom_errorbar(position=position_dodge(width = 0.9))+theme_bw()
a=a+labs(x="",y="Coefficient Value (Dependent Variable = Log SCC ($ per ton CO2))",col="Variable Type",pch="Model")
a=a+coord_flip()+geom_hline(yintercept = 0)
a
#####------------Random Forest Analysis ----------

#TO DO: Need to actually resample from original dataset, not just down-sample distribution with equal paper weighting

#define sample type
#1. random sample - rand
#2. draw each structural change with equal probability -struc
#3. draw each two-way interaction of structural changes with equal probability - struc2

source("src/analysis/find_distribution.R")
#if not equal-weighting by paper, then need to resample from original dataset
#this code fits a distribution to each row in the dataset
set.seed(12345)
all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

#start by generating distributions for each row
dists=list()
for (ii in 1:nrow(dat)) {
  print(ii)
  all.as <- t(dat[ii, all.as.cols])
  qs <- all.qs[!is.na(all.as)]
  as <- all.as[!is.na(all.as)]
  mu <- dat$`Central Value ($ per ton CO2)`[ii]
  if (is.na(mu) && length(qs) == 0) {
    next
  }
  
  dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
}
source("src/data_cleaining_scripts/cleaning_master.R")
dat=dat%>%
  mutate(across("Carbon Cycle":"Learning",~replace_na(.x, 0)))%>%
  mutate(across("Carbon Cycle":"Learning",~fct_collapse(.x,No=c("-1.0","0"),Yes=c("1.0","Calibrated"))))
cols=c(which(colnames(dat)=="Carbon Cycle"):which(colnames(dat)=="Learning"))
colnames(dat)[cols]=paste0(colnames(dat)[cols],"_struc")

dat$Earth_system_struc=factor(ifelse(dat$"Carbon Cycle_struc"=="Yes","Yes",ifelse(dat$"Climate Model_struc"=="Yes","Yes","No")))
#dat$Tipping_points_struc=factor(ifelse(dat$"Tipping Points_struc"=="Yes","Yes",ifelse(dat$"Tipping Points2_struc"=="Yes","Yes","No")))

type="struc" #possible values - rand, struc, struc2
samp=1e6 #down-sample full 10e6 distribution to fit random forests

if(type=="rand"){
  distrf=distreg[which(distreg$draw>0),]
  distrf$Earth_system_struc=ifelse(distrf$Carbon.Cycle_struc=="Yes","Yes",ifelse(distrf$Climate.Model_struc=="Yes","Yes","No"))
  distrf$Tipping_points_struc=ifelse(distrf$Tipping.Points_struc=="Yes","Yes",ifelse(distrf$Tipping.Points2_struc=="Yes","Yes","No"))
  distrf=distrf[sample(1:nrow(distrf),size=samp,replace=FALSE),]
} 
if(type=="struc"){
  struc=grep("_struc",colnames(dat))[-c(1:4)]
  
  weights_struc=matrix(nrow=nrow(dat),ncol=length(struc))
  for(i in 1:length(struc)){
    #equally-weight rows with and without structural change
    struc_yes=which(dat[,struc[i]]=="Yes");struc_no=which(dat[,struc[i]]=="No")
    
    #assign equal total weighting to the set of obserations with and without the structural change, with uniform sampling within each group
    weights_struc[struc_yes,i]=1/length(struc_yes);weights_struc[struc_no,i]=1/length(struc_no)
  }
  #sum probability weights across rows and normalize to get probability weight for each row
  weights=rowSums(weights_struc)/sum(rowSums(weights_struc))
  
  distrf=matrix(nrow=samp,ncol=2)
  struc_samp=sample(1:nrow(dat),size=samp,replace=TRUE,prob=weights)
  for(i in 1:samp){
    if(i%%10000==0) print(i)
    
    distrf[i,1]=sample(dists[[struc_samp[i]]],1)
    distrf[i,2]=struc_samp[i]
  }
  colnames(distrf)=c("draw","row")
  
  #bind in covariates
  covars=cbind(dat[,struc],param,dicemodel,fundmodel,pagemodel,backstop,sccyear_from2020,declining,discountrate)
  distreg=cbind(distrf,covars[as.matrix(distrf)[,2],])
  
  distreg=distreg[-which(distreg$draw<quantile(distreg$draw,0.01)|distreg$draw>quantile(distreg$draw,0.99)),]
  distreg=distreg[complete.cases(distreg),]

  #fix column names
  colnames(distreg) <- gsub(" ", ".", colnames(distreg))
  colnames(distreg) <- gsub("/", ".", colnames(distreg))
  colnames(distreg) <- gsub("-", ".", colnames(distreg))
  colnames(distreg) <- gsub("\\(" ,".", colnames(distreg))
  colnames(distreg) <- gsub(")", ".", colnames(distreg))
  
  distrf=distreg
  fwrite(distrf,file="outputs/distribution_structuralchangeweighted_withcovars.csv")
  
}

distrf$y=log(distrf$draw)
distrf=as.data.frame(distrf)

#drop some variables
distrf=distrf[,-which(colnames(distrf)%in%c("dicemodel","fundmodel","pagemodel"))]
#limit to pre-2100 - vast majority of observations
distrf=distrf[-which(dat$`SCC Year`[distrf$row]>2100),]
#remove 2.6% of distribution with values <=0 that can't be logged
distrf=distrf[-which(is.na(distrf$y)|is.infinite(distrf$y)),]

rfmod=ranger(y~.,data=distrf%>%select(-c(draw,row)),num.trees=500,min.node.size=200,max.depth=12,verbose=TRUE,importance="permutation")

rfmod_explained=explain(rfmod,data=distrf%>%select(-c(draw,row,y)),y=distrf$y)
rfmod_diag=model_diagnostics(rfmod_explained)
rfmod_mod=model_parts(rfmod_explained);plot(rfmod_mod)
rfmod_prof=model_profile(rfmod_explained,N=500,variable="sccyear_from2020",groups="Persistent...Growth.Damages_struc")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="discountrate",k=3)
