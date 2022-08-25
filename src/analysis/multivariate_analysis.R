## setwd("~/research/scciams/scc_structural/")

library(data.table)
library(DALEX)
library(ranger)
library(fixest)
library(modelsummary)
library(forcats)
library(lfe)
library(patchwork)
library(MetBrewer)

dist=fread(file="outputs/distribution_v2.csv")
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

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
failure=numeric(length=nrow(dat));failure[which(!is.na(dat$`Other Market Failure?`))]=1
sccyear_from2020=as.numeric(dat$`SCC Year`)-2020
marketonly=numeric(length=nrow(dat));marketonly[which(dat$`Market Only Damages`=="1.0")]=1
declining=numeric(length=nrow(dat));declining[which(dat$`Declining Discounting?` =="1.0")]=1
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

#simplest linear regression to start
mod_ols=feols(fml=I(log(draw))~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+sccyear_from2020+I(sccyear_from2020^2)+discountrate+I(discountrate^2)+declining+dicemodel+fundmodel+pagemodel+backstop+failure+log.scc.synth + missing.scc.synth,cluster="row",data=distreg[which(distreg$draw>0),])

mod_paperfe=feols(fml=I(log(draw))~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+sccyear_from2020+I(sccyear_from2020^2)+discountrate+I(discountrate^2)+declining+dicemodel+fundmodel+pagemodel+backstop+failure+log.scc.synth + missing.scc.synth|paper,cluster="row",data=distreg[which(distreg$draw>0),])

library(lfe)
mod_paperfe=felm(log(draw)~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+sccyear_from2020+I(sccyear_from2020^2)+discountrate+I(discountrate^2)+declining+dicemodel+fundmodel+pagemodel+backstop+failure+log.scc.synth + missing.scc.synth|paper | 0 | row,data=distreg[which(distreg$draw>0),])

varnames=c("sccyear_from2020" ="SCC Year","I(I(sccyear_from2020^2))" ="SCC Year^2","discountrate" ="Discount Rate","I(I(discountrate^2))"="Discount Rate^2","declining"="Declining DR","Earth_system_strucYes"="Earth System", "Tipping.Points_strucYes" ="Climate Tipping Points","Tipping.Points2_strucYes"="Damages Tipping Points", "Persistent...Growth.Damages_strucYes"="Growth Damages","Epstein.Zin_strucYes"="Epstein Zin","Ambiguity.Model.Uncertainty_strucYes"="Ambiguity","Limitedly.Substitutable.Goods_strucYes" ="Limited-Substitutability","Inequality.Aversion_strucYes"="Inequality Aversion","Learning_strucYes"="Learning","TFP.Growth_paramYes"="TFP Growth","Population.Growth_paramYes"="Pop Growth","Emissions.Growth_paramYes"="Emissions Growth", "Transient.Climate.Response_paramYes" ="Trans. Climate Resp.","Carbon.Cycle2_paramYes" ="Carbon Cycle (Param)","Equilibrium.Climate.Sensitivity_paramYes"="Eqm. Climate Sens.", "Tipping.Point.Magnitude_paramYes"="Tipping Point Size","Damage.Function_paramYes"="Damage Function","Adaptation.Rates_paramYes" ="Adaptation Rates","Income.Elasticity_paramYes"="Income Elasticity","Constant.Discount.Rate_paramYes" ="Const. Discount Rate","EMUC2_paramYes" ="EMUC", "PRTP2_paramYes"="PRTP","Risk.Aversion..EZ.Utility._paramYes" ="Risk Aversion", "dicemodel"="DICE","fundmodel"="FUND","pagemodel" ="PAGE","backstop"="Backstop Price","failure"="Other Market Failure","log.scc.synth"="Damage-based SCC")

modelsummary(models=list(mod_ols,mod_paperfe),coef_rename = varnames,output="outputs\\sccolsanalysis.docx",stars=TRUE)


#set up regression of difference from "baseSCC"
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

df <- get.all.scc(dat)

df$log.scc <- log(df$scc)
df$log.scc[!is.finite(df$log.scc)] <- NA

allcols <- names(dat)[c(1, 8, 10, 12:13,15:16, 18:23, 28:36)]
allcols[grep("Alternative ethical approaches", allcols)] <- "Alternative ethical approaches"

df <- multivar.prep(df)

#collapse Carbon Cycle and Climate Model into single Earth System category
df$Earth_system=factor(ifelse(df$"Carbon Cycle"=="1.0",1,ifelse(df$"Climate Model"=="1.0",1.0,0)))

allcols=append(allcols,"Earth_system")

#collapse "Calibrated" into 1
df=df%>%
  mutate(across(c("Carbon Cycle":"Alternative ethical approaches","Earth_system"),~fct_collapse(.x,"0"=c("-1.0","0"),"1"=c("1.0","Calibrated"))))

form <- as.formula(paste0("log.scc ~ `", paste(allcols[!(allcols %in% c("IAM Calibrated To (if applicable)", basemodelcols))], collapse="` + `"), "` + modified |  basecode|0|basecode"))
mod_basescc <- felm(form, data=df) #2167 observations deleted due to missingness?

varnames_base=c("Earth_system1"="Earth System","`Tipping Points`1" ="Climate Tipping Points","`Tipping Points2`1"="Damages Tipping Points", "`Persistent / Growth Damages`1"="Growth Damages","`Epstein-Zin`1"="Epstein Zin","`Ambiguity/Model Uncertainty`1"="Ambiguity","`Limitedly-Substitutable Goods`1" ="Limited-Substitutability","`Inequality Aversion`1"="Inequality Aversion","Learning1"="Learning","`Backstop Price?`1.0"="Backstop","`Other Market Failure?`1.0"="Other Market Failure","modifiedTRUE"="Other Modification","`Alternative ethical approaches`1"="Other Ethical Approaches")
modelsummary(models=list(mod_basescc),coef_rename = varnames_base,output="outputs\\sccolsanalysis_basescc.docx",stars=TRUE)



#impossible to show properly in table to make a plot of coefficient values instead
coefs=data.frame(coefs=c(mod_ols$coefficients,mod_paperfe$coefficients),sd=c(mod_ols$se,mod_paperfe$se),names=c(names(mod_ols$coefficients),names(mod_paperfe$coefficients)),mod=c(rep("No Fixed-Effects",length(mod_ols$coefficients)),rep("Paper Fixed-Effects",length(mod_paperfe$coefficients))))
coefs$type="Structural";coefs$type[grep("param",coefs$names)]="Parametric";coefs$type[c(grep("Intercept",coefs$names),grep("sccyear",coefs$names),grep("declining",coefs$names),grep("discountrate",coefs$names),grep("dicemodel",coefs$names),grep("fundmodel",coefs$names),grep("pagemodel",coefs$names),grep("backstop",coefs$names),grep("market",coefs$names),grep("failure",coefs$names))]="Other"
coefs$names=fct_recode(coefs$names,'SCC Year'="sccyear_from2020", "SCC Year^2"="I(I(sccyear_from2020^2))","Discount Rate"="discountrate","Discount Rate^2"="I(I(discountrate^2))","Declining DR"="declining","Earth System"="Earth_system_strucYes", "Climate Tipping Points"="Tipping.Points_strucYes","Damages Tipping Points"="Tipping.Points2_strucYes","Growth Damages"= "Persistent...Growth.Damages_strucYes","Epstein Zin"="Epstein.Zin_strucYes","Ambiguity"="Ambiguity.Model.Uncertainty_strucYes","Limited-Substitutability"="Limitedly.Substitutable.Goods_strucYes" ,"Inequality Aversion"="Inequality.Aversion_strucYes","Learning"="Learning_strucYes","TFP Growth"="TFP.Growth_paramYes","Pop Growth"="Population.Growth_paramYes","Emissions Growth"="Emissions.Growth_paramYes", "Trans. Climate Resp."="Transient.Climate.Response_paramYes" ,"Carbon Cycle"="Carbon.Cycle2_paramYes" ,"Eqm. Climate Sens."="Equilibrium.Climate.Sensitivity_paramYes", "Tipping Point Size"="Tipping.Point.Magnitude_paramYes","Damage Function"="Damage.Function_paramYes","Adaptation Rates"="Adaptation.Rates_paramYes" ,"Income Elasticity"="Income.Elasticity_paramYes","Const. Discount Rate"="Constant.Discount.Rate_paramYes" ,"EMUC"="EMUC2_paramYes", "PRTP"="PRTP2_paramYes","Risk Aversion"="Risk.Aversion..EZ.Utility._paramYes" , "DICE"="dicemodel","FUND"="fundmodel","PAGE"="pagemodel" ,"Backstop Price"="backstop","Market Only Damages"="marketonly","Other Market Failure"="failure")


#add base scc regression coefficients
coefs_base=data.frame(coefs=mod_basescc$coefficients,sd=mod_basescc$se,names=names(coefficients(mod_basescc)),mod="Base SCC Comparison")
coefs_base$type=c(rep("Other",2),rep("Structural",10),"Other")
coefs_base$names=fct_recode(coefs_base$names,"Earth System"="Earth_system1", "Climate Tipping Points"="`Tipping Points`1","Damages Tipping Points"="`Tipping Points2`1","Growth Damages"= "`Persistent / Growth Damages`1","Epstein Zin"="`Epstein-Zin`1","Ambiguity"="`Ambiguity/Model Uncertainty`1","Limited-Substitutability"="`Limitedly-Substitutable Goods`1" ,"Inequality Aversion"="`Inequality Aversion`1","Learning"="Learning1","Backstop Price"="`Backstop Price?`1.0","Other Market Failure"="`Other Market Failure?`1.0")

#drop some coefficients
coefs_base=coefs_base[-c(grep("modified",coefs_base$names),grep("Alternative",coefs_base$names)),]
colnames(coefs_base)=colnames(coefs)
coefs=rbind(coefs,coefs_base)

coefs$names=ordered(coefs$names,levels=rev(c("Earth System","Climate Tipping Points","Damages Tipping Points","Limited-Substitutability","Growth Damages","Epstein Zin","Inequality Aversion","Learning","Ambiguity","TFP Growth","Pop Growth","Emissions Growth","Trans. Climate Resp.","Carbon Cycle","Eqm. Climate Sens.","Tipping Point Size","Damage Function","Adaptation Rates","Income Elasticity","Const. Discount Rate","EMUC","PRTP","Risk Aversion","DICE","PAGE","FUND","Discount Rate","Discount Rate^2","Declining DR","SCC Year","SCC Year^2","Backstop Price","Market Only Damages","Other Market Failure","(Intercept)")))
coefs$mod=ordered(coefs$mod,levels=rev(c("Base SCC Comparison","Paper Fixed-Effects","No Fixed-Effects")))

coeftypes=c("Structural","Parametric","Other")
titles=c("Structural Changes","Parametric Variation","Other Parameters")
cols=met.brewer("Navajo",3,type="discrete")

coefs$title <- fct_recode(coefs$type, "Structural Changes"="Structural", "Parametric Variation"="Parametric",
                          "Other Parameters"="Other")
coefs$title <- factor(coefs$title, levels=titles)

ggplot(subset(coefs, !is.na(names) & !names%in%c("(Intercept)","SCC Year","SCC Year^2")), aes(coefs, names, col=mod)) +
    facet_grid(title ~ "Standardized coefficient values, by group", scales="free_y", space="free_y") +
    geom_vline(xintercept=0) +
    geom_point(position=position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin=coefs-1.67*sd,xmax=coefs+1.67*sd), position=position_dodge(width = 0.5)) +
    theme_bw() + labs(y=NULL,x="Coefficient Value\n(Dependent Variable = Log SCC)",col="Model") +
    geom_hline(yintercept = 0) + scale_color_manual(values=cols)+ theme(strip.background =element_rect(fill="white"))


plots=list()
count=1
for(i in coeftypes){
  plots[[count]]=ggplot(coefs%>%filter(!names%in%c("(Intercept)","SCC Year","SCC Year^2"))%>%filter(type==i),aes(x=names,y=coefs,ymin=coefs-1.67*sd,ymax=coefs+1.67*sd,col=mod))+geom_point(position=position_dodge(width = 0.9))+geom_errorbar(position=position_dodge(width = 0.9))+theme_bw()
  plots[[count]]=plots[[count]]+labs(x="",y="Coefficient Value\nDependent Variable = Log SCC",col="Model",title=titles[count])
  plots[[count]]=plots[[count]]+coord_flip()+geom_hline(yintercept = 0)+ylim(-2,3)
  plots[[count]]=plots[[count]]+scale_color_manual(values=cols,labels=c("Model 1: No Fixed Effects", "Model 2: Paper Fixed-Effects", "Model 3: Base SCC"))
  if(count!=3) plots[[count]]=plots[[count]]+theme(legend.position="none")
  plots[[count]]=plots[[count]]+geom_vline(xintercept = 1:13+0.5)
  count=count+1
}
x11()
plots[[1]]/plots[[2]]/plots[[3]]

save(coefs,file="outputs/regression_coefficients.RDat")


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
