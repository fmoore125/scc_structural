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

#simplest linear regression to start
mod=feols(fml=I(log(draw))~sccyear_from2020+I(sccyear_from2020^2)+discountrate+declining+Carbon.Cycle_struc+Climate.Model_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+dicemodel+fundmodel+pagemodel+backstop,cluster="row",data=distreg[which(distreg$draw>0),])

mod_paperfe=feols(fml=I(log(draw))~sccyear_from2020+I(sccyear_from2020^2)+discountrate+declining+Carbon.Cycle_struc+Climate.Model_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+dicemodel+fundmodel+pagemodel+backstop|paper,cluster="row",data=distreg[which(distreg$draw>0),])

# #add in "base SCC values" to get additional power for within-paper variation analysis
# basevals=dat[-which(is.na(dat$`Reported Base Model SCC (if applicable)`)),]
# basedat=matrix(nrow=nrow(basevals),ncol=ncol(distreg));colnames(basedat)=colnames(distreg);basedat=as.data.frame(basedat)
# basedat$draw=basevals$`Reported Base Model SCC (if applicable)`;basedat$row=rownames(basevals);basedat$row=as.numeric(basedat$row)
# basedat[,3:26]="No"
# basedat$dicemodel=0;basedat$dicemodel[which(basedat$row%in%c(grep("DICE",dat$`Base IAM (if applicable)`),grep("DICE",dat$`IAM Calibrated To (if applicable)`)))]=1
# basedat$fundmodel=0;basedat$fundmodel[which(basedat$row%in%c(grep("FUND",dat$`Base IAM (if applicable)`),grep("FUND",dat$`IAM Calibrated To (if applicable)`)))]=1
# basedat$pagemodel=0;basedat$pagemodel[which(basedat$row%in%c(grep("PAGE",dat$`Base IAM (if applicable)`),grep("PAGE",dat$`IAM Calibrated To (if applicable)`)))]=1
# basedat$backstop=dat$`Backstop Price?`[basedat$row];basedat$backstop[which(is.na(basedat$backstop))]=0
# basedat$declining=dat$`Declining Discounting?`[basedat$row];basedat$declining[which(is.na(basedat$declining))]=0
# basedat$sccyear_from2020=as.numeric(dat$`SCC Year`[basedat$row])-2020;basedat$discountrate=dat$discountrate[basedat$row];basedat$paper=dat$ID_number[basedat$row]
# 
# distreg_mod=rbind(distreg,basedat)
# 
# mod_paperfe_withbase=feols(fml=I(log(draw))~sccyear_from2020+I(sccyear_from2020^2)+discountrate+declining+Carbon.Cycle_struc+Climate.Model_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
#                     Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
#                     Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
#                     Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
#                     Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+dicemodel+fundmodel+pagemodel+backstop|paper,cluster="row",data=distreg_mod[which(distreg_mod$draw>0),])
# 

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
#downsample distribution for random forest
distrf=distreg[which(distreg$draw>0),]
samp=1e6
distrf=distrf[sample(1:nrow(distrf),size=samp,replace=FALSE),]
distrf$y=log(distrf$draw)
distrf=as.data.frame(distrf)

#drop some variables
distrf=distrf[,-which(colnames(distrf)%in%c("dicemodel","fundmodel","pagemodel"))]
#limit to pre-2100 - vast majority of observations
distrf=distrf[-which(dat$`SCC Year`[distrf$row]>2100),]

rfmod=ranger(num.trees=200,min.node.size=200,max.depth=12,y=distrf$y,x=distrf[,-c(1:2,31)],verbose=TRUE,importance="permutation")

rfmod_explained=explain(rfmod,data=distrf[,-c(1:2,31)],y=distrf$y)
rfmod_diag=model_diagnostics(rfmod_explained)
rfmod_mod=model_parts(rfmod_explained);plot(rfmod_mod)
rfmod_prof=model_profile(rfmod_explained,N=500,variable="sccyear_from2020",groups="Persistent / Growth Damages_struc")
rfmod_prof=model_profile(rfmod_explained,N=500,variable="discountrate",k=3)
