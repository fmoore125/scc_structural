#combine base SCC estimates with expert input distributions to construct alternative to Synthetic SCC in paper

#get base SCC estimates (code copied from figure_two.R)

library(lfe)
library(data.table)
library(patchwork)
library(MASS)

#Base SCC regression
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

df <- get.all.scc(dat)

df$log.scc <- log(df$scc)
df$log.scc[!is.finite(df$log.scc)] <- NA

allcols <- names(dat)[c(1, 8, 10, 12:13,15:16, 18:23, 28:36)]
allcols[grep("Alternative ethical approaches", allcols)] <- "Alternative ethical approaches"

df <- multivar.prep(df)

allcols=append(allcols,"Earth_system")

#collapse "Calibrated" into 1
df=df%>%
  mutate(across(c("Carbon Cycle":"Alternative ethical approaches"),~fct_collapse(.x,"0"=c("-1.0","0","-1"),"1"=c("1.0","Calibrated"))))
#collapse Carbon Cycle and Climate Model into single Earth System category
df$Earth_system=factor(ifelse(df$"Carbon Cycle"==1,1,ifelse(df$"Climate Model"=="1",1,0)))

#drop Nordhaus outlier row
todrop=which(df$log.scc>log(70000))
df=df[-todrop,]

form <- as.formula(paste0("log.scc ~ `", paste(allcols[!(allcols %in% c("IAM Calibrated To (if applicable)", basemodelcols))], collapse="` + `"), "` + modified |  basecode|0|basecode"))
mod_basescc <- felm(form, data=df) #2167 observations deleted due to missingness?

varnames_base=c("Earth_system1"="Earth System","`Tipping Points`1" ="Climate Tipping Points","`Tipping Points2`1"="Damages Tipping Points", "`Persistent / Growth Damages`1"="Growth Damages","`Epstein-Zin`1"="Epstein Zin","`Ambiguity/Model Uncertainty`1"="Ambiguity","`Limitedly-Substitutable Goods`1" ="Limited-Substitutability","`Inequality Aversion`1"="Inequality Aversion","Learning1"="Learning","`Backstop Price?`1.0"="Backstop","`Other Market Failure?`1.0"="Other Market Failure","modifiedTRUE"="Other Modification","`Alternative ethical approaches`1"="Other Ethical Approaches")

#plot base scc regression coefficients for structural changes

coefs=c("Earth_system1","`Tipping Points`1","`Tipping Points2`1","`Limitedly-Substitutable Goods`1","`Persistent / Growth Damages`1","`Inequality Aversion`1","`Epstein-Zin`1","Learning1","`Ambiguity/Model Uncertainty`1")
coeflocs=which(rownames(mod_basescc$coefficients)==coefs[1])
for(i in 2:length(coefs)) coeflocs=append(coeflocs,which(rownames(mod_basescc$coefficients)==coefs[i]))

regdat=data.frame(coef=mod_basescc$coefficients[coeflocs],se=mod_basescc$se[coeflocs])
rownames(regdat)=c("Earth System","Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Distributional Weights","Epstein-Zin","Learning","Ambiguity/Model Uncertainty")
regdat$var=rownames(regdat);regdat$var=factor(regdat$var,levels=c("Earth System","Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Distributional Weights","Epstein-Zin","Learning","Ambiguity/Model Uncertainty"))

varcovar=mod_basescc$vcv[coeflocs,coeflocs]
colnames(varcovar)=rownames(regdat); rownames(varcovar)=rownames(regdat)
  
#get multivariate Bayesian distributions over inputs 
bayespost=read.csv("data/expert_survey/meta-analysis-distribution.csv")
bayespost$question[which(bayespost$question=="Inequality Aversion")]="Distributional Weights"

#reorder question to match regression order
bayespost$question=factor(bayespost$question,levels=rownames(regdat))
bayespost=bayespost[order(bayespost$question),]

#1000 draws over multivariate joint posterior distribution

#--------------Alternate Synthetic SCC Sampler ----------------------

n=4000

coefdraws=mvrnorm(n=n,mu=regdat$coef,Sigma=varcovar)

#start with DICE model value of $43 in 2020 dollars based on Nordhaus 2018
scc_start=43
log_scc_end=numeric(length=n)

for(i in 1:n){
  #draw structural model from expert distribution and covert to binary based on Bernoulli distribution
  struc_draw=bayespost$prob[which(bayespost$iterations==i)]
  draw=sapply(struc_draw,FUN=function(x) rbinom(1,1,x))
  #calculate new log SCC based on draws
  log_scc_end[i]=log(scc_start)+coefdraws[i,]%*%draw
}

#add Duan smearing estimate to calculate modified SCCs
smear=mean(exp(mod_basescc$residuals))

scc_end=exp(log_scc_end)*smear

#get distribution stats
relprobs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)

truncmean=function(data,trim=0.001){
  highquant=quantile(data,1-trim)
  lowquant=quantile(data,trim)
  return(mean(data[which(data>lowquant&data<highquant)]))
}

scc_end_quants=quantile(scc_end,relprobs)

scc_end_truncmean=truncmean(scc_end)
