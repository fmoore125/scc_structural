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

#library(lfe)
# mod_paperfe=felm(log(draw)~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
#             Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
#             Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
#             Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
#             Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+sccyear_from2020+I(sccyear_from2020^2)+discountrate+I(discountrate^2)+declining+dicemodel+fundmodel+pagemodel+backstop+failure+log.scc.synth + missing.scc.synth|paper | 0 | row,data=distreg[which(distreg$draw>0),])

distreg$whichmodel <- 'other'
distreg$whichmodel[distreg$dicemodel == 1] <- 'dice'
distreg$whichmodel[distreg$fundmodel == 1] <- 'fund'
distreg$whichmodel[distreg$pagemodel == 1] <- 'page'

labels <- data.frame(pred=c('Transient.Climate.Response_param', 'Carbon.Cycle2_param', 'Equilibrium.Climate.Sensitivity_param',
                            'Earth_system_struc',
                            'Tipping.Points_struc', 'Tipping.Point.Magnitude_param', 'Tipping.Points2_struc',
                            'Persistent...Growth.Damages_struc', 'Limitedly.Substitutable.Goods_struc', 'Inequality.Aversion_struc',
                            'Risk.Aversion..EZ.Utility._param', 'Epstein.Zin_struc',
                            'whichmodel', 'Ambiguity.Model.Uncertainty_struc',
                            'Learning_struc', 'failure',
                            'TFP.Growth_param', 'Population.Growth_param', 'Emissions.Growth_param',
                            'Constant.Discount.Rate_param', 'discountrate', 'I(discountrate^2)',
                            'PRTP2_param', 'EMUC2_param',
                            'log.scc.synth', 'missing.scc.synth',
                            'Damage.Function_param', 'Adaptation.Rates_param', 'Income.Elasticity_param'),
                     type=c(rep('param', 3),
                             'struc',
                             'struc', 'param', 'struc',
                             'struc', 'struc', 'struc',
                             'param', 'struc',
                             'struc', 'param', # NOTE: treating model unc. as parametric
                             'struc', 'struc',
                             rep('param', 3),
                             'param', 'value', 'value',
                             'param', 'param',
                             'value', 'value',
                             rep('param', 3)),
                     label=c(rep('Earth System', 3),
                             'Earth System',
                             'Climate Tipping Points', 'Damages Tipping Points', 'Damages Tipping Points',
                             'Growth Damages', 'Limited Substitutability', 'Distributional Weights',
                             'Epstein-Zin', 'Epstein-Zin',
                             'Model', 'Model', # Label as model uncertainty
                             'Learning', 'Market Failures',
                             rep('Socioeconomic Uncertainty', 3),
                             rep('Discounting', 3),
                             'Pure Time Preference', 'Elasticity of Marginal Utility',
                             'Damage Function Parameters', 'Damage Function Parameters',
                             rep('Damage Function Parameters', 3)),
                     color=c(rep('#a6611a', 3),
                             '#a6611a',
                             '#dfc27d', '#fed976', '#fed976',
                             '#feb24c', '#a1d99b', '#74c476',
                             '#31a354', '#31a354',
                             rep('#006d2c', 2),
                             '#9ecae1', '#3182bd',
                             rep('#08519c', 3),
                             rep('#de8682', 3),
                             '#de2d26', '#a50f15',
                             '#54278f', '#54278f',
                             rep('#54278f', 3)))

mod_paperfe2 <- lm(log(draw)~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
            Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
            Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
            Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
            Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+discountrate+I(discountrate^2)+declining+whichmodel+backstop+failure+log.scc.synth + missing.scc.synth + factor(paper),data=distreg[which(distreg$draw>0 & distreg$sccyear_from2020 >= -10 & distreg$sccyear_from2020 <= 10),])

anv <- anova(mod_paperfe2)
anvdf <- as.data.frame(anv)
anvdf$pred <- rownames(anvdf)

anvdf <- anvdf[!(anvdf$pred %in% c("Residuals", "factor(paper)")),]

stopifnot(nrow(anvdf %>% left_join(labels)) == nrow(anvdf))
anvdf2 <- labels %>% left_join(anvdf)
anvdf2$pred <- factor(anvdf2$pred, levels=anvdf2$pred)
anvdf2$cumss <- 1 - cumsum(anvdf2$`Sum Sq`) / sum(anvdf2$`Sum Sq`)
anvdf2.label <- anvdf2 %>% group_by(label) %>% summarize(cumss.lo=min(cumss)) %>% arrange(cumss.lo)
anvdf2.label$cumss.hi <- c(anvdf2.label$cumss.lo[-1], 1)

anvdf2.label$`Sum Sq` <- (anvdf2.label$cumss.hi + anvdf2.label$cumss.lo) / 2
anvdf2.label$width <- anvdf2.label$cumss.hi - anvdf2.label$cumss.lo

ggplot(anvdf2, aes(x=1, y=`Sum Sq`)) +
    geom_bar(aes(group=pred, fill=color, linetype=type), colour='#000000',
             linewidth=.25, stat="identity", position="fill", width=1) +
    geom_text(data=subset(anvdf2.label, width > .002), aes(label=label, hjust='left'), nudge_x=.65, size=2) +
    geom_label(data=subset(anvdf2.label, width > .002 & label != "Epstein-Zin"), aes(label=round(width * 100, 0)), size=2) +
    geom_label(data=subset(anvdf2.label, label == "Epstein-Zin"), aes(label=round(width * 100, 0)), nudge_y=-.02, size=2) +
    scale_fill_identity() + scale_linetype_manual("Source of\nVariation\n(outline):",
                                                  breaks=c('struc', 'param', 'value'), values=c('solid', 'blank', 'dashed'),
                                                  labels=c("Structural", "Parametric", "Specified")) +
    scale_y_continuous("Percent of explained variance", labels=scales::percent) + xlim(0.5, 6) + xlab(NULL) +
    theme_bw() + theme(legend.justification=c(1,1), legend.position=c(.995,.995),
                            axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("figures/anova.pdf", width=3.5, height=5)

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
df$Earth_system=factor(ifelse(df$"Carbon Cycle"==1,1,ifelse(df$"Climate Model"=="1",1,0)))

allcols=append(allcols,"Earth_system")

#collapse "Calibrated" into 1
df=df%>%
  mutate(across(c("Carbon Cycle":"Alternative ethical approaches","Earth_system"),~fct_collapse(.x,"0"=c("-1.0","0","-1"),"1"=c("1.0","Calibrated"))))

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

