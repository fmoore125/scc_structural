#multi variate analysis figure - Base SCC comparison and anova
library(lfe)
library(data.table)
library(patchwork)

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


d=ggplot(regdat,aes(x=var,y=coef))+geom_point(col="black")+theme_classic()+
  geom_errorbar(aes(ymin=coef-1.96*se,ymax=coef+1.96*se),col="black",width=0.2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_hline(yintercept = 0)+labs(x="",y="Regression Coefficient on Log SCC")+
  #scale_color_manual(values=c("Earth System"="#a6611a","Tipping Points: Climate"="#dfc27d","Tipping Points: Damages"="#fed976","Limited Substitutability"="#a1d99b","Persistent / Growth Damages"="#feb24c","Distributional Weights"="#99d8c9","Epstein-Zin"="#41ae76","Learning"="#9ecae1","Ambiguity/Model Uncertainty"="#005824"))+
  theme(legend.position = "none",text=element_text(size=16))

#------ANOVA Plot----Figure 2b

#multivariate analysis - explain scc variance as a function of structural changes, parametric variance, SCC Year, discount rate
dist=fread(file="outputs/distribution_v2_Jan2024.csv")
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/damage_funcs_lib.R")

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

#drop outlier Nordhaus row
todrop=which(dat$`Central Value ($ per ton CO2)`>70000)
distreg=distreg[-which(distreg$row==todrop)]

#concatenate carbon cycle and climate model changes into a single "Earth System Change"
distreg$Earth_system_struc=factor(ifelse(distreg$"Carbon.Cycle_struc"=="Yes","Yes",ifelse(distreg$"Climate.Model_struc"=="Yes","Yes","No")))

distreg=distreg[-which(distreg$draw<quantile(distreg$draw,0.01)|distreg$draw>quantile(distreg$draw,0.99)),]
distreg=distreg[complete.cases(distreg),]

distreg$whichmodel <- 'other'
distreg$whichmodel[distreg$dicemodel == 1] <- 'dice'
distreg$whichmodel[distreg$fundmodel == 1] <- 'fund'
distreg$whichmodel[distreg$pagemodel == 1] <- 'page'

colnames(distreg) <- gsub(" ", ".", colnames(distreg))

colnames(distreg)[grep("Persistent",colnames(distreg))]="Persistent...Growth.Damages_struc"

#regression of log SCC with paper fixed effects
df_paperfe2 <- distreg[which(distreg$draw>0 & distreg$sccyear_from2020 >= -10 & distreg$sccyear_from2020 <= 10),]
mod_paperfe2 <- lm(log(draw)~Earth_system_struc+Tipping.Points_struc+Tipping.Points2_struc+Persistent...Growth.Damages_struc+
                     Epstein.Zin_struc+Ambiguity.Model.Uncertainty_struc+Limitedly.Substitutable.Goods_struc+Inequality.Aversion_struc+
                     Learning_struc+TFP.Growth_param+Population.Growth_param+Emissions.Growth_param+Transient.Climate.Response_param+Carbon.Cycle2_param+
                     Equilibrium.Climate.Sensitivity_param+Tipping.Point.Magnitude_param+Damage.Function_param+Adaptation.Rates_param+Income.Elasticity_param+
                     Constant.Discount.Rate_param+EMUC2_param+PRTP2_param+Risk.Aversion..EZ.Utility._param+discountrate+I(discountrate^2)+declining+whichmodel+backstop+failure+log.scc.synth + missing.scc.synth + factor(paper),data=df_paperfe2)
anv <- anova(mod_paperfe2)
anvdf <- as.data.frame(anv)
anvdf$pred <- rownames(anvdf)
anvdf <- anvdf[!(anvdf$pred %in% c("Residuals", "factor(paper)")),]

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
                             'Learning', 'Other Market Failures',
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

stopifnot(nrow(anvdf %>% left_join(labels)) == nrow(anvdf))
anvdf2 <- labels %>% left_join(anvdf)
anvdf2$pred <- factor(anvdf2$pred, levels=anvdf2$pred)
anvdf2$cumss <- 1 - cumsum(anvdf2$`Sum Sq`) / sum(anvdf2$`Sum Sq`)
anvdf2.label <- anvdf2 %>% group_by(label) %>% summarize(cumss.lo=min(cumss)) %>% arrange(cumss.lo)
anvdf2.label$cumss.hi <- c(anvdf2.label$cumss.lo[-1], 1)

anvdf2.label$`Sum Sq` <- (anvdf2.label$cumss.hi + anvdf2.label$cumss.lo) / 2
anvdf2.label$width <- anvdf2.label$cumss.hi - anvdf2.label$cumss.lo

e=ggplot(anvdf2, aes(x=1, y=`Sum Sq`)) +
  geom_bar(aes(group=pred, fill=color),
           linewidth=.25, stat="identity", position="fill", width=1) +
  geom_text(data=subset(anvdf2.label, width > .002), aes(label=paste0(label," (",round(width*100,0),"%)"), hjust='left'), nudge_x=.65, size=4)+
  scale_fill_identity()+scale_y_continuous("Percent of explained variance", labels=scales::percent) + xlim(0.5, 4) + xlab(NULL) +
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text=element_text(size=10),axis.line.x = element_blank(),axis.ticks.length.x=unit(0,'lines'),plot.margin = margin(5.5,5.5,1,5.5) )


layout="
######BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
AAAAA#BBBB
######BBBB
"
patchwork=d+e+plot_layout(design=layout)

ggsave("figures/Science Revision/figure2.pdf",plot=patchwork,width=11,height=7,units="in")



