#show expert assessment of discount rate and model structure compared to literature
library(tidyverse)
library(patchwork)
library(MetBrewer)


#discount rate
#Expert discount rates from Drupp and Hansel 2018 AEJ paper
discountsurvey=read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")

#get literature discount rates
source("src/data_cleaining_scripts/cleaning_master.R")

fig4adat=data.frame(dr=discountsurvey$SDR,type="Expert")
fig4adat=rbind(fig4adat,data.frame(dr=dat$discountrate,type="Literature"))

fig4a=ggplot(fig4adat)+geom_density(aes(x=dr,group=type,lty=type,col=type),lwd=0.75)+
  theme_classic()+labs(x="Discount Rate",y="")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+
  scale_x_continuous(expand = c(0, 0), limits = c(0,12.5))+scale_color_manual(values=c("deepskyblue4","firebrick2"))+
  annotate("text",x=c(6.8,5.5),y=c(0.18,0.29),col=c("firebrick2","deepskyblue4"),label=c("Literature","Expert Assessment"))

#this is output from bayesian meta-analysis over joint probability of structural change inclusion (in bayes.R)
bayespost=read.csv("data/expert_survey/meta-analysis-distribution.csv")

#get fraction of dataset containing structural changes
dat$'Earth System'=ifelse(dat$`Climate Model`=="1.0"|dat$`Carbon Cycle`=="1.0","1.0",NA)

strucfrac=dat%>%
  #dplyr::select("Tipping Points":"Alternative ethical approaches (not Discounted Utilitarianism)","Earth System")%>%
  dplyr::summarise(across(c("Tipping Points":"Learning","Earth System"), ~ sum(.x=="1.0",na.rm=T)/n()))%>%
  pivot_longer(everything())

strucfrac$name[c(1:2,6)]=c("Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability")
colnames(strucfrac)[1]="question"

bayespost$question=factor(bayespost$question,levels=c("Earth System","Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Inequality Aversion","Epstein-Zin","Learning","Ambiguity/Model Uncertainty"))

fig4b=ggplot(bayespost)+geom_density(aes(x=prob,group=question),col="deepskyblue4")+facet_wrap(~as.factor(question),ncol=3)+
  theme_bw()+geom_vline(data=strucfrac,aes(xintercept=value),lty=2,lwd=1,color="firebrick2")+labs(x="Probability of Inclusion, Expert Assessment",y="")+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.background=element_rect(fill="white"))+xlim(-0.02,1)+
  geom_text(data=strucfrac,aes(x=value-0.025,y=10),label="Literature Frequency",angle=90,color="firebrick2")

fig4a+fig4b

#add variable importance plot
vip=read.csv("outputs/randomforest_plots/variable_importance_plot_Jan2024.csv")
vip$column=fct_reorder(vip$column,vip$import)

fig4c=ggplot(vip,aes(x=0,xend=import,y=column,yend=column,col=Type))+geom_segment(linewidth=4)+
  scale_color_manual(values=met.brewer(name="Hokusai1", n=4, type="discrete"),labels=c("Damage\nFunction","Other","Parametric\nUncertainty","Model Structure"))+
  theme_bw()+labs(x="Varaible Importance",y="",color="")


