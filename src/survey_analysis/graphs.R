library(tidyverse)
library(forcats)
library(rriskDistributions)
library(EnvStats)
library(data.table)

dat=read.csv("C:/Users/fmoore/Dropbox/SCC_expert_survey/Author_survey/Data/data_SCC-expert-survey_final_Fran.csv")

#first figure - distributions of literature and all things considered values

fig1dat=dat%>%
  select(contains(c("SCC_Lit","SCC_True_","Interview.number")))%>%
  pivot_longer(!"Interview.number..ongoing.")%>%
  mutate(type=substr(name,5,7))%>%
  mutate(quantile=as.factor(substr(name,9,11)))%>%
  mutate(quantile=fct_collapse(quantile,lower=c("2p5","_2p"),upper=c("7p5","_97"),central=c("Cen","_Ce")))%>%
  filter(!is.na(value))%>%
  pivot_wider(id_cols=c("Interview.number..ongoing.",type),values_from=value,names_from=quantile)
colnames(fig1dat)[1]="id"

#create jittered x variable
fig1dat$xj=jitter(as.numeric(as.factor(fig1dat$type)),amount=0.25)

#sample true and literature distributions over authors and add to plot, with 2010-2030 SCC distribution from meta-analysis
#assume triangular distribution - for simplicity use given 2.5% and 97.5% bounds

nsamp=10000

litdist=numeric(length=nsamp);truedist=numeric(length=nsamp)

#randomly draw experts uniformly
lit_samp=sample(unique(fig1dat$id[which(fig1dat$type=="Lit")]),nsamp,replace=TRUE);true_samp=sample(unique(fig1dat$id[which(fig1dat$type=="Tru")]),nsamp,replace=TRUE)

for(i in 1:nsamp){
  j=which(fig1dat$id==lit_samp[i]&fig1dat$type=="Lit")
  litdist[i]=ifelse(is.na(fig1dat$lower[j])|is.na(fig1dat$upper[j]),fig1dat$central[j],rtri(1,min=fig1dat$lower[j],max=fig1dat$upper[j],mode=fig1dat$central[j]))
  k=which(fig1dat$id==true_samp[i]&fig1dat$type=="Tru")
  truedist[i]=ifelse(is.na(fig1dat$lower[k])|is.na(fig1dat$upper[k]),fig1dat$central[k],rtri(1,min=fig1dat$lower[k],max=fig1dat$upper[k],mode=fig1dat$central[k]))
  if(i%%1000==0) print(i)
}

distdat=data.frame(dist=c(litdist,truedist),type=c(rep(1,nsamp),rep(2,nsamp)))

#read in meta-analysis distribution, subset to 2010-2030
metadist=fread("outputs/distribution.csv")
source("src/data_cleaining_scripts/cleaning_master.R")

metadist$SCCyear=dat$`SCC Year`[metadist$row]

metadist=metadist%>%filter(SCCyear>2010&SCCyear<2031)
metadist$type=3

#downsample distribution for ease of plotting
metadist=metadist[sample(1:nrow(metadist),0.2*nrow(metadist),replace=FALSE),]


#add boxplots to graph
a=ggplot(fig1dat,aes(x=xj,y=central,ymin=lower,ymax=upper,col=type))+geom_line(aes(group=id,x=xj,y=central),col="lightgrey")+geom_pointrange(lty=2)
a=a+coord_cartesian(ylim=c(-100,1100))+theme_bw()+theme(text=element_text(size=14))+labs(x="",y="2020 SCC ($ per ton CO2)")+scale_color_manual(values=c("#af2436","#51a154"),guide=FALSE)
a=a+scale_x_continuous(breaks=c(1,2), labels=c("Literature", "All Things Considered"), limits=c(0.5, 4.5)) 
a=a+geom_boxplot(data=distdat%>%filter(type==1),aes(y=dist,group=type,x=type),inherit.aes = FALSE,position=position_nudge(x=1.75),fill="#af2436",width=0.25)
a=a+geom_boxplot(data=distdat%>%filter(type==2),aes(y=dist,group=type,x=type),inherit.aes = FALSE,position=position_nudge(x=1.25),fill="#51a154",width=0.25)
a=a+geom_boxplot(data=metadist,aes(y=draw,group=type,x=type),inherit.aes = FALSE,position=position_nudge(x=0.75),fill="slateblue3",width=0.25)
a


