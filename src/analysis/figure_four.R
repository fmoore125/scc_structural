#load random forest
library(data.table)
library(tidyverse)

load("outputs\\randomforest_plots\\rfdistsmodel-final.RData")
source("src/data_cleaining_scripts/cleaning_master.R")
dist=fread("outputs/distribution_v2_Jan2024.csv")

#get bar plots for of distributions for example trees

examples=c(149,350)

#example draw:
#discount rate 3%
#included structural elements: Damages TP, Growth Damages, Limited Subs, Learning
#excluded structural elements: Earth system, climate TP, EZ, Ambiguity, Distributional Weighting
#SCC Year: 2020
#Pub Year: 2020
#All Param Uncertainty, declining DR, no other market failures, no backstop price, all damages


path_149=c("FALSE",'0',"FALSE","TRUE",'1',"FALSE")

path_350=c("FALSE",'1',"TRUE",'0','0')

#get rows for each split
rows_149=list()
rows_149[[1]]=forest[[149]]$rows;rows_149[[2]]=forest[[149]]$children$`FALSE`$rows;rows_149[[3]]=forest[[149]]$children$`FALSE`$children$'0'$rows
rows_149[[4]]=forest[[149]]$children$`FALSE`$children$'0'$children$`FALSE`$rows
rows_149[[5]]=forest[[149]]$children$`FALSE`$children$'0'$children$`FALSE`$children$`TRUE`$rows
rows_149[[6]]=forest[[149]]$children$`FALSE`$children$'0'$children$`FALSE`$children$`TRUE`$children$'1'$rows
rows_149[[7]]=forest[[149]]$children$`FALSE`$children$'0'$children$`FALSE`$children$`TRUE`$children$'1'$children$`FALSE`$rows

rows_350=list()
rows_350[[1]]=forest[[350]]$rows;rows_350[[2]]=forest[[350]]$children$`FALSE`$rows;rows_350[[3]]=forest[[350]]$children$`FALSE`$children$'1'$rows
rows_350[[4]]=forest[[350]]$children$`FALSE`$children$'1'$children$`TRUE`$rows
rows_350[[5]]=forest[[350]]$children$`FALSE`$children$'1'$children$`TRUE`$children$'0'$rows
rows_350[[6]]=forest[[350]]$children$`FALSE`$children$'1'$children$`TRUE`$children$'0'$children$'0'$rows

#loop through examples and generate quantiles for each 

relprobs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)

truncmean=function(data,trim=0.001){
  highquant=quantile(data,1-trim)
  lowquant=quantile(data,trim)
  return(mean(data[which(data>lowquant&data<highquant)]))
}
pfuns=c(map(relprobs,~partial(quantile,probs=.x,na.rm=T)), truncmean)

fulldist=dist%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(fulldist)=c("lowest","min","lower","middle","upper","max","highest", "mu")

dists_149=data.frame(fulldist)
for(i in 1:length(rows_149)){
  tempdist=dist%>%
    filter(row%in%rows_149[[i]])%>%
    summarize_at(vars(draw),funs(!!!pfuns))
  colnames(tempdist)=c("lowest","min","lower","middle","upper","max","highest", "mu")
  dists_149=rbind(dists_149,tempdist)
}

dists_149$name=c("Full",paste("Split",1:length(rows_149)))
dists_149$y=rev(1:8)

dists_350=data.frame(fulldist)
for(i in 1:length(rows_350)){
  tempdist=dist%>%
    filter(row%in%rows_350[[i]])%>%
    summarize_at(vars(draw),funs(!!!pfuns))
  colnames(tempdist)=c("lowest","min","lower","middle","upper","max","highest", "mu")
  dists_350=rbind(dists_350,tempdist)
}

dists_350$name=c("Full",paste("Split",1:length(rows_350)))
dists_350$y=c(8,7,6,5,4,3,1)+0.5

boxplots_149=ggplot(dists_149)+coord_flip(ylim=c(-5,1200))+theme_classic()+
  geom_boxplot(aes(x=y,min=min,lower=lower,middle=middle,upper=upper,max=max,group=name),inherit.aes=FALSE,stat="identity",width=0.3)+
  geom_segment(aes(x=y,xend=y,y=lowest,yend=min,group=name),lty=2)+geom_segment(aes(x=y,xend=y,y=max,yend=highest,group=name),lty=2)+
  geom_point(aes(x=y,y=mu,group=name))+theme(text=element_text(size=16))+labs(y="SCC ($ per ton CO2)",x="")

boxplots_350=ggplot(dists_350)+coord_flip(ylim=c(-5,1200))+theme_classic()+
  geom_boxplot(aes(x=y,min=min,lower=lower,middle=middle,upper=upper,max=max,group=name),inherit.aes=FALSE,stat="identity",width=0.3)+
  geom_segment(aes(x=y,xend=y,y=lowest,yend=min,group=name),lty=2)+geom_segment(aes(x=y,xend=y,y=max,yend=highest,group=name),lty=2)+
  geom_point(aes(x=y,y=mu,group=name))+theme(text=element_text(size=16))+labs(y="SCC ($ per ton CO2)",x="")

#combine into one plot
dists_149$tree=149;dists_350$tree=350
dists_both=rbind(dists_149,dists_350)

boxplots_both=ggplot(dists_both,aes(x=y,min=min,lower=lower,middle=middle,upper=upper,max=max,group=interaction(name,tree),col=as.factor(tree)))+coord_flip(ylim=c(-5,1200))+theme_classic()+
  geom_boxplot(stat="identity")+geom_segment(aes(xend=y,y=lowest,yend=min),lty=2)+geom_segment(aes(xend=y,y=max,yend=highest),lty=2)+
  geom_point(aes(y=mu))+theme(text=element_text(size=16))+labs(y="SCC ($ per ton CO2)",x="")+scale_color_manual(values=c("#41b3a9","#870cb3"),guide="none")

#get prediction from random forest for example input
source("src/analysis/randomforest_dists_load.R")

sampdat=dat[1,]
##set obvious ones
sampdat[,c(which(colnames(sampdat)=="TFP Growth"):which(colnames(sampdat)=="Risk Aversion (EZ Utility)"))]=factor("Yes",levels=c('No', 'Yes'))

sampdat$backstop=factor("No", levels=c('No', 'Yes'));sampdat$declining=factor("Yes", levels=c('No', 'Yes'));sampdat$marketonly=factor("No", levels=c('No', 'Yes'))
sampdat$failure=factor("No", levels=c('No', 'Yes'));sampdat$missing.scc.synth=FALSE
#for synthetic scc - draw from literature values
sampdat$log.scc.synth=sample(dat$log.scc.synth[-which(is.na(dat$log.scc.synth))],samppred,replace=TRUE)
sampdat$PublicationYear=2020
sampdat$discountrate=3

#structural elements for sample: Excluded elements
sampdat$Earth_system_struc=factor("No",levels=c("Yes","No"))
sampdat$Tipping.Points_struc=factor("No",levels=c("Yes","No"))
sampdat$Epstein.Zin_struc=factor("No",levels=c("Yes","No"))
sampdat$Ambiguity.Model.Uncertainty_struc=factor("No",levels=c("Yes","No"))
sampdat$Inequality.Aversion_struc=factor("No",levels=c("Yes","No"))
#structural elements for sample: Included elements
sampdat$Limitedly.Substitutable.Goods_struc=factor("Yes",levels=c("Yes","No"))
sampdat$Persistent...Growth.Damages_struc=factor("Yes",levels=c("Yes","No"))
sampdat$Tipping.Points2_struc=factor("Yes",levels=c("Yes","No"))
sampdat$Learning_struc=factor("Yes",levels=c("Yes","No"))
sampdat$sccyearformerge=2020

#set irrelevant columns to NA
sampdat$ID_number=NA;sampdat$`Bibtex Name`=NA;sampdat$Reference=NA;sampdat$`Reported Base Model SCC (if applicable)`=NA;sampdat$Authors=NA;sampdat$`Empirical Improvement or Sensitvity Analysis?`=NA;sampdat$`Added By`=NA
sampdat$`SCC Dollar Year`=NA;sampdat$`IAM Calibrated To (if applicable)`=NA;
sampdat$log.scc.2020usd=NA;sampdat$`SCC Year`=NA;sampdat$discountrate2=NA
sampdat$`PAPER LOCATION`=NA;sampdat$deflator=NA;sampdat$`Declining Discounting?`=NA
sampdat$temp.2100=NA;sampdat$temp.2100.source=NA
sampdat[,c(which(colnames(sampdat)=="Carbon Cycle"):which(colnames(sampdat)=="Max"))]=NA
sampdat$`Earth system`=NA;sampdat$EMUC=NA;sampdat$PRTP=NA;sampdat$`Market Only Damages`=NA;sampdat$`Socio-Economic Scenario`=NA
sampdat$cons_growth_percap=NA;sampdat$`Central Value ($ per ton CO2)`=NA;sampdat$`Base IAM (if applicable)`=NA

dist=predict.forest(forest,sampdat,dist=NULL)

#empirical cdf

