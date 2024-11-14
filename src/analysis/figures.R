## Generates extra figures on residual variance, official SCC comparisons, and the distribution of structural changes.

library(ggplot2)
library(data.table)
library(ggridges)
library(fixest)
library(plyr)
library(MetBrewer)
library(lfe)
library(forcats)
library(zoo)
library(patchwork)

dist=fread(file="outputs/distribution_v2_Jan2024.csv")
dist_weighted=fread(file="outputs/distribution_coauthorweighted_v2.csv")
dist_weighted_citations=fread(file="outputs/distribution_citationweighted_v2.csv")

source("src/data_cleaining_scripts/cleaning_master.R")
#drop outlier Nordhaus row
todrop=which(dat$`Central Value ($ per ton CO2)`>70000)
if(is.finite(todrop)) dist=dist[-which(dist$row==todrop),]

#add year information
dist$year=as.numeric(dat$`SCC Year`[dist$row])
dist$yeargroup=cut(dist$year,breaks=c(1990,2010,2030,2070,2100,2400))
dist$yeargroup=fct_recode(dist$yeargroup,'<2010'="(1.99e+03,2.01e+03]",'2010-2030'="(2.01e+03,2.03e+03]",'2030-2070'="(2.03e+03,2.07e+03]",'2070-2100'="(2.07e+03,2.1e+03]",'>2100'="(2.1e+03,2.4e+03]")

#different figure looking at publication date
dist$pubyear=dat$Year[dist$row]
distshort=dist%>%filter(yeargroup%in%c('2010-2030',"2030-2070"))%>%
  group_by(yeargroup,pubyear)%>%dplyr::summarize(mean=mean(draw,na.rm=T),median=quantile(draw,0.5,na.rm=T),lower_25=quantile(draw,0.25,na.rm=T),upper_25=quantile(draw,0.75,na.rm=T),lower_5=quantile(draw,0.05,na.rm=T),upper_5=quantile(draw,0.95,na.rm=T))

distshort=pivot_longer(distshort,cols=3:8,names_to="Stats",values_to="SCC")
distshort$Stats=factor(distshort$Stats,levels=c("mean","lower_5","lower_25","median","upper_25","upper_5"))

distshort=distshort%>%group_by(Stats,yeargroup)%>%dplyr::arrange(pubyear)%>%mutate(rollingSCC=zoo::rollmean(SCC,5,fill="extend"))

#just for purposes of plotting - replace negative values with 1
distshort$rollingSCC[which(distshort$rollingSCC<0)]=1

breaks_ln=c(0,2,4,6)
b=ggplot(distshort,aes(x=pubyear,y=log(rollingSCC)))+geom_line(aes(lty=Stats),data=distshort%>%filter(Stats!="mean"))+facet_wrap(~yeargroup,ncol=1)
b=b+scale_linetype_manual(values=c(3,2,1,2,3),guide=FALSE)+theme_bw()+labs(x="Publication Year",y="SCC ($ per ton CO2)")
b=b+scale_y_continuous(limits=c(-1,7),breaks=breaks_ln,labels=c(round_any(exp(breaks_ln),10)))
b=b+geom_line(data=distshort%>%filter(Stats=="mean"),lwd=1.25,col="#f0aa3d")
b=b+theme(strip.background =element_rect(fill="white"))
b

#downsample distribution so avoid crashing computer
library(ggridges)

dist$pubyear=dat$Year[dist$row]
distplot=dist[which(dist$year%in%2010:2030),]
distplot=distplot[-which(distplot$draw<quantile(distplot$draw,0.05)|distplot$draw>quantile(distplot$draw,0.95)),]
distplotmeans=distplot%>%group_by(pubyear)%>%summarize(meanval=mean(draw))

a=ggplot(distplot,aes(x=draw,y=as.factor(pubyear),group=pubyear,fill=pubyear))
a=a+stat_density_ridges(geom="density_ridges_gradient",calc_ecdf = TRUE,bandwidth=2)
a=a+theme_ridges()+theme(text=element_text(size=18))+theme(legend.position = "none")+labs(x="2010-2030 SCC Distribution ($ per ton CO2), Central 90% of Distribution",y="Publication Year")
a=a+geom_point(data=distplotmeans,aes(x=meanval),col="magenta",size=3)
a

#look at differences in distribution for "Sensitivity Analysis" vs "Empirical Improvement / Framework Expansion
dat$`Empirical Improvement or Sensitvity Analysis?`=fct_collapse(dat$`Empirical Improvement or Sensitvity Analysis?`,"Empirical Improvement"=c("Empirical improvement","Empirical Improvement","Empriical Improvement"),"Sensitvity Analysis"=c("Sensitivity analysis","Sensitivity Analysis"))

distplot=dist[-which(dist$draw<quantile(dist$draw,0.005)|dist$draw>quantile(dist$draw,0.995)|dist$year<=2010|dist$year>2100|is.na(dist$year)),]
distplot$type=dat$`Empirical Improvement or Sensitvity Analysis?`[distplot$row]
distplot=distplot%>%filter(type!="Other"&year>2009)

cols=met.brewer("Isfahan2",3,type="discrete")

a=ggplot(distplot,aes(x=draw,y=type,group=type,fill=type))+geom_boxplot()+facet_grid(yeargroup~.)
a=a+theme_bw()+theme(text=element_text(size=18),strip.background =element_rect(fill="white"),legend.position="none")
a=a+labs(x="SCC ($ per ton CO2)",y="")+scale_fill_manual(values=cols)
a

#distributions of strucutral changes - with and without structural change

#concatenate Carbon Cycle and Climate System changes into Earth System changes
dat$'Earth System'=as.character(ifelse(dat$`Carbon Cycle`==1,1.0,ifelse(dat$`Climate Model`==1,1.0,0)))

struc=dat%>%
  dplyr::select("Earth System","Tipping Points":"Learning")%>%
  replace(is.na(.),0)

diststruc=cbind(dist,struc[dist$row,])

diststruc=diststruc%>%dplyr::rename("Tipping Points: Climate"="Tipping Points","Tipping Points: Damages"="Tipping Points2","Limited Substitutability"="Limitedly-Substitutable Goods","Distributional Weighting"="Inequality Aversion")

diststruc=diststruc%>%
  pivot_longer(cols="Earth System":"Learning",names_to="StructuralChange",values_to="Presence")

diststruc$Presence=fct_collapse(diststruc$Presence,No=c("-1.0","0"),Yes=c("1.0","Calibrated",1))

diststrucdensities=diststruc%>%
  filter(yeargroup=="2010-2030")%>%
  group_by(StructuralChange)%>%
  filter(Presence=="Yes"&draw>0)%>%
  mutate(logscc=log(draw))%>%
  group_modify(~ ggplot2:::compute_density(.x$logscc, NULL,bw=0.4))%>%
  dplyr::rename(logscc=x)

#find rows with no structural changes at all for "reference" density
reference=which(apply(struc,MARGIN=1,FUN=function(x) sum(x=="0"))==9)
referencedensity=diststruc%>%
  filter(yeargroup=="2010-2030"&row%in%reference&draw>0)%>%
  mutate(logscc=log(draw))%>%
  group_modify(~ ggplot2:::compute_density(.x$logscc, NULL,bw=0.4))%>%
  dplyr::rename(logscc=x)

changes=unique(diststruc$StructuralChange);referencedensity$StructuralChange=changes[1]
for(i in 2:length(changes)){temp=referencedensity[,1:6];temp$StructuralChange=changes[i];referencedensity=rbind(referencedensity,temp)}
referencedensity=referencedensity[,c(7,1:6)]

diststrucdensities$type="Changed";referencedensity$type="Reference"
diststrucdensities=dplyr::bind_rows(diststrucdensities,referencedensity)

diststrucdensities$StructuralChange=ordered(diststrucdensities$StructuralChange,levels=rev(c("Earth System" , "Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Inequality Aversion","Epstein-Zin","Learning" ,"Ambiguity/Model Uncertainty")))

#add number of papers and number of observations
diststruc$paper=dat$ID_number[diststruc$row]
papers=diststruc%>%
  group_by(StructuralChange)%>%
  filter(Presence=="Yes"&yeargroup=="2010-2030")%>%
  dplyr::summarise(npapers=length(unique(paper)),n=length(unique(row)))

diststrucdensities$type=as.factor(diststrucdensities$type)

breaks_ln=c(0,2.5,5,7.5,10)
a=ggplot(diststrucdensities,aes(x=logscc,y=StructuralChange,height=density,group=interaction(type,StructuralChange),fill=as.factor(type)))+geom_density_ridges(stat = "identity",scale=0.92,lwd=1)
a=a+theme_ridges()+theme_bw()+theme(text=element_text(size=12),legend.position="right")+labs(x="2010-2030 SCC (2020 $ per ton CO2)",y="",fill="")
a=a+scale_x_continuous(limits=c(0,10),breaks=breaks_ln,labels=c(round_any(exp(breaks_ln[1:3]),10),round_any(exp(breaks_ln[4:5]),100)))+scale_fill_manual(values=c("steelblue4",NA),na.value=NA)
a=a+geom_text(data=papers,aes(label=paste0("n=",npapers," (",n,")"),y=StructuralChange,x=9.2),inherit.aes=FALSE,size=4,nudge_y=0.35)

#add in figure 2 data
load(file="outputs/expert_survey_data_products/fig2surveydata.Rdat")

# #use reference value given in survey of $41.50 to convert % changes into SCC values
# refval=41.5
# helper=function(x) log((1+x/100)*refval)
# for(i in 2:10) fig2dat_vals[,i]=helper(fig2dat_vals[,i])
# fig2dat_vals=pivot_longer(fig2dat_vals,cols=2:10,names_to="StructuralChange",values_to="logscc")
#
# a=a+geom_density_ridges2(data=fig2dat_vals,aes(x=logscc,y=StructuralChange,group=StructuralChange),fill="darkgoldenrod",inherit.aes=FALSE,stat = "binline", binwidth = 0.1,scale=0.88)
# a=a+geom_vline(xintercept=log(refval),lty=3,col="#36c687",lwd=1)
# a
#
#add graph to the side showing histograms of assessed quality
fig2dat_qual=fig2dat_qual%>%
  pivot_longer(2:10,names_to="StructuralChange",values_to="ans")%>%
  filter(!is.na(ans))
fig2dat_qual$ans=as.factor(fig2dat_qual$ans)
fig2dat_qual$ans=fct_relevel(fig2dat_qual$ans,"Strongly Disagree","Disagree","Neither Agree nor Disagree","Agree","Strongly Agree")
fig2dat_qual$ans=fct_recode(fig2dat_qual$ans,Neutral="Neither Agree nor Disagree")
fig2dat_qual=fig2dat_qual%>%
  group_by(StructuralChange,ans)%>%
  dplyr::summarise(tot=n())
fig2dat_qual$StructuralChange=fct_relevel(fig2dat_qual$StructuralChange,levels(diststrucdensities$StructuralChange))

b=ggplot(fig2dat_qual,aes(y=StructuralChange,x=tot,fill=ans,group=StructuralChange))
b=b+geom_bar(position="fill",stat="identity")+theme_bw()+scale_fill_manual(values=c('#7b3294','#c2a5cf','#f7f7f7','#a6dba0','#008837'))
b=b+theme(text=element_text(size=12),legend.position="right")
b=b+labs(x="To What Extent Do You Agree with the Statement \"Papers that Include This Structural Change in the Current Literature Provide a Better SCC than Those That Exclude It\"",y="",fill="")
#
# b=b+geom_density_ridges2(fill="#f7e057",col="black",stat = "binline", binwidth = 1,scale=0.95)
# b=b+theme_bw()+theme(text=element_text(size=12),axis.text.y=element_blank(),axis.text.x = element_text(angle = 90,vjust = 1, hjust=0.5))+labs(y="",x="")
# b=b+scale_x_discrete(labels = function(x) str_wrap(x, width = 8))

a/b+plot_annotation(tag_levels="A",theme=theme(plot.tag=element_text(size=16)))


#calculate means
means=diststrucdensities%>%
  group_by(type,StructuralChange)%>%
  dplyr::summarize(mean=exp(sum(logscc*density)/sum(density)))

#calculate 90th percentile
upper=diststruc%>%
  group_by(StructuralChange)%>%
  dplyr::filter(Presence=="Yes")%>%
  dplyr::filter(draw>0)%>%
  dplyr::summarize(upper75=exp(quantile(log(draw),0.75)))

referenceupper=dist%>%
  filter(row%in%reference&draw>0)%>%
  summarize(upper75=exp(quantile(log(draw),0.75)))


#"balance table" for interpreting structural scc distributions

rows=diststruc%>%
  group_by(StructuralChange)%>%
  filter(Presence==1)%>%
  dplyr::summarise(rows=unique(row))

for(i in 1:length(changes)){
  relrows=rows%>%filter(StructuralChange==changes[i])%>%select(rows)
  relrows=as.numeric(relrows$rows)

  temp=dat[relrows,]
  summary=c(round(mean(as.numeric(temp$`SCC Year`,na.rm=T))),round(sd(as.numeric(temp$`SCC Year`,na.rm=T)),1),round(mean(temp$discountrate,na.rm=T),2),round(sd(temp$discountrate,na.rm=T),2))

  if(i==1) balance=summary; if(i>1) balance=rbind(balance,summary)
}
temp=dat[reference,]
summary_ref=c(round(mean(as.numeric(temp$`SCC Year`),na.rm=T)),round(sd(as.numeric(temp$`SCC Year`),na.rm=T),1),round(mean(temp$discountrate,na.rm=T),2),round(sd(temp$discountrate,na.rm=T),2))
balance=rbind(balance,summary_ref)
balance=as.data.frame(balance);balance$StructuralChange=c(changes,"Reference")
colnames(balance)[1:4]=c("SCC Year","SCC Year SD","Discount Rate","Discount Rate SD");rownames(balance)=1:(length(changes)+1)
write.csv(balance,"outputs/distribution_balancetable.csv")

###----------variance due to within paper vs between paper ------------
dist$paper=as.factor(dat$ID_number[dist$row])
mod=feols(I(log(draw))~1|paper,data=dist[-which(dist$draw<=0),])
mod2=feols(I(log(draw))~1|row,data=dist[-which(dist$draw<=0),])
r2(mod);r2(mod2)
