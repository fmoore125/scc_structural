library(ggplot2)
library(data.table)
library(ggridges)
library(fixest)
library(plyr)
library(MetBrewer)
library(lfe)
library(forcats)
library(zoo)

dist=fread(file="outputs/distribution_v2.csv")
dist_weighted=fread(file="outputs/distribution_coauthorweighted_v2.csv")
dist_weighted_citations=fread(file="outputs/distribution_citationweighted_v2.csv")

source("src/data_cleaining_scripts/cleaning_master.R")

#make plot of distribution taking out different sources of varition
#1. take out means
dist_demeaned=dist$draw-mean(dist$draw)
#2 take out paper fixed-effects
dist$paper=dat$ID_number[dist$row]
femod=felm(draw~1|paper,data=dist)
dist_paperfe=femod$residuals
#3. take out row fixed-effects
femod_row=felm(draw~1|row,data=dist)
dist_rowfe=femod_row$residuals

#plot the quantiles of the distributions
relprobs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)

quants=data.frame(Demeaned=quantile(dist_demeaned,relprobs),PaperFE=quantile(dist_paperfe,relprobs),RowFE=quantile(dist_rowfe,relprobs))
quants$probs=relprobs
quants=pivot_longer(quants,cols=1:3,names_to="Model",values_to="Quantiles")
quants$Model=fct_recode(quants$Model,'De-Meaned'="Demeaned",'Paper Means Removed'="PaperFE","Observables Removed"="RowFE")

a=ggplot(quants,aes(x=Model,y=Quantiles))
a=a+geom_hline(yintercept = 0,lty=3)
temp=quants%>%filter(probs%in%c(0.05,0.95))%>%pivot_wider(id_cols = Model,names_from=probs,values_from=Quantiles)
colnames(temp)=c("Model","Start","End")
a=a+geom_segment(aes(y=Start,yend=End,xend=Model),data=temp,lwd=1.25)
temp=quants%>%filter(probs%in%c(0.25,0.75))%>%pivot_wider(id_cols = Model,names_from=probs,values_from=Quantiles)
colnames(temp)=c("Model","Start","End")
a=a+geom_rect(aes(ymin=Start,ymax=End),data=temp,xmin=c(0.75,1.75,2.75),xmax=c(1.25,2.25,3.25),inherit.aes=FALSE,fill="#3d427b",color=NA)
a=a+geom_segment(data=quants%>%filter(probs==0.5),x=c(0.75,1.75,2.75),xend=c(1.25,2.25,3.25),aes(yend=Quantiles),col="#edc123",lwd=1.5)
a=a+theme_bw()+labs(x="",y="Residual SCC Distribution ($ per ton CO2)")+theme(text=element_text(size=18))

#full distribution
dist$year=as.numeric(dat$`SCC Year`[dist$row])
dist$yeargroup=cut(dist$year,breaks=c(1990,2010,2030,2070,2100,2400))
dist$yeargroup=fct_recode(dist$yeargroup,'<2010'="(1.99e+03,2.01e+03]",'2010-2030'="(2.01e+03,2.03e+03]",'2030-2070'="(2.03e+03,2.07e+03]",'2070-2100'="(2.07e+03,2.1e+03]",'>2100'="(2.1e+03,2.4e+03]")

distplot=dist[-which(dist$draw<quantile(dist$draw,0.005)|dist$draw>quantile(dist$draw,0.995)|dist$year<=2010|is.na(dist$year)|dist$year>2100),]
distplot$y=ifelse(distplot$yeargroup%in%c('2010-2030','2030-2070'),-0.004,-0.008)

a=ggplot(distplot,aes(x=draw,fill=yeargroup,y=y))+geom_boxplot(width=0.006)+geom_density(aes(x=draw,fill=yeargroup),inherit.aes=FALSE,adjust=3)+facet_grid(yeargroup~.)
a=a+theme_bw()+labs(x="SCC ($ per ton CO2)",y="")+scale_fill_discrete(guide="none")+theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),text=element_text(size=18),strip.background =element_rect(fill="white"))
a=a+geom_hline(yintercept = 0)+scale_x_continuous(breaks=c(-100,0,100,200,300,400,500,1000,1500,2000))

#add IWG values
iwg=read.csv("C:/Users/fmoore/Documents/GitHub/scc_structural/outputs/iwgruns.csv",row.names=1)
iwg=iwg[,-(1+grep("PAGE.59",colnames(iwg)):dim(iwg)[2])]
#just keep central years from ranges - 2020 and 2050
iwg=iwg[,which(iwg[1,]==2020|iwg[1,]==2050)]
iwgdist=data.frame(early=as.numeric(as.character(unlist((iwg[5:dim(iwg)[1],which(iwg[1,]==2020)])))))
iwgdist$late=as.numeric(as.character(unlist((iwg[5:dim(iwg)[1],which(iwg[1,]==2050)]))))
iwgdist=pivot_longer(iwgdist,cols=1:2,names_to=c("yeargroup"))
extra=data.frame(yeargroup=c('2070-2100'),value=NA)
iwgdist=rbind(iwgdist,extra)
iwgdist$yeargroup=fct_recode(iwgdist$yeargroup,"2010-2030"="early","2030-2070"="late")
iwgdist=iwgdist[-which(iwgdist$value<quantile(iwgdist$value,0.01,na.rm=T)|iwgdist$value>quantile(iwgdist$value,0.99,na.rm=T)),]

a=a+geom_boxplot(aes(x=value,alpha=0,y=-0.01),width=0.0035,data=iwgdist,fill="white",inherit.aes=FALSE)+scale_alpha_continuous(guide="none")
ann_text=data.frame(text=c("-- IWG 2020","-- IWG 2050"),yeargroup=factor(c("2010-2030","2030-2070"),levels=levels(iwgdist$yeargroup)),x=c(690,690))
a=a+geom_text(data=ann_text,aes(label=text,y=-0.01,x=x))

x11()
a

# #different figure looking at publication date
# 
# #downsample distribution so avoid crashing computer
# library(ggridges)
# 
# dist$pubyear=dat$Year[dist$row]
# distplot=dist[which(dist$year%in%2010:2030),]
# distplot=distplot[-which(distplot$draw<quantile(distplot$draw,0.05)|distplot$draw>quantile(distplot$draw,0.95)),]
# distplotmeans=distplot%>%group_by(pubyear)%>%summarize(meanval=mean(draw))
# 
# a=ggplot(distplot,aes(x=draw,y=as.factor(pubyear),group=pubyear,fill=pubyear))
# a=a+stat_density_ridges(geom="density_ridges_gradient",calc_ecdf = TRUE,bandwidth=2)
# a=a+theme_ridges()+theme(text=element_text(size=18))+theme(legend.position = "none")+labs(x="2010-2030 SCC Distribution ($ per ton CO2), Central 90% of Distribution",y="Publication Year")
# a=a+geom_point(data=distplotmeans,aes(x=meanval),col="magenta",size=3)
# a

#look at differences in distribution for "Sensitivity Analysis" vs "Empirical Improvement / Framework Expansion
dat$`Empirical Improvement or Sensitvity Analysis?`=fct_collapse(dat$`Empirical Improvement or Sensitvity Analysis?`,"Empirical Improvement"=c("Empirical improvement","Empirical Improvement","Empriical Improvement"),"Sensitvity Analysis"=c("Sensitivity analysis","Sensitivity Analysis"))

distplot=dist[-which(dist$draw<quantile(dist$draw,0.005)|dist$draw>quantile(dist$draw,0.995)|dist$year<=2010|dist$year>2100|is.na(dist$year)),]
distplot$type=dat$`Empirical Improvement or Sensitvity Analysis?`[distplot$row]
distplot=distplot%>%filter(type!="Other"&year>2009)

a=ggplot(distplot,aes(x=draw,y=type,group=type,fill=type))+geom_boxplot()+facet_grid(yeargroup~.)
a=a+theme_bw()+theme(text=element_text(size=18),strip.background =element_rect(fill="white"),legend.position="none")
a=a+labs(x="SCC ($ per ton CO2)",y="")+scale_fill_manual(values=c("#c56423","#e8ab32","#305962"))
a

#distributions of strucutral changes - with and without structural change
struc=dat%>%
  select("Carbon Cycle":"Learning")%>%
  replace(is.na(.),0)%>%
  mutate_all(funs(replace(.,which(.=="Calibrated"),1.0)))%>%
  mutate_all(as.numeric)

diststruc=cbind(dist,struc[dist$row,])

colnames(diststruc)[c(7:8,12)]=c("Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability")
diststruc$'Earth System'=ifelse(diststruc$`Carbon Cycle`==1.0,1.0,ifelse(diststruc$`Climate Model`==1.0,1.0,0))

diststruc=diststruc%>%
  select(draw,row,"Tipping Points: Climate":"Earth System")%>%
  pivot_longer(cols="Tipping Points: Climate":"Earth System",names_to="StructuralChange",values_to="Presence")

diststruc$Presence=replace(diststruc$Presence,which(diststruc$Presence==-1),0)

diststrucdensities=diststruc%>%
  group_by(StructuralChange)%>%
  filter(Presence==1&draw>0)%>%
  mutate(logscc=log(draw))%>%
  group_modify(~ ggplot2:::compute_density(.x$logscc, NULL,bw=0.4))%>%
  rename(logscc=x)

#find rows with no structural changes at all for "reference" density
reference=which(apply(struc,MARGIN=1,FUN=function(x) sum(x=="0"))==10)
referencedensity=diststruc%>%
  filter(row%in%reference&draw>0)%>%
  mutate(logscc=log(draw))%>%
  group_modify(~ ggplot2:::compute_density(.x$logscc, NULL,bw=0.4))%>%
  rename(logscc=x)

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
  filter(Presence==1)%>%
  dplyr::summarise(npapers=length(unique(paper)),n=length(unique(row)))

diststrucdensities$type=as.factor(diststrucdensities$type)

breaks_ln=c(0,2.5,5,7.5,10)
a=ggplot(diststrucdensities,aes(x=logscc,y=StructuralChange,height=density,group=interaction(type,StructuralChange),fill=as.factor(type)))+geom_density_ridges(stat = "identity",scale=0.92,lwd=1)
a=a+theme_ridges()+theme_bw()+theme(text=element_text(size=18))+labs(x="SCC (2020 $ per ton CO2)",y="",fill="")
a=a+scale_x_continuous(limits=c(0,10),breaks=breaks_ln,labels=c(round_any(exp(breaks_ln[1:3]),10),round_any(exp(breaks_ln[4:5]),100)))+scale_fill_manual(values=c("steelblue4",NA))
a=a+geom_text(data=papers,aes(label=paste0("n=",npapers," (",n,")"),y=StructuralChange,x=9.2),inherit.aes=FALSE,size=6,nudge_y=0.35)
a

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
