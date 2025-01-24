## Performs validation checks on the process of translating quantile information into Monte Carlo distribution draws.

## setwd("~/research/scciams/scc_structural")

library(data.table)
library(ggridges)

source("src/data_cleaining_scripts/cleaning_master.R")

dat$ID_number=as.integer(dat$ID_number)
papers=unique(dat$ID_number)

dist <- read.csv("outputs/distribution_v2.csv")

## make some figures analyzing variance in distribution
mod=numeric(length=nrow(dist))
mod[c(grep("DICE",dat$`Base IAM (if applicable)`[dist$row]),grep("DICE",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="DICE"
mod[c(grep("FUND",dat$`Base IAM (if applicable)`[dist$row]),grep("FUND",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="FUND"
mod[c(grep("PAGE",dat$`Base IAM (if applicable)`[dist$row]),grep("PAGE",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="PAGE"
mod[which(mod==0)]="Other"

dist$mod=ordered(mod,levels=c("Other","PAGE","FUND","DICE"))

a=ggplot(dist[which(dist$draw>quantile(dist$draw,0.01)&dist$draw<quantile(dist$draw,0.99)),],aes(x=draw,y=mod,group=mod,fill=factor(stat(quantile))))
a=a+stat_density_ridges(geom="density_ridges_gradient",calc_ecdf = TRUE,quantiles=4,quantile_lines = TRUE,bandwidth=2,scale=0.9)
a=a+theme_bw()+labs(x="SCC ($ per ton CO2)",y="Base or Calibration Model")+ scale_fill_viridis_d(name = "Quartiles")
a=a+theme(text=element_text(size=20))

## look at differences in distribution by structural parameters, year, discount rate etc
dist$year=as.numeric(dat$`SCC Year`[dist$row])
dist$yeargroup=cut(dist$year,breaks=c(1990,2010,2050,2090,2400))
dist$yeargroup=fct_recode(dist$yeargroup,'<2010'="(1.99e+03,2.01e+03]",'2010-2050'="(2.01e+03,2.05e+03]",'2050-2090'="(2.05e+03,2.09e+03]",'>2090'="(2.09e+03,2.4e+03]")

dist$yeargroup=ordered(dist$yeargroup,levels=c("<2010","2010-2050","2050-2090",">2090"))

a=ggplot(dist[which(dist$draw>quantile(dist$draw,0.01)&dist$draw<quantile(dist$draw,0.99)),]%>%filter(!is.na(yeargroup)),aes(x=draw,y=yeargroup,group=yeargroup,fill=factor(stat(quantile))))
a=a+stat_density_ridges(geom="density_ridges_gradient",calc_ecdf = TRUE,quantiles=4,quantile_lines = TRUE,bandwidth=2,scale=0.9)
a=a+theme_bw()+labs(x="SCC ($ per ton CO2)",y="SCC Year")+ scale_fill_viridis_d(name = "Quartiles")
a=a+theme(text=element_text(size=20))
x11()
a

## plot of structural parameters

struc=dat%>%select('Carbon Cycle':'Learning')%>% #note- this is dropping one of the strucutral columns that have very few rows in them
    replace(is.na(.),"No") %>%
    replace(.=="1.0","Yes")
colnames(struc)[which(colnames(struc)=="Tipping Points")]="Climate Tipping Points"
colnames(struc)[which(colnames(struc)=="Tipping Points2")]="Damage Tipping Points"

dist=cbind(dist,struc[dist$row,])
dist_struc=pivot_longer(dist,cols='Carbon Cycle':'Learning',names_to="StructuralChange",values_to = "Changed")%>%
    filter(Changed%in%c("Yes","No"))#filter out some odd PAGE values where tipping points were removed from model

a=ggplot(dist_struc[which(dist_struc$draw>quantile(dist_struc$draw,0.01)&dist_struc$draw<quantile(dist_struc$draw,0.99)),],aes(x=draw,y=Changed))+
  geom_density_ridges(bandwidth=2)+facet_wrap(~StructuralChange)+theme_bw()+labs(x="SCC ($ per ton CO2)",y="Change Present?")

## identify high-leverage papers - how much does mean change if paper is dropped? For 2010-2030 period

'%notin%'=Negate('%in%')

meanchange=numeric(length=length(papers))
meanval=dist%>%filter(year%in%2010:2030)%>%dplyr::summarize(mean(draw))
for(i in 1:length(meanchange)){
  print(i)
  rows=which(dat$ID_number==papers[i])
  meanchange[i]=dist%>%filter(year%in%2010:2030&row%notin%rows)%>%summarize(mean(draw)-meanval)
}
meanchange=data.frame(ID_number=papers,change=unlist(meanchange))

bibs=dat%>%
  select(ID_number,Reference)%>%
  distinct()

meanchange=merge(meanchange,bibs)
meanchange=meanchange%>%arrange(desc(abs(change)))
write.csv(meanchange,file="outputs/meanleverage.csv")

## make graph showing high-leverage papers

mc=read.csv("outputs/meanleverage.csv")
mc=mc[1:15,]

mc$Reference=ordered(mc$Reference,levels=as.character(mc$Reference))

a=ggplot(mc,aes(y=change,x=Reference))+geom_bar(stat="identity")
a=a+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(x="",y="Change in Mean SCC After Dropping Paper ($)",fill="Paper Focus")
#a=a+scale_fill_manual(values=c("#335361","#4ec1a2","#c46692","#f57d51"))
a

###----compare distributions with and without co-author weighting ---------------

dist=fread(file="outputs/distribution_v2.csv")
dist_weighted=fread(file="outputs/distribution_coauthorweighted_v2.csv")
dist_weighted_citations=fread(file="outputs/distribution_citationweighted_v2.csv")

#drop outlier from Nordhaus 2019
dist=dist[-which(dist$row==1251),];dist_weighted=dist_weighted[-which(dist_weighted$row==1251),];dist_weighted_citations=dist_weighted_citations[-which(dist_weighted_citations$row==1251),]

#compare quantiles
breaks=c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
quantiles=data.frame(quantilesplits=breaks,unweighted=quantile(dist$draw,breaks),coauthorweights=quantile(dist_weighted$draw,breaks,na.rm=T),citationweights=quantile(dist_weighted_citations$draw,breaks,na.rm=T))
write.csv(quantiles,file="outputs/quantiles_differentweightings.csv")

###--------discount rate over time ----------------

a=ggplot(dat%>%group_by(ID_number)%>%summarise(Year=Year[1],discountrate=mean(discountrate,na.rm=T)),aes(x=Year,y=discountrate))+geom_point()+geom_smooth(method="lm")
a=a+theme_bw()+labs(x="Publication Year",y="Discount Rate")
a

