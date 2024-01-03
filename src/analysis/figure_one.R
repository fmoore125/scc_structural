library(tidyverse)
library(ggridges)
library(patchwork)

dist=fread(file="outputs/distribution_v2_Jan2024.csv")

source("src/data_cleaining_scripts/cleaning_master.R")

#drop outlier Nordhaus row
todrop=which(dat$`Central Value ($ per ton CO2)`>70000)
if(is.finite(todrop)) dist=dist[-which(dist$row==todrop),]

#add year, discount rate, pub year, and damages information
dist$year=as.numeric(dat$`SCC Year`[dist$row])
dist$dr=as.numeric(dat$discountrate[dist$row])
dist$pubyear=as.numeric(dat$Year[dist$row])
dist$damages=dat$`Damage Function Info: Model, Commonly-Used Function, or Function`[dist$row]
dist$damages[which(is.na(dist$damages))]=as.character(dat$`Base IAM (if applicable)`[dist$row[which(is.na(dist$damages))]])

#visualize distributions for 2010-2030, 2030-2070, 2070-2100
dist$yeargroup=cut(dist$year,breaks=c(2009,2030,2070,2101),labels=c("2010-2030","2031-2070","2071-2100"))

a1=ggplot(dist%>%filter(complete.cases(dist$yeargroup)&dist$yeargroup=="2010-2030"))+geom_density(aes(group=yeargroup,x=draw),adjust=3,lwd=0.75)
a1=a1+scale_x_continuous(breaks=c(-100,0,100,200,300,400,500,1000,1500),minor_breaks=c(-50,seq(0,175,by=25),seq(250,450,by=50),seq(600, 1100, by=100)), limits=c(-100,1100), expand=c(0, 0))
a1=a1+theme_bw()
a1=a1+theme(legend.position =c(0.8,0.8),text=element_text(size=16),strip.background =element_rect(fill="white"),axis.text = element_blank(),axis.ticks=element_blank(), plot.margin = unit(c(1,1,0,1), "cm")) 
a1=a1+labs(y="",x="")

#retain 2010-2030 values as 2020 equivalent SCC
distplot=dist[which(dist$year%in%c(2010:2030)),]

#create DR and pub year groups
distplot$drgroup=cut(distplot$dr,breaks=c(0,2.5,12),na.rm=T,labels=c("<2.5",">=2.5"))
distplot$pubyeargroup=cut(distplot$pubyear,breaks=c(2000,2009,2016,2022),labels=c("2000-2009","2010-2015","2016-2021"))
#integrate some missing damage function info from Base IAM column
#simplify damages
distplot$damages=as.factor(distplot$damages)
distplot$damages=fct_collapse(distplot$damages,DICE=levels(distplot$damages)[c(grep("DICE",levels(distplot$damages)),grep("RICE",levels(distplot$damages)))],FUND=levels(distplot$damages)[grep("FUND",levels(distplot$damages))],PAGE=levels(distplot$damages)[grep("PAGE",levels(distplot$damages))])
distplot$damages=fct_recode(distplot$damages,HowardSterner="HowardSterner (0.007438*T^2)")

#plot the quantiles of the distributions
relprobs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)

truncmean=function(data,trim=0.001){
  highquant=quantile(data,1-trim)
  lowquant=quantile(data,trim)
  return(mean(data[which(data>lowquant&data<highquant)]))
}
pfuns=c(map(relprobs,~partial(quantile,probs=.x,na.rm=T)), truncmean)

summarytot=distplot%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(summarytot)=c("lowest","min","lower","middle","upper","max","highest", "mu")
summarytot=cbind(data.frame(group="Full Distribution"),summarytot)

summarydist_dr=distplot%>%
  drop_na(drgroup)%>%
  group_by(drgroup)%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(summarydist_dr)=c("group","lowest","min","lower","middle","upper","max","highest", "mu")
#summarydist_dr$y=c(-0.01-0.0035*c(0:1),NA)

summarydist_pubyear=distplot%>%
  drop_na(pubyeargroup)%>%
  group_by(pubyeargroup)%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(summarydist_pubyear)=c("group","lowest","min","lower","middle","upper","max","highest", "mu")
#summarydist_pubyear$y=c(-0.02-0.0035*c(0:2),NA)

summarydist_dam=distplot%>%
  filter(damages%in%c("DICE","FUND","PAGE","HowardSterner","Weitzman"))%>%
  group_by(damages)%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(summarydist_dam)=c("group","lowest","min","lower","middle","upper","max","highest", "mu")

#get variation by model structure and reference distribution
dat$'Earth System'=as.character(ifelse(dat$`Carbon Cycle`%in%c("1.0",1),1.0,ifelse(dat$`Climate Model`%in%c("1.0",1),1.0,0)))

struc=dat%>%
  select("Earth System","Tipping Points":"Learning")%>%
  replace(is.na(.),0)

diststruc=cbind(dist,struc[dist$row,])

diststruc=diststruc%>%dplyr::rename("Tipping Points: Climate"="Tipping Points","Tipping Points: Damages"="Tipping Points2","Limited Substitutability"="Limitedly-Substitutable Goods")

diststruc=diststruc%>%
  pivot_longer(cols="Earth System":"Learning",names_to="StructuralChange",values_to="Presence")

diststruc$Presence=fct_collapse(diststruc$Presence,No=c("-1.0","0","-1"),Yes=c("1.0","Calibrated",1))

diststrucdensities=diststruc%>%
  group_by(StructuralChange)%>%
  filter(year%in%c(2010:2030)&Presence=="Yes")%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(diststrucdensities)=c("group","lowest","min","lower","middle","upper","max","highest", "mu")

#find rows with no structural changes at all for "reference" density
reference=which(apply(struc,MARGIN=1,FUN=function(x) sum(x=="0"))==9)
referencedensity=diststruc%>%
  filter(year%in%c(2010:2030)&row%in%reference)%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(referencedensity)=c("lowest","min","lower","middle","upper","max","highest", "mu")
referencedensity=cbind(data.frame(group="Reference"),referencedensity)

densities=rbind(summarytot,summarydist_dr,summarydist_pubyear,summarydist_dam,diststrucdensities,referencedensity)
densities$group=factor(densities$group, levels=rev(c("Full Distribution","<2.5",">=2.5","2000-2009","2010-2015","2016-2021","DICE","FUND","PAGE","HowardSterner","Weitzman","Reference","Earth System","Tipping Points: Climate","Tipping Points: Damages","Limited Substitutability","Persistent / Growth Damages","Inequality Aversion","Epstein-Zin","Learning","Ambiguity/Model Uncertainty")))
ymin=-100;ymax=1100
densities=densities%>%mutate(across(lowest:mu, function(x) ifelse(x<ymin,ymin,ifelse(x>ymax,ymax,x))))

#number of papers and number of observations for each structural change
diststruc$paper=dat$ID_number[diststruc$row]
papers=diststruc%>%
  group_by(StructuralChange)%>%
  filter(Presence=="Yes"&year%in%c(2010:2030))%>%
  dplyr::summarise(npapers=length(unique(paper)),n=length(unique(row)))
colnames(papers)[1]=c("group")

distplot$paper=dat$ID_number[distplot$row]

papers_dr=distplot%>%
  group_by(drgroup)%>%
  dplyr::summarize(npapers=length(unique(paper)),n=length(unique(row)))
colnames(papers_dr)[1]=c("group")

papers_pubyear=distplot%>%
  group_by(pubyeargroup)%>%
  dplyr::summarize(npapers=length(unique(paper)),n=length(unique(row)))
colnames(papers_pubyear)[1]=c("group")

papers_damagegroup=distplot%>%
  filter(damages%in%c("DICE","FUND","PAGE","HowardSterner","Weitzman"))%>%
  group_by(damages)%>%
  dplyr::summarize(npapers=length(unique(paper)),n=length(unique(row)))
colnames(papers_damagegroup)[1]=c("group")

papers_ref=data.frame(group="Reference",npapers=length(unique(distplot$paper[which(distplot$row%in%reference)])),n=length(unique(distplot$row[which(distplot$row%in%reference)])))

paperstot=data.frame(group="Full Distribution",npapers=length(unique(distplot$paper)),n=length(unique(distplot$row)))

papersfull=rbind(paperstot, papers_dr[1:2,],papers_pubyear[1:3,],papers_damagegroup,papers_ref,papers)


#Figure 1b

a=ggplot(densities)+coord_flip()+theme_bw()
a=a+geom_boxplot(aes(x=group,min=min,lower=lower,middle=middle,upper=upper,max=max,col=group),inherit.aes=FALSE,stat="identity")
a=a+scale_y_continuous(breaks=c(-100,0,100,200,300,400,500,1000,1500),minor_breaks=c(-50,seq(0,175,by=25),seq(250,450,by=50),seq(600, 1100, by=100)), limits=c(-100,1100), expand=c(0, 0))
a=a+geom_segment(aes(x=group,xend=group,y=lowest,yend=min,col=group),lty=2)+geom_segment(aes(x=group,xend=group,y=max,yend=highest,col=group),lty=2)
a=a+geom_point(aes(x=group,y=mu,col=group))
a=a+annotate("text",x=c(19.4,17.5,13,2.5),y=1000,label=c("Discount\nRate","Publication\nYear","Damages","Model\nStructure"))
a=a+geom_vline(xintercept=c(20.5,18.5,15.5,10.5))
a=a+scale_color_manual(values=c("Full Distribution"="black","<2.5"="#253494", ">=2.5"="#41b6c4","2000-2009"="#fbb4b9","2010-2015"="#f768a1","2016-2021"="#7a0177","DICE"='#fed976',"FUND"='#feb24c',"PAGE"='#fd8d3c',"HowardSterner"='#f03b20',"Weitzman"='#bd0026',"Reference"="grey50","Earth System"="#00B7A7","Tipping Points: Climate"="#554258","Tipping Points: Damages"="#943D67","Limited Substitutability"="#C97B72","Persistent / Growth Damages"="#FFCD12","Inequality Aversion"="#3F9127","Epstein-Zin"="#0A5755","Learning"="#39245D","Ambiguity/Model Uncertainty"="#A40000"))
a=a+theme(legend.position = "none",text=element_text(size=16),strip.background =element_rect(fill="white"),plot.margin = unit(c(0,1,1,1), "cm"))
a=a+labs(x="",y="2010-2030 SCC ($ per ton CO2)")
a=a+geom_text(data=papersfull,aes(label=paste0("n=",npapers," (",n,")"),x=group,col=group),y=700,position=position_nudge(x=0.25))

pdf(file="figures/Science Revision/figure1_full.pdf",width=11,height=8.5)
a1+plot_spacer()+a+plot_layout(nrow=3,heights=c(1,-0.5,5))+ plot_annotation(theme = theme(plot.margin = margin()))
dev.off()

#Figure 1b - Variance Decomposition


