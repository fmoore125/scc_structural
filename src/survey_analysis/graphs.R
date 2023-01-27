library(tidyverse)
library(forcats)
library(rriskDistributions)
library(EnvStats)
library(data.table)
library(vctrs)

source("src/analysis/find_distribution.R")

dat=read.csv("data/expert_survey/data_SCC-expert-survey_final_anonymous.csv")

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

#fit distributions to mean, 2.5th and 97.5th responses for all provided distributions
distwrapper=function(central,lower,upper){
  #creates a wrapper to use pmap on the fig1dat dataframe to fit distribution and take 1000 samples
    case=ifelse(is.na(lower),ifelse(is.na(upper),1,2),ifelse(is.na(upper),3,4))

  if(case==2){generate.exactmu(mu=central,qs=0.975,as=upper,1000)}
  else if(case==1){rep(central,1000)}
  else if(case==3){generate.exactmu(mu=central,qs=0.025,as=lower,1000)}
  else if(case==4){generate.exactmu(mu=central,qs=c(0.025,0.975),as=c(lower,upper),1000)}
}

## errors <- data.frame()
## for (ii in 1:nrow(fig1dat)) {
##     vals <- distwrapper(fig1dat$central[ii], fig1dat$lower[ii], fig1dat$upper[ii])
##     if (is.null(vals)) {
##         if (is.na(fig1dat$lower[ii]) && is.na(fig1dat$upper[ii]))
##             vals <- rep(fig1dat$central[ii], 1000)
##         else
##             vals <- rnorm(1000, fig1dat$central[ii], 0.5 * (fig1dat$upper[ii] - fig1dat$lower[ii]) / 1.96)
##     }

##     errors <- rbind(errors, data.frame(central=fig1dat$central[ii], mu=mean(vals),
##                                        lower=fig1dat$lower[ii], q025=quantile(vals, .025),
##                                        upper=fig1dat$upper[ii], q975=quantile(vals, .975)))
## }
## errors$type <- fig1dat$type

twodists=fig1dat%>%
    mutate(samp=pmap(list(central,lower,upper),function(central,lower,upper){
        vals <- distwrapper(central,lower,upper)
        if (is.null(vals))
            vals <- rnorm(1000, central, 0.5 * (upper - lower) / 1.96)
        vals
    }))%>%
  group_by(type)%>%
  dplyr::summarise(dist=vec_c(samp))%>%
  unnest(cols = c(dist))

fwrite(twodists,file="outputs/expert_survey_data_products/question1_distributions.csv")

#read in meta-analysis distribution, subset to 2010-2030
metadist=fread("outputs/distribution_v2.csv")
source("src/data_cleaining_scripts/cleaning_master.R")

metadist$SCCyear=dat$`SCC Year`[metadist$row]

metadist=metadist%>%filter(SCCyear>2010&SCCyear<2031)
metadist$type=3

#downsample distribution for ease of plotting
metadist=metadist[sample(1:nrow(metadist),0.2*nrow(metadist),replace=FALSE),]

#add boxplots to graph
a=ggplot(fig1dat,aes(x=xj,y=central,ymin=lower,ymax=upper,col=type))+geom_line(aes(group=id,x=xj,y=central),col="lightgrey")+geom_pointrange(lty=2)
a=a+coord_cartesian(ylim=c(-100,1500))+theme_bw()+theme(text=element_text(size=14))+labs(x="",y="2020 SCC ($ per ton CO2)")+scale_color_manual(values=c("#af2436","#51a154"),guide=FALSE)
a=a+scale_x_continuous(breaks=c(1,2), labels=c("Literature", "All Things Considered"), limits=c(0.5, 4.25))
a=a+geom_boxplot(data=twodists%>%filter(type=="Lit"&dist>quantile(dist,0.01)&dist<quantile(dist,0.99)),aes(y=dist,group=type,x=1),inherit.aes = FALSE,position=position_nudge(x=1.75),fill="#af2436",width=0.25)
a=a+geom_boxplot(data=twodists%>%filter(type=="Tru"&dist>quantile(dist,0.01)&dist<quantile(dist,0.99)),aes(y=dist,group=type,x=2),inherit.aes = FALSE,position=position_nudge(x=1.25),fill="#51a154",width=0.25)
a=a+geom_boxplot(data=metadist%>%filter(draw>quantile(draw,0.01)&draw<quantile(draw,0.99)),aes(y=draw,group=type,x=3),inherit.aes = FALSE,position=position_nudge(x=0.75),fill="slateblue3",width=0.25)
a

#process data for integration of question into Figure 2

codingvals=data.frame(bin=1:13,meanval=c(-200,-150,-75,-40,-20,-5,0,5,20,40,75,150,200))

convertbins=function(x,lookup=codingvals){
  return(ifelse(is.na(x),NA,lookup$meanval[x]))
}

fig2dat_vals=dat%>%
  select("Interview.number..ongoing.","Q2_MeanSCC_earth_system":"Q2_MeanSCC_distributional_weighting")%>%
  rename("ID"="Interview.number..ongoing.")%>%
  mutate(across(.cols="Q2_MeanSCC_earth_system":"Q2_MeanSCC_distributional_weighting",.fns=~convertbins(.)))

#simplify names of structural changes
strucs=c("Earth System","Tipping Points: Climate","Tipping Points: Damages","Persistent / Growth Damages","Limited Substitutability","Epstein-Zin","Ambiguity/Model Uncertainty","Learning","Inequality Aversion")
colnames(fig2dat_vals)[2:10]=strucs

codingqual=data.frame(bin=1:5,meanval=c("Strongly Disagree","Disagree","Neither Agree nor Disagree","Agree","Strongly Agree"))

fig2dat_qual=dat%>%
  select("Interview.number..ongoing.","Q2_Include_earth_system":"Q2_Include_distributional_weighting")%>%
  rename("ID"="Interview.number..ongoing.")%>%
  mutate(across(.cols="Q2_Include_earth_system":"Q2_Include_distributional_weighting",.fns=~convertbins(.,lookup=codingqual)))
colnames(fig2dat_qual)[2:10]=strucs

save(fig2dat_vals,fig2dat_qual,file="outputs/expert_survey_data_products/fig2surveydata.Rdat")
