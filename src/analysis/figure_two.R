library(tidyverse)
library(see)

dat=read.csv("data/expert_survey/data_SCC-expert-survey_final_anonymous.csv")

fig2dat=dat%>%
  select(contains(c("SCC_Lit","SCC_True_","Interview.number")))%>%
  pivot_longer(!"Interview.number..ongoing.")%>%
  mutate(type=substr(name,5,7))%>%
  mutate(quantile=as.factor(substr(name,9,11)))%>%
  mutate(quantile=fct_collapse(quantile,lower=c("2p5","_2p"),upper=c("7p5","_97"),central=c("Cen","_Ce")))%>%
  filter(!is.na(value))%>%
  pivot_wider(id_cols=c("Interview.number..ongoing.",type),values_from=value,names_from=quantile)
colnames(fig2dat)[1]="id"

#create jittered x variable
fig2dat$xj=jitter(as.numeric(as.factor(fig2dat$type)),amount=0.1)

b=ggplot(fig2dat)
b=b+geom_violinhalf(data=fig2dat%>%filter(type=="Tru"),aes(x=type,y=central),position=position_nudge(0.17),fill="#FF495C")
b=b+geom_violinhalf(data=fig2dat%>%filter(type=="Lit"),aes(x=type,y=central),position=position_nudge(-0.17),fill="#E5E059",flip=TRUE)
b=b+theme_bw()+scale_x_discrete(labels=c("Lit"="Literature","Tru"="Comprehensive"))+theme(text=element_text(size=16))
b=b+labs(x="",y="SCC ($ per ton CO2)")+scale_fill_manual(values=c("#E5E059","#FF495C"),guide="none")
b=b+geom_segment(data=pivot_wider(fig2dat,id_cols=id,names_from=type,values_from = c(central,xj)),aes(x=xj_Lit,xend=xj_Tru,y=central_Lit,yend=central_Tru),col="grey50")
b=b+geom_errorbar(aes(x=xj,ymin=lower,ymax=upper),lwd=0.09,alpha=0.3,col="grey85")
b=b+scale_y_continuous(limits=c(-30,1200), expand=c(0, 0))
b=b+geom_point(aes(x=xj,y=central,fill=type),size=2,pch=21,col="black")

pdf(file="figures/Science Revision/figure2a.pdf",width=11,height=8.5)
a1+a+plot_layout(nrow=2,heights=c(1,5))
dev.off()