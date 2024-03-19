library(tidyverse)
library(data.table)

rfdist_dir="C:\\Users\\fmoore\\Documents\\GitHub2\\scc_structural\\scc_structural\\outputs\\synthetic_scc\\Structural SCC RF Experiments\\"

#load best prediction from random forest
load(paste0(rfdist_dir,"RFD_best.Rdata"))
rfdist=data.frame(dist=allsamp)

#plot the quantiles of the distributions
relprobs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)

truncmean=function(data,trim=0.001){
  highquant=quantile(data,1-trim)
  lowquant=quantile(data,trim)
  return(mean(data[which(data>lowquant&data<highquant)]))
}
pfuns=c(map(relprobs,~partial(quantile,probs=.x,na.rm=T)), truncmean)

summarytot=rfdist%>%
  summarize_at(vars(dist),funs(!!!pfuns))
colnames(summarytot)=c("lowest","min","lower","middle","upper","max","highest", "mu")

#load IWG and EPA SCC distributions for comparison

#add EPA values
eparuns=list.files("outputs/epa_scc",full.names = TRUE) #2020 distribution of co2 for 3 damage functions
epa=fread(eparuns[2]);epa=filter(epa,sector=="total")$scghg
epa=append(epa,fread(eparuns[1])$scghg);epa=append(epa,fread(eparuns[3])$scghg)
epa=data.frame(draw=epa)
epa_summary=epa%>%summarize_at(vars(draw),funs(!!!pfuns))
colnames(epa_summary)=c("lowest","min","lower","middle","upper","max","highest", "mu")

#add IWG values
iwg=read.csv("outputs/iwgruns.csv",row.names=1)
iwg=iwg[,-(1+grep("PAGE.59",colnames(iwg)):dim(iwg)[2])]
#just keep 2020 SCC year
iwg=iwg[,which(iwg[1,]==2020)]
iwgdist=data.frame(draw=as.numeric(as.character(unlist((iwg[5:dim(iwg)[1],which(iwg[1,]==2020)])))))
iwgdist_summary=iwgdist%>%
  summarize_at(vars(draw),funs(!!!pfuns))
colnames(iwgdist_summary)=c("lowest","min","lower","middle","upper","max","highest", "mu")

dists=rbind(summarytot,epa_summary,iwgdist_summary)
dists$y=c(3,2.5,2.25)
dists$group=c("Random Forest","EPA","IWG")

a=ggplot(dists)+coord_flip()+theme_bw()+
  geom_boxplot(aes(x=y,min=min,lower=lower,middle=middle,upper=upper,max=max,col=group),inherit.aes=FALSE,stat="identity",width=c(0.4,0.2,0.2),lwd=0.4)+
  scale_y_continuous(breaks=seq(0,1200,by=100),minor_breaks=c(-50,seq(0,400,by=25),seq(400,1200,by=50)), limits=c(-50,1300), expand=c(0,0))+
  geom_segment(aes(x=y,xend=y,y=lowest,yend=min,col=group),lty=2,lwd=0.4)+geom_segment(aes(x=y,xend=y,y=max,yend=highest,col=group),lty=2,lwd=0.4)+
  scale_color_manual(values=c("#233350","#3e778a","black"),guide=NULL)+labs(x="",y="2020 SCC ($ per ton CO2)")+geom_point(aes(x=y,y=mu,col=group),size=2)+
  theme(text=element_text(size=16),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.grid.minor.x = element_line(linewidth = 0.25), panel.grid.major.x = element_line(linewidth = 0.5,color="#959595"), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+
  annotate("text",x=c(3.04,2.55,2.3),y=c(650,400,150),label=c("Synthetic SCC","EPA","IWG"),col=c("black","#233350","#3e778a"),size=6)+
  geom_vline(xintercept=2.67,lty=4)

pdf(file="figures/Science Revision/figure5_a.pdf",width=9,height=4)
a
dev.off()
