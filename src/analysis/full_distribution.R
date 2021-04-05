library(data.table)
library(readxl)
library(tidyverse)
library(ggridges)

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

#start by generating distributions for each row
dists=list()
for (ii in 1:nrow(dat)) {
  print(ii)
  all.as <- t(dat[ii, all.as.cols])
  qs <- all.qs[!is.na(all.as)]
  as <- all.as[!is.na(all.as)]
  mu <- dat$`Central Value ($ per ton CO2)`[ii]
  if (is.na(mu) && length(qs) == 0) {
    next
  }
  
  dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
}

papers=unique(dat$ID_number)

nsamp=1e7
dist=matrix(nrow=nsamp,ncol=2)

for(i in 1:nsamp){
  if(i%%10000==0) print(i)
  
  #first uniform draw from papers
  paper=sample(papers,1)
  
  #draw from rows for each paper
  rows=which(dat$ID_number==paper)
  row=ifelse(length(rows)==1,rows,sample(rows,1))

  draw=sample(dists[[row]],1)
  dist[i,]=c(draw,row)
}

colnames(dist)=c("draw","row")
fwrite(dist,file="outputs/distribution.csv")

#make some figures analyzing variance in distribution
mod=numeric(length=nrow(dist))
mod[c(grep("DICE",dat$`Base IAM (if applicable)`[dist$row]),grep("DICE",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="DICE"
mod[c(grep("FUND",dat$`Base IAM (if applicable)`[dist$row]),grep("FUND",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="FUND"
mod[c(grep("PAGE",dat$`Base IAM (if applicable)`[dist$row]),grep("PAGE",dat$`IAM Calibrated To (if applicable)`[dist$row]))]="PAGE"
mod[which(mod==0)]="Other"

dist$mod=ordered(mod,levels=c("DICE","PAGE","FUND","Other"))

a=ggplot(dist[which(dist$draw>quantile(dist$draw,0.01)&dist$draw<quantile(dist$draw,0.99)),],aes(x=draw,y=mod,group=mod))+geom_density_ridges(scale=0.9)
a=a+theme_bw()+labs(x="SCC ($ per ton CO2)",y="Base or Calibration Model")
