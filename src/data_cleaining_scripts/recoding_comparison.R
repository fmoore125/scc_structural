library(readxl)
library(reshape2)
library(tidyverse)

# generate distributions for initial coding and recoding

source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

orig=dat

#read in recoded data and adjust dollar values for comparison

dat=as.data.frame(read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Recoded Data"))
colnames(dat)=dat[2,];dat=dat[-c(1:2),]

source("src/data_cleaining_scripts/cleaning_modelnames.R")
source("src/data_cleaining_scripts/cleaning_standardizedollaryears.R")

recode=dat

#generate pdfs for recoded data and compare

all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(orig) == 'Min'):which(names(orig) == 'Max')
all.as.cols_recode=which(names(recode)=="Min"):which(names(recode)=="Max")

papers=unique(recode$ID_number)
orig_fulldist=list();recode_fulldist=list()
nsamps=1e6

for(i in 1:length(papers)){
  print(i)
  rows_orig=which(orig$ID_number==papers[i])
  rowdist_orig=list()
  
  for (ii in rows_orig) {
    all.as <- t(orig[ii, all.as.cols])
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- orig$`Central Value ($ per ton CO2)`[ii]
    if (is.na(mu) && length(qs) == 0) {
      next
    }
    
    rowdist_orig[[ii]] <- generate.pdf(mu, qs, as, 1e6)
  }
  orig_fulldist[[i]]=numeric(length=nsamps)
  for(j in 1:nsamps){
    #chose a row and sample from distribution in that row
    orig_fulldist[[i]][j]=sample(rowdist_orig[[ifelse(length(rows_orig)==1,rows_orig,sample(rows_orig,1))]],1)
  }
  
  #repeat for recoded data
  
  rows_recode=which(recode$ID_number==papers[i])
  if(155%in%rows_recode) rows_recode=rows_recode[-which(rows_recode%in%c(155,156))] #more problemeatic rows
  rowdist_recode=list()
  
  for (ii in rows_recode) {
    all.as <- t(recode[ii, all.as.cols_recode])
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- recode$`Central Value ($ per ton CO2)`[ii]
    if (is.na(mu) && length(qs) == 0) {
      next
    }
    rowdist_recode[[ii]]=generate.pdf(mu,qs,as,1e6)
  }
  recode_fulldist[[i]]=numeric(length=nsamps)
  for(j in 1:nsamps){
    #chose a row and sample from distribution in that row
    recode_fulldist[[i]][j]=sample(rowdist_recode[[ifelse(length(rows_recode)==1,rows_recode,sample(rows_recode,1))]],1)
  }
  
}

#compare original and recoded distributions
x11()
par(mfrow=c(3,3),mar=c(2,2,2,2))
for(i in 1:9){
  fulldat=data.frame(orig=orig_fulldist[[i]],recode=recode_fulldist[[i]])
  fulldat=pivot_longer(fulldat,cols=c(orig,recode))
  fulldat$name=as.factor(fulldat$name)
  
  boxplot(value~name,data=fulldat,names=c("Original","Recoded"),main=paste("Paper", papers[i]))

}

x11()
par(mfrow=c(3,3),mar=c(2,2,2,2))
for(i in 10:length(papers)){
  fulldat=data.frame(orig=orig_fulldist[[i]],recode=recode_fulldist[[i]])
  fulldat=pivot_longer(fulldat,cols=c(orig,recode))
  fulldat$name=as.factor(fulldat$name)
  
  boxplot(value~name,data=fulldat,names=c("Original","Recoded"),main=paste("Paper", papers[i]))
  
}
