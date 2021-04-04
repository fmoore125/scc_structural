source("src/analysis/find_distribution.R")
source("src/data_cleaining_scripts/cleaning_master.R")

library(data.table)

all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

#start by generating distributions for each row
dists=list()
for (ii in 1:nrow(dat)) {
  print(ii)
  if(dat$ID_number[ii]=="1005.0") next #skip this problematic paper causing errors in generate. for the time being 
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

#for the time being, skip paper 1005 which is causing an error
papers=papers[-which(papers=="1005.0")]

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
fwrite(dist,file="C:/Users/fmoore/Box/Davis Stuff/SCC Structural Uncertainty/distribution.csv")
