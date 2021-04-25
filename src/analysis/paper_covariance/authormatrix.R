library(readxl)
library(tidyverse)

source("src/data_cleaining_scripts/cleaning_master.R")

authors=read.csv("data/authorlist.csv")
for(i in 2:dim(authors)[2]) levels(authors[,i])[levels(authors[,i])==""]=NA

sharedauthors=matrix(nrow=nrow(authors),ncol=nrow(authors))

for(i in 1:nrow(authors)){
  print(i)
  authors1=authors[i,2:10]
  authors1=authors1[!is.na(authors1)]
  for(j in 1:nrow(authors)){
    if(i==j){sharedauthors[i,j]=0;next}
    authors2=authors[j,2:10]
    authors2=authors2[!is.na(authors2)]
    #what fraction of authors from paper i are also on paper j?
    sharedauthors[i,j]=sum(authors1%in%authors2)/length(authors1)
  }
}

#calculate null value of shared authorship by reshuffling authors, keeping number of authors for each paper the same
nreps=250

distinctauthors=unique(unlist(authors[,2:10]))
distinctauthors=distinctauthors[-is.na(distinctauthors)]

nullboots=matrix(nrow=dim(authors)[1],ncol=nreps)

for(k in 1:nreps){
  print(k)
  temp=matrix(nrow=nrow(authors),ncol=nrow(authors))
  for(i in 1:nrow(authors)){
    authors1=authors[i,2:10]
    authors1=authors1[!is.na(authors1)]
    authors1_rand=sample(distinctauthors,length(authors1),replace=FALSE)
    for(j in 1:nrow(authors)){
      if(i==j){temp[i,j]=0;next}
      authors2=authors[j,2:10]
      authors2=authors2[!is.na(authors2)]
      authors2_rand=sample(distinctauthors,length(authors2),replace=FALSE)
      #what fraction of authors from paper i are also on paper j?
      temp[i,j]=sum(authors1_rand%in%authors2_rand)/length(authors1_rand)
    }
  }
  nullboots[,k]=rowSums(temp)
}

#get average null value
nullval=mean(colMeans(nullboots))
save(nullval,file="src/analysis/paper_covariance/nullvalue.Rdat")

#weighting function based on shared authorship across papers
weights=data.frame(ID_number=authors$ID_number,weight=ifelse(rowSums(sharedauthors)<nullval,1,1/(rowSums(sharedauthors)-nullval+1)))

write.csv(weights,file="src/analysis/paper_covariance/paperweightings.csv")
