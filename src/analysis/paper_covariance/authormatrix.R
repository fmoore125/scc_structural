library(readxl)

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