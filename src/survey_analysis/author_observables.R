library(tidyverse)

source("src/data_cleaining_scripts/cleaning_master.R")
dat$ID_number=as.integer(dat$ID_number)

authorlist=read.csv("data/authorlist.csv")

authors=as.character(authorlist$Author1)
for(i in 3:10) authors=append(authors,as.character(authorlist[,i]))
authors=unique(authors)

outdat=matrix(nrow=length(authors),ncol=16)


for(i in 1:length(authors)){
  #find paper ids with that author listed
  relrows=authorlist%>%
    mutate(across(contains("Author"),~case_when(str_detect(.,authors[i])~1),.names = 'new_{col}'))%>%
    unite(papers, starts_with('new'), na.rm = TRUE, sep = ' ')%>%
    filter(papers==1)%>%
    select(ID_number)
  
  #find relevant rows of data-frame based on authored papers
  reldat=dat%>%
    filter(ID_number%in%as.numeric(relrows$ID_number))
  
  #average 2020 SCC and discount rate
  
  scc=reldat%>%
    filter(as.numeric(`SCC Year`)%in%c(2010:2030))%>%
    dplyr::summarize(SCC=quantile(`Central Value ($ per ton CO2)`,0.5,na.rm=T),dr=quantile(discountrate,0.5,na.rm=T))
  
  #Weighting for each structural change
  
  strucs=reldat%>%
    mutate_at(vars("Carbon Cycle":"Alternative ethical approaches (not Discounted Utilitarianism)"),~replace_na(.,0))%>%
    mutate_at(vars("Carbon Cycle":"Alternative ethical approaches (not Discounted Utilitarianism)"),function(x) replace(x, which(x=="Calibrated"),1))%>%
    mutate_at(vars("Carbon Cycle":"Alternative ethical approaches (not Discounted Utilitarianism)"),~as.numeric(.))%>%
    group_by(ID_number)%>%
    dplyr::summarise(across("Carbon Cycle":"Alternative ethical approaches (not Discounted Utilitarianism)",~max(.)))%>%
    mutate_at(vars("Carbon Cycle":"Alternative ethical approaches (not Discounted Utilitarianism)"),function(x) replace(x, which(x<0),0))
  
  #fraction empirical improvement and framework expansion
  
  type=reldat%>%
    group_by(ID_number)%>%
    dplyr::summarise(type=.data[["Empirical Improvement or Sensitvity Analysis?"]][1])
  
  npapers=reldat%>%
    mutate(npapers=length(unique(ID_number)))%>%
    select(npapers)
  
  outdat[i,]=as.numeric(c(scc,colMeans(strucs[,-1]),sum(type$type=="Framework Expansion")/dim(type)[1],sum(type$type=="Empirical Improvement")/dim(type)[1],npapers$npapers[1]))
  print(i)
  
}

outdat=as.data.frame(outdat);colnames(outdat)=c("MedianSCC_2010-2030","MedianDR_2010-2030",colnames(dat)[26:36],"Framework Expansion","Empirical Improvement","npapers")
outdat$Author=authors

write.csv(outdat,file="outputs/author_observables.csv")


