library(tidyverse)
library(compare)

dat=as.data.frame(read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Data Entry"))
colnames(dat)=dat[2,];dat=dat[-c(1:2),]
names(dat)[duplicated(names(dat))] <- paste0(names(dat)[duplicated(names(dat))], '2')
#check for rows that are identical except for SCC values

#keep bibliographic, strucutral and parametric information and find unique values
uniquevals=dat%>%
  dplyr::select('ID_number','Base IAM (if applicable)','IAM Calibrated To (if applicable)','SCC Year','Scenario (e.g. Optimal, BAU)','Constant Discount Rate':'Alternative ethical approaches (not Discounted Utilitarianism)','TFP Growth':'PRTP2')

uniquevals_dist=distinct(uniquevals)
uniquevals_dist$rowid=1:nrow(uniquevals_dist)

test=merge(uniquevals,uniquevals_dist,all=T,sort=FALSE)

probs=names(which(table(test$rowid)>1))

duplicatedrows=dat[which(test$rowid%in%probs),which(colnames(dat)%in%c("ID_number","Added By"))]
duplicatedrows=distinct(duplicatedrows)
write.csv(duplicatedrows,"outputs/duplicated.csv")

#second check - ranges reported with no parametric uncertainty?

probs=numeric(length=nrow(dat))
for(i in 1:length(probs)){
  paramuncertainty=dat[i,which(colnames(dat)=="Min"):which(colnames(dat)=="Max")]
  if(sum(is.na(paramuncertainty))==length(paramuncertainty)) next #no parametric uncertainty
  reportedvar=dat[i,which(colnames(dat)=="TFP Growth"):which(colnames(dat)=="Risk Aversion (EZ Utility)")]
  if(sum(is.na(reportedvar))==length(reportedvar)) probs[i]=1
}

nouncertainty=distinct(dat[which(probs==1),which(colnames(dat)%in%c("ID_number","Added By"))])
