dat=as.data.frame(read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Data Entry"))
colnames(dat)=dat[2,];dat=dat[-c(1:2),]

source("src/data_cleaining_scripts/cleaning_modelnames.R")
source("src/data_cleaining_scripts/cleaning_standardizedollaryears.R")
