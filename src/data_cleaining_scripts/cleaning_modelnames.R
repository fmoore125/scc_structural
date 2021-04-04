library(readxl)
library(forcats)

#code to standardize model names

#dat=read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Data Entry")
#colnames(dat)=dat[2,];dat=dat[-c(1:2),]

#remove spaces and hyphens from model names
dat$`Base IAM (if applicable)`=str_replace_all(dat$`Base IAM (if applicable)`, fixed(" "), "")
dat$`Base IAM (if applicable)`=str_replace_all(dat$`Base IAM (if applicable)`, fixed("-"), "")
dat$`IAM Calibrated To (if applicable)`=str_replace_all(dat$`IAM Calibrated To (if applicable)`, fixed(" "), "")
dat$`IAM Calibrated To (if applicable)`=str_replace_all(dat$`IAM Calibrated To (if applicable)`, fixed("-"), "")

dat$`Base IAM (if applicable)`=as.factor(dat$`Base IAM (if applicable)`)
dat$`IAM Calibrated To (if applicable)`=as.factor(dat$`IAM Calibrated To (if applicable)`)

#manually fix remainders
dat$`Base IAM (if applicable)`=fct_collapse(dat$`Base IAM (if applicable)`,FUNDIAWG=c("FUNDIAWG","IAWGFUND"),PAGE2009=c("PAGE2009","PAGE09"))
dat$`Base IAM (if applicable)`=fct_recode(dat$`Base IAM (if applicable)`,DICE1994="DICE94",DICE1998="DICE98",DICE1999="DICE99") 
levels(dat$`Base IAM (if applicable)`)[levels(dat$`Base IAM (if applicable)`)=='n/a'] <- NA
levels(dat$`IAM Calibrated To (if applicable)`)[levels(dat$`IAM Calibrated To (if applicable)`)=='n/a'] <- NA

       