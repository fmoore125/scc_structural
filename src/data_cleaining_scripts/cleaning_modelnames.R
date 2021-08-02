library(readxl)
library(forcats)
library(stringr)

#code to standardize model names

#dat=read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Data Entry")
#colnames(dat)=dat[2,];dat=dat[-c(1:2),]

#remove spaces and hyphens from model names
dat$`Base IAM (if applicable)`=str_replace_all(dat$`Base IAM (if applicable)`, fixed(" "), "")
dat$`Base IAM (if applicable)`=str_replace_all(dat$`Base IAM (if applicable)`, fixed("-"), "")
dat$`IAM Calibrated To (if applicable)`=str_replace_all(dat$`IAM Calibrated To (if applicable)`, fixed(" "), "")
dat$`IAM Calibrated To (if applicable)`=str_replace_all(dat$`IAM Calibrated To (if applicable)`, fixed("-"), "")
dat$`Socio-Economic Scenario`=str_replace_all(dat$`Socio-Economic Scenario`, fixed(" "), "")
dat$`Socio-Economic Scenario`=str_replace_all(dat$`Socio-Economic Scenario`, fixed("-"), "")


dat$`Base IAM (if applicable)`=as.factor(dat$`Base IAM (if applicable)`)
dat$`IAM Calibrated To (if applicable)`=as.factor(dat$`IAM Calibrated To (if applicable)`)
dat$`Socio-Economic Scenario`=as.factor(dat$`Socio-Economic Scenario`)

#manually fix remainders
dat$`Base IAM (if applicable)`=fct_collapse(dat$`Base IAM (if applicable)`,FUNDIAWG=c("FUNDIAWG","IAWGFUND"),PAGE2009=c("PAGE2009","PAGE09"))
dat$`Base IAM (if applicable)`=fct_recode(dat$`Base IAM (if applicable)`,DICE1994="DICE94",DICE1998="DICE98",DICE1999="DICE99")
dat$`IAM Calibrated To (if applicable)`=fct_collapse(dat$`IAM Calibrated To (if applicable)`,PAGE2009=c("PAGE09"))
dat$`IAM Calibrated To (if applicable)`=fct_recode(dat$`IAM Calibrated To (if applicable)`,DICE1994="DICE94")

#concatenate some socio-economic scenarios
dat$`Socio-Economic Scenario`=fct_collapse(dat$`Socio-Economic Scenario`,A1=c("A1","SRESA1","A1FIoutput;WorldBankpopulation","A1B","SRESA1b","SRESA1B"),A2=c("A2","SRESA2"),B1=c("B1","SRESB1"),B2=c("B2","SRESB2"),DICE2016R=c("DICE2016R","DICE2016R2,extended"),SSPAverage=c("AcrossSSPs","SSPAverage"),IAWG_BAU=c("IWG2010BAU","IWG2015BAU","Baseline(AverageoverIAWGScenarios)"),PAGE2009=c("PAGE2009","PAGE09"))
dat$`Socio-Economic Scenario`=fct_recode(dat$`Socio-Economic Scenario`,DICE1994="DICE94",DICE1999="DICE99")

levels(dat$`Base IAM (if applicable)`)[levels(dat$`Base IAM (if applicable)`)=='n/a'] <- NA
levels(dat$`IAM Calibrated To (if applicable)`)[levels(dat$`IAM Calibrated To (if applicable)`)=='n/a'] <- NA

clean.modelnames <- function(names) {
    oldw <- getOption("warn")
    options(warn=-1)

    names <- str_replace_all(names, fixed(" "), "")
    names <- str_replace_all(names, fixed("-"), "")
    names <- as.factor(names)
    names <- fct_collapse(names, FUNDIAWG=c("FUNDIAWG","IAWGFUND"), PAGE2009=c("PAGE2009","PAGE09"))
    names <- fct_recode(names, DICE1994="DICE94", DICE1998="DICE98", DICE1999="DICE99")
    levels(names)[levels(names)=='n/a'] <- NA

    options(warn=oldw)

    names
}
