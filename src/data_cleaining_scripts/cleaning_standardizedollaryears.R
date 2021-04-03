library(readxl)
library(lubridate)
library(tidyverse)

#code to adjust values to common dollar years

dat=as.data.frame(read_excel("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx",sheet="Data Entry"))
colnames(dat)=dat[2,];dat=dat[-c(1:2),]

#convert values in Euros to dollars
eurovals=c(grep("EU",dat$`Central Value ($ per ton CO2)`),grep("EU",dat$`Reported Base Model SCC (if applicable)`))
exchange=read.csv("data/exchangerates.csv")

#find columns that need to be transformed
distcols=which(names(dat) == 'Min'):which(names(dat) == 'Max')
relcols=c(which(colnames(dat)%in%c("Central Value ($ per ton CO2)","Reported Base Model SCC (if applicable)")),distcols)

for(j in eurovals){
  year=as.numeric(dat$Year)[j]
  rate=exchange$EurotoDol[which(exchange$Year==year)]
  dat[j,relcols]=substr(dat[j,relcols],4,5)
  dat[j,relcols]=as.numeric(dat[j,relcols])*rate
}

#convert columns to numerics
dat$`SCC Dollar Year`=as.numeric(dat$`SCC Dollar Year`)
dat$Year=as.numeric(dat$Year)
dat$`Central Value ($ per ton CO2)`=as.numeric(dat$`Central Value ($ per ton CO2)`)
dat$`Reported Base Model SCC (if applicable)`=as.numeric(dat$`Reported Base Model SCC (if applicable)`)
for(i in distcols) dat[,i]=as.numeric(dat[,i])

#run modelnames cleaningscript (cleaning_modelnames.R) to standardize model names

#for entries with missing dollar years try and infer based on calibration model where available
lookup=read.csv("data/model_dollaryear_lookup.csv") #matching of models to dollar year based on other data entered in spreadsheet

missing=which(is.na(dat$`SCC Dollar Year`))

for(i in 1:length(missing)){
  mod=ifelse(is.na(dat$`Base IAM (if applicable)`[missing[i]]),as.character(dat$`IAM Calibrated To (if applicable)`[missing[i]]),as.character(dat$`Base IAM (if applicable)`[missing[i]]))
  if(mod%in%lookup$Model){
    dat$`SCC Dollar Year`[missing[i]]=lookup$Dollar.Year[which(lookup$Model==mod)]
    next
  }
  else{
    dat$`SCC Dollar Year`[missing[i]]=dat$Year[missing[i]]-5 #if no dollar year value and no model match then use 5 years before the publication year
  }
}

#read in deflator time series from St Louis Fed: https://fred.stlouisfed.org/series/GDPDEF
#Reference year = 2012

def=read.csv("data/usgdp_deflator.csv")
def$DATE=mdy(def$DATE)
def$year=year(def$DATE)
#get annual average
def=def%>%
  group_by(year)%>%
  dplyr::summarise(deflator=mean(GDPDEF))

#adjust to 2020 index
def$deflator=def$deflator/def$deflator[which(def$year==2020)]*100

#merge into dataset and adjust all dollar values (Central SCC, Base SCC, and distribution ranges) to 2020 values
names(dat)[duplicated(names(dat))] <- paste0(names(dat)[duplicated(names(dat))], '2')
dat = dat %>% left_join(def, by=c("SCC Dollar Year"="year"))

relcols=c(which(colnames(dat)%in%c("Central Value ($ per ton CO2)","Reported Base Model SCC (if applicable)")),which(colnames(dat)=="Min"):which(colnames(dat)=="Max"))
for(i in relcols) dat[,i]=dat[,i]/(dat$deflator/100)

