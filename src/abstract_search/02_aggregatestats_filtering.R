dat=read.csv("C:/Users/fmoore/Box/Davis Stuff/SCC Structural Uncertainty/Abstract Search/compiledpapers_20200930_finalreview.csv")

#first code
firstcode=dat$Abstract.Filter
firstcode[which(is.na(firstcode))]=dat$Second.NA.Filter[which(is.na(firstcode))]
firstcode[which(is.na(firstcode))]=dat$Third.Abstract.Filter[which(is.na(firstcode))]

#second code
secondcode=dat$Second.NA.Filter
#remove values where this is actually the first code
secondcode[which(is.na(dat$Abstract.Filter))]=NA
#fill in values where the third code is the second code
secondcode[which(is.na(dat$Abstract.Filter)&is.finite(dat$Second.NA.Filter))]=dat$Third.Abstract.Filter[which(is.na(dat$Abstract.Filter)&is.finite(dat$Second.NA.Filter))]

agree=length(c(which(firstcode==0&secondcode==0),which(firstcode==1&secondcode==1)))
#fraction agree
agree/sum(is.finite(secondcode))

#third code
thirdcode=dat$Third.Abstract.Filter
thirdcode[which(is.na(dat$Abstract.Filter)|is.na(dat$Second.NA.Filter))]=NA

agree_third=length(which(secondcode==0&thirdcode==0|secondcode==1&thirdcode==1))
