## setwd("~/research/scciams/scc_structural")
## source("src/data_cleaining_scripts/cleaning_master.R")

library(reshape2)
library(tidyr)

## Clean up inputted data

dat$`Constant Discount Rate (%)` <- as.numeric(dat$`Constant Discount Rate (%)`)
dat$`PRTP` <- as.numeric(dat$`PRTP`)
dat$`EMUC` <- as.numeric(dat$`EMUC`)
dat$`RRA` <- as.numeric(dat$`RRA`)
dat$`IES` <- as.numeric(dat$`IES`)

# find entries using endogenous discounting
constant <- which(is.finite(dat$`Constant Discount Rate (%)`))

## Prepare scenario data

# read in data on consumption growth scenarios from different models etc
scens <- as.data.frame(read_excel("data/scenarios/consumption_temperature_trajectories.xlsx",sheet="Consumption Growth Rates"))
colnames(scens) <- scens[1,]
scens <- scens[-1,]
scens <- scens%>%pivot_longer(cols="2010":"2200",names_to="year",values_to="cons_growth_percap")
scens$`Short Name for Merging`=as.factor(scens$`Short Name for Merging`)
colnames(scens)[2] <- "shortname"
# interpolate values for each socio-economic scenario
scens <- scens%>%
  group_by(Model)%>%
  tidyr::fill(cons_growth_percap,.direction="downup")%>%
  ungroup()%>%
  group_by(shortname,year)%>%
  dplyr::summarise(cons_growth_percap=mean(cons_growth_percap,na.rm=T)*100)
scens$year <- as.numeric(scens$year)

scensyears <- scens%>%group_by(year)%>%summarise(cons_growth_percap=mean(cons_growth_percap,na.rm=T))

## Calculate average discount rates

## Year for discount rate merging
mround <- function(x,base){
  base*round(x/base)
}
dat$sccyearformerge <- mround(as.numeric(dat$`SCC Year`),5)
dat$sccyearformerge[which(dat$sccyearformerge<2010)] <- 2010 #match early years to most recent year available

## Version 1:

ptm <- proc.time()

# for ramsey discounting, add per-capita consumption growth rates based on socio-economic scenarios
given <- which(is.finite(as.numeric(as.character(dat$`Socio-Economic Scenario`))))

# merge in by year and model
dat <- left_join(dat,scens,by=c("sccyearformerge"="year","Socio-Economic Scenario"="shortname"), suffix=c('.old', ''), keep=FALSE)

# for scc years greater than 2200, set consumption growth to average of 2200 scenarios
dat$cons_growth_percap[which(dat$`SCC Year`>2200)] <- mean(scens$cons_growth_percap[which(scens$year==2200)],na.rm=T)

# deal with a couple weird cases manually
slowa1b <- which(dat$`Socio-Economic Scenario`=="A1B2%GDPgrowth")
fasta1b=which(dat$`Socio-Economic Scenario`=="A1B+3%GDPgrowth")
for(i in slowa1b) dat$cons_growth_percap[i] <- scens$cons_growth_percap[which(scens$shortname=="A1"&scens$year==dat$sccyearformerge[i])]-2
for(i in fasta1b) dat$cons_growth_percap[i] <- scens$cons_growth_percap[which(scens$shortname=="A1"&scens$year==dat$sccyearformerge[i])]+3

# fill in remaining with average consumption growth over all scenarios for that year
for(i in 1:dim(dat)[1]){
  if(i %in% constant) next
  if(i %in% given) next
  if(!is.na(dat$cons_growth_percap[i])) next
  dat$cons_growth_percap[i] <- scensyears$cons_growth_percap[which(scensyears$year==dat$sccyearformerge[i])]
}

# generate discount rate column
discount <- rep(NA,dim(dat)[1])
discount[constant] <- dat$`Constant Discount Rate (%)`[constant]
for(i in 1:length(discount)){
  if(!is.na(discount[i])) next
  prtp <- dat$PRTP[i]
  emuc <- ifelse(!is.na(dat$EMUC[i]),dat$EMUC[i],1/dat$IES[i])
  discount[i] <- ifelse(i%in%given,prtp+emuc*as.numeric(as.character(dat$`Socio-Economic Scenario`[i])),prtp+emuc*dat$cons_growth_percap[i])
}

dat$discountrate <- discount #just 42 missing discount rate entries, mostly due to unreported PRTP

print("Version 1:")
print(proc.time() - ptm)

## Version 2:

ptm <- proc.time()

## for scc years greater than 2200, set consumption growth to average of 2200 scenarios
dat$sccyearformerge[which(dat$sccyearformerge > 2200)] <- 2200

get.cons.growth.percap <- function(scenario) {
    if (is.finite(suppressWarnings(as.numeric(as.character(scenario)))))
        return(data.frame(year=seq(min(scens$year), max(scens$year), by=5),
                          cons_growth_percap=as.numeric(as.character(scenario))))

    if (scenario %in% scens$shortname) {
        df <- subset(scens, shortname == scenario)
        if (any(!is.nan(df$cons_growth_percap))) # not true for IAWG_BAU
            return(df)
    }

    if (!is.na(scenario)) {
        if (scenario == "A1B2%GDPgrowth") {
            df <- subset(scens, shortname == 'A1')
            df$cons_growth_percap <- df$cons_growth_percap - 2
            return(df)
        }
        if (scenario == "A1B+3%GDPgrowth") {
            df <- subset(scens, shortname == 'A1')
            df$cons_growth_percap <- df$cons_growth_percap + 3
            return(df)
        }
    }

    return(scensyears)
}

dat$discountrate2 <- NA
dat$discountrate2[constant] <- dat$`Constant Discount Rate (%)`[constant]

emuc <- ifelse(!is.na(dat$EMUC), dat$EMUC, 1/dat$IES)

for (scenario in unique(dat$`Socio-Economic Scenario`)) {
    if (is.na(scenario))
        rows <- which(is.na(dat$`Socio-Economic Scenario`) & is.na(dat$discountrate2))
    else
        rows <- which(dat$`Socio-Economic Scenario` == scenario & is.na(dat$discountrate2))

    scendf <- get.cons.growth.percap(scenario)

    scendf.torow <- data.frame(year=dat$sccyearformerge[rows]) %>% left_join(scendf, by='year')

    dat$discountrate2[rows] <- dat$PRTP[rows] + emuc[rows]*scendf.torow$cons_growth_percap
}

print("Version 2:")
print(proc.time() - ptm)

#stopifnot(sum((dat$discountrate != dat$discountrate2 & dat$`SCC Year` < 2200) | is.na(dat$discountrate) != is.na(dat$discountrate2)) == 0)

#
# ## Develop scenarios list
#
# ## May be a single average growth rate, or a data.frame of year, rate
# scenarios.best <- list()
# scenarios.second <- list()
#
# ## Load SSP data
#
# ssp.gdp <- read.csv("data/scenarios/ssps/gdp.csv")
# ssp.gdp.long <- melt(ssp.gdp, c('MODEL', 'SCENARIO', 'REGION', 'UNIT', 'VAR1', 'VAR2', 'VAR3', 'VAR4'))
#
# ssp.pop <- read.csv("data/scenarios/ssps/total_population.csv")
# ssp.pop.long <- melt(ssp.pop, c('MODEL', 'SCENARIO', 'REGION', 'UNIT', 'VAR1', 'VAR2', 'VAR3', 'VAR4'))
#
# sspdf <- ssp.gdp.long %>% left_join(ssp.pop.long, by=c("MODEL", "SCENARIO", "REGION", 'variable'), suffix=c('.gdp', '.pop'))
# sspdf$year <- as.numeric(substring(sspdf$variable, 2, 5))
#
# ## Calculate growth rates
#
# sspdf2 <- sspdf %>% group_by(MODEL, SCENARIO, REGION) %>% summarize(year, value.gdp, value.pop,
#                                                                     gdppc2010=(value.gdp / value.pop)[year == 2010],
#                                                                     dlog=c(NA, log(value.gdp / value.pop)[-1] - log(value.gdp / value.pop)[-length(value.gdp)]) / 5)
# sspdf2$ssp <- substring(sspdf2$SCENARIO, 1, 4)
#
# ## Convert to Weitzman-Gollier discounting
# ## Do this first within regions, since different models produce different regions
#
# sspdf3 <- sspdf2 %>% group_by(ssp, REGION, year) %>% summarize(mindlog=min(dlog, na.rm=T), maxdlog=max(dlog, na.rm=T),
#                                                                rate=log(mean(exp(dlog), na.rm=T)),
#                                                                pop=mean(value.pop, na.rm=T),
#                                                                gdppc2010=exp(mean(log(gdppc2010), na.rm=T)))
#
# ## Combine across all regions
#
# sspdf3.withgdppc <- sspdf3 %>% group_by(ssp, REGION) %>% summarize(growth=c(1, exp(rate[year >= 2010] * 5)),
#                                                                    gdppc=gdppc2010[1] * cumprod(growth),
#                                                                    pop=c(pop[year >= 2010], NA),
#                                                                    year=c(year[year >= 2010], max(year) + 5))
# sspdf3.withgdppc$valid <- is.finite(sspdf3.withgdppc$gdppc) & is.finite(sspdf3.withgdppc$pop)
#
# sspdf4 <- sspdf3.withgdppc %>% group_by(ssp, year) %>% summarize(gdppc=sum((gdppc * pop)[valid]) / sum(pop[valid]))
# sspdf4.withrate <- sspdf4 %>% group_by(ssp) %>% summarize(year, gdppc,
#                                                           dlog=c(NA, log(gdppc)[-1] - log(gdppc)[-length(gdppc)]) / 5)
#
# ## For informational purposes
# sspdf4.withrate %>% group_by(ssp) %>% summarize(rate=mean(dlog, na.rm=T))
#
# ## Record these
# for (ssp in unique(sspdf4.withrate$ssp))
#     scenarios.best[[ssp]] <- data.frame(year=sspdf4.withrate$year[sspdf4.withrate$ssp == ssp], rate=sspdf4.withrate$dlog[sspdf4.withrate$ssp == ssp])
#
# ## Use these for other standard scenarios
#
# scenarios.second[['BAU']] <- scenarios.best[['SSP5']]
# scenarios.second[['Baseline']] <- scenarios.best[['SSP5']]
# scenarios.second[['Base']] <- scenarios.best[['SSP5']]
# scenarios.second[['Optimal']] <- scenarios.best[['SSP2']]
# scenarios.second[['Stabilization']] <- scenarios.best[['SSP1']]
#
# ## Record DICE average growth rates
#
# scenarios.best[['Baseline:DICE1992']] <- (22272 / 11293) ** (1 / (2100 - 2015)) - 1
# scenarios.best[['Baseline:DICE2016R']] <- (73367 / 14183) ** (1 / (2100 - 2015)) - 1
#
# for (dicemodel in c('DICE1994', 'DICE1998', 'DICE1999'))
#     scenarios.second[[paste0('Baseline:', dicemodel)]] <- scenarios.best[['Baseline:DICE1992']]
# for (dicemodel in c('DICE2007', 'DICE2010', 'DICE2013', 'DICE2016R2'))
#     scenarios.second[[paste0('Baseline:', dicemodel)]] <- scenarios.best[['Baseline:DICE2016R']]
#
# ## Construct effective discount rate
#
# default.rate <- median(dat$`Constant Discount Rate (%)`, na.rm=T) / 100
#
# ## Check for inconsistency
#
# dat$discount.problems <- NA
# dat$discount.problems[is.na(dat$PRTP) & is.na(dat$`Constant Discount Rate (%)`)] <- "no discounting"
# dat$discount.problems[!is.na(dat$PRTP) & !is.na(dat$`Constant Discount Rate (%)`)] <- "PRTP and constant"
# dat$discount.problems[!is.na(dat$EMUC) & !is.na(dat$RRA)] <- "EMUC and EZ"
# dat$discount.problems[!is.na(dat$PRTP) & (is.na(dat$EMUC) & is.na(dat$RRA))] <- "no elasticity"
# dat$discount.problems[(is.na(dat$RRA) & !is.na(dat$IES)) | (!is.na(dat$RRA) & is.na(dat$IES))] <- "incomplete EZ"
#
# dat$effective.discount.rate <- dat$`Constant Discount Rate (%)` / 100
#
# for (ii in which(!is.na(dat$PRTP))) {
#     scenario <- dat$`Scenario (e.g. Optimal, BAU)`[ii]
#     scendata <- NULL
#
#     ## First-best options
#     if (scenario %in% names(scenarios.best)) {
#         scendata <- scenarios.best[[scenario]]
#     }
#     ## Model:Scenario entries
#     if (paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.best)) {
#        scendata <- scenarios.best[[paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':')]]
#     }
#     if (paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.best)) {
#         scendata <- scenarios.best[[paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':')]]
#     }
#     if (sum(sapply(names(scenarios.best), function(known) length(grep(known, scenario)))) > 0) {
#         found <- names(scenarios.best)[sapply(names(scenarios.best), function(known) length(grep(known, scenario))) > 0]
#         scendata <- scenarios.best[[found]]
#         dat$discount.problems[ii] <- "Ambiguous"
#     }
#
#     ### Second-best options
#     if (scenario %in% names(scenarios.second)) {
#         scendata <- scenarios.second[[scenario]]
#         dat$discount.problems[ii] <- "Second-best"
#     }
#     ## Model:Scenario entries
#     if (paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.second)) {
#         scendata <- scenarios.second[[paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':')]]
#     }
#     if (paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.second)) {
#         scendata <- scenarios.second[[paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':')]]
#     }
#     if (sum(sapply(names(scenarios.second), function(known) length(grep(known, scenario)))) > 0) {
#         found <- names(scenarios.second)[sapply(names(scenarios.second), function(known) length(grep(known, scenario))) > 0]
#         scendata <- scenarios.second[[found]]
#         dat$discount.problems[ii] <- "Ambiguous second-best"
#     }
#
#     if (is.null(scendata)) {
#         dat$discount.problems[ii] <- "No scenario"
#         next
#     }
#
#     sccyear <- dat$`SCC Year`[ii]
#     if (is.data.frame(scendata)) {
#         growthrate <- mean(scendata$rate[scendata$year >= sccyear])
#     } else {
#         growthrate <- scendata
#     }
#
#     if (!is.na(dat$IES[ii])) {
#         dat$effective.discount.rate[ii] <- dat$PRTP[ii] + growthrate / dat$IES[ii]
#     } else {
#         dat$effective.discount.rate[ii] <- dat$PRTP[ii] + dat$EMUC[ii] * growthrate
#     }
# }
#
#
# dat$discount.problems[is.na(dat$discount.problems) & is.na(dat$effective.discount.rate)] <- "Unknown problem"
#
# if (any(!is.na(dat$discount.problems))) {
#     print("Discounting problems:")
#     print(table(dat$discount.problems))
# }
#
# dat$effective.discount.rate[is.na(dat$effective.discount.rate)] <- default.rate
# dat$effective.discount.rate.percent <- dat$effective.discount.rate * 100
