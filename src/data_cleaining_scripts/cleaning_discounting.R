## setwd("~/research/scciams/scc_structural")
## source("src/data_cleaining_scripts/cleaning_master.R")

library(reshape2)

## Clean up inputted data

dat$`Constant Discount Rate (%)` <- as.numeric(dat$`Constant Discount Rate (%)`)
dat$`PRTP` <- as.numeric(dat$`PRTP`)
dat$`EMUC` <- as.numeric(dat$`EMUC`)
dat$`RRA` <- as.numeric(dat$`RRA`)
dat$`IES` <- as.numeric(dat$`IES`)

## Develop scenarios list

## May be a single average growth rate, or a data.frame of year, rate
scenarios.best <- list()
scenarios.second <- list()

## Load SSP data

ssp.gdp <- read.csv("data/scenarios/ssps/gdp.csv")
ssp.gdp.long <- melt(ssp.gdp, c('MODEL', 'SCENARIO', 'REGION', 'UNIT', 'VAR1', 'VAR2', 'VAR3', 'VAR4'))

ssp.pop <- read.csv("data/scenarios/ssps/total_population.csv")
ssp.pop.long <- melt(ssp.pop, c('MODEL', 'SCENARIO', 'REGION', 'UNIT', 'VAR1', 'VAR2', 'VAR3', 'VAR4'))

sspdf <- ssp.gdp.long %>% left_join(ssp.pop.long, by=c("MODEL", "SCENARIO", "REGION", 'variable'), suffix=c('.gdp', '.pop'))
sspdf$year <- as.numeric(substring(sspdf$variable, 2, 5))

## Calculate growth rates

sspdf2 <- sspdf %>% group_by(MODEL, SCENARIO, REGION) %>% summarize(year, value.gdp, value.pop,
                                                                    gdppc2010=(value.gdp / value.pop)[year == 2010],
                                                                    dlog=c(NA, log(value.gdp / value.pop)[-1] - log(value.gdp / value.pop)[-length(value.gdp)]) / 5)
sspdf2$ssp <- substring(sspdf2$SCENARIO, 1, 4)

## Convert to Weitzman-Gollier discounting
## Do this first within regions, since different models produce different regions

sspdf3 <- sspdf2 %>% group_by(ssp, REGION, year) %>% summarize(mindlog=min(dlog, na.rm=T), maxdlog=max(dlog, na.rm=T),
                                                               rate=log(mean(exp(dlog), na.rm=T)),
                                                               pop=mean(value.pop, na.rm=T),
                                                               gdppc2010=exp(mean(log(gdppc2010), na.rm=T)))

## Combine across all regions

sspdf3.withgdppc <- sspdf3 %>% group_by(ssp, REGION) %>% summarize(growth=c(1, exp(rate[year >= 2010] * 5)),
                                                                   gdppc=gdppc2010[1] * cumprod(growth),
                                                                   pop=c(pop[year >= 2010], NA),
                                                                   year=c(year[year >= 2010], max(year) + 5))
sspdf3.withgdppc$valid <- is.finite(sspdf3.withgdppc$gdppc) & is.finite(sspdf3.withgdppc$pop)

sspdf4 <- sspdf3.withgdppc %>% group_by(ssp, year) %>% summarize(gdppc=sum((gdppc * pop)[valid]) / sum(pop[valid]))
sspdf4.withrate <- sspdf4 %>% group_by(ssp) %>% summarize(year, gdppc,
                                                          dlog=c(NA, log(gdppc)[-1] - log(gdppc)[-length(gdppc)]) / 5)

## For informational purposes
sspdf4.withrate %>% group_by(ssp) %>% summarize(rate=mean(dlog, na.rm=T))

## Record these
for (ssp in unique(sspdf4.withrate$ssp))
    scenarios.best[[ssp]] <- data.frame(year=sspdf4.withrate$year[sspdf4.withrate$ssp == ssp], rate=sspdf4.withrate$dlog[sspdf4.withrate$ssp == ssp])

## Use these for other standard scenarios

scenarios.second[['BAU']] <- scenarios.best[['SSP5']]
scenarios.second[['Baseline']] <- scenarios.best[['SSP5']]
scenarios.second[['Base']] <- scenarios.best[['SSP5']]
scenarios.second[['Optimal']] <- scenarios.best[['SSP2']]
scenarios.second[['Stabilization']] <- scenarios.best[['SSP1']]

## Record DICE average growth rates

scenarios.best[['Baseline:DICE1992']] <- (22272 / 11293) ** (1 / (2100 - 2015)) - 1
scenarios.best[['Baseline:DICE2016R']] <- (73367 / 14183) ** (1 / (2100 - 2015)) - 1

for (dicemodel in c('DICE1994', 'DICE1998', 'DICE1999'))
    scenarios.second[[paste0('Baseline:', dicemodel)]] <- scenarios.best[['Baseline:DICE1992']]
for (dicemodel in c('DICE2007', 'DICE2010', 'DICE2013', 'DICE2016R2'))
    scenarios.second[[paste0('Baseline:', dicemodel)]] <- scenarios.best[['Baseline:DICE2016R']]

## Construct effective discount rate

default.rate <- median(dat$`Constant Discount Rate (%)`, na.rm=T) / 100

## Check for inconsistency

dat$discount.problems <- NA
dat$discount.problems[is.na(dat$PRTP) & is.na(dat$`Constant Discount Rate (%)`)] <- "no discounting"
dat$discount.problems[!is.na(dat$PRTP) & !is.na(dat$`Constant Discount Rate (%)`)] <- "PRTP and constant"
dat$discount.problems[!is.na(dat$EMUC) & !is.na(dat$RRA)] <- "EMUC and EZ"
dat$discount.problems[!is.na(dat$PRTP) & (is.na(dat$EMUC) & is.na(dat$RRA))] <- "no elasticity"
dat$discount.problems[(is.na(dat$RRA) & !is.na(dat$IES)) | (!is.na(dat$RRA) & is.na(dat$IES))] <- "incomplete EZ"

dat$effective.discount.rate <- dat$`Constant Discount Rate (%)` / 100

for (ii in which(!is.na(dat$PRTP))) {
    scenario <- dat$`Scenario (e.g. Optimal, BAU)`[ii]
    scendata <- NULL

    ## First-best options
    if (scenario %in% names(scenarios.best)) {
        scendata <- scenarios.best[[scenario]]
    }
    ## Model:Scenario entries
    if (paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.best)) {
       scendata <- scenarios.best[[paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':')]]
    }
    if (paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.best)) {
        scendata <- scenarios.best[[paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':')]]
    }
    if (sum(sapply(names(scenarios.best), function(known) length(grep(known, scenario)))) > 0) {
        found <- names(scenarios.best)[sapply(names(scenarios.best), function(known) length(grep(known, scenario))) > 0]
        scendata <- scenarios.best[[found]]
        dat$discount.problems[ii] <- "Ambiguous"
    }

    ### Second-best options
    if (scenario %in% names(scenarios.second)) {
        scendata <- scenarios.second[[scenario]]
        dat$discount.problems[ii] <- "Second-best"
    }
    ## Model:Scenario entries
    if (paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.second)) {
        scendata <- scenarios.second[[paste(dat$`Base IAM (if applicable)`[ii], scenario, sep=':')]]
    }
    if (paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':') %in% names(scenarios.second)) {
        scendata <- scenarios.second[[paste(dat$`IAM Calibrated To (if applicable)`[ii], scenario, sep=':')]]
    }
    if (sum(sapply(names(scenarios.second), function(known) length(grep(known, scenario)))) > 0) {
        found <- names(scenarios.second)[sapply(names(scenarios.second), function(known) length(grep(known, scenario))) > 0]
        scendata <- scenarios.second[[found]]
        dat$discount.problems[ii] <- "Ambiguous second-best"
    }

    if (is.null(scendata)) {
        dat$discount.problems[ii] <- "No scenario"
        next
    }
    
    sccyear <- dat$`SCC Year`[ii]
    if (is.data.frame(scendata)) {
        growthrate <- mean(scendata$rate[scendata$year >= sccyear])
    } else {
        growthrate <- scendata
    }
    
    if (!is.na(dat$IES[ii])) {
        dat$effective.discount.rate[ii] <- dat$PRTP[ii] + growthrate / dat$IES[ii]
    } else {
        dat$effective.discount.rate[ii] <- dat$PRTP[ii] + dat$EMUC[ii] * growthrate
    }
}


dat$discount.problems[is.na(dat$discount.problems) & is.na(dat$effective.discount.rate)] <- "Unknown problem"

if (any(!is.na(dat$discount.problems))) {
    print("Discounting problems:")
    print(table(dat$discount.problems))
}

dat$effective.discount.rate[is.na(dat$effective.discount.rate)] <- default.rate
dat$effective.discount.rate.percent <- dat$effective.discount.rate * 100
