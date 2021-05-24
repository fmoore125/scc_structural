library(reshape2)

## Clean up inputted data

dat$`Constant Discount Rate (%)` <- as.numeric(dat$`Constant Discount Rate (%)`)
dat$`PRTP` <- as.numeric(dat$`PRTP`)
dat$`EMUC` <- as.numeric(dat$`EMUC`)
dat$`RRA` <- as.numeric(dat$`RRA`)
dat$`IES` <- as.numeric(dat$`IES`)

##### Construct "effective discount rate"

scenarios <- list()

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
    scenarios[[ssp]] <- data.frame(year=sspdf4.withrate$year[sspdf4.withrate$ssp == ssp], rate=sspdf4.withrate$dlog[sspdf4.withrate$ssp == ssp])

## Use these for other standard scenarios

scenarios[['BAU']] <- scenarios[['SSP5']]
scenarios[['Baseline']] <- scenarios[['SSP5']]
scenarios[['Base']] <- scenarios[['SSP5']]
scenarios[['Optimal']] <- scenarios[['SSP2']]
scenarios[['Stabilization']] <- scenarios[['SSP1']]

## Record DICE average growth rates

scenarios[['Baseline:DICE-1992']] <- log(22272 / 11293) / (2100 - 2015)
scenarios[['Baseline:DICE-2016R']] <- log(73367 / 14183) / (2100 - 2015)

