## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

dat$`Ambiguity/Model Uncertainty`[is.na(dat$`Ambiguity/Model Uncertainty`)] <- "0"
dat$`Carbon Cycle`[is.na(dat$`Carbon Cycle`)] <- "0"
dat$`Climate Model`[is.na(dat$`Climate Model`)] <- "0"
dat$`Tipping Points`[is.na(dat$`Tipping Points`)] <- "0"
dat$`Epstein-Zin`[is.na(dat$`Epstein-Zin`)] <- "0"
dat$`Inequality Aversion`[is.na(dat$`Inequality Aversion`)] <- "0"
dat$`Learning`[is.na(dat$`Learning`)] <- "0"
dat$`Limitedly-Substitutable Goods`[is.na(dat$`Limitedly-Substitutable Goods`)] <- "0"
dat$`Persistent / Growth Damages`[is.na(dat$`Persistent / Growth Damages`)] <- "0"

emitdf <- get.emits("RCP 8.5")
dmgfunc <- function(T) 0.0023888 * T^2
year0 <- 2020
gdp0 <- 84.54e12
discountrate <- 3

source("src/analysis/damage_funcs_lib.R")

library(rpart)

mod <- rpart(log.scc.2020usd ~ `SCC Year` + `Year` + discountrate + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + scc.synth, data=dat)

library(rpart.plot)

rpart.plot(mod)

