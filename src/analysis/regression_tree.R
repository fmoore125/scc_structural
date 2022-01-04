## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

dat <- multivar.prep(dat)

## emitdf <- get.emits("RCP 8.5")
## dmgfunc <- function(T) 0.0023888 * T^2
## year0 <- 2020
## gdp0 <- 84.54e12
## discountrate <- 3

source("src/analysis/damage_funcs_lib.R")

dat$log.scc.synth <- log(dat$scc.synth)

library(rpart)

mod <- rpart(log.scc.2020usd ~ `SCC Year` + `Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + scc.synth, data=dat)

mod <- rpart(log.scc.2020usd ~ `SCC Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=dat)

mod <- rpart(log.scc.2020usd ~ `SCC Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=dat)

library(rpart.plot)

rpart.plot(mod)

get.allvar <- function(framevar) {
    head <- framevar[1]
    if (framevar[2] == "<leaf>")
        allvar <- head
    else
        allvar <- c(head, get.allvar(framevar[-1]))
    if (framevar[length(allvar)+2] == "<leaf>")
        allvar <- c(allvar, head)
    else
        allvar <- c(allvar, head, get.allvar(framevar[-1:-(length(allvar)+1)]))
    allvar
}

my.rpart.plot <- function(mod) {
    split.fun <- function(x, labs, digits, varlen, faclen)
    {
        allvar <- c("root", get.allvar(x$frame$var))
        lts <- grep("<", labs)
        newlts <- labs[lts]
        newlts[allvar[lts] == "Synthetic SCC"] <- paste("<", round(exp(as.numeric(sapply(str_split(newlts[allvar[lts] == "Synthetic SCC"], " "), function(x) x[3])))))
        newlts[allvar[lts] != "Synthetic SCC"] <- paste("<", sapply(str_split(newlts[allvar[lts] != "Synthetic SCC"], " "), function(x) x[3]))
        labs[lts] <- newlts
        gts <- grep(">", labs)
        newgts <- labs[gts]
        newgts[allvar[gts] == "Synthetic SCC"] <- paste("≥", round(exp(as.numeric(sapply(str_split(newgts[allvar[gts] == "Synthetic SCC"], " "), function(x) x[3])))))
        newgts[allvar[gts] != "Synthetic SCC"] <- paste("≥", sapply(str_split(newgts[allvar[gts] != "Synthetic SCC"], " "), function(x) x[3]))
        labs[gts] <- newgts
        labs[labs == "0"] <- "No"
        labs[labs == "1.0"] <- "Yes"
        labs[labs == "1.0,Calibrated"] <- "Yes OR Calibrated"
        labs
    }

    mod$frame$yval <- exp(mod$frame$yval)
    mod$frame$var[mod$frame$var == 'log.scc.synth'] <- "Synthetic SCC"
    mod$frame$var[mod$frame$var == 'discountrate'] <- "Discount Rate (%)"

    rpart.plot(mod, type=5, split.fun=split.fun)
}

my.rpart.plot(mod)

## First, remove effect of discount rate and SCC Year

summary(lm(log.scc.2020usd ~ `SCC Year` + discountrate, data=dat[!is.na(dat$`SCC Year`) & !is.na(dat$discountrate),]))
summary(lm(log.scc.2020usd ~ poly(`SCC Year`, 2) + poly(discountrate, 2), data=dat[!is.na(dat$`SCC Year`) & !is.na(dat$discountrate),]))

factout <- lm(log.scc.2020usd ~ poly(`SCC Year`, 2) + poly(discountrate, 2), data=dat[!is.na(dat$`SCC Year`) & !is.na(dat$discountrate),])
summary(factout)

dat$resid <- NA
dat$resid[-factout$na.action] <- factout$resid

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=dat)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=dat)

my.rpart.plot(mod)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat)

my.rpart.plot(mod)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat[dat$`Persistent / Growth Damages` == "0",])

my.rpart.plot(mod)

## Try manipulating it

## Also include comparison

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")
source("src/analysis/damage_funcs_lib.R")

df <- get.all.scc(dat)
df$log.scc <- log(df$scc)
df$log.scc[!is.finite(df$log.scc)] <- NA
df$`SCC Year` <- as.numeric(df$`SCC Year`)
df$log.scc.synth <- log(df$scc.synth)

factout <- lm(log.scc ~ `SCC Year` + discountrate, data=df)
summary(factout)

df$resid <- NA
df$resid[-factout$na.action] <- factout$resid

library(rpart)

df <- multivar.prep(df)

mod <- rpart(log.scc ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=df)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=df)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Tipping Points2` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=df)

library(rpart.plot)

split.fun <- function(x, labs, digits, varlen, faclen)
{
    lts <- grep("<", labs)
    newlts <- paste("<", round(exp(as.numeric(sapply(str_split(labs[lts], " "), function(x) x[3])))))
    labs[lts] <- newlts
    gts <- grep(">", labs)
    newgts <- paste("≥", round(exp(as.numeric(sapply(str_split(labs[gts], " "), function(x) x[3])))))
    labs[gts] <- newgts
    labs[labs == "0"] <- "No"
    labs[labs == "1.0"] <- "Yes"
    labs[labs == "1.0,Calibrated"] <- "Yes OR Calibrated"
    labs
}

mod$frame$yval <- exp(mod$frame$yval)
mod$frame$var[mod$frame$var == 'log.scc.synth'] <- "log(Synth. SCC)"

rpart.plot(mod, type=5, split.fun=split.fun)
