setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

if (.Platform$OS.type != "unix")
    pdffunc <- grDevices::cairo_pdf
else
    pdffunc <- pdf


dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

## emitdf <- get.emits("RCP 8.5")
## dmgfunc <- function(T) 0.0023888 * T^2
## year0 <- 2020
## gdp0 <- 84.54e12
## discountrate <- 3

source("src/analysis/damage_funcs_lib.R")

## Redo log.scc.synth so I get NAs
dat$log.scc.synth <- log(dat$scc.synth)

datall <- get.all.scc(dat)

dat <- multivar.prep(dat)
datall <- multivar.prep(datall)

library(rpart)
library(splines)
library(lfe)

dat$`Earth System` <- "0"
dat$`Earth System`[paste(dat$`Carbon Cycle`, dat$`Climate Model`) != "0 0"] <- "1.0"

datall$`Earth System` <- "0"
datall$`Earth System`[paste(datall$`Carbon Cycle`, datall$`Climate Model`) != "0 0"] <- "1.0"

struccols <- sapply(c("Backstop Price?", "Other Market Failure?", "Market Only Damages", "Ambiguity/Model Uncertainty", "Earth System", "Tipping Points", "Tipping Points2", "Epstein-Zin", "Inequality Aversion", "Learning", "Limitedly-Substitutable Goods", "Persistent / Growth Damages", "Alternative ethical approaches"), function(col) which(names(dat) == col))

weights <- get.struct.weights(dat, struccols)

struccolsall <- sapply(c("Backstop Price?", "Other Market Failure?", "Market Only Damages", "Ambiguity/Model Uncertainty", "Earth System", "Tipping Points", "Tipping Points2", "Epstein-Zin", "Inequality Aversion", "Learning", "Limitedly-Substitutable Goods", "Persistent / Growth Damages", "Alternative ethical approaches"), function(col) which(names(datall) == col))

weightsall <- get.struct.weights(datall, struccolsall)

names(dat)[names(dat) == 'Tipping Points'] <- "Climate Tipping Point"
names(dat)[names(dat) == 'Tipping Points2'] <- "Damages Tipping Point"
names(datall)[names(datall) == 'Tipping Points'] <- "Climate Tipping Point"
names(datall)[names(datall) == 'Tipping Points2'] <- "Damages Tipping Point"

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

for (period in c('2010-2030', '2030-2070', '2070-2100', '2010-2100')) {
    if (period == '2010-2030')
        in.period <- dat$`SCC Year` > 2010 & dat$`SCC Year` <= 2030
    else if (period == '2030-2070')
        in.period <- dat$`SCC Year` > 2030 & dat$`SCC Year` <= 2070
    else if (period == '2070-2100')
        in.period <- dat$`SCC Year` > 2070 & dat$`SCC Year` <= 2100
    else if (period == '2010-2100')
        in.period <- dat$`SCC Year` > 2010 & dat$`SCC Year` <= 2100
    else {
        in.period <- F
        print("Unknown period")
    }

    ## mod <- rpart(log.scc.2020usd ~ `SCC Year` + `Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + scc.synth, data=dat[in.period,], weights=weights[in.period])

    mod <- rpart(log.scc.2020usd ~ `Backstop Price?` + `Other Market Failure?` + discountrate + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat[in.period,], weights=weights[in.period])

    library(rpart.plot)

    ## rpart.plot(mod)

    pdffunc(paste0("outputs/figures/levels-", period, ".pdf"), width=6.5, height=4)
    my.rpart.plot(mod)
    dev.off()

    ## First, remove effect of discount rate and SCC Year

    ## summary(lm(log.scc.2020usd ~ `SCC Year` + discountrate, data=dat[in.period & !is.na(dat$`SCC Year`) & !is.na(dat$discountrate),]))
    ## summary(lm(log.scc.2020usd ~ poly(`SCC Year`, 2) + poly(discountrate, 2), data=dat[in.period & !is.na(dat$`SCC Year`) & !is.na(dat$discountrate),]))

    ## factout <- lm(log.scc.2020usd ~ poly(`SCC Year`, 2) + poly(discountrate, 2), data=dat[in.period & !is.na(dat$`SCC Year`) & !is.na(dat$discountrate),])
    ## summary(factout)

    ## top 3 in sort(table(dat$discountrate)[table(dat$discountrate) > 10])
    ##
    splines <- as.data.frame(cbind(ns(dat$discountrate, knots=c(1.5, 3, 5)),
                                   ns(dat$`SCC Year`, knots=c(2010, 2050, 2100))))
    names(splines) <- c('dr1', 'dr2', 'dr3', 'dr4', 'sy1', 'sy2', 'sy3', 'sy4')
    dat2 <- cbind(dat, splines)
    factout <- lm(log.scc.2020usd ~ dr1 + dr2 + dr3 + dr4 + sy1 + sy2 + sy3 + sy4, data=dat2[in.period & !is.na(dat2$`SCC Year`) & !is.na(dat2$discountrate),])

    subdat <- dat[in.period & !is.na(dat$`SCC Year`) & !is.na(dat$discountrate),]
    subdat$resid <- NA
    subdat$resid[-factout$na.action] <- factout$resid
    subwgt <- weights[in.period & !is.na(dat$`SCC Year`) & !is.na(dat$discountrate)]

    mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=subdat, weights=subwgt)

    pdffunc(paste0("outputs/figures/ratios-", period, ".pdf"), width=6.5, height=4)
    my.rpart.plot(mod)
    dev.off()

    ## Try with base SCCs

    factout <- felm(log.scc.2020usd ~ 1 | basecode, data=datall[in.period,])

    subdat <- datall[in.period,]
    subdat$resid <- NA
    subdat$resid[-factout$na.action] <- factout$resid
    subwgt <- weightsall[in.period]

    mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=subdat, weights=subwgt)

    pdffunc(paste0("outputs/figures/ratios-all-", period, ".pdf"), width=6.5, height=4)
    my.rpart.plot(mod)
    dev.off()
}

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat)

my.rpart.plot(mod)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=dat[dat$`Persistent / Growth Damages` == "0",])

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

mod <- rpart(log.scc ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=df)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth, data=df)

mod <- rpart(resid ~ `Backstop Price?` + `Other Market Failure?` + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches`, data=df)

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
