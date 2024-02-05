## setwd("~/research/scciams/scc_structural")

source("src/analysis/randomforest_dists_load.R")
library(ggplot2)

load("outputs/rfdistsmodel-final-precompute.RData")

## Construct synthetic SCC

predict.forest.all <- function(expdat, incrows, outfile) {
    pdf <- data.frame()
    for (ii in sample(which(incrows))) {
        if (ii %in% pdf$row)
            next
        quants.pred <- predict.forest(forest, expdat[ii,], NULL)
        pdf <- rbind(pdf, data.frame(row=ii, quant=all.qs, pred=quants.pred))
    }

    ## Construct Monte Carlo draws
    allsamp <- c()
    for (row in unique(pdf$row)) {
        inv.cdf <- approx(pdf$quant[pdf$row == row], pdf$pred[pdf$row == row])
        uu <- runif(1000)
        samples <- inv.cdf$y[findInterval(uu, inv.cdf$x)]
        allsamp <- c(allsamp, samples)
    }

    print(quantile(allsamp, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.)))
    print(mean(allsamp))

    save(pdf, allsamp, file=outfile)
}

idealdat <- dat # include all rows

## set obvious ones
for (cc in which(names(dat) == 'TFP Growth'):which(names(dat) == 'Risk Aversion (EZ Utility)'))
    idealdat[, cc] <- "1"
idealdat$`Backstop Price?` <- "0"
idealdat$`Declining Discounting?` <- "1"
idealdat$`Market Only Damages` <- "0"
idealdat$`Other Market Failure?` <- "0"
## for damage-function-based scc - draw from literature values
idealdat$log.scc.synth[idealdat$missing.scc.synth] <- sample(idealdat$log.scc.synth[!idealdat$missing.scc.synth],
                                                            sum(idealdat$missing.scc.synth),replace=TRUE)
idealdat$Year <- 2020
idealdat$sccyearformerge <- 2020

bayespost <- read.csv("data/expert_survey/meta-analysis-distribution.csv")
bayespost2 <- dcast(bayespost, iterations ~ question, value.var='prob')
bayespost3 <- bayespost2[sample(nrow(bayespost2), nrow(dat), replace=T),]
bayespost3[, -1] <- ifelse(matrix(runif((ncol(bayespost3)-1)*nrow(bayespost3)), nrow(bayespost3), ncol(bayespost3)-1) < bayespost3[, -1], "1", "0")

idealdat <- cbind(idealdat[-c(which(names(dat) == 'Carbon Cycle'):which(names(dat) == 'Learning'), which(names(dat) == 'Earth system'))], bayespost3)
names(idealdat)[names(idealdat) == 'Tipping Points: Climate'] <- 'Tipping Points'
names(idealdat)[names(idealdat) == 'Tipping Points: Damages'] <- 'Tipping Points2'
names(idealdat)[names(idealdat) == 'Limited Substitutability'] <- 'Limitedly-Substitutable Goods'
names(idealdat)[names(idealdat) == 'Earth System'] <- 'Earth system'

discountsurvey <- read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
idealdat$discountrate <- sample(discountsurvey$SDR[!is.na(discountsurvey$SDR)], nrow(idealdat), replace=TRUE)

predict.forest.all(idealdat, incrows, "outputs/rf_experiments/RFD_best.RData")

##----A. No structural changes (classic DICE assumptions)-----

## set obvious ones
dicedat <- idealdat
for (cc in which(names(dicedat) == 'TFP Growth'):which(names(dicedat) == 'Risk Aversion (EZ Utility)'))
    dicedat[, cc] <- "0"
dicedat$`Backstop Price?` <- "0"
dicedat$`Declining Discounting?` <- "0"
dicedat$`Market Only Damages` <- "0"
dicedat$`Other Market Failure?` <- "0"

## for damage-function-based scc - draw from literature values
rel <- dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("DICE-2007","DICE-2013R","DICE-2016R2","DICE 2007","DICE 2010","DICE 2013","DICE 2013R","DICE 2016","DICE2007","DICE2010","DICE2013","DICE2016R")
dicedat$log.scc.synth <- sample(dat$log.scc.synth[rel & !dat$missing.scc.synth],
                                nrow(dat),replace=TRUE)

for (cc in which(names(dicedat) == 'Ambiguity/Model Uncertainty'):which(names(dicedat) == 'Tipping Points2'))
    dicedat[, cc] <- "0"

dicedat$discountrate <- 4.6

predict.forest.all(dicedat, incrows, "outputs/rf_experiments/RFD_A_dice.RData")

##----B. EPA assumptions -----
epadat <- dicedat

## structural changes to Earth System
epadat$`Earth system` <- "1"
## parametric uncertainty in tfp growth, pop growth, earth system and damage functions
epadat$`TFP Growth` <- "1"
epadat$`Population Growth` <- "1"
epadat$`Emissions Growth` <- "1"
epadat$`Transient Climate Response` <- "1"
epadat$`Carbon Cycle2` <- "1"
epadat$`Equilibrium Climate Sensitivity` <- "1"
epadat$`Damage Function` <- "1"
## central discount rate of 2%
epadat$discountrate <- 2
## damage function from Howard and Sterner
rel <- dat$`Damage Function Info: Model, Commonly-Used Function, or Function`%in%c("HowardSterner","HowardSterner (0.007438*T^2)")
epadat$log.scc.synth <- sample(dat$log.scc.synth[rel & !dat$missing.scc.synth],
                               nrow(dat),replace=TRUE)

predict.forest.all(epadat, incrows, "outputs/rf_experiments/RFD_B_epa.RData")

##---- C. All structural changes and no structural changes
alldat <- idealdat
for (cc in which(names(alldat) == 'Ambiguity/Model Uncertainty'):which(names(alldat) == 'Tipping Points2'))
    alldat[, cc] <- "1"

predict.forest.all(alldat, incrows, "outputs/rf_experiments/RFD_C_all.RData")

##---- E. Multiple constant discount rates (1, 1.5, 2.5, 3, 5)
discs <- c(1,1.5,2,2.5,3,5)

for (disc in discs) {
    print(disc)
    discdat <- idealdat
    discdat$discountrate <- disc
    predict.forest.all(discdat, incrows,
                       paste0("outputs/rf_experiments/RFD_E_", disc, ".RData"))
}

##---- F. Multiple publication years (2000, 2010, 2020)
pubyears <- c(2000, 2010, 2020)

for (pubyear in pubyears) {
    print(pubyear)
    pubyeardat <- idealdat
    pubyeardat$Year <- pubyear
    predict.forest.all(pubyeardat, incrows, paste0("outputs/rf_experiments/RFD_F_", pubyear, ".RData"))
}

##---- G. Multiple SCC pulse years (2020, 2050)
sccyears <- c(2020, 2050, 2100)

for (sccyear in sccyears) {
    print(sccyear)
    sccyeardat <- idealdat
    sccyeardat$sccyearformerge <- sccyear
    predict.forest.all(sccyeardat, incrows, paste0("outputs/rf_experiments/RFD_G_", sccyear, ".RData"))
}

##---- H. Multiple synthetic SCCs (DICE, FUND, PAGE, Howard & Sterner)

damages <- c(11.72858,6.337801,14.02219,38.60881)
damagemods <- c("DICE2016r2","FUND38","PAGE2009","HowardSterner")

#From James: DICE 2016r2: 11.72858
#FUND 3.8: 6.337801
#PAGE 2009: 14.02219
#Howard & Sterner: 38.60881

for (ii in 1:length(damages)) {
    print(damagemods[ii])
    dmgdat <- idealdat
    dmgdat$log.scc.synth <- log(damages[ii])
    predict.forest.all(dmgdat, incrows, paste0("outputs/rf_experiments/RFD_H_", damagemods[ii], ".RData"))
}

## I. Build-up from DICE

structural.uncertainty <- list("Ambiguity/Model Uncertainty"=c(), "Earth system"=c("Carbon Cycle2"),
                               "Epstein-Zin"=c("Risk Aversion (EZ Utility)"),
                               "Inequality Aversion"=c(), "Learning"=c(), "Limitedly-Substitutable Goods"=c(),
                               "Persistent / Growth Damages"=c(), "Tipping Points"=c("Tipping Point Magnitude"),
                               "Tipping Points2"=c("Tipping Point Magnitude"))

## I1. DICE + damages
i1dat <- dicedat
i1dat$log.scc.synth <- idealdat$log.scc.synth
predict.forest.all(i1dat, incrows, "outputs/rf_experiments/RFD_I1.RData")

## I2. I1 + discounting
i2dat <- i1dat
i2dat$discountrate <- idealdat$discountrate
i2dat$`Declining Discounting?` <- "1"
predict.forest.all(i2dat, incrows, "outputs/rf_experiments/RFD_I2.RData")

## I3. I2 + structural
i3dat <- i2dat
for (cc in which(names(dicedat) == 'Ambiguity/Model Uncertainty'):which(names(dicedat) == 'Tipping Points2')) {
    i3dat[, cc] <- idealdat[, cc]
    predict.forest.all(i3dat, incrows, paste0("outputs/rf_experiments/RFD_I3_", cc, ".RData"))
}

## I4. I3 + uncertainty
i4dat <- i3dat
for (cc in which(names(dicedat) == 'TFP Growth'):which(names(dicedat) == 'Risk Aversion (EZ Utility)')) {
    i4dat[, cc] <- "1"
    predict.forest.all(i4dat, incrows, paste0("outputs/rf_experiments/RFD_I4_", cc, ".RData"))
}

## I5. I4 - damages
i5dat <- i4dat
i5dat$log.scc.synth <- dicedat$log.scc.synth
predict.forest.all(i5dat, incrows, "outputs/rf_experiments/RFD_I5.RData")

## I6. I5 - discounting
i6dat <- i5dat
i6dat$discountrate <- 4.6
i6dat$`Declining Discounting?` <- "0"
predict.forest.all(i6dat, incrows, "outputs/rf_experiments/RFD_I6.RData")

## I7. I6 - structural
removed.uncertainties <- c()
i7dat <- i6dat
for (cc in which(names(dicedat) == 'Ambiguity/Model Uncertainty'):which(names(dicedat) == 'Tipping Points2')) {
    i7dat[, cc] <- "0"
    for (name in structural.uncertainty[[names(dicedat)[cc]]]) {
        i7dat[, name] <- "0"
        removed.uncertainties <- c(removed.uncertainties, name)
    }
    predict.forest.all(i7dat, incrows, paste0("outputs/rf_experiments/RFD_I7_", cc, ".RData"))
}

## I8. I7 - uncertainty
i8dat <- i7dat
for (cc in which(names(dicedat) == 'TFP Growth'):which(names(dicedat) == 'Risk Aversion (EZ Utility)')) {
    if (names(dicedat)[cc] %in% removed.uncertainties)
        next # don't save it
    i8dat[, cc] <- "0"
    predict.forest.all(i8dat, incrows, paste0("outputs/rf_experiments/RFD_I8_", cc, ".RData"))
}
