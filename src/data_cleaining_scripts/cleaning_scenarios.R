library(readxl)

## Clean up emissions names

old <- dat$`Emissions Scenario`
dat$`Emissions Scenario` <- gsub("^RCP\\s*([0-9.]+)$", "RCP \\1", old)

## Functions to translate scenarios to temperatures

temptrajdf <- read_excel("data/scenarios/consumption_temperature_trajectories.xlsx",sheet="Temperature Trajectories", skip=1)
temptrajdf$Model <- clean.modelnames(temptrajdf$Model)

## Interpret emissions
get.emits <- function(emitscen) {
    if (emitscen %in% c('BAU', 'RCP 8.5')) {
        df <- read_excel("data/scenarios/rcps/RCP85_EMISSIONS.xls",sheet="RCP85_EMISSIONS", skip=37)
        return(data.frame(year=df$`v YEARS/GAS >`, gtcemit=df$FossilCO2 + df$OtherCO2))
    }
    if (emitscen %in% c('RCP 6.0')) {
        df <- read_excel("data/scenarios/rcps/RCP6_EMISSIONS.xls",sheet="RCP6_EMISSIONS", skip=37)
        return(data.frame(year=df$`v YEARS/GAS >`, gtcemit=df$FossilCO2 + df$OtherCO2))
    }
    if (emitscen %in% c('RCP 4.5')) {
        df <- read_excel("data/scenarios/rcps/RCP45_EMISSIONS.xls",sheet="RCP45_EMISSIONS", skip=37)
        return(data.frame(year=df$`v YEARS/GAS >`, gtcemit=df$FossilCO2 + df$OtherCO2))
    }
    if (emitscen %in% c('RCP 2.6')) {
        df <- read_excel("data/scenarios/rcps/RCP3PD_EMISSIONS.xls",sheet="RCP3PD_EMISSIONS", skip=37)
        return(data.frame(year=df$`v YEARS/GAS >`, gtcemit=df$FossilCO2 + df$OtherCO2))
    }

    NULL
}

## get.emits("RCP 8.5")
## get.emits("Wierdo")

## Model emissions to temperatures
dice.2013 <- function(emit.gtc.peryr.by5, lastyear=2305) {
    ## Parameters
    MAT <- 830.4
    MU <- 1527
    ML <- 10010
    TATM <- 0.8
    TOCEAN <- 0.0068

    b12 <- 0.088
    b23 <- 0.0025
    b11 <- 0.912
    b21 <- 0.038328889
    b22 <- 0.959171111
    b32 <- 0.0003375
    b33 <- 0.9996625

    t2xco2 <- 2.9
    c1 <- 0.098
    c3 <- 0.088
    c4 <- 0.025

    forcoth <- c(0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.025, 1.05, 1.075, 1.1, 1.125, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1, 1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675, 1.7, 1.725)
    fco22x <- 3.8
    eqmat <- 588

    ## Modeling
    results <- data.frame(year=2010, dtemp.1900=TATM)
    for (year in seq(2015, lastyear, by=5)) {
        tt <- (year - 2010) / 5 + 1
        MAT_next <- MAT * b11 + MU * b21 + (emit.gtc.peryr.by5[tt - 1] * (5 / 3.666))
        MU_next <- MAT * b12 + MU * b22 + ML * b32
        ML_next <- ML * b33 + MU * b23

        FORC <- fco22x * (log((MAT_next / eqmat)) / log(2)) + forcoth[tt]

        TATM_next <- TATM + c1 * ((FORC - (fco22x / t2xco2) * TATM) - (c3 * (TATM - TOCEAN)))
        TOCEAN_next <- TOCEAN + c4 * (TATM - TOCEAN)

        MAT <- MAT_next
        MU <- MU_next
        ML <- ML_next
        TATM <- TATM_next
        TOCEAN <- TOCEAN_next

        results <- rbind(results, data.frame(year, dtemp.1900=TATM))
    }

    results
}

## emitdf <- get.emits("RCP 8.5")
## dice.2013(approx(emitdf$year, emitdf$gtcemit, seq(2010, 2300, by=5))$y)

## Get direct temperatures or model emissions
get.temps <- function(emitscen, lastyear=2300) {
    emitscen.asmodel <- as.character(clean.modelnames(emitscen))
    if (emitscen.asmodel %in% temptrajdf$Model) {
        row <- which(temptrajdf$Model == emitscen.asmodel)[1]
        if (any(!is.na(temptrajdf[row, -1]))) {
            temps <- temptrajdf[row, -1]
            valid <- !is.na(temps)
            return(data.frame(year=names(temps)[valid], dtemp.1900=temps[valid]))
        }
    }

    emitdf <- get.emits(emitscen)
    if (is.null(emitdf))
        return(NULL)
    return(dice.2013(approx(emitdf$year, emitdf$gtcemit, seq(2010, lastyear, by=5))$y, lastyear+5))
}

## get.temps("RCP 8.5")
## get.temps("Weirdo")
## get.temps("DICE2013R")

## Guess 2100 for every emissions scenario
dat$temp.2100 <- NA
for (emitscen in unique(dat$`Emissions Scenario`)) {
    temps <- get.temps(emitscen, 2100)
    if (!is.null(temps)) {
        dt.2100 <- temps$dtemp.1900[temps$year == 2100]
        if (length(dt.2100) == 0)
            dt.2100 <- spline(temps$year, temps$dtemp.1900, method='natural', xout=2100)$y
        dat$temp.2100[dat$`Emissions Scenario` == emitscen] <- dt.2100
    }
}

if (F) {
    ## Check which scenarios we can't model yet
    for (emitscen in unique(dat$`Emissions Scenario`)) {
        if (is.null(get.temps(emitscen)))
            print(emitscen)
    }
}
