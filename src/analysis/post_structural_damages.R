setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

df <- multivar.prep(dat)

damagecols <- names(df)[c(1, 7:8, 12:13, 24:36)] # paper, model & damage func info, structural changes
df$damagecode <- sapply(1:nrow(df), function(ii) paste(df[ii, damagecols], collapse=', '))

colabbr <- c('BP', 'OF', 'MO', 'CC', 'CM', 'TC', 'TD', 'PD', 'EZ', 'MU', 'LS', 'IA', 'Le', 'AE')

structcols <- names(df)[c(12:13, 24, 26:36)] # paper, model & damage func info, structural changes
df$structcode <- sapply(1:nrow(df), function(ii) paste(colabbr[df[ii, structcols] != 0], collapse=''))

## Calculate linear and quadratic predictors for each row
df$dmg1 <- NA
df$dmg2 <- NA
for (emitscen in unique(df$`Emissions Scenario`)) {
    temps <- get.temps(emitscen, 2200)
    if (!is.null(temps)) {
        for (ii in which(df$`Emissions Scenario` == emitscen)) {
            discountrate  <- df$discountrate[ii]
            years <- (as.numeric(df$`SCC Year`[ii])+1):2200
            discfactors <- (1 + discountrate / 100)^-(years - as.numeric(df$`SCC Year`[ii]))
            df$dmg1[ii] <- sum(discfactors)
            df$dmg2[ii] <- sum(2 * discfactors * spline(temps$year, temps$dtemp.1900, method='natural', xout=years)$y)
        }
    }
}

summary(lm(`Central Value ($ per ton CO2)` ~ 0 + dmg1 + dmg2, data=df))
summary(lm(`Central Value ($ per ton CO2)` ~ 0 + damagecode : dmg1 + damagecode : dmg2, data=df))
summary(lm(`Central Value ($ per ton CO2)` ~ 0 + structcode : dmg1 + structcode : dmg2, data=df))
