setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

damagecols <- names(dat)[c(1, 7:8, 12:13, 24:36)] # paper, model & damage func info, structural changes
damagecols[damagecols == "Alternative ethical approaches (not Discounted Utilitarianism)"] <- "Alternative ethical approaches"
structcols <- names(dat)[c(12:13, 24, 26:36)] # paper, model & damage func info, structural changes
structcols[structcols == "Alternative ethical approaches (not Discounted Utilitarianism)"] <- "Alternative ethical approaches"

do.include.basescc <- T

if (do.include.basescc) {
    dat2 <- get.all.scc(dat)
    df <- multivar.prep(dat2)
} else {
    df <- multivar.prep(dat)
}

df$`SCC Year` <- as.numeric(as.character(df$`SCC Year`))
gamma <- .45 / 1000e9
year0 <- 1990
cons0 <- 3972 # average global
pop0 <- 5.28e9

## Example
## calc.scc(

df$damagecode <- sapply(1:nrow(df), function(ii) paste(df[ii, damagecols], collapse=', '))

colabbr <- c('BP', 'OF', 'MO', 'CC', 'CM', 'TC', 'TD', 'PD', 'EZ', 'MU', 'LS', 'IA', 'Le', 'AE')
df$structcode <- sapply(1:nrow(df), function(ii) paste(colabbr[df[ii, structcols] != 0], collapse='-'))

## Calculate linear and quadratic predictors for each row
fullscenarios <- paste(df$`Emissions Scenario`, df$`Socio-Economic Scenario`)
df$dmg1 <- NA
df$dmg2 <- NA
for (fullscen in unique(fullscenarios)) {
    temps <- get.temps(df$`Emissions Scenario`[fullscenarios == fullscen][1], 2200)
    if (is.null(temps))
        temps <- get.temps(df$`Socio-Economic Scenario`[fullscenarios == fullscen][1], 2200)

    if (!is.null(temps)) {
        grows <- get.cons.growth.percap(as.character(df$`Socio-Economic Scenario`[fullscenarios == fullscen][1]))
        if (is.null(grows))
            next

        allgrows <- spline(grows$year, grows$cons_growth_percap, method='natural', xout=year0:2200)$y
        allcons <- cumprod(c(cons0, (1 + allgrows / 100)))

        for (ii in which(fullscenarios == fullscen)) {
            if (is.na(df$discountrate[ii]) || is.na(df$`SCC Year`[ii]))
                next
            discountrate  <- df$discountrate[ii]
            years <- (as.numeric(df$`SCC Year`[ii])+1):2200
            discfactors <- (1 + discountrate / 100)^-(years - as.numeric(df$`SCC Year`[ii]))
            disccons <- discfactors * allcons[years - year0 + 1]
            df$dmg1[ii] <- sum(disccons)
            df$dmg2[ii] <- sum(2 * disccons * spline(temps$year, temps$dtemp.1900, method='natural', xout=years)$y)
        }
    }
}

if (!do.include.basescc)
    df$scc <- df$`Central Value ($ per ton CO2)`

summary(lm(scc ~ 0 + dmg1 + dmg2, data=df))
summary(lm(scc ~ 0 + modified:dmg1 + modified:dmg2, data=df))

plotdfs <- function(mod, prefix, siglim=1.64, addci=F, reqcodes=c()) {
    ses <- sqrt(diag(vcov(mod)))
    coeffs <- coef(mod)

    signif <- abs(coeffs / ses) > siglim
    codes <- sapply(names(coeffs)[signif], function(str) substring(str, nchar(prefix) + 1, nchar(str) - 5))
    validcodes <- unique(c(reqcodes, names(table(codes))[table(codes) >= 1]))

    TT <- seq(0, 6, length.out=100)
    pdf <- data.frame()
    for (code in validcodes) {
        alpha <- coeffs[paste0(prefix, code, ":dmg1")]
        beta <- coeffs[paste0(prefix, code, ":dmg2")]
        codepdf <- data.frame(code, TT, dmg=(alpha * TT + beta * TT^2) / (gamma * pop0))
        if (addci) {
            cols <- grep(paste0(prefix, code, ":"), colnames(vcov(mod)))
            se.fit <- sqrt(diag(cbind(TT, TT^2) %*% vcov(mod)[cols, cols] %*% rbind(TT, TT^2)))
            codepdf$se <- se.fit / (gamma * pop0)
        }
        pdf <- rbind(pdf, codepdf)
    }

    library(ggplot2)
    gp <- ggplot(pdf, aes(TT, dmg, colour=code)) +
        geom_line() + scale_x_continuous("Temperature change from 1990", expand=c(0, 0)) + theme_bw() +
        scale_y_continuous("Damages (% GDP)", labels=scales::percent)
    if (addci)
        gp <- gp + geom_ribbon(aes(ymin=dmg - 1.96*se, ymax=dmg + 1.96*se, group=code), alpha=.5, colour=NA)
    gp
}

plotdfs(lm(scc ~ 0 + modified:dmg1 + modified:dmg2, data=df), "modified", reqcodes='FALSE')

summary(lm(scc ~ 0 + damagecode : dmg1 + damagecode : dmg2, data=df))
summary(lm(scc ~ 0 + structcode : dmg1 + structcode : dmg2, data=df))

mod <- lm(scc ~ 0 + structcode : dmg1 + structcode : dmg2, data=df)
plotdfs(mod, "structcode")

summary(lm(scc ~ 0 + structcode : dmg1, data=df))

subset(df, structcode == "MO-CC-CM-PD")$scc

library(quantreg)
allowed <- names(table(df$structcode))[table(df$structcode) > 18] # below this get singular design matrix
mod <- rq(scc ~ 0 + structcode : dmg1 + structcode : dmg2, data=df[df$structcode %in% allowed, ])
coeffs <- coef(mod)
signif <- sign(summary(mod, alpha=.2)$coefficients[, 2]) == sign(summary(mod, alpha=.2)$coefficients[, 3])
codes <- sapply(names(coeffs)[signif], function(str) substring(str, 11, nchar(str) - 5))
validcodes <- names(table(codes))[table(codes) == 2]

TT <- seq(0, 6, length.out=100)
pdf <- data.frame()
for (code in validcodes) {
    alpha <- coeffs[paste0("structcode", code, ":dmg1")]
    beta <- coeffs[paste0("structcode", code, ":dmg2")]
    pdf <- rbind(pdf, data.frame(code, TT, dmg=alpha * TT + beta * TT^2))
}

ggplot(pdf, aes(TT, dmg, colour=code)) +
    geom_line() + scale_x_continuous("Temperature change from 1990", expand=c(0, 0)) + theme_bw() +
    scale_y_continuous("Damages (% GDP)", labels=scales::percent)

ggplot(pdf, aes(TT, dmg, colour=code)) +
    coord_cartesian(ylim=c(0, .1)) +
    geom_line() + scale_x_continuous("Temperature change from 1990", expand=c(0, 0)) + theme_bw() +
    scale_y_continuous("Damages (% GDP)", labels=scales::percent)
