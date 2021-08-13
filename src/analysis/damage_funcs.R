## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

emitdf <- get.emits("RCP 8.5")
dmgfunc <- function(T) 0.0023888 * T^2
year0 <- 2020
gdp0 <- 84.54e12
discountrate <- 3

## calc.welfare <- function(emuc, cons) {
##     if (emuc == 1)
##         log(cons)
##     else
##         (cons^(1 - emuc) - 1) / (1 - emuc)
## }

calc.discounted.damages <- function(year0, years, damages, discountrate) {
    annualdamage <- approx(years, damages, year0:(year0 + 300), rule=2)$y
    sum(annualdamage / ((1 + discountrate / 100)^(0:300)))
}

calc.scc <- function(dmgfunc, emitdf, year0, gdp0, discountrate) {
    emit.gtc.peryr.by5 <- approx(emitdf$year, emitdf$gtcemit, seq(year0, year0 + 300, by=5), rule=2)$y
    tempdf <- dice.2013(emit.gtc.peryr.by5)
    tempdf$damage <- sapply(tempdf$dtemp.1900, dmgfunc)
    loss0 <- calc.discounted.damages(year0, tempdf$year, tempdf$damage, discountrate)

    emit.gtc.peryr.by5[1] <- emit.gtc.peryr.by5[1] + 1 / 5
    tempdf <- dice.2013(emit.gtc.peryr.by5)
    tempdf$damage <- sapply(tempdf$dtemp.1900, dmgfunc)
    loss1 <- calc.discounted.damages(year0, tempdf$year, tempdf$damage, discountrate)

    (loss1 - loss0) * gdp0 / 1e9 # fraction / Gt -> $/t
}

fund.3.8.xy <- list(x=c(0.3910223144784639, 1.032304099636741, 2.3195381421899324, 3.2903476907109495, 4.518681888946549, 5.677088738972496, 7.761805915931499),
                    y=c(-0.44337750120701674, -0.5218335926184218, -0.04707365484684298, 0.6980580440963468, 2.005659567619763, 3.5884609194785684, 7.331017649267743))

dice.2010.xy <- list(x=c(0.3901141670991178, 1.0132330046704723, 2.22106901920083, 3.190062272963155, 4.372340425531915, 5.748313440581214, 7.753373118837571),
                     y=c(0.301351858805858, 0.6594335067861166, 1.7199989270961857, 2.9660425942814226, 4.888418003326001, 7.551901722010621, 12.27415374711657))

page.09.xy <- list(x=c(0.39244940321743643, 1.0132330046704723, 2.232615464452517, 3.190062272963155, 4.3816813700051895, 5.727036844836533, 7.775428126621692),
                   y=c(0.4272839439944209, 0.6594335067861166, 1.6994796416501259, 2.9660425942814226, 5.919210342792768, 9.638431414623678, 15.789791320208142))

dice.2007 <- function(T) 0.0023888 * T^2
fund.3.5 <- function(T) -0.007457*T + 0.002685*T^2 - 0.000100*T^3 + 0.000001*T^4
dice.2013r <- function(T) 0.00267 * T^2
fund.3.8 <- approxfun(fund.3.8.xy$x, fund.3.8.xy$y)
dice.2010 <- approxfun(dice.2010.xy$x, dice.2010.xy$y)
page.09 <- approxfun(page.09.xy$x, page.09.xy$y)

hardcodeddfs <- list("DICE-2007"=dice.2007, "DICE 2007"=dice.2007, "DICE2007"=dice.2007,
                     "DICE2010"=dice.2010, "DICE 2010"=dice.2010,
                     "DICE-2013R"=dice.2013r, "DICE 2013R"=dice.2013r, "DICE2013"=dice.2013r, "DICE 2013"=dice.2013r,
                     "DICE2007, cubic damages (T^3)"=function(T) 0.0023888 * T^3,
                     "DICE-2007 + 0.00644 * (T / 4)^3"=function(T) dice.2007(T) + 0.00644 * (T / 4)^3,
                     "PAGE2009"=page.09,
                     "FUND 3.10"=fund.3.8, "FUND 3.6"=fund.3.5, "FUND 3.9"=fund.3.8, "FUND"=fund.3.8,
                     "FUND 3.7"=fund.3.8, "FUND 3.5"=fund.3.5, "FUND 3.4"=fund.3.5,
                     "Weitzman"=function(T) 0.0023888 * T^2 + 0.0000051 * T^6.754,
                     "HowardSterner"=function(T) 0.01145 * T^2, "HowardSterner (0.007438*T^2)"=function(T) 0.007438 * T^2)
hardcodeddfs.source <- list("DICE-2007"="DICE", "DICE 2007"="DICE", "DICE2007"="DICE",
                            "DICE2010"="DICE", "DICE 2010"="DICE",
                            "DICE-2013R"="DICE", "DICE 2013R"="DICE", "DICE2013"="DICE", "DICE 2013"="DICE",
                            "DICE2007, cubic damages (T^3)"="DICE+",
                            "DICE-2007 + 0.00644 * (T / 4)^3"="DICE+",
                            "PAGE2009"="PAGE",
                            "FUND 3.10"="FUND", "FUND 3.6"="FUND", "FUND 3.9"="FUND", "FUND"="FUND",
                            "FUND 3.7"="FUND", "FUND 3.5"="FUND", "FUND 3.4"="FUND",
                            "Weitzman"="Weitzman",
                            "HowardSterner"="HowardSterner", "HowardSterner (0.007438*T^2)"="HowardSterner")

## QUESTIONS:
## HowardSterner is with catstrophic and/or with productivity?
## TODO: "DietzStern"

simpleenv <- new.env(parent=baseenv())
assign("T", 1, envir=simpleenv)

unknowns <- c()
evaleddfs <- list()
for (dmgfunc in unique(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`)) {
    if (dmgfunc %in% names(hardcodedfs)) {
        dmg <- hardcodedfs[[dmgfunc]](1)
    } else {
        origdmgfunc <- dmgfunc
        dmgfunc <- gsub("−", "-", dmgfunc)

        dmg <- tryCatch({
            eval(str2lang(dmgfunc), simpleenv)
        }, error=function(e) {
            "error"
        })

        if (!is.na(dmg) && dmg == 'error') {
            dmgfunc <- gsub("(\\d)T", "\\1*T", dmgfunc)
            dmgfunc <- gsub("\\)\\(", ")*(", dmgfunc)
            dmgfunc <- gsub("(\\d)\\(", "\\1*(", dmgfunc)
            dmgfunc <- gsub("\\]", ")", gsub("\\[", "(", dmgfunc))
            dmgfunc <- gsub("\\bC\\b", "T", dmgfunc)
            dmgfunc <- gsub("T\\(t\\)", "T", dmgfunc)

            dmg <- tryCatch({
                eval(str2lang(dmgfunc), simpleenv)
            }, error=function(e) {
                "error"
            })
        }

        if (!is.na(dmg)) {
            if (dmg != 'error')
                evaleddfs[[origdmgfunc]] <- eval(str2lang(paste0("function(T) ", dmgfunc)), baseenv())
            else
                unknowns <- c(unknowns, origdmgfunc)
        }
    }
}

dat$scc.synth <- NA
dat$scc.source <- NA
for (dmgfunc in unique(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`)) {
    founddf <- NULL
    if (dmgfunc %in% names(hardcodedfs)) {
        founddf <- hardcodeddfs[[dmgfunc]]
        scc.source <- hardcodeddfs.source[[dmgfunc]]
    } else if (dmgfunc %in% names(evaleddfs)) {
        founddf <- evaleddfs[[dmgfunc]]
        scc.source <- "Explicit"
    }

    if (!is.null(founddf)) {
        dmg <- calc.scc(founddf, emitdf, year0, gdp0, discountrate)
        dat$scc.synth[dat$`Damage Function Info: Model, Commonly-Used Function, or Function` == dmgfunc] <- dmg
        dat$scc.source[dat$`Damage Function Info: Model, Commonly-Used Function, or Function` == dmgfunc] <- scc.source
    }
}

library(ggplot2)

ggplot(dat, aes(scc.synth, `Central Value ($ per ton CO2)`)) +
    geom_point(aes(colour=scc.source)) +
    geom_smooth(method='lm') +
    geom_abline(slope=1) +
    scale_x_log10() + scale_y_log10() + theme_bw() + xlab("Synthetic SCC") +
    scale_colour_discrete(name="Damage function")

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

dat$log.scc.synth <- log(dat$scc.synth)
dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

mod <- lm(log.scc.2020usd ~ `SCC Year` + log.scc.synth + discountrate, data=dat)
mod <- lm(log.scc.2020usd ~ `SCC Year` + scc.source + log.scc.synth + discountrate, data=dat)
