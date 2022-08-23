## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

ggplot(subset(dat, !is.na(temp.2100.source)), aes(temp.2100, `Central Value ($ per ton CO2)`)) +
    geom_point(aes(colour=temp.2100.source)) +
    geom_smooth(method='lm') +
    scale_y_log10() + scale_colour_discrete("Emissions scenario") +
    theme_bw() + xlab("Temperature in 2100")

source("src/analysis/damage_funcs_lib.R")

library(ggplot2)

ggplot(plotdfs, aes(T, dmg, group=dmgfunc, colour=scc.source)) +
    geom_line() +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_y_continuous(labels=scales::percent) +
    scale_colour_manual(name=NULL, breaks=c('FUND', 'DICE', 'DICE+', 'PAGE', 'HowardSterner', 'Weitzman', 'DietzStern', 'Explicit'),
                        values=c('#1b9e77','#e6ab02','#a6761d','#d95f02','#7570b3','#e7298a','#66a61e', '#808080')) +
xlab(NULL) + ylab("Consumption damages")

## ggplot(subset(plotdfs, scc.source == 'Explicit'), aes(T, dmg, group=dmgfunc, colour=added.by)) +
##     geom_line() +
##     coord_cartesian(ylim=c(-.1, 1)) +
##     scale_y_continuous(labels=scales::percent) +
##     scale_colour_manual(name="Added By", values=c('#ffff33', '#4daf4a', '#377eb8', '#e41a1c','#984ea3','#ff7f00')) +
## xlab(NULL) + ylab("Consumption damages")

## ggplot(subset(plotdfs, scc.source == 'Explicit' & added.by == 'Fran, Gernot'), aes(T, dmg, group=dmgfunc, colour=substring(dmgfunc, 1, 20))) +
##     geom_line() +
##     coord_cartesian(ylim=c(-.1, 1)) +
##     scale_y_continuous(labels=scales::percent) +
##     scale_colour_discrete(name="Damage function:") +
## xlab(NULL) + ylab("Consumption damages")

## T <- seq(0, 8, length.out=100)
## plot(T, dice.2010(T))

ggplot(subset(dat, scc.synth > .01 & scc.synth < 1e5 & `Central Value ($ per ton CO2)` > .01), aes(scc.synth, `Central Value ($ per ton CO2)`)) +
    geom_point(aes(colour=scc.source)) +
    geom_smooth(method='lm') +
    geom_abline(slope=1) +
    scale_x_log10() + scale_y_log10() +
    theme_bw() + xlab("Synthetic SCC") +
    scale_colour_discrete(name="Damage function")

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

dat$log.scc.synth <- log(dat$scc.synth)
dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

unique(dat$`Emissions Scenario`[dat$temp.2100])

dat2 <- subset(dat, scc.synth > .01)
dat2$temp.2100.known <- ifelse(is.na(dat2$temp.2100), 'NO', 'YES')
dat2$temp.2100[is.na(dat2$temp.2100)] <- 1

summary(lm(log.scc.2020usd ~ log.scc.synth, data=dat2))
summary(lm(log.scc.2020usd ~ `SCC Year` + discountrate + log.scc.synth, data=dat2))
summary(lm(log.scc.2020usd ~ `SCC Year` + discountrate + log.scc.synth + temp.2100 : temp.2100.known, data=dat2))
summary(lm(log.scc.2020usd ~ `SCC Year` + discountrate + temp.2100, data=dat2))
