## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

ggplot(subset(dat, !is.na(temp.2100.source)), aes(temp.2100, `Central Value ($ per ton CO2)`)) +
    geom_point(aes(colour=temp.2100.source)) +
    geom_smooth(method='lm') +
    scale_y_log10() + scale_colour_discrete("Emissions scenario") +
    theme_bw() + xlab("Temperature in 2100")

source("src/analysis/damage_funcs_lib.R")

## Materials for Synthetic SCC section

mean(!is.na(dat$scc.synth))
table(dat$scc.source)
for (scc.source in unique(dat$scc.source)) {
    print(c(scc.source, sum(dat$scc.source == scc.source, na.rm=T), length(unique(dat$`Damage Function Info: Model, Commonly-Used Function, or Function`[!is.na(dat$scc.source) & dat$scc.source == scc.source]))))
}

library(ggplot2)

ggplot(plotdfs, aes(T, dmg, group=dmgfunc, colour=scc.source)) +
    geom_line() +
    coord_cartesian(ylim=c(-.1, 1)) +
    scale_y_continuous(labels=scales::percent) + scale_x_continuous(expand=c(0, 0)) +
    scale_colour_manual(name=NULL, breaks=c('FUND', 'DICE', 'DICE+', 'PAGE', 'HowardSterner', 'Weitzman', 'DietzStern', 'Explicit'),
                        values=c('#1b9e77','#e6ab02','#a6761d','#d95f02','#7570b3','#e7298a','#66a61e', '#808080')) +
    xlab(NULL) + ylab("Consumption damages") + theme_bw()
ggsave("outputs/figures/dmgcurves.pdf", width=5, height=4)

## Make heat map
plotdfs$dbin <- cut(plotdfs$dmg, c(-Inf, seq(0, 1, by=.02)))
plotdfs$tbin <- factor(floor(plotdfs$T * 5) / 5 + .05)
grid <- matrix(0, length(levels(plotdfs$tbin)), length(levels(plotdfs$dbin)))
for (ii in 1:nrow(plotdfs))
    grid[as.numeric(plotdfs$tbin[ii]), as.numeric(plotdfs$dbin[ii])] <- grid[as.numeric(plotdfs$tbin[ii]), as.numeric(plotdfs$dbin[ii])] + plotdfs$count[ii]
rownames(grid) <- levels(plotdfs$tbin)
colnames(grid) <- levels(plotdfs$dbin)

griddf <- melt(grid, c('T', 'dbin'))
griddf$dmg <- seq(0, 1, by=.02)[as.numeric(griddf$dbin)]

tops <- c('DICE-2016R2', 'FUND 3.7', 'HowardSterner (0.007438*T^2)', 'Weitzman', 'PAGE2009')

ggplot(griddf, aes(T, dmg)) +
    geom_raster(aes(fill=value)) +
    geom_line(data=subset(plotdfs, dmgfunc %in% tops), aes(group=dmgfunc, colour=scc.source)) +
    scale_y_continuous(labels=scales::percent, expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    scale_fill_gradient("Estimates:", trans='log10', na.value='#ffffff', low='#fff7ec', high='#7f0000') +
    scale_colour_manual(name="Examples:", breaks=c('FUND', 'DICE', 'PAGE', 'HowardSterner', 'Weitzman'),
                        values=c('#1b9e77','#e6ab02','#66a61e','#7570b3','#e7298a')) +
    xlab("Warming level (C)") + ylab("Consumption damages")
ggsave("outputs/figures/dmgheat.png", width=6.5, height=4)

## Plot distribution
bydmgfunc <- dat %>% group_by(ID_number) %>% summarize(dmgfunc=`Damage Function Info: Model, Commonly-Used Function, or Function`, citweight=1/length(Reference)) %>% group_by(dmgfunc) %>% summarize(citweight=sum(citweight))
plotdfs2 <- plotdfs %>% left_join(bydmgfunc)

ggplot(subset(plotdfs2, T %in% c(1, 2, 3, 4, 6) & dmg > -1),
       aes(pmin(.5, dmg), colour=scc.source, fill=scc.source)) +
    facet_grid(paste(T, '°C') ~ .) +
    geom_histogram(aes(weight=citweight)) +
    scale_y_continuous("Citations represented", expand=c(0, 0)) +
    scale_x_continuous("Consumption damages (%)", expand=c(0, 0), labels=scales::percent) +
    scale_colour_manual(name=NULL, breaks=c('FUND', 'DICE', 'DICE+', 'PAGE', 'HowardSterner', 'Weitzman', 'DietzStern', 'Explicit'),
                        values=c('#1b9e77','#e6ab02','#a6761d','#d95f02','#7570b3','#e7298a','#66a61e', '#808080')) +
    scale_fill_manual(name=NULL, breaks=c('FUND', 'DICE', 'DICE+', 'PAGE', 'HowardSterner', 'Weitzman', 'DietzStern', 'Explicit'),
                        values=c('#1b9e77','#e6ab02','#a6761d','#d95f02','#7570b3','#e7298a','#66a61e', '#808080')) +
    theme_bw()
ggsave("outputs/figures/dmgdist.pdf", width=5, height=4)

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
ggsave("outputs/figures/sccsynth.pdf", width=6.5, height=4)

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
