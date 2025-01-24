## setwd("~/research/scciams/scc_structural")

library(ggplot2)
library(reshape2)
library(dplyr)

add.bigtbl.row <- function(label, mcs) {
    bigtbl <<- rbind(bigtbl, data.frame(label, `2.5%`=quantile(mcs, .025), `5%`=quantile(mcs, .05), `10%`=quantile(mcs, .1), `25%`=quantile(mcs, .25),
                                        `50%`=quantile(mcs, .5), `75%`=quantile(mcs, .75), `90%`=quantile(mcs, .9), `95%`=quantile(mcs, .95),
                                        `97.5%`=quantile(mcs, .975), Mean=mean(mcs)))
}

experiments <- c("best"="Synthetic SCC", "A_dice"="DICE", "B_epa"="EPA", "C_none"="No Structural Changes", "C_all"="All Structural Changes",
                 "E_1"="1% Discount Rate", "E_1.5"="1.5% Discount Rate", "E_2"="2% Discount Rate",
                 "E_2.5"="2.5% Discount Rate", "E_3"="3% Discount Rate", "E_5"="5% Discount Rate",
                 "F_2000"="Publication in 2000", "F_2010"="Publication in 2010", "F_2020"="Publication in 2020",
                 "G_2020"="SCC Year 2020", "G_2050"="SCC Year 2050", "G_2100"="SCC Year 2100",
                 "H_DICE2016r2"="DICE 2016R2 Damages", "H_FUND38"="FUND 3.8 Damages",
                 "H_HowardSterner"="Howard & Sterner Damages", "H_PAGE2009"="PAGE 2009 Damages")

bigtbl <- data.frame() # label c(0.025,0.05,0.25,0.5,0.75,0.95,0.975) mean

for (name in names(experiments)) {
    load(paste0("outputs/rf_experiments/RFD_", name, ".RData"))
    add.bigtbl.row(experiments[[name]], allsamp)
}

## Make the build-up graph

load("outputs/idealdat.RData")

get.bar2 <- function(base, modified, label, colour=NA) {
    data.frame(ymin=mean(base), ymax=mean(modified), ci25=quantile(sort(modified) - sort(base), .25) + mean(base),
               ci75=quantile(sort(modified) - sort(base), .75) + mean(base), label, colour)
}

load("outputs/rf_experiments/RFD_A_dice.RData")
barpdf <- get.bar2(0, allsamp, 'DICE')
allsamp.last <- allsamp

load("outputs/rf_experiments/RFD_I1.RData")
barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '+ Discounting'))
allsamp.last <- allsamp

load("outputs/rf_experiments/RFD_I2.RData")
barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '+ Damages'))
allsamp.last <- allsamp

for (cc in 75:83) {
    load(paste0("outputs/rf_experiments/RFD_I3_", cc, ".RData"))
    barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '+ Structural', names(idealdat)[cc]))
    allsamp.last <- allsamp
}

for (cc in 44:57) {
    load(paste0("outputs/rf_experiments/RFD_I4_", cc, ".RData"))
    barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '+ Uncertainty', names(idealdat)[cc]))
    allsamp.last <- allsamp
}

load(paste0("outputs/rf_experiments/RFD_best.RData"))
barpdf <- rbind(barpdf, get.bar2(0, allsamp, '= Synthetic\nSCC'))

load("outputs/rf_experiments/RFD_I5.RData")
barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '- Discounting'))
allsamp.last <- allsamp

load("outputs/rf_experiments/RFD_I6.RData")
barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '- Damages'))
allsamp.last <- allsamp

for (cc in 75:83) {
    load(paste0("outputs/rf_experiments/RFD_I7_", cc, ".RData"))
    barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '- Structural', names(idealdat)[cc]))
    allsamp.last <- allsamp
}

for (cc in c(44:47, 49, 51:56)) {
    load(paste0("outputs/rf_experiments/RFD_I8_", cc, ".RData"))
    barpdf <- rbind(barpdf, get.bar2(allsamp.last, allsamp, '- Uncertainty', names(idealdat)[cc]))
    allsamp.last <- allsamp
}

barpdf$label <- factor(barpdf$label, levels=c('DICE', '+ Discounting', '+ Damages', '+ Structural', '+ Uncertainty',
                                              '= Synthetic\nSCC',
                                              '- Discounting', '- Damages', '- Structural', '- Uncertainty'))

forboth <- barpdf %>% group_by(colour) %>% summarize(diffs=sum(abs(ymax - ymin)))
major <- forboth$colour[forboth$diffs > 10] # actually, min is 20

barpdf$colour[!(barpdf$colour %in% major) & !is.na(barpdf$colour)] <- "Minor"
barpdf$xmin <- as.numeric(barpdf$label) - .5
barpdf$xmax <- as.numeric(barpdf$label) + .5
barpdf$xmin[barpdf$label == '+ Structural'] <- 3.5 + seq(0, by=1/9, length.out=9)
barpdf$xmax[barpdf$label == '+ Structural'] <- 3.5 + seq(1/9, by=1/9, length.out=9)
barpdf$xmin[barpdf$label == '+ Uncertainty'] <- 4.5 + seq(0, by=1/14, length.out=14)
barpdf$xmax[barpdf$label == '+ Uncertainty'] <- 4.5 + seq(1/12, by=1/14, length.out=14)
barpdf$xmin[barpdf$label == '- Structural'] <- 8.5 + seq(0, by=1/9, length.out=9)
barpdf$xmax[barpdf$label == '- Structural'] <- 8.5 + seq(1/9, by=1/9, length.out=9)
barpdf$xmin[barpdf$label == '- Uncertainty'] <- 9.5 + seq(0, by=1/14, length.out=11)
barpdf$xmax[barpdf$label == '- Uncertainty'] <- 9.5 + seq(1/12, by=1/14, length.out=11)

ggplot(barpdf[1:27,]) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, fill=colour, colour=colour)) +
    geom_errorbar(aes(x=(xmin + xmax) / 2, ymin=ci25, ymax=ci75), width=.1) +
    geom_hline(yintercept=0) +
    scale_fill_manual("Components:", breaks=c(NA,
                                              "Earth system", "Epstein-Zin", "Learning", "Persistent / Growth Damages",
                                              "Minor",
                                              "Adaptation Rates", "Damage Function", "EMUC2",
                                              "Equilibrium Climate Sensitivity", "PRTP2", "TFP Growth"),
                      values=c("#808080", '#e31a1c', '#33a02c', '#6a3d9a', '#b15928',
                               '#a6cee3', '#8dd3c7', '#b2df8a',
                               '#cab2d6', '#ffff99', '#fb9a99',
                               '#fdbf6f'),
                      labels=c('Aggregate Bars', "Earth system", "Epstein-Zin", "Learning", "Persistent/Growth Damages",
                               '\nMinor (< 10 $/t change,\nstructural or uncertainty)\n',
                               "Adaptation Rates Uncertainty", "Damage Function Uncertainty", "EMUC Uncertainty",
                               "ECS Uncertainty", "PRTP Uncertainty", "TFP Growth Uncertainty")) +
    scale_colour_manual("Components:", breaks=c(NA,
                                                "Earth system", "Epstein-Zin", "Learning", "Persistent / Growth Damages",
                                                "Minor",
                                                "Adaptation Rates", "Damage Function", "EMUC2",
                                                "Equilibrium Climate Sensitivity", "PRTP2", "TFP Growth"),
                        values=c("#808080", '#e31a1c', '#33a02c', '#6a3d9a', '#b15928',
                                 '#a6cee3', '#8dd3c7', '#b2df8a',
                                 '#cab2d6', '#ffff99', '#fb9a99',
                                 '#fdbf6f'),
                        labels=c('Aggregate Bars', "Earth system", "Epstein-Zin", "Learning", "Persistent/Growth Damages",
                                 '\nMinor (< 10 $/t change,\nstructural or uncertainty)\n',
                                 "Adaptation Rates Uncertainty", "Damage Function Uncertainty", "EMUC Uncertainty",
                                 "ECS Uncertainty", "PRTP Uncertainty", "TFP Growth Uncertainty")) +
    scale_x_continuous(NULL, breaks=1:length(levels(barpdf$label)), labels=levels(barpdf$label)) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")

## Perform Monte Carlo setup

## Other bar plots

ggplot.compare <- function(experiments, outfile) {
    load("outputs/rf_experiments/RFD_best.RData")
    mu <- mean(allsamp)

    allpdf <- data.frame()
    for (ii in 1:length(experiments)) {
        load(paste0("outputs/rf_experiments/RFD_", names(experiments)[ii], ".RData"))
        allpdf <- rbind(allpdf, data.frame(mu=mean(allsamp), q50=median(allsamp),
                                           ci25=quantile(allsamp, .25), ci75=quantile(allsamp, .75), label=experiments[ii]))
    }
    allpdf$label <- factor(allpdf$label, levels=experiments)

    ggplot(allpdf, aes(label, mu)) +
        geom_hline(data=data.frame(mu), aes(yintercept=mu)) +
        geom_hline(yintercept=0) +
        geom_col() + geom_point(aes(y=q50)) +
        geom_errorbar(aes(ymin=ci25, ymax=ci75), width=.5) +
        theme_bw() + scale_y_continuous(NULL, expand=c(0, 0), limits=c(0, max(allpdf$ci75)*1.1)) +
        xlab("SCCO2 (USD/tCO2)")
    ggsave(file.path("figures/", outfile), width=3.25, height=2.5)
}

ggplot.compare(c("C_none"="None", "best"="Synthetic\nSCC", "C_all"="All Structural\nChanges"), "rfdists-compare-struct.pdf") # "B_epa"="EPA",
ggplot.compare(c("best"="Synthetic\nSCC", "E_1"="1%", "E_2"="2%", "E_5"="5%"), "rfdists-compare-drate.pdf")
ggplot.compare(c("best"="Synthetic\nSCC", "F_2000"="2000", "F_2010"="2010", "F_2020"="2020"), "rfdists-compare-pubyear.pdf")
ggplot.compare(c("best"="Synthetic\nSCC", "H_DICE2016r2"="DICE", "H_FUND38"="FUND",
                 "H_HowardSterner"="H & S", "H_PAGE2009"="PAGE"), "rfdists-compare-dmg.pdf")
ggplot.compare(c("best"="Synthetic\nSCC", "G_2020"="2020", "G_2050"="2050", "G_2100"="2100"),
               "rfdists-compare-pulseyear.pdf")
