## setwd("~/research/scciams/scc_structural")

library(ggplot2)
library(dplyr)

load("outputs/idealdat.RData")
rfdist_dir <- "outputs/Structural SCC RF Experiments"

get.bar2 <- function(base, modified, label, colour=NA) {
    data.frame(ymin=mean(base), ymax=mean(modified), ci25=quantile(sort(modified) - sort(base), .25) + mean(base),
               ci75=quantile(sort(modified) - sort(base), .75) + mean(base), label, colour)
}

has.structural.uncertainty <- c("Earth system", "Epstein-Zin", "Tipping Points", "Tipping Points2")

seqs <- c('forward', 1:29)

todo <- c("Discounting", "Damages", names(idealdat)[which(names(idealdat) == 'Ambiguity/Model Uncertainty'):which(names(idealdat) == 'Tipping Points2')], names(idealdat)[which(names(idealdat) == 'TFP Growth'):which(names(idealdat) == 'Risk Aversion (EZ Utility)')])

allpdf <- data.frame()

for (seq in seqs) {
    print(seq)
    if (seq == 'forward') {
        files <- c('RFD_I1', 'RFD_I2', paste0('RFD_I3_', 75:83), paste0('RFD_I4_', 44:57))
        effects <- c('Discounting', 'Damages', names(idealdat)[75:83], names(idealdat)[44:57])
    } else if (seq == 'backward') {
        files <- c(paste0('RFD_I8_', rev(c(44:47, 49, 51:56))), paste0('RFD_I7_', rev(75:83)), 'RFD_I6', 'RFD_I5', 'RFD_best')[-1] # drop first one (= DICE)
        effects <- c(NA, names(idealdat)[rev(c(44:47, 49, 51:56))], names(idealdat)[rev(75:83)], 'Damages', 'Discounting')[-1]
        effects[effects %in% has.structural.uncertainty] <- NA # don't measure these
    } else {
        load(paste0(rfdist_dir, "/RFD_J", seq, "_sequence.RData"))
        files <- c(paste0('RFD_J', seq, '_', 1:length(sequence)), 'RFD_best')
        effects <- c(sequence, todo[!(todo %in% sequence)])
    }

    load(file.path(rfdist_dir, "RFD_A_dice.RData"))
    allsamp.last <- allsamp

    for (ii in 1:length(files)) {
        load(paste0(rfdist_dir, "/", files[ii], ".RData"))
        allpdf <- rbind(allpdf, get.bar2(allsamp.last, allsamp, effects[ii]))
        allsamp.last <- allsamp
    }
}

allpdf2 <- allpdf %>% group_by(label) %>% summarize(ydiff=mean(ymax - ymin), ci25=mean(ci25 - ymin), ci75=mean(ci75 - ymin))

myorder <- c('Discounting', 'Damages', names(idealdat)[75:83], names(idealdat)[44:57])
label <- c('Discounting', 'Damages', rep('Structural', length(75:83)), rep('Uncertainty', length(44:57)))
colour <- c(NA, NA, names(idealdat)[75:83], names(idealdat)[44:57])

load(file.path(rfdist_dir, "/RFD_A_dice.RData"))
finpdf <- get.bar2(0, allsamp, 'DICE')

for (ii in 1:length(myorder)) {
    row <- allpdf2[!is.na(allpdf2$label) & allpdf2$label == myorder[ii],]
    ymin <- finpdf$ymax[nrow(finpdf)]
    finpdf <- rbind(finpdf, data.frame(ymin=ymin, ymax=ymin + row$ydiff,
                                       ci25=ymin + row$ci25, ci75=ymin + row$ci75, label=paste("+", label[ii]), colour=colour[ii]))
}

load(file.path(rfdist_dir, "/RFD_best.RData"))
finpdf <- rbind(finpdf, get.bar2(0, allsamp, '= Synthetic\nSCC'))

finpdf$label <- factor(finpdf$label, levels=c('DICE', '+ Discounting', '+ Damages', '+ Structural', '+ Uncertainty',
                                              '= Synthetic\nSCC'))

finpdf$diffs <- finpdf$ymax - finpdf$ymin
major <- finpdf$colour[abs(finpdf$diffs) > 10]

finpdf$colour[!(finpdf$colour %in% major) & !is.na(finpdf$colour)] <- "Minor"
finpdf$xmin <- as.numeric(finpdf$label) - .5
finpdf$xmax <- as.numeric(finpdf$label) + .5
finpdf$xmin[finpdf$label == '+ Structural'] <- 3.5 + seq(0, by=1/9, length.out=9)
finpdf$xmax[finpdf$label == '+ Structural'] <- 3.5 + seq(1/9, by=1/9, length.out=9)
finpdf$xmin[finpdf$label == '+ Uncertainty'] <- 4.5 + seq(0, by=1/14, length.out=14)
finpdf$xmax[finpdf$label == '+ Uncertainty'] <- 4.5 + seq(1/12, by=1/14, length.out=14)
finpdf$xmin[finpdf$label == '- Structural'] <- 8.5 + seq(0, by=1/9, length.out=9)
finpdf$xmax[finpdf$label == '- Structural'] <- 8.5 + seq(1/9, by=1/9, length.out=9)
finpdf$xmin[finpdf$label == '- Uncertainty'] <- 9.5 + seq(0, by=1/14, length.out=11)
finpdf$xmax[finpdf$label == '- Uncertainty'] <- 9.5 + seq(1/12, by=1/14, length.out=11)

gp <- ggplot(finpdf) +
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
    scale_x_continuous(NULL, breaks=1:length(levels(finpdf$label)), labels=levels(finpdf$label)) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")
ggsave("figures/figure5b.pdf", width=8, height=3)
