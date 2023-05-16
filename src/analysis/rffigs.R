## setwd("~/research/scciams/scc_structural")

library(ggplot2)
library(reshape2)
library(dplyr)

bigtbl <- data.frame() # label c(0.025,0.05,0.25,0.5,0.75,0.95,0.975) mean

add.bigtbl.row <- function(label, mcs) {
    bigtbl <<- rbind(bigtbl, data.frame(label, `2.5%`=quantile(mcs, .025), `5%`=quantile(mcs, .05), `10%`=quantile(mcs, .1), `25%`=quantile(mcs, .25),
                                        `50%`=quantile(mcs, .5), `75%`=quantile(mcs, .75), `90%`=quantile(mcs, .9), `95%`=quantile(mcs, .95),
                                        `97.5%`=quantile(mcs, .975), Mean=mean(mcs)))
}

distrf=read.csv("outputs/distribution_structuralchangeweighted_withcovars_v2.csv")
distrf=distrf%>%filter(sccyear_from2020<=80)
load(file="outputs/randomforestmodel.Rdat")
distrf$yhat <- rfmod$predictions

## Question: How do residuals differ with param values
distrf$paramscore <- (distrf$TFP.Growth_param == "Yes") + (distrf$Population.Growth_param == "Yes") + (distrf$Emissions.Growth_param == "Yes") + (distrf$Transient.Climate.Response_param == "Yes") + (distrf$Carbon.Cycle2_param == "Yes") + (distrf$Equilibrium.Climate.Sensitivity_param == "Yes") + (distrf$Tipping.Point.Magnitude_param == "Yes") + (distrf$Damage.Function_param == "Yes") + (distrf$Adaptation.Rates_param == "Yes") + (distrf$Income.Elasticity_param == "Yes") + (distrf$Constant.Discount.Rate_param == "Yes") + (distrf$EMUC2_param == "Yes") + (distrf$PRTP2_param == "Yes") + (distrf$Risk.Aversion..EZ.Utility._param == "Yes")
distrf$resid <- rfmod_explained$residuals

distrf.sum <- distrf %>% group_by(row) %>% summarize(mu=mean(yhat), sd=sd(yhat), sd.draw=sd(draw), prob=length(draw) / nrow(distrf),
                                                     paramscore=mean(paramscore), resid=mean(resid)) %>% filter(sd.draw > 0)

## Who is getting low yhats?
lowrows <- distrf.sum$row[distrf.sum$mu < 40]
unique(subset(distrf, row %in% lowrows)[, 3:(ncol(distrf)-3)])

quantile(distrf.sum$resid[distrf.sum$paramscore == 0])
quantile(distrf.sum$resid[distrf.sum$paramscore > 1])

get.resids <- function(preds, times=1) {
  resids <- rep(NA, length(preds) * times)
  roundeddraws <- rep(round(preds), times)
  for (draw in unique(roundeddraws)) {
      print(draw)
      probs <- dnorm(draw, distrf.sum$mu, distrf.sum$sd) * distrf.sum$prob
      rows <- sample(distrf.sum$row, sum(roundeddraws == draw), prob=probs, replace=T)
      for (row in unique(rows)) {
          resids[roundeddraws == draw][rows == row] <- sample(rfmod_explained$residuals[distrf$row == row], sum(rows == row), replace=T)
      }
  }
  resids
}

df <- read.csv("outputs/rf_experiments/2020_2050_2100.csv")
names(df) <- c('X2020', 'X2050', 'X2100')

get.density <- function(xx, group) {
    dens <- density(xx, adjust=3)
    data.frame(x=dens$x, y=dens$y, group=group)
}

add.bigtbl.row('2020', df$X2020)
add.bigtbl.row('2050', df$X2050)
add.bigtbl.row('2100', df$X2100)

df2 <- rbind(get.density(df$X2020, '2020'), get.density(df$X2050, '2050'), get.density(df$X2100, '2100'))
df3 <- melt(df) %>% left_join(data.frame(variable=c('X2020', 'X2050', 'X2100'), group=c('2020', '2050', '2100')))
df3.sum <- df3 %>% group_by(group) %>% summarize(x=mean(value), ci25=quantile(value, .25), ci75=quantile(value, .75))
df3.sum$y <- sapply(1:nrow(df3.sum), function(ii) {
    approxfun(df2$x[df2$group == df3.sum$group[ii]], df2$y[df2$group == df3.sum$group[ii]])(df3.sum$x[df3.sum$group == df3.sum$group[ii]])})

ggplot(df2, aes(x, y, colour=group)) +
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, .0017)) +
    geom_line(lwd=0.75) + geom_segment(data=df3.sum, aes(x=x, xend=x, y=0, yend=y, colour=group)) +
    theme_bw() + scale_y_continuous(NULL, expand=c(0, 0)) + xlab("SCC (USD/tCO2)") +
    scale_colour_manual("Year:", breaks=c('2020', '2050', '2100'), values=c('purple2', '#a6cee3', '#b2df8a')) +
    theme(#legend.justification=c(1,1), legend.position=c(1,1),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave("figures/figure4-year.pdf", width=3.25, height=2.5)

read.and.resid <- function(filepath, group, times=1) {
    subdf <- read.csv(filepath)
    names(subdf) <- 'yhat'
    resids <- get.resids(subdf$yhat, times=times)
    if (times > 1)
        subdf <- data.frame(yhat=rep(subdf$yhat, times))
    subdf$value <- subdf$yhat + resids
    subdf$group <- group
    subdf
}

get.densities <- function(dflong) {
    df <- data.frame()
    for (group in unique(dflong$group)) {
        subdf <- get.density(dflong$value[dflong$group == group], group)
        df <- rbind(df, subdf)
    }
    df
}

ggplot.draws <- function(df.x, levels, ltitle, outfile) {
    for (group in unique(df.x$group))
        add.bigtbl.row(group, df.x$value[df.x$group == group])

    df2.x <- get.densities(df.x)
    df3.x.sum <- df.x %>% group_by(group) %>% summarize(x=mean(value))
    df3.x.sum$y <- sapply(1:nrow(df3.x.sum), function(ii) {
        approxfun(df2.x$x[df2.x$group == df3.x.sum$group[ii]], df2.x$y[df2.x$group == df3.x.sum$group[ii]])(df3.x.sum$x[df3.x.sum$group == df3.x.sum$group[ii]])
    })

    df2.x$group <- factor(df2.x$group, levels=levels)

    colours <- c('purple2', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f')

    ggplot(df2.x, aes(x, y, colour=group)) +
        coord_cartesian(xlim=c(0, 1000), ylim=c(0, .0017)) +
        geom_line(lwd=0.75) + geom_segment(data=df3.x.sum, aes(x=x, xend=x, y=0, yend=y, colour=group)) +
        theme_bw() + scale_y_continuous(NULL, expand=c(0, 0)) + xlab("SCCO2 (USD/tCO2)") +
        scale_colour_manual(ltitle, breaks=levels, values=colours[1:length(levels)]) +
        theme(#legend.justification=c(1,1), legend.position=c(1,1),
            axis.text.y=element_blank(), axis.ticks.y=element_blank())
    ggsave(file.path("figures/", outfile), width=3.25, height=2.5)
}

df.xstruct <- rbind(read.and.resid("outputs/rf_experiments/C_nostruc.csv", 'None'),
                    data.frame(yhat=NA, value=df$X2020, group='Expert + Lit.'),
                    read.and.resid("outputs/rf_experiments/C_allstruc.csv", 'All'))

ggplot.draws(df.xstruct, c('Expert + Lit.', 'None', 'All'), "Structural Set:", "figure4-struct.pdf")

df.xstruct2 <- rbind(read.and.resid("outputs/rf_experiments/C_nostruc.csv", 'None'),
                     read.and.resid("outputs/rf_experiments/A_DICE.csv", 'DICE'),
                    read.and.resid("outputs/rf_experiments/B_EPA.csv", 'EPA'),
                    data.frame(yhat=NA, value=df$X2020, group='Expert + Lit.'),
                    read.and.resid("outputs/rf_experiments/C_allstruc.csv", 'All'))

ggplot.draws(df.xstruct2, c('Expert + Lit.', 'None', 'DICE', 'EPA', 'All'), "Structural Set:", "figure4-struct2.pdf")

df.xdrate <- rbind(read.and.resid("outputs/rf_experiments/E_discountrates_1.csv", '1%'),
                   ##read.and.resid("outputs/rf_experiments/E_discountrates_1.5.csv", '1.5%'),
                   read.and.resid("outputs/rf_experiments/E_discountrates_2.csv", '2%'),
                   ##read.and.resid("outputs/rf_experiments/E_discountrates_2.5.csv", '2.5%'),
                   ##read.and.resid("outputs/rf_experiments/E_discountrates_3.csv", '3%'),
                   read.and.resid("outputs/rf_experiments/E_discountrates_5.csv", '5%'),
                   data.frame(yhat=NA, value=df$X2020, group='Expert + Lit.'))

ggplot.draws(df.xdrate, c('Expert + Lit.', '1%', '2%', '5%'), "Discount Rate:", "figure4-drate.pdf")

df.xpub <- rbind(read.and.resid("outputs/rf_experiments/F_pubyears_2000.csv", '2000'),
                 read.and.resid("outputs/rf_experiments/F_pubyears_2010.csv", '2010'),
                 read.and.resid("outputs/rf_experiments/F_pubyears_2020.csv", '2020'),
                 data.frame(yhat=NA, value=df$X2020, group='Expert + Lit.'))

ggplot.draws(df.xpub, c('Expert + Lit.', '2000', '2010', '2020'), "Pub. Year:", "figure4-pub.pdf")

df.xdmg <- rbind(read.and.resid("outputs/rf_experiments/H_damages_DICE2016r2.csv", 'DICE'),
                 read.and.resid("outputs/rf_experiments/H_damages_FUND38.csv", 'FUND'),
                 read.and.resid("outputs/rf_experiments/H_damages_HowardSterner.csv", 'H. & S.'),
                 read.and.resid("outputs/rf_experiments/H_damages_PAGE2009.csv", 'PAGE'),
                 data.frame(yhat=NA, value=df$X2020, group='Expert + Lit.'))

ggplot.draws(df.xdmg, c('Expert + Lit.', 'DICE', 'FUND', 'H. & S.', 'PAGE'), "Damage func.:", "figure4-dmg.pdf")

library(xtable)
print(xtable(bigtbl), include.rownames=F)

## setwd("~/research/scciams/scc_structural")

library(ggplot2)
library(dplyr)

distrf=read.csv("outputs/distribution_structuralchangeweighted_withcovars_v2.csv")
distrf=distrf%>%filter(sccyear_from2020<=80)
load(file="outputs/randomforestmodel.Rdat")
distrf$yhat <- rfmod$predictions

distrf.sum <- distrf %>% group_by(row) %>% summarize(mu=mean(yhat), sd=sd(yhat), sd.draw=sd(draw), prob=length(draw) / nrow(distrf)) %>%
    filter(sd.draw > 0)

get.resids <- function(preds, times=1) {
  resids <- rep(NA, length(preds) * times)
  roundeddraws <- rep(round(preds), times)
  for (draw in unique(roundeddraws)) {
      print(draw)
      probs <- dnorm(draw, distrf.sum$mu, distrf.sum$sd) * distrf.sum$prob
      rows <- sample(distrf.sum$row, sum(roundeddraws == draw), prob=probs, replace=T)
      for (row in unique(rows)) {
          resids[roundeddraws == draw][rows == row] <- sample(rfmod_explained$residuals[distrf$row == row], sum(rows == row), replace=T)
      }
  }
  resids
}

read.and.resid <- function(filepath, times=1) {
    subdf <- read.csv(filepath)
    names(subdf) <- 'yhat'
    resids <- get.resids(subdf$yhat, times=times)
    if (times > 1)
        subdf <- data.frame(yhat=rep(subdf$yhat, times))
    subdf$value <- subdf$yhat + resids
    subdf
}

df0 <- read.and.resid("outputs/rf_experiments/A_DICE.csv")
df1 <- read.and.resid("outputs/rf_experiments/I_1.csv")
df2 <- read.and.resid("outputs/rf_experiments/I_2.csv")
df3 <- read.and.resid("outputs/rf_experiments/I_3.csv")
df4 <- read.and.resid("outputs/rf_experiments/I_4.csv")

get.bar <- function(base, modified, label) {
    data.frame(ymin=mean(base), ymax=mean(modified), ci25=quantile(sort(modified) - sort(base), .25) + mean(base),
               ci75=quantile(sort(modified) - sort(base), .75) + mean(base), label)
}

pdf <- rbind(get.bar(0, df0$value, 'DICE'),
             get.bar(df0$value, df1$value, '+ Structural'),
             get.bar(df1$value, df2$value, '+ Uncertainty'),
             get.bar(df2$value, df3$value, '+ Discounting'),
             get.bar(df3$value, df4$value, '+ Damages'),
             get.bar(0, df4$value, '= Expert + Lit.'))
pdf$label <- factor(pdf$label, levels=c('DICE', '+ Structural', '+ Uncertainty', '+ Discounting', '+ Damages', '= Expert + Lit.'))

ggplot(pdf) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=as.numeric(label) - .5, xmax=as.numeric(label) + .5)) +
    geom_errorbar(aes(x=as.numeric(label), ymin=ci25, ymax=ci75), width=.5) +
    geom_hline(yintercept=0) +
    scale_x_continuous(NULL, breaks=as.numeric(pdf$label), labels=pdf$label) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")

df5 <- read.and.resid("outputs/rf_experiments/I_5.csv")
df6 <- read.and.resid("outputs/rf_experiments/I_6.csv")
df7 <- read.and.resid("outputs/rf_experiments/I_7.csv")

pdf <- rbind(get.bar(0, df0$value, 'DICE'),
             get.bar(df0$value, df5$value, '+ Damages'),
             get.bar(df5$value, df6$value, '+ Discounting'),
             get.bar(df6$value, df7$value, '+ Structural'),
             get.bar(df7$value, df4$value, '+ Uncertainty'),
             get.bar(0, df4$value, '= Expert + Lit.'))
pdf$label <- factor(pdf$label, levels=c('DICE', '+ Damages', '+ Discounting', '+ Structural', '+ Uncertainty', '= Expert + Lit.'))

ggplot(pdf) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=as.numeric(label) - .5, xmax=as.numeric(label) + .5)) +
    geom_errorbar(aes(x=as.numeric(label), ymin=ci25, ymax=ci75), width=.5) +
    geom_hline(yintercept=0) +
    scale_x_continuous(NULL, breaks=as.numeric(pdf$label), labels=pdf$label) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")

df8 <- read.and.resid("outputs/rf_experiments/I_8.csv")
df9 <- read.and.resid("outputs/rf_experiments/I_9.csv")

pdf <- rbind(get.bar(0, df0$value, 'DICE'),
             get.bar(df0$value, df1$value, '+ Structural'),
             get.bar(df1$value, df8$value, '+ Damages'),
             get.bar(df8$value, df9$value, '+ Discounting'),
             get.bar(df9$value, df4$value, '+ Uncertainty'),
             get.bar(0, df4$value, '= Expert + Lit.'))
pdf$label <- factor(pdf$label, levels=c('DICE', '+ Structural', '+ Damages', '+ Discounting', '+ Uncertainty', '= Expert + Lit.'))

ggplot(pdf) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=as.numeric(label) - .5, xmax=as.numeric(label) + .5)) +
    geom_errorbar(aes(x=as.numeric(label), ymin=ci25, ymax=ci75), width=.5) +
    geom_hline(yintercept=0) +
    scale_x_continuous(NULL, breaks=as.numeric(pdf$label), labels=pdf$label) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")

get.bar2 <- function(base, modified, label, colour=NA) {
    data.frame(ymin=mean(base), ymax=mean(modified), ci25=quantile(sort(modified) - sort(base), .25) + mean(base),
               ci75=quantile(sort(modified) - sort(base), .75) + mean(base), label, colour)
}

library(data.table)
sampdat=fread(file="outputs/rf_experiments/basepredictiondata.csv")
sampdat <- as.data.frame(sampdat)

pdf <- get.bar2(0, df0$value, 'DICE')

lastdf <- df0
for (cc in grep("struc",names(sampdat))) {
    nextdf <- read.and.resid(paste0("outputs/rf_experiments/J_", cc, ".csv"))
    subpdf <- get.bar2(lastdf$value, nextdf$value, '+ Structural', names(sampdat)[cc])
    pdf <- rbind(pdf, subpdf)
    lastdf <- nextdf
}

pdf <- rbind(pdf, get.bar2(df1$value, df8$value, '+ Damages'),
             get.bar2(df8$value, df9$value, '+ Discounting'))

lastdf <- df9
for (cc in grep("param",names(sampdat))) {
    nextdf <- read.and.resid(paste0("outputs/rf_experiments/J_", cc, ".csv"))
    subpdf <- get.bar2(lastdf$value, nextdf$value, '+ Uncertainty', names(sampdat)[cc])
    pdf <- rbind(pdf, subpdf)
    lastdf <- nextdf
}

pdf <- rbind(pdf, get.bar2(0, df4$value, '= Expert + Lit.'))

pdf$label <- factor(pdf$label, levels=c('DICE', '+ Structural', '+ Damages', '+ Discounting', '+ Uncertainty', '= Expert + Lit.'))

pdf2 <- pdf
## pdf2$sign <- sign(pdf2$ymax - pdf2$ymin)
## pdf2$label2
## pdf2 %>% group_by(label, ) %>%
##     summarize(ymin=ymin[1], ymax=ymax[length(ymax)], ci25=min(ci25), ci75=max(ci75), colour=colour[1]
pdf2$colour[abs(pdf2$ymax - pdf2$ymin) < 10] <- "Minor"
pdf2$xmin <- as.numeric(pdf2$label) - .5
pdf2$xmax <- as.numeric(pdf2$label) + .5
pdf2$xmin[pdf2$label == '+ Structural'] <- 1.5 + seq(0, by=1/9, length.out=9)
pdf2$xmax[pdf2$label == '+ Structural'] <- 1.5 + seq(1/9, by=1/9, length.out=9)
pdf2$xmin[pdf2$label == '+ Uncertainty'] <- 4.5 + seq(0, by=1/14, length.out=14)
pdf2$xmax[pdf2$label == '+ Uncertainty'] <- 4.5 + seq(1/12, by=1/14, length.out=14)

ggplot(pdf2) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, fill=colour)) +
    geom_errorbar(aes(x=(xmin + xmax) / 2, ymin=ci25, ymax=ci75), width=.1, position='dodge') +
    geom_hline(yintercept=0) +
    scale_x_continuous(NULL, breaks=as.numeric(pdf$label), labels=pdf$label) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")

df3.sum.x <- df3.sum
names(df3.sum.x)[2] <- 'mu'

pdfextra <- rbind(cbind(df.xdmg %>% group_by(group) %>% summarize(mu=mean(value), ci25=quantile(value, .25), ci75=quantile(value, .75)), label="Alternative\nDamages"),
                  cbind(df.xdrate %>% group_by(group) %>% summarize(mu=mean(value), ci25=quantile(value, .25), ci75=quantile(value, .75)), label="Alternative\nDiscounting"),
                  cbind(df3.sum.x[, 1:4], label="Alternative\nPulse Years"))

pdf2$label <- factor(pdf2$label, levels=c('DICE', '+ Structural', '+ Damages', '+ Discounting', '+ Uncertainty', '= Expert + Lit.', 'Alternative\nDamages', 'Alternative\nDiscounting', 'Alternative\nPulse Years'))
pdfextra$label <- factor(pdfextra$label, levels=c('DICE', '+ Structural', '+ Damages', '+ Discounting', '+ Uncertainty', '= Expert + Lit.', 'Alternative\nDamages', 'Alternative\nDiscounting', 'Alternative\nPulse Years'))
pdfextra <- subset(pdfextra, !(group %in% c("Expert + Lit.", "DICE", "FUND", "2020")))

ggplot(pdf2) +
    geom_rect(aes(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, fill=colour, colour=colour)) +
    geom_errorbar(aes(x=(xmin + xmax) / 2, ymin=ci25, ymax=ci75), width=.1) +
    geom_segment(aes(x=pdf2$xmax[nrow(pdf2)] + .05, y=pdf2$ymax[nrow(pdf2)],
                     xend=as.numeric(pdfextra$label[nrow(pdfextra)]) + .5, yend=pdf2$ymax[nrow(pdf2)])) +
    geom_label(data=pdfextra, aes(x=as.numeric(label), y=mu, label=group)) +
    geom_hline(yintercept=0) +
    scale_fill_manual("Components:", breaks=c(NA, "Persistent...Growth.Damages_struc", "Learning_struc", "Earth_system_struc", "Minor",
                                                "TFP.Growth_param", "Population.Growth_param", "Emissions.Growth_param",
                                                "Equilibrium.Climate.Sensitivity_param", "Tipping.Point.Magnitude_param", "Damage.Function_param",
                                                "Risk.Aversion..EZ.Utility._param"),
                        values=c("#808080", '#e31a1c', '#33a02c', '#6a3d9a', '#b15928',
                                 '#a6cee3', '#8dd3c7', '#b2df8a',
                                 '#cab2d6', '#ffff99', '#fb9a99',
                                 '#fdbf6f'),
                      labels=c('Aggregate Bars', 'Persistence/Growth Damages', 'Learning', 'Earth System', '\nMinor (< 10 $/t change,\nstructural or uncertainty)\n',
                               "TFP Growth", "Population Growth", "Emissions Growth",
                               "Equilibrium Climate Sensitivity", "Tipping Point Magnitude", "Damage Function",
                               "Risk Aversion (EZ Utility)")) +
    scale_colour_manual("Components:", breaks=c(NA, "Persistent...Growth.Damages_struc", "Learning_struc", "Earth_system_struc", "Minor",
                                                "TFP.Growth_param", "Population.Growth_param", "Emissions.Growth_param",
                                                "Equilibrium.Climate.Sensitivity_param", "Tipping.Point.Magnitude_param", "Damage.Function_param",
                                                "Risk.Aversion..EZ.Utility._param"),
                        values=c("#808080", '#e31a1c', '#33a02c', '#6a3d9a', '#b15928',
                                 '#a6cee3', '#8dd3c7', '#b2df8a',
                                 '#cab2d6', '#ffff99', '#fb9a99',
                                 '#fdbf6f'),
                      labels=c('Aggregate Bars', 'Persistence/Growth Damages', 'Learning', 'Earth System', '\nMinor (< 10 $/t change,\nstructural or uncertainty)\n',
                               "TFP Growth", "Population Growth", "Emissions Growth",
                               "Equilibrium Climate Sensitivity", "Tipping Point Magnitude", "Damage Function",
                               "Risk Aversion (EZ Utility)")) +
    scale_x_continuous(NULL, breaks=1:length(levels(pdf2$label)), labels=levels(pdf2$label)) +
    theme_bw() + ylab("Social cost of carbon (USD/tCO2)")
