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
distrf.sum <- distrf %>% group_by(row) %>% summarize(mu=mean(yhat), sd=sd(yhat), sd.draw=sd(draw), prob=length(draw) / nrow(distrf)) %>% filter(sd.draw > 0)

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
df3.sum <- df3 %>% group_by(group) %>% summarize(x=mean(value))
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
                    ##read.and.resid("outputs/rf_experiments/A_DICE.csv", 'None'),
                    ##read.and.resid("outputs/rf_experiments/B_EPA.csv", 'EPA'),
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
