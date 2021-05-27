## setwd("~/research/scciams/scc_structural")

## TODOs:
## From Fran:
## 
##  - Get DICE growth rates for more model versions
##  - Get some average PAGE and FUND growth rates too - I think these are exogenous right, so we should just be able to pull from model documentation?
## 
## Match DICE models, PAGE and FUND and for any remainders using Ramsey (not too many I think?) match to a reasonable SSP scenario

library(readxl)
library(dplyr)
library(ggplot2)

## Prepare data
source("src/data_cleaining_scripts/cleaning_master.R")

perref <- data.frame(table(dat$Reference))
names(perref) <- c("Reference", "rows")

df2 <- dat %>% left_join(perref)
df2$`SCC Year` <- as.numeric(df2$`SCC Year`)

df2$log.scc.2020usd <- log(df2$`Central Value ($ per ton CO2)`)
df2$log.scc.2020usd[!is.finite(df2$log.scc.2020usd)] <- NA

## Basic log model

mod <- lm(log.scc.2020usd ~ `SCC Year` + `Year` + effective.discount.rate.percent, data=df2, weights=1 / df2$rows)

## Try including
df2$`Ambiguity/Model Uncertainty`[is.na(df2$`Ambiguity/Model Uncertainty`)] <- "0"
df2$`Carbon Cycle`[is.na(df2$`Carbon Cycle`)] <- "0"
df2$`Climate Model`[is.na(df2$`Climate Model`)] <- "0"
df2$`Tipping Points`[is.na(df2$`Tipping Points`)] <- "0"
df2$`Epstein-Zin`[is.na(df2$`Epstein-Zin`)] <- "0"
df2$`Inequality Aversion`[is.na(df2$`Inequality Aversion`)] <- "0"
df2$`Learning`[is.na(df2$`Learning`)] <- "0"
df2$`Non-Substitutable Goods`[is.na(df2$`Non-Substitutable Goods`)] <- "0"
df2$`Persistent / Growth Damages`[is.na(df2$`Persistent / Growth Damages`)] <- "0"

mod2 <- lm(log.scc.2020usd ~ `SCC Year` + `Year` + effective.discount.rate.percent + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Tipping Points` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Non-Substitutable Goods` + `Persistent / Growth Damages`, data=df2, weights=1 / df2$rows)

library(stargazer)
stargazer(list(mod, mod2), single.row=T)

## Predict into the future
projdf <- tibble(Year=rep(2000:2100, 3), `SCC Year`=rep(c(2020, 2050, 2100), each=101))
projdf <- cbind(projdf, predict(mod, projdf, interval='confidence'))
projdf$scc.fit <- exp(projdf$fit + var(mod$resid)/2)
projdf$scc.upr <- exp(projdf$upr + var(mod$resid)/2)
projdf$scc.lwr <- exp(projdf$lwr + var(mod$resid)/2)

ggplot(projdf, aes(Year, scc.fit, colour=factor(`SCC Year`))) +
    geom_ribbon(data=subset(projdf, `SCC Year` == 2020), aes(ymin=scc.lwr, ymax=scc.upr), alpha=.5) +
    geom_line(lwd=1) + scale_x_continuous(expand=c(0, 0)) +
    scale_y_log10() + xlab("Year of publication") + ylab("SCC (2020 USD / tCO2)") +
    theme_bw() + scale_colour_discrete(name="Pulse Year:")

ggplot(projdf, aes(Year, scc.fit, colour=factor(`SCC Year`))) +
    geom_line() + coord_cartesian(ylim=c(0, 10000)) +
    xlab("Year of publication") + ylab("SCC (2020 USD / tCO2)") +
    theme_bw()

## Year of SCC surpassing VSL
vsl <- 7.4e6 * gdp.deflator$factor[gdp.deflator$year == 2020] / gdp.deflator$factor[gdp.deflator$year == 2006]
year.vsl <- log(vsl / exp(mod$coeff[1] + var(mod$resid)/2)) / sum(mod$coeff[2:3])

## Try cubic polynomials

df3 <- subset(df2, !is.na(`SCC Year`) & !is.na(`Year`))
mod <- lm(`Central Value ($ per ton CO2)` ~ poly(`SCC Year`, 3) + poly(`Year`, 3), data=df3, weights=1 / df3$rows)

