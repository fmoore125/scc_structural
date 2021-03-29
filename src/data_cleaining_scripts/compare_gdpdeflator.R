## setwd("~/research/scciams/scc_structural")

## James's code

library(readxl)
library(dplyr)
library(ggplot2)

deflator.growth <- read.csv("data/old/gdp-deflator-growth.csv")

gdp.deflator <- data.frame(year=1930:2020, growth=as.numeric(deflator.growth[1, -1:-2]))
gdp.deflator$factor <- 1
gdp.deflator$factor[1] <- 1
for (ii in 2:nrow(gdp.deflator))
    gdp.deflator$factor[ii] <- gdp.deflator$factor[ii-1] * (1 + gdp.deflator$growth[ii-1] / 100)
factor.2020 <- gdp.deflator$factor[gdp.deflator$year == 2020]

df <- read_xlsx("data/data_collection/SCC Meta-Analysis Data Template_Revised.xlsx", skip=2, sheet=1)
df$`SCC Dollar Year` <- as.numeric(df$`SCC Dollar Year`)
df$`SCC Year` <- as.numeric(df$`SCC Year`)
df$`Central Value ($ per ton CO2)` <- as.numeric(df$`Central Value ($ per ton CO2)`)

df2 <- df %>% left_join(gdp.deflator, by=c(`SCC Dollar Year`='year'))
df2$scc.2020usd <- df2$`Central Value ($ per ton CO2)` * factor.2020 / df2$factor

## Fran's code

source("src/data_cleaining_scripts/cleaning_standardizedollaryears.R")

## Comparison

compare <- data.frame(james=df2$scc.2020usd, fran=dat$`Central Value ($ per ton CO2)`, dat$`SCC Dollar Year`)
sum(is.na(compare$james) & !is.na(compare$fran))
sum(!is.na(compare$james) & is.na(compare$fran))

ggplot(compare, aes(james, fran)) +
    geom_point() + scale_x_log10() + scale_y_log10()
