## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

## Approach:
## Create a full dataset with normal and baseline values, along with all change columns
## Then run this as a fixed effect model, with a unique FE for each base value conditions.

## Then I can handle both base-to-changed case and the case where
##   Row i has modifications x, y
##   Row j has modifications x, y, z

## Also, only do comparisons within a paper.

source("src/analysis/all_scc_lib.R")

allcols <- names(dat)[c(1, 8, 10, 12:16, 18:24, 26:36)]

df <- get.all.scc(dat)

names(df)[37] <- "Alternative ethical approaches"
allcols[grep("Alternative ethical approaches", allcols)] <- "Alternative ethical approaches"

library(lfe)

summary(felm(scc ~ modified | basecode, data=df))

df$log.scc <- log(df$scc)
df$log.scc[!is.finite(df$log.scc)] <- NA

df$`IAM Calibrated To (if applicable)`[is.na(df$`IAM Calibrated To (if applicable)`)] <- "None"
df$`Backstop Price?`[is.na(df$`Backstop Price?`)] <- "0"
df$`Other Market Failure?`[is.na(df$`Other Market Failure?`)] <- "0"
df$`Other Market Failure?`[df$`Other Market Failure?` != "0"] <- "1.0"
df$`Market Only Damages`[is.na(df$`Market Only Damages`)] <- "0"
df$`Carbon Cycle`[is.na(df$`Carbon Cycle`)] <- "0"
df$`Climate Model`[is.na(df$`Climate Model`)] <- "0"
df$`Tipping Points`[is.na(df$`Tipping Points`)] <- "0"
df$`Tipping Points2`[is.na(df$`Tipping Points2`)] <- "0"
df$`Persistent / Growth Damages`[is.na(df$`Persistent / Growth Damages`)] <- "0"
df$`Epstein-Zin`[is.na(df$`Epstein-Zin`)] <- "0"
df$`Ambiguity/Model Uncertainty`[is.na(df$`Ambiguity/Model Uncertainty`)] <- "0"
df$`Limitedly-Substitutable Goods`[is.na(df$`Limitedly-Substitutable Goods`)] <- "0"
df$`Inequality Aversion`[is.na(df$`Inequality Aversion`)] <- "0"
df$`Learning`[is.na(df$`Learning`)] <- "0"
df$`Alternative ethical approaches`[is.na(df$`Alternative ethical approaches`)] <- "0"

form <- as.formula(paste0("scc ~ `", paste(allcols[!(allcols %in% basemodelcols)], collapse="` + `"), "` + modified | basecode"))
summary(felm(form, data=df))

form <- as.formula(paste0("log.scc ~ `", paste(allcols[!(allcols %in% basemodelcols)], collapse="` + `"), "` + modified | basecode"))
mod <- felm(form, data=df)

form <- as.formula(paste0("log.scc ~ `", paste(allcols[!(allcols %in% c("IAM Calibrated To (if applicable)", basemodelcols))], collapse="` + `"), "` + modified | basecode"))
mod <- felm(form, data=df)

form <- as.formula(paste0("log.scc ~ `", paste(allcols[!(allcols %in% c("IAM Calibrated To (if applicable)", basemodelcols))], collapse="` + `"), "` + modified | `IAM Calibrated To (if applicable)` + basecode"))
mod <- felm(form, data=df)

pdf <- data.frame(name=names(coef(mod)), coeff=coef(mod), se=mod$se)
pdf <- pdf[!is.na(pdf$coeff),]
pdf$name <- factor(pdf$name, levels=rev(pdf$name))

library(ggplot2)
library(scales)

ggplot(pdf, aes(name, coeff, colour=pmin(4, abs(coeff/se)))) +
    coord_flip(ylim=c(-3, 2)) + theme_bw() +
    geom_hline(yintercept=0) + geom_point() +
    geom_linerange(aes(ymin=coeff - 1.96*se, ymax=coeff + 1.96*se)) +
    scale_colour_gradient2(name="|t value|", low=rainbow(4)[1], mid=rainbow(4)[2], high=rainbow(4)[3], midpoint=1.96) +
    xlab(NULL) + scale_y_continuous("Relative change in SCC", labels=scales::percent)

## TODO: Apply this to every pair of rows in each paper. Find all non-inval results.

compare.rows <- function(row1, row2) {
    sames <- data.frame()
    diffs <- data.frame(row=1:2)
    valid <- 'either'
    for (col in allcols) {
        if ((is.na(row1[1, col]) && is.na(row2[1, col])) || row1[1, col] == row2[2, col])
            sames[, col] <- row1[1, col]
        else {
            diffs[, col] <- c(row1[1, col], row2[2, col])
            if (valid == 'either') {
                if (is.na(row1[1, col]) && !is.na(row2[1, col]))
                    valid <- 'base1'
                if (!is.na(row1[1, col]) && is.na(row2[1, col]))
                    valid <- 'base2'
                if (!is.na(row1[1, col]) && !is.na(row2[1, col]))
                    valid <- '1diff'
            } else if (valid == 'base1') {
                if (!is.na(row1[1, col]))
                    valid <- 'inval'
            } else if (valid == 'base2') {
                if (!is.na(row2[1, col]))
                    valid <- 'inval'
            } else if (valid == '1diff') {
                if (!is.na(row1[1, col]) || !is.na(row2[1, col]))
                    valid <- 'inval'
            }
        }
    }

    list(valid=valid, sames=sames, diffs=diffs)
}
