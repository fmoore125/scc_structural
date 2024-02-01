## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")
library(ggplot2)

dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

source("src/analysis/damage_funcs_lib.R")
dat <- multivar.prep(dat)

dist <- read.csv(file="outputs/distribution_v2.csv")

#drop outlier Nordhaus row
todrop=which(dat$`Central Value ($ per ton CO2)`>70000)
if(is.finite(todrop)) dist=dist[-which(dist$row==todrop),]

dat$row <- 1:nrow(dat)
dat$`Earth system` <- ifelse(dat$`Carbon Cycle` == "1.0" | dat$`Carbon Cycle` == "1", "1",
                      ifelse(dat$`Climate Model` == "1.0" | dat$`Climate Model` == "1", "1", "0"))
dat$`Inequality Aversion`[dat$`Inequality Aversion` == "Calibrated"] <- "1.0"
dat$`Inequality Aversion`[dat$`Inequality Aversion` == "1.0"] <- "1"
dat$`Persistent / Growth Damages`[dat$`Persistent / Growth Damages` == "Calibrated"] <- "1.0"
dat$`Persistent / Growth Damages`[dat$`Persistent / Growth Damages` == "1.0"] <- "1"
dat$`Tipping Points2`[dat$`Tipping Points2` == "-1.0" | dat$`Tipping Points2` == "-1"] <- "0"
dat$`TFP Growth`[is.na(dat$`TFP Growth`)] <- "0"
dat$`Population Growth`[is.na(dat$`Population Growth`)] <- "0"
dat$`Emissions Growth`[is.na(dat$`Emissions Growth`)] <- "0"
dat$`Transient Climate Response`[is.na(dat$`Transient Climate Response`)] <- "0"
dat$`Carbon Cycle2`[is.na(dat$`Carbon Cycle2`)] <- "0"
dat$`Equilibrium Climate Sensitivity`[is.na(dat$`Equilibrium Climate Sensitivity`)] <- "0"
dat$`Tipping Point Magnitude`[is.na(dat$`Tipping Point Magnitude`)] <- "0"
dat$`Damage Function`[is.na(dat$`Damage Function`)] <- "0"
dat$`Adaptation Rates`[is.na(dat$`Adaptation Rates`)] <- "0"
dat$`Income Elasticity`[is.na(dat$`Income Elasticity`)] <- "0"
dat$`Constant Discount Rate`[is.na(dat$`Constant Discount Rate`)] <- "0"
dat$`EMUC2`[is.na(dat$`EMUC2`)] <- "0"
dat$`PRTP2`[is.na(dat$`PRTP2`)] <- "0"
dat$`Risk Aversion (EZ Utility)`[is.na(dat$`Risk Aversion (EZ Utility)`)] <- "0"
dat$`Declining Discounting?`[is.na(dat$`Declining Discounting?`)] <- "0"
dat$log.scc.synth[dat$missing.scc.synth] <- NA # We can handle this

incrows <- !is.na(dat$sccyearformerge) & dat$sccyearformerge <= 2100 &
    dat$row %in% dist$row # Some have no distribution

get.numeric.branches <- function(datcol, splitpt) {
    if (is.na(splitpt))
        as.character(is.na(datcol)) # TRUE or FALSE
    else
        ifelse(is.na(datcol), "NA", as.character(datcol < splitpt)) # TRUE, FALSE, NA
}

cols <- c("Tipping Points", "Tipping Points2", "Persistent / Growth Damages", "Epstein-Zin",
          "Ambiguity/Model Uncertainty", "Limitedly-Substitutable Goods", "Inequality Aversion",
          "Learning", "Earth system", "TFP Growth", "Population Growth", "Emissions Growth",
          "Transient Climate Response", "Carbon Cycle2", "Equilibrium Climate Sensitivity",
          "Tipping Point Magnitude", "Damage Function", "Adaptation Rates", "Income Elasticity",
          "Constant Discount Rate", "EMUC2", "PRTP2", "Risk Aversion (EZ Utility)",
          "Backstop Price?", "Declining Discounting?", "Market Only Damages", "Other Market Failure?",
          "sccyearformerge", "discountrate", "log.scc.synth", "Year")
if (F) {
    for (col in cols)
        print(c(col, unique(dat[, col])[1:min(3, length(unique(dat[, col])))], "NAs:", sum(is.na(dat[, col]))))
}

predict.tree <- function(tree, datpred, dist, ndraw=1e3) {
    if (nrow(datpred) > 1) {
        draws <- c()
        for (ii in 1:nrow(datpred))
            draws <- c(draws, predict.tree(tree, datpred[ii,], dist, ndraw=ndraw))
        return(draws)
    }

    if (!is.na(tree$split) && tree$split == 'terminal')
        return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))

    if (!is.na(tree$split) && tree$split == 'categorical')
        branch <- as.character(datpred[, tree$col])
    else if (is.numeric(tree$split) || is.na(tree$split))
        branch <- get.numeric.branches(datpred[, tree$col], tree$split)
    else
        return(NULL) # unknown split

    if (is.null(tree$children[[branch]])) # branch not available in training data
        return(sample(dist$draw[dist$row %in% tree$rows], ndraw, replace=T))
    return(predict.tree(tree$children[[branch]], datpred, dist, ndraw=ndraw))
}

all.qs <- c(0, 0.001, 0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99, 0.999, 1)

load("outputs/rfdistsmodel-final.RData")

predict.forest <- function(forest, datpred, dist, ndraw=1e5, quants=all.qs) {
    qmat <- matrix(NA, length(forest), length(quants))
    for (ii in 1:length(forest)) {
        if (!is.list(forest[[ii]]))
            next
        draws <- predict.tree(forest[[ii]], datpred, dist, ndraw=ndraw)
        qvals <- quantile(draws, quants)
        qmat[ii,] <- qvals
    }

    colMeans(qmat, na.rm=T)
}

calc.varimport.tree <- function(tree, baseimport=1) {
    results <- data.frame(column=tree$column, import=baseimport)
    for (child in names(tree$children)) {
        childtree <- tree$children[[child]]
        if (is.na(childtree$split) || childtree$split != 'terminal') {
            rows <- calc.varimport.tree(childtree, baseimport=baseimport * length(childtree$rows) / length(tree$rows))
            results <- rbind(results, rows)
        }
    }
    results
}

calc.varimport.forest <- function(forest) {
    alltrees <- data.frame()
    count <- 0
    for (ii in 1:length(forest)) {
        if (!is.list(forest[[ii]]))
            next
        results <- calc.varimport.tree(forest[[ii]])
        alltrees <- rbind(alltrees, calc.varimport.tree(forest[[ii]]) %>% group_by(column) %>% summarize(import=max(import)))
        count <- count + 1
    }

    alltrees %>% group_by(column) %>% summarize(import=sum(import) / count)
}

varimport <- calc.varimport.forest(forest)
label.chg <- c('Other Market Failure?'='Feature: Other Market Failure',
               'Income Elasticity'='Uncertainty: Income Elasticity',
               'Equilibrium Climate Sensitivity'='Uncertainty: Equilibrium Climate Sensitivity',
               'Population Growth'='Uncertainty: Population Growth',
               'Tipping Point Magnitude'='Uncertainty: Tipping Point Magnitude',
               'Inequality Aversion'='Structural: Inequality Aversion',
               'Earth system'='Structural: Earth System',
               'Tipping Points'='Structural: Climate Tipping Points',
               'TFP Growth'='Uncertainty: TFP Growth',
               'Declining Discounting?'='Feature: Declining Discounting',
               'Ambiguity/Model Uncertainty'='Structural: Ambiguity/Model Uncertainty',
               'Carbon Cycle2'='Uncertainty: Carbon Cycle',
               'Adaptation Rates'='Uncertainty: Adaptation Rates',
               'Backstop Price?'='Feature: Backstop Price',
               'Tipping Points2'='Structural: Damage Tipping Points',
               'EMUC2'='Uncertainty: EMUC',
               'PRTP2'='Uncertainty: PRTP',
               'Limitedly-Substitutable Goods'='Structural: Limitedly-Substitutable Goods',
               'Damage Function'='Uncertainty: Damage Function',
               'log.scc.synth'='Quantity: Damage Function-based SCC',
               'Year'='Quantity: Publication Year',
               'discountrate'='Quantity: Standardized Discount Rate',
               'Learning'='Structural: Learning',
               'Persistent / Growth Damages'='Structural: Persistent/Growth Damages',
               'sccyearformerge'='Quantity: SCC Year',
               'Market Only Damages'='Feature: Market Only Damages',
               'Constant Discount Rate'='Uncertainty: Constant Discount Rate',
               'Emissions Growth'='Uncertainty: Emissions Growth',
               'Epstein-Zin'='Structural: Epstein-Zin Preferences',
               'Risk Aversion (EZ Utility)'='Uncertainty: EZ Risk Aversion',
               'Transient Climate Response'='Uncertainty: Transient Climate Response')
varimport$label <- label.chg[varimport$column]
varimport$label <- factor(varimport$label, levels=varimport$label[order(varimport$import)])

ggplot(varimport, aes(label, import)) +
    coord_flip() +
    geom_col() + theme_bw() +
    scale_y_continuous("Variable Importance", labels=scales::percent) + xlab(NULL)
ggsave("figures/rfdists-varimport.pdf", width=6.5, height=5)
write.csv(varimport, "outputs/rfdists-varimport.csv", row.names=F)

## Construct synthetic SCC

idealdat <- dat # include all rows

## set obvious ones
for (cc in which(names(dat) == 'TFP Growth'):which(names(dat) == 'Risk Aversion (EZ Utility)'))
    idealdat[, cc] <- "1"
idealdat$`Backstop Price?` <- "0"
idealdat$`Declining Discounting?` <- "1"
idealdat$`Market Only Damages` <- "0"
idealdat$`Other Market Failure?` <- "0"
## for damage-function-based scc - draw from literature values
idealdat$log.scc.synth[idealdat$missing.scc.synth] <- sample(idealdat$log.scc.synth[!idealdat$missing.scc.synth],
                                                            sum(idealdat$missing.scc.synth),replace=TRUE)
idealdat$Year <- 2020
idealdat$sccyearformerge <- 2020

bayespost <- read.csv("data/expert_survey/meta-analysis-distribution.csv")
bayespost2 <- dcast(bayespost, iterations ~ question, value.var='prob')
bayespost3 <- bayespost2[sample(nrow(bayespost2), nrow(dat), replace=T),]
bayespost3[, -1] <- ifelse(matrix(runif((ncol(bayespost3)-1)*nrow(bayespost3)), nrow(bayespost3), ncol(bayespost3)-1) < bayespost3[, -1], "1", "0")

idealdat <- cbind(idealdat[-c(which(names(dat) == 'Carbon Cycle'):which(names(dat) == 'Learning'), which(names(dat) == 'Earth system'))], bayespost3)
names(idealdat)[names(idealdat) == 'Tipping Points: Climate'] <- 'Tipping Points'
names(idealdat)[names(idealdat) == 'Tipping Points: Damages'] <- 'Tipping Points2'
names(idealdat)[names(idealdat) == 'Limited Substitutability'] <- 'Limitedly-Substitutable Goods'
names(idealdat)[names(idealdat) == 'Earth System'] <- 'Earth system'

discountsurvey <- read.csv("data/Drupp_et_al_2018_AEJ_Constant_SDR.csv")
idealdat$discountrate <- sample(discountsurvey$SDR[!is.na(discountsurvey$SDR)], nrow(idealdat), replace=TRUE)

predict.forest.all <- function(expdat, incrows, outfile) {
    pdf <- data.frame()
    for (ii in sample(which(incrows))) {
        if (ii %in% pdf$row)
            next
        print(ii)
        quants.pred <- predict.forest(forest, expdat[ii,], dist)
        pdf <- rbind(pdf, data.frame(row=ii, quant=all.qs, pred=quants.pred))
    }

    ## Construct Monte Carlo draws
    allsamp <- c()
    for (row in unique(pdf$row)) {
        inv.cdf <- approx(pdf$quant[pdf$row == row], pdf$pred[pdf$row == row])
        uu <- runif(1000)
        samples <- inv.cdf$y[findInterval(uu, inv.cdf$x)]
        allsamp <- c(allsamp, samples)
    }

    print(quantile(allsamp, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.)))
    print(mean(allsamp))

    save(pdf, allsamp, file=outfile)
}

predict.forest.all(idealdat, incrows, "outputs/rf_experiments/RFD_best.RData")

##---- E. Multiple constant discount rates (1, 1.5, 2.5, 3, 5)
discs <- c(1,1.5,2,2.5,3,5)

for (disc in discs) {
    print(disc)
    discdat <- idealdat
    discdat$discountrate <- disc
    predict.forest.all(discdat, 1:nrow(dat) %in% sample(which(incrows), 100),
                       paste0("outputs/rf_experiments/RFD_E_", disc, ".RData"))
}

#---- F. Multiple publication years (2000, 2010, 2020)
pubyears=c(2000,2010,2020)

for (pubyear in pubyears) {
    print(pubyear)
    pubyeardat <- idealdat
    pubyeardat$Year <- pubyear
    predict.forest.all(pubyeardat, 1:nrow(dat) %in% sample(which(incrows), 100),
                       paste0("outputs/rf_experiments/RFD_F_", pubyear, ".RData"))
}
