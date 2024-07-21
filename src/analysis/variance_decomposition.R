setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")

dat$`SCC Year` <- as.numeric(dat$`SCC Year`)

dat$log.scc.2020usd <- log(dat$`Central Value ($ per ton CO2)`)
dat$log.scc.2020usd[!is.finite(dat$log.scc.2020usd)] <- NA

dat <- multivar.prep(dat)

source("src/analysis/damage_funcs_lib.R")

## Use MICE to infer remaining SCC synths
library(mice)
micedf <- dat[, c('scc.synth', 'Base IAM (if applicable)', 'IAM Calibrated To (if applicable)', 'Empirical Improvement or Sensitvity Analysis?',
                  'SCC Year', 'Backstop Price?', 'Other Market Failure?', 'Emissions Scenario', 'Socio-Economic Scenario',
                  'Reported Base Model SCC (if applicable)', 'Constant Discount Rate (%)', 'PRTP', 'Declining Discounting?', 'EMUC', 'RRA', 'IES',
                  'Market Only Damages', 'Damage Function Info: Model, Commonly-Used Function, or Function', 'Carbon Cycle', 'Climate Model',
                  'Tipping Points', 'Tipping Points2', 'Persistent / Growth Damages', 'Epstein-Zin', 'Ambiguity/Model Uncertainty',
                  'Limitedly-Substitutable Goods', 'Inequality Aversion', 'Learning',
                  'Alternative ethical approaches', 'log.scc.2020usd')]

micedf$`Base IAM (if applicable)` <- as.character(micedf$`Base IAM (if applicable)`)
micedf$`Base IAM (if applicable)`[grep("IWG", micedf$`Base IAM (if applicable)`)] <- "IWG"
micedf$`Base IAM (if applicable)`[grep("Combined", micedf$`Base IAM (if applicable)`)] <- "Combined"

micedf$`Emissions Scenario` <- as.character(micedf$`Emissions Scenario`)
micedf$`Emissions Scenario`[grep("Optimal", micedf$`Emissions Scenario`)] <- "Optimal"
micedf$`Emissions Scenario`[grep("BAU", micedf$`Emissions Scenario`)] <- "BAU"
micedf$`Emissions Scenario`[grep("Base", micedf$`Emissions Scenario`)] <- "BAU"
micedf$`Emissions Scenario`[grep("Reference", micedf$`Emissions Scenario`)] <- "BAU"
micedf$`Emissions Scenario`[grep("IWG|IAWG", micedf$`Emissions Scenario`)] <- "IWG"
micedf$`Emissions Scenario`[micedf$`Emissions Scenario` %in% c('Pulse', "Constraint - p(T > 2.5) = .5")] <- NA
micedf$`Emissions Scenario`[micedf$`Emissions Scenario` == 'A1b'] <- 'A1B'
micedf$`Emissions Scenario` <- sub("^SRES ", "", micedf$`Emissions Scenario`)
micedf$`Emissions Scenario` <- sub("^Constraint - ", "", micedf$`Emissions Scenario`)
micedf$`Emissions Scenario` <- sub(" stabilization$", "", micedf$`Emissions Scenario`)
micedf$`Emissions Scenario` <- sub(" CO2$", "", micedf$`Emissions Scenario`)

micedf$`Socio-Economic Scenario` <- as.character(micedf$`Socio-Economic Scenario`)
micedf$`Socio-Economic Scenario`[micedf$`Socio-Economic Scenario` %in% c('1.5', '2.0', '1.3', 'n/a', '1.0', 'N/A')] <- NA
micedf$`Socio-Economic Scenario`[grep("Optimal", micedf$`Socio-Economic Scenario`)] <- "Optimal"
micedf$`Socio-Economic Scenario`[grep("IWG|IAWG", micedf$`Socio-Economic Scenario`)] <- "IWG"

sapply(names(micedf), function(col) ifelse(is.numeric(micedf[, col]), NA, length(unique(micedf[!is.na(micedf[, col]), col]))))

names(micedf) <- paste0('X', 1:ncol(micedf))
imputed <- mice(micedf)
withall <- complete(imputed)

dat$scc.synth.imputed <- withall$X1
dat$log.scc.synth.imputed <- log(dat$scc.synth.imputed)

dat$`Earth System` <- "0"
dat$`Earth System`[paste(dat$`Carbon Cycle`, dat$`Climate Model`) != "0 0"] <- "1.0"

struccols <- sapply(c("Backstop Price?", "Other Market Failure?", "Market Only Damages", "Ambiguity/Model Uncertainty", "Earth System", "Tipping Points", "Tipping Points2", "Epstein-Zin", "Inequality Aversion", "Learning", "Limitedly-Substitutable Goods", "Persistent / Growth Damages", "Alternative ethical approaches"), function(col) which(names(dat) == col))

weights <- get.struct.weights(dat, struccols)

names(dat)[names(dat) == 'Tipping Points'] <- "Climate Tipping Point"
names(dat)[names(dat) == 'Tipping Points2'] <- "Damages Tipping Point"

get.anvdf <- function(mod) {
    anv <- anova(mod)
    anvdf <- as.data.frame(anv)
    anvdf$pred <- rownames(anvdf)
    anvdf
}

plot.anova <- function(mod, anvdf=NULL, minlabel=.02) {
    if (is.null(anvdf))
        anvdf <- get.anvdf(mod)
    anvdf <- anvdf[anvdf$pred != "Residuals",]

    anvdf2 <- anvdf %>%
        arrange(desc(pred)) %>%
        mutate(ypos = cumsum(`Sum Sq`)- 0.5*`Sum Sq` )

    ggplot(anvdf2, aes(x="", y=`Sum Sq`, fill=pred)) +
        geom_bar(stat="identity", width=1, color='white') +
        coord_polar("y", start=0) +
        theme_void() + theme(legend.position="none") +
        geom_text(data=subset(anvdf2, `Sum Sq` > minlabel), aes(y = ypos, label = pred), color = "white", size=3)
}

mod <- lm(log.scc.2020usd ~ `SCC Year` + `Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth.imputed, data=dat, weights=weights)
plot.anova(mod)

mod <- lm(log.scc.2020usd ~ `SCC Year` + `Year` + `Backstop Price?` + `Other Market Failure?` + discountrate + `Market Only Damages` + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth + missing.scc.synth, data=dat, weights=weights)
plot.anova(mod)  # sum is less than imputed

mod0 <- lm(log.scc.2020usd ~ `Year` + `Other Market Failure?` + discountrate + `Ambiguity/Model Uncertainty` + `Carbon Cycle` + `Climate Model` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth.imputed, data=dat[dat$`SCC Year` %in% 2010:2030 & dat$`Market Only Damages` == "0",], weights=weights[dat$`SCC Year` %in% 2010:2030 & dat$`Market Only Damages` == "0"])
anvdf0 <- get.anvdf(mod0)

mod1 <- lm(log.scc.2020usd ~ `Year` + discountrate + `Carbon Cycle` + `Climate Model` + `Persistent / Growth Damages` + log.scc.synth.imputed, data=dat[dat$`SCC Year` %in% 2010:2030 & dat$`Market Only Damages` == "1",], weights=weights[dat$`SCC Year` %in% 2010:2030 & dat$`Market Only Damages` == "1"])
anvdf1 <- get.anvdf(mod1)

anvdf <- anvdf0 %>% left_join(anvdf1, by='pred', suffix=c('.0', '.1'))
anvdf$`Sum Sq` <- anvdf$`Sum Sq.0`
both <- !is.na(anvdf$`Sum Sq.1`)
anvdf$`Sum Sq`[both] <- anvdf$`Sum Sq.0`[both] * sum(weights[dat$`Market Only Damages` == "0"]) / sum(weights) +
    anvdf$`Sum Sq.1`[both] * sum(weights[dat$`Market Only Damages` == "1"]) / sum(weights)
plot.anova(NULL, anvdf)

anvdf$`Sum Sq`[anvdf$pred == 'discountrate'] / sum(anvdf$`Sum Sq`)
anvdf$`Sum Sq`[anvdf$pred == 'log.scc.synth.imputed'] / sum(anvdf$`Sum Sq`)

unique(dat[dat$`SCC Year` %in% 2010:2030, struccols])

mod <- lm(log.scc.2020usd ~ `Other Market Failure?` + discountrate + `Ambiguity/Model Uncertainty` + `Earth System` + `Climate Tipping Point` + `Damages Tipping Point` + `Epstein-Zin` + `Inequality Aversion` + `Learning` + `Limitedly-Substitutable Goods` + `Persistent / Growth Damages` + `Alternative ethical approaches` + log.scc.synth.imputed, data=dat[dat$`SCC Year` %in% 2010:2030,], weights=weights[dat$`SCC Year` %in% 2010:2030])
plot.anova(mod, minlabel=.01)

anvdf <- get.anvdf(mod)
anvdf <- anvdf[anvdf$pred != "Residuals",]

anvdf$`Sum Sq`[anvdf$pred == 'discountrate'] / sum(anvdf$`Sum Sq`)
anvdf$`Sum Sq`[anvdf$pred == 'log.scc.synth.imputed'] / sum(anvdf$`Sum Sq`)

## Try dropping every other component
preds <- c("`Other Market Failure?`", "discountrate", "`Ambiguity/Model Uncertainty`", "`Earth System`", "`Climate Tipping Point`", "`Damages Tipping Point`", "`Epstein-Zin`", "`Inequality Aversion`", "`Learning`", "`Limitedly-Substitutable Goods`", "`Persistent / Growth Damages`", "`Alternative ethical approaches`", "log.scc.synth.imputed")

results <- data.frame()
for (ii in 1:length(preds)) {
    if (preds[ii] %in% c("`Earth System`", "`Persistent / Growth Damages`"))
        next
    print(preds[ii])
    usepreds <- preds[-ii]
    mod <- lm(as.formula(paste("log.scc.2020usd ~", paste(usepreds, collapse=" + "))),
              data=dat[dat$`SCC Year` %in% 2010:2030,], weights=weights[dat$`SCC Year` %in% 2010:2030])

    anvdf <- get.anvdf(mod)
    residpart <- anvdf$`Sum Sq`[anvdf$pred == "Residuals"] / sum(anvdf$`Sum Sq`)
    anvdf <- anvdf[anvdf$pred != "Residuals",]

    results <- rbind(results, data.frame(drop=preds[ii], pred=anvdf$pred, frac=anvdf$`Sum Sq` / sum(anvdf$`Sum Sq`), residpart))
}

ressum <- results %>% group_by(pred) %>% summarize(mu=mean(frac), min=min(frac), max=max(frac))

basemod <- lm(as.formula(paste("log.scc.2020usd ~", paste(preds, collapse=" + "))),
              data=dat[dat$`SCC Year` %in% 2010:2030,], weights=weights[dat$`SCC Year` %in% 2010:2030])
anvdf <- get.anvdf(basemod)
anvdf <- anvdf[anvdf$pred != "Residuals",]

anvdf$base <- anvdf$`Sum Sq` / sum(anvdf$`Sum Sq`)
ressum2 <- ressum %>% left_join(anvdf[, c('pred', 'base')])
ressum2$mu <- ressum2$mu / sum(ressum2$mu)

labels <- list("log.scc.synth.imputed"="Damage function", "discountrate"="Discount rate",
               "`Other Market Failure?`"="Other Market Failures")
ressum2$label <- sapply(ressum2$pred, function(pred) ifelse(pred %in% names(labels), labels[[pred]], gsub("`", "", pred)))

ressum3 <- ressum2 %>%
    arrange(desc(label)) %>%
    mutate(ypos = cumsum(mu)- 0.5*mu )

ggplot(ressum3, aes(x="", y=mu, fill=label)) +
    geom_bar(stat="identity", width=1, color='grey') +
    coord_polar("y", start=pi/2) +
    theme_void() + theme(legend.position="none") +
    geom_text(data=subset(ressum3, mu > .01), aes(y = ypos, label = label), color = "white", size=3)
ggsave("outputs/figures/anova-pie.pdf", width=5, height=5)

library(xtable)
ressum2$perc <- paste0(format(100 * ressum2$base, digits=1), "%")
ressum2$range <- paste0('[', round(1000 * ressum2$min) / 10, " - ", round(1000 * ressum2$max) / 10, '%]')
print(xtable(ressum2[, c('label', 'perc', 'range')]), include.rownames=F)

