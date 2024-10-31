## Generates tables and figures describing the content of the literature review.

## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

authors <- read.csv("data/authorlist.csv")
authnum <- length(unique(as.vector(as.matrix(authors[, -1])))) - 1
authors$count <- rowSums(authors[, -1] != "")
authors$included <- authors$ID_number %in% dat$ID_number

dat$`SCC Year` <- as.numeric(as.character(dat$`SCC Year`))
dat$uncinfo.ext <- ifelse(!is.na(dat$Min) | !is.na(dat$`0.1th`) | !is.na(dat$`99.9th`) | !is.na(dat$`Max`), 1, NA)
dat$uncinfo.tail <- ifelse(!is.na(dat$`1th`) | !is.na(dat$`2.5th`) | !is.na(dat$`5th`) | !is.na(dat$`95th`) | !is.na(dat$`97.5th`) | !is.na(dat$`99th`), 1, NA)
dat$uncinfo.central <- ifelse(!is.na(dat$`10th`) | !is.na(dat$`17th`) | !is.na(dat$`25th`) | !is.na(dat$`75th`) | !is.na(dat$`83rd`) | !is.na(dat$`90th`) | (!is.na(dat$`Central Value ($ per ton CO2)`) & !is.na(dat$`50th`) & dat$`Central Value ($ per ton CO2)` != dat$`50th`), 1, NA)
dat$tail.level <- ifelse(!is.na(dat$`0.1th`), 99.9, ifelse(!is.na(dat$`99th`), 99, ifelse(!is.na(dat$`97.5th`), 97.5, ifelse(!is.na(dat$`95th`), 95, ifelse(!is.na(dat$`90th`), 90, ifelse(!is.na(dat$`83rd`), 83, ifelse(!is.na(dat$`75th`), 75, ifelse((!is.na(dat$`Central Value ($ per ton CO2)`) & !is.na(dat$`50th`) & dat$`Central Value ($ per ton CO2)` != dat$`50th`), 50, NA))))))))

dat.paper <- dat %>% group_by(ID_number) %>% summarize(numest=length(`Central Value ($ per ton CO2)`))

tbl.unique <- data.frame(label=c("Papers", "Estimates", "Authors"), N=nrow(dat), unique=as.integer(c(length(unique(dat$ID_number)), length(unique(dat$`Central Value ($ per ton CO2)`)), authnum)))
tbl.value <- rbind(data.frame(label="Authors per paper", N=sum(dat$ID_number %in% authors$ID_number), mean=mean(authors$count), median=median(authors$count), min=min(authors$count), max=max(authors$count)),
                   data.frame(label="Estimates per paper", N=nrow(dat), mean=mean(dat.paper$numest), median=median(dat.paper$numest), min=min(dat.paper$numest), max=max(dat.paper$numest)))
tbl.present <- data.frame(label=c(), N=c(), present=c())
for (col in c('SCC Year', 'Central Value ($ per ton CO2)', 'Reported Base Model SCC (if applicable)', 'Constant Discount Rate (%)', 'PRTP', 'EMUC', 'RRA', 'IES', 'tail.level'))
    tbl.value <- rbind(tbl.value, data.frame(label=col, N=sum(!is.na(dat[, col])), mean=mean(dat[, col], na.rm=T), median=median(dat[, col], na.rm=T), min=min(dat[, col], na.rm=T), max=max(dat[, col], na.rm=T)))
for (col in c('Backstop Price?', 'Other Market Failure?', 'Declining Discounting?', 'Market Only Damages', 'Carbon Cycle', 'Climate Model', 'Tipping Points', 'Tipping Points2', 'Persistent / Growth Damages', 'Epstein-Zin', 'Ambiguity/Model Uncertainty', 'Limitedly-Substitutable Goods', 'Inequality Aversion', 'Learning', 'Alternative ethical approaches (not Discounted Utilitarianism)', 'uncinfo.ext', 'uncinfo.tail', 'uncinfo.central', 'TFP Growth', 'Population Growth', 'Emissions Growth', 'Transient Climate Response', 'Carbon Cycle2', 'Equilibrium Climate Sensitivity', 'Tipping Point Magnitude', 'Damage Function', 'Adaptation Rates', 'Income Elasticity', 'Constant Discount Rate', 'EMUC2', 'PRTP2', 'Risk Aversion (EZ Utility)'))
    tbl.present <- rbind(tbl.present, data.frame(label=col, N=nrow(dat), present=sum(!is.na(dat[, col]))))
for (col in c('Emissions Scenario', 'Socio-Economic Scenario', 'Damage Function Info: Model, Commonly-Used Function, or Function'))
    tbl.unique <- rbind(tbl.unique, data.frame(label=col, N=sum(!is.na(dat[, col])), unique=length(unique(dat[, col]))))

library(xtable)

print(xtable(tbl.unique), include.rownames=F)
print(xtable(tbl.present), include.rownames=F)
print(xtable(tbl.value), include.rownames=F)

library(ggplot2)

dat2 <- dat %>% group_by(ID_number) %>% summarize(Year=Year[1])
ggplot(dat2, aes(Year)) +
    geom_bar() + theme_bw() + scale_x_continuous("Publication year", expand=c(0, 0)) +
    scale_y_continuous("Number of publications", expand=c(0, 0))

dat$scc.year <- as.numeric(dat$`SCC Year`)
ggplot(dat, aes(scc.year)) +
    geom_histogram() + theme_bw() + scale_x_continuous("Year of SCC Pulse", expand=c(0, 0)) +
    scale_y_continuous("Number of estimates", expand=c(0, 0))

dat$discountshape <- ifelse(!is.na(dat$`Declining Discounting?`), "Declining", "Static")
dat$discounttype <- ifelse(!is.na(dat$`Constant Discount Rate (%)`), "Constant",
                    ifelse(!is.na(dat$PRTP) & !is.na(dat$EMUC), "Ramsey",
                    ifelse(!is.na(dat$PRTP) & !is.na(dat$RRA) & !is.na(dat$IES), "Epstein-Zin", "Other")))

ggplot(dat, aes(Year, discountrate)) +
    ##geom_jitter(aes(colour=discounttype, shape=discountshape), width=.1) +
    geom_point(aes(colour=discounttype, shape=discountshape), size=1) +
    geom_smooth() +
    scale_shape_manual("Parameters:", breaks=c('Static', 'Declining'), values=c(16, 25)) +
    scale_colour_discrete("Approach:", breaks=c("Constant", "Ramsey", "Epstein-Zin", "Other")) +
    theme_bw() + xlab("Publication Year") + ylab("SCC estimate discount rate (%)")

labels <- list('3IAW'='IWG Combination', 'Anal'='Analytic IAM', 'Broc'='Other', 'Comb'='Other', 'COME'='Other',
               'DSIC'='DSICE', 'FAIR'='Other', 'GHKT'='Other', 'Golo'='Other', 'IWG2'='IWG Combination',
               'matD'='matDICE', 'RESP'='RESPONSE', 'Reza'='Other', 'SCM4'='Other', 'SCOR'='Other',
               'Tol('='Other', 'WIAG'='Other', 'WITC'='Other', 'NICE'='Other')
dat$modelbase <- sapply(dat$`Base IAM (if applicable)`, function(ss) ifelse(substring(ss, 1, 4) %in% names(labels), labels[[substring(ss, 1, 4)]], substring(ss, 1, 4)))
dat$modelbase[is.na(dat$modelbase)] <- "N/A"

ggplot(dat, aes(Year)) +
    geom_bar(aes(fill=modelbase)) + theme_bw() + scale_x_continuous("Publication year", expand=c(0, 0)) +
    scale_y_continuous("Number of estimates", expand=c(0, 0)) +
    scale_fill_manual("Model group:", breaks=c(unique(dat$modelbase)[!(unique(dat$modelbase) %in% c('Other', 'N/A'))], 'Other', 'N/A'),
                      values=c('#984ea3', '#377eb8', '#4daf4a', '#e41a1c', '#ff7f00', '#ffff33', '#cab2d6',
                               '#a65628', '#f781bf', '#999999', '#333333'))
