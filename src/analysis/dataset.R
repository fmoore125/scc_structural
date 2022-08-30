## setwd("~/research/scciams/scc_structural")

source("src/data_cleaining_scripts/cleaning_master.R")

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
