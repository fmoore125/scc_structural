source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/all_scc_lib.R")
source("src/analysis/randomforest_dists_lib.R")

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
