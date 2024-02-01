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

experiments <- c("best"="Synthetic SCC", "A_dice"="DICE", "B_epa"="EPA", "C_all"="All Structural Changes",
                 "E_1"="1% Discount Rate", "E_1.5"="1.5% Discount Rate", "E_2"="2% Discount Rate",
                 "E_2.5"="2.5% Discount Rate", "E_3"="3% Discount Rate") #, "E_5"="5% Discount Rate")

for (name in names(experiments)) {
    load(paste0("outputs/rf_experiments/RFD_", name, ".RData"))
    add.bigtbl.row(experiments[[name]], allsamp)
}

