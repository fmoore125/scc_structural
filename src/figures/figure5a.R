## setwd("~/research/scciams/scc_structural")

library(ggplot2)
library(dplyr)

do.pool.distdisc <- T
do.plot.log <- F

quants <- c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)

## Load the RF result
rfdist_dir <- "outputs/Structural SCC RF Experiments"
load(file.path(rfdist_dir, "RFD_best.Rdata"))
rfdist <- data.frame(draw=allsamp) %>% summarize(lowest=quantile(draw, quants[1]), min=quantile(draw, quants[2]),
                                                 lower=quantile(draw, quants[3]), middle=quantile(draw, quants[4]),
                                                 upper=quantile(draw, quants[5]), max=quantile(draw, quants[6]),
                                                 highest=quantile(draw, quants[7]), mu=mean(draw), discount='Pooled')

## EPA 2023 results
eparuns <- list.files("data/epa_scc", full.names=T) #2020 distribution of co2 for 3 damage functions
epa1 <- read.csv(eparuns[1])[, c('scghg', 'discount_rate')]
epa2 <- subset(read.csv(eparuns[2]), sector=="total")[, c('scghg', 'discount_rate')]
epa3 <- read.csv(eparuns[3])[, c('scghg', 'discount_rate')]
if (do.pool.distdisc) {
    epadist <- rbind(epa1, epa2, epa3) %>% summarize(lowest=quantile(scghg, quants[1]), min=quantile(scghg, quants[2]),
                                                     lower=quantile(scghg, quants[3]), middle=quantile(scghg, quants[4]),
                                                     upper=quantile(scghg, quants[5]), max=quantile(scghg, quants[6]),
                                                     highest=quantile(scghg, quants[7]), mu=mean(scghg))
    epadist$discount <- 'Pooled'
} else {
    epadist <- rbind(epa1, epa2, epa3) %>% group_by(discount_rate) %>% summarize(lowest=quantile(scghg, quants[1]), min=quantile(scghg, quants[2]),
                                                                                 lower=quantile(scghg, quants[3]), middle=quantile(scghg, quants[4]),
                                                                                 upper=quantile(scghg, quants[5]), max=quantile(scghg, quants[6]),
                                                                                 highest=quantile(scghg, quants[7]), mu=mean(scghg))
    epadist$discount <- paste0(substring(epadist$discount_rate, 1, 3), "%")
}

## IWG results
iwg <- read.csv("data/iwgruns.csv", row.names=1)
iwgvals <- data.frame()
for (rate in unique(unlist(iwg[3,]))) {
    if (is.na(rate))
        next
    values <- iwg[5:nrow(iwg), which(iwg[1,] == 2020 & iwg[3,] == rate & iwg[4,] == "CO2")]
    iwgvals <- rbind(iwgvals, data.frame(discount_rate=rate, scghg=as.numeric(as.character(unlist(values)))))
}
if (do.pool.distdisc) {
    iwgdist <- iwgvals %>% summarize(lowest=quantile(scghg, quants[1]), min=quantile(scghg, quants[2]),
                                     lower=quantile(scghg, quants[3]), middle=quantile(scghg, quants[4]),
                                     upper=quantile(scghg, quants[5]), max=quantile(scghg, quants[6]),
                                     highest=quantile(scghg, quants[7]), mu=mean(scghg))
    iwgdist$discount <- "Pooled"
} else {
    iwgdist <- iwgvals %>% group_by(discount_rate) %>% summarize(lowest=quantile(scghg, quants[1]), min=quantile(scghg, quants[2]),
                                                                 lower=quantile(scghg, quants[3]), middle=quantile(scghg, quants[4]),
                                                                 upper=quantile(scghg, quants[5]), max=quantile(scghg, quants[6]),
                                                                 highest=quantile(scghg, quants[7]), mu=mean(scghg))
    iwgdist$discount <- paste0(substring(iwgdist$discount_rate, 1, 3), "%")
}

## Other official SCCs
df <- rbind(data.frame(region="New York", discount=c(0.03, 0.02, 0.01), year=2023, sccyear=2020, value=c(53, 130, 420), reference="https://dec.ny.gov/regulatory/guidance-and-policy-documents/climate-change-guidance-documents"), # 0, 2200
            ## data.frame(region="US EPA", discount=0.02, year=2023, sccyear=2020, value=190, reference="https://www.epa.gov/system/files/documents/2023-12/epa_scghg_2023_report_final.pdf"),
            ## data.frame(region="US IWG", discount=0.025, year=2010, sccyear=2020, value=mean(iwgdist), reference="iwgruns.csv"),
            data.frame(region="Germany", discount=c(0.015, 0.025), year=2019, sccyear=2020, value=c(777, 223), reference="Fran's Calculations"), # Used 0 PRTP to match with https://media.rff.org/documents/Newell_Pizer_Prest_21-16.pdf
            data.frame(region="Canada", discount=0.02, year=2023, sccyear=2020, value=247 * 0.7978 * 105.366 / 110.185, reference="https://www.canada.ca/en/environment-climate-change/services/climate-change/science-research-data/social-cost-ghg.html"), # https://www.exchangerates.org.uk/CAD-USD-spot-exchange-rates-history-2021.html, https://fred.stlouisfed.org/series/GDPDEF
            data.frame(region="California", discount=c(0.05, 0.03, 0.025, 0.03), year=2017, sccyear=2020, value=c(14, 50, 74, 147) * 105.366 / 102.288, reference="https://www.gao.gov/assets/gao-20-254.pdf"),
            data.frame(region="Minnesota", discount=c(0.05, 0.03), year=2018, sccyear=2020, value=c(10, 45) * 105.366 / 102.288, reference="https://www.gao.gov/assets/gao-20-254.pdf"))
df$panel <- "US States"
df$panel[df$region %in% c('Canada', 'Germany')] <- "Other"
df$panel <- factor(df$panel, levels=c('Syn.', 'Distrib.', 'US States', 'Other'))
df$discount <- paste0(df$discount * 100, "%")

## df2 <- rbind(subset(df, !(region %in% c("US EPA", "Canada"))),
##              data.frame(region="US EPA & Canada", discount=0.02, year=2023, sccyear=2020, value=df$value[df$region == "US EPA"], reference=NA))

## Produce combined figure

dfd <- rbind(cbind(region="US EPA", year=2023, sccyear=2020, epadist, panel="Distrib."),
             cbind(region="US IWG", year=2010, sccyear=2020, iwgdist, panel="Distrib."),
             cbind(region="This paper", year=2024, sccyear=2020, rfdist, panel="Syn."))
dfd$panel <- factor(dfd$panel, levels=c('Syn.', 'Distrib.', 'US States', 'Other'))

if (do.plot.log) {
    ggplot(df, aes(region, shape=factor(discount), group=paste(region, discount))) +
        facet_grid(panel ~ ., scales="free_y", space="free_y") +
        coord_flip() +
        geom_boxplot(data=dfd, aes(min=pmax(1, min), lower=pmax(1, lower), middle=middle, upper=upper, max=max), stat="identity", position=position_dodge(width=0.75), width=.75, show.legend=F) + # width=c(0.4,0.2,0.2), lwd=0.4,
        geom_point(data=dfd, aes(y=mu), position=position_dodge(width=0.75)) +
        geom_point(aes(y=value)) +
        theme_bw() +
        scale_y_log10("Central 2020 SCC ($2020 per ton CO2)", expand=expansion(mult=c(0, .05))) + xlab(NULL) +
        scale_colour_discrete("Region:") +
        scale_shape_manual("Near-term Discount Rate:", breaks=c("Pooled", 0.01, 0.015, 0.02, 0.025, 0.03, 0.05), values=c(1, 15, 16, 17, 3, 4, 6)) +
        theme(legend.position="bottom") +
        guides(shape=guide_legend(nrow=2))
} else {
    gp <- ggplot(df, aes(region, shape=factor(discount), group=paste(region, discount), colour=region == 'This paper')) +
        facet_grid(panel ~ ., scales="free_y", space="free_y") +
        coord_flip() +
        geom_boxplot(data=dfd, aes(min=min, lower=lower, middle=middle, upper=upper, max=max), stat="identity", position=position_dodge(width=0.75), width=.75, show.legend=F) + # width=c(0.4,0.2,0.2), lwd=0.4,
        geom_point(data=dfd, aes(y=mu), position=position_dodge(width=0.75)) +
        geom_segment(data=dfd, aes(xend=region, y=lowest, yend=min), position=position_dodge(width=0.75), lty=2, lwd=0.4) +
        geom_segment(data=dfd, aes(xend=region, y=highest, yend=max), position=position_dodge(width=0.75), lty=2, lwd=0.4) +
        geom_point(aes(y=value)) +
        theme_bw() +
        scale_y_continuous("Central 2020 SCC ($2020 per ton CO2)") +
        scale_x_discrete(NULL) +
        scale_colour_manual(NULL, breaks=c(F, T), values=c('#000000', '#d95f02')) +
        scale_shape_manual("Near-term Discount Rate:", breaks=c("Pooled", "1%", "1.5%", "2%", "2.5%", "3%", "5%"), values=c(1, 15, 16, 17, 3, 4, 10)) +
        theme(legend.position="bottom") +
        guides(shape=guide_legend(nrow=2), colour='none')
    ggsave("figures/figure5a.pdf", width=8, height=4)
}
