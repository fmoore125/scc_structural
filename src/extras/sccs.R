library(ggplot2)

iwg=read.csv("outputs/iwgruns.csv",row.names=1)
values <- iwg[5:nrow(iwg), which(iwg[1,] == 2020 & iwg[3,] == "2.50%" & iwg[4,] == "CO2")]
iwgdist <- as.numeric(as.character(unlist(values)))

df <- rbind(data.frame(region="New York", discount=0.02, year=2023, sccyear=2020, value=130, reference="https://dec.ny.gov/regulatory/guidance-and-policy-documents/climate-change-guidance-documents"),
            data.frame(region="US EPA", discount=0.02, year=2023, sccyear=2020, value=190, reference="https://www.epa.gov/system/files/documents/2023-12/epa_scghg_2023_report_final.pdf"),
            data.frame(region="US IWG", discount=0.025, year=2010, sccyear=2020, value=mean(iwgdist), reference="iwgruns.csv"),
            data.frame(region="Germany", discount=0.015, year=2019, sccyear=2020, value=777, reference="Fran's Calculations"), # Used 0 PRTP to match with https://media.rff.org/documents/Newell_Pizer_Prest_21-16.pdf
            data.frame(region="Canada", discount=0.02, year=2023, sccyear=2020, value=247 * 0.7978 * 105.366 / 110.185, reference="https://www.canada.ca/en/environment-climate-change/services/climate-change/science-research-data/social-cost-ghg.html"), # https://www.exchangerates.org.uk/CAD-USD-spot-exchange-rates-history-2021.html, https://fred.stlouisfed.org/series/GDPDEF
            data.frame(region="California", discount=0.03, year=2017, sccyear=2020, value=50 * 105.366 / 102.288, reference="https://www.gao.gov/assets/gao-20-254.pdf"),
            data.frame(region="Minnesota", discount=0.03, year=2018, sccyear=2020, value=45 * 105.366 / 102.288, reference="https://www.gao.gov/assets/gao-20-254.pdf"))

df2 <- rbind(subset(df, !(region %in% c("US EPA", "Canada"))),
             data.frame(region="US EPA & Canada", discount=0.02, year=2023, sccyear=2020, value=df$value[df$region == "US EPA"], reference=NA))

ggplot(df2, aes(x="Estimates", y=value, colour=region, shape=factor(discount))) +
    coord_flip() +
    geom_point() + theme_bw() +
    scale_y_log10("Central 2020 SCC ($2020 per ton CO2)", breaks=c(50, 70, 100, 200, 300, 500, 700)) + xlab(NULL) +
    scale_colour_discrete("Region:") + scale_shape_discrete("Near-term Discount Rate:") +
    theme(legend.position="bottom") +
    guides(shape=guide_legend(nrow=2))

