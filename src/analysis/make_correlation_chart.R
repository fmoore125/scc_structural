# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table, janitor, magrittr, fixest,
  broom, tidyverse, tidylog, corrplot
)
options("tidylog.display" = NULL)
`%notin%` <- Negate(`%in%`)

# runs data cleaning to generate main dataframe as `dat`
source("src/data_cleaining_scripts/cleaning_master.R")

names_vec = c("PRTP", "EMUC", "Discount Rate", "Carbon Cycle", "Climate Model", "Climate Tipping", "Damage Tipping", "Persist/Growth Damage",
              "Epstein-Zin", "Ambiguity", "Limited Substitutability", "Inequality Aversion", "Learning", "Alt. Ethical Approaches")

# read in data for correlations
# NAs are 0s so replace them
M <- dat %>%
  select(PRTP, EMUC, discountrate, `Carbon Cycle`:`Alternative ethical approaches (not Discounted Utilitarianism)`) %>%
  as_tibble() %>%
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  replace(is.na(.), 0) %>%
  cor()

rownames(M) = names_vec
colnames(M) = names_vec

# write the corr chart
png("outputs/corr_chart.png", width = 15, height = 15, units = "in", res = 500)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M,
  method = "color", col = col(100),
  type = "lower", order = "hclust",
  addCoef.col = "black", # Add coefficient of correlation
  tl.col = "black",
  cl.cex = 1.5,
  tl.cex = 1,
  # Combine with significance
  # hide correlation coefficient on the principal diagonal
  diag = T
)

dev.off()

dat2 <- dat %>%
    select(`Carbon Cycle`:`Alternative ethical approaches (not Discounted Utilitarianism)`)

Ndf <- data.frame()
for (cc1 in 1:ncol(dat2))
    for (cc2 in cc1:ncol(dat2))
        Ndf <- rbind(Ndf, data.frame(col1=names_vec[3+cc1], col2=names_vec[3+cc2], count=sum(!is.na(dat2[, cc1]) & !is.na(dat2[, cc2]))))
Ndf$col1 <- factor(Ndf$col1, names_vec[4:length(names_vec)])
Ndf$col2 <- factor(Ndf$col2, rev(names_vec[4:length(names_vec)]))

ggplot(Ndf, aes(col1, col2, fill=count)) +
    geom_raster() + geom_label(aes(label=count)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_distiller("Estimates", palette='BuPu', direction=1, trans='log10') +
    xlab(NULL) + ylab(NULL)
ggsave("outputs/num_chart.png", width=6.5, height=5)
