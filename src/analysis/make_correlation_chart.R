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
