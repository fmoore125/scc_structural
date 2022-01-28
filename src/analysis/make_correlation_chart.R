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

# read in data for correlations
# NAs are 0s so replace them
M <- dat %>%
  select(`Carbon Cycle`:`Alternative ethical approaches (not Discounted Utilitarianism)`) %>%
  as_tibble() %>%
  mutate(across(everything(), ~ as.numeric(.x))) %>%
  replace(is.na(.), 0) %>%
  cor()

# write the corr chart
png("outputs/corr_chart.png", width = 10, height = 10, units = "in", res = 500)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M,
  method = "color", col = col(100),
  type = "upper", order = "hclust",
  addCoef.col = "black", # Add coefficient of correlation
  tl.col = "black",
  # Combine with significance
  # hide correlation coefficient on the principal diagonal
  diag = FALSE
)

dev.off()
