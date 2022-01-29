# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table, janitor, magrittr, fixest, 
  broom, tidyverse, tidylog, mddelsummary
)
options("tidylog.display" = NULL)
`%notin%` = Negate(`%in%`)

# load main dataset
source("src/data_cleaining_scripts/cleaning_master.R")

dat = dat %>% 
  clean_names() %>% 
  rename(damage_tipping = tipping_points2,
         climate_tipping = tipping_points)

# get row ids for resampling
row_samples = fread("outputs/distribution.csv")

# resample the dataset
resampled_df = dat %>% 
  slice(row_samples$row) %>% 
  select(central_value_per_ton_co2) %>% 
  drop_na()

#############################
#### COMPUTE MEAN EXCESS
#############################

# test for fat tails: plot mean excess function, positive slope = fat tail
mef_df = dat$central_value_per_ton_co2
mef_df = mef_df[!is.na(mef_df)]

# compute mean excess function for each unique SCC value
# MEF(x) = E[SCC - x | SCC > x]
# slope > 0 -> shape parameter of GPD > 1 and fat tails
# Why? Expectation is increasing faster than the threshold if slope > 1
# Conte and Kelly (2018) pg 685, Ghosh (2010)
# Need a lower and upper threshold, pick multiple values to test
mef_df = map_df(unique(mef_df),
       function(x) {
         tibble(
           x = x,
           me = mean(mef_df[mef_df > x] - x, na.rm = T)
           )
       })


# set up etable
dict = c("(Intercept)" = "Constant", 
         x = "Threshold (x)",
         me = "Mean Excess",
         note1 = dsb("*Notes*: This is a note that illustrates how to access notes ",
                     "from the dictionary."))

setFixest_etable(digits = 1, 
                 keep = "x",
                 dict = dict,
                 fitstat = c('n'),
                 style.tex = style.tex("aer", model.format = "(i)")
                 )

# run MEF regressions
results = list()
results[[1]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .3, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .7, na.rm = T))
  )

results[[2]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .3, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .9, na.rm = T))
)

results[[3]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .7, na.rm = T))
)

results[[4]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .9, na.rm = T))
)

results[[5]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .99, na.rm = T))
)

results[[6]] = feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .7, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .99, na.rm = T))
)

etable(results, 
       se.row = F,
       extralines = list(
         "Lower Percentile Threshold" = list(".3" = 2, ".5" = 3, ".7"),
         "Upper Percentile Threshold" = list(".7", ".9", ".7", ".9", ".99", ".99")
       ))

# make binscattered MEF plot
mef_df %>% 
  filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) %>% 
  filter(x < quantile(resampled_df$central_value_per_ton_co2, .99, na.rm = T)) %>% 
  mutate(bins = ntile(x, 100)) %>% 
  group_by(bins) %>% 
  summarise(me = mean(me, na.rm = T),
            x = mean(x, na.rm = T)) %>% 
  ggplot() +
  geom_point(aes(x = x, y = me), size = 3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    title = element_text(size = 24),
    axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24),
    panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black"), axis.ticks = element_line()
  ) +
  labs(
    y = "",
    x = "Threshold SCC (x)",
    subtitle = "Mean Excess Function (MEF(x) = E[SCC - x | SCC > x])"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1700)) +
  scale_y_continuous(breaks = scales::pretty_breaks())

ggsave("outputs/mean_excess_function.png", height = 10, width = 10)


#############################
#### Pr(tail | characteristic)
#############################

tail_prob = .99

#### Early period

tail_df = dat %>% 
  slice(row_samples$row) %>% 
  filter(as.numeric(scc_year) < 2025 & as.numeric(scc_year) > 0) %>% 
  mutate(in_tail = central_value_per_ton_co2 > quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) %>% 
  mutate(across(carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism, ~as.numeric(.x)))

tail_df[is.na(tail_df)] = 0
tail_df$scc_year = as.numeric(tail_df$scc_year)

feols(
  as.formula(paste("in_tail ~ discountrate + scc_year + ", paste0(names(dat)[26:36], collapse = "+")))
             , tail_df) %>% 
  tidy() %>% 
  arrange(desc(estimate))

#### Mid period

tail_df = dat %>% 
  slice(row_samples$row) %>% 
  filter(as.numeric(scc_year) < 2070 & as.numeric(scc_year) > 2030) %>% 
  mutate(in_tail = central_value_per_ton_co2 > quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) %>% 
  mutate(across(carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism, ~as.numeric(.x)))

tail_df[is.na(tail_df)] = 0
tail_df$scc_year = as.numeric(tail_df$scc_year)

feols(
  as.formula(paste("in_tail ~ discountrate + scc_year + ", paste0(names(dat)[26:36], collapse = "+")))
  , tail_df) %>% 
  tidy() %>% 
  arrange(desc(estimate))

#### Late period

tail_df = dat %>% 
  slice(row_samples$row) %>% 
  filter(as.numeric(scc_year) < 2120 & as.numeric(scc_year) > 2080) %>% 
  mutate(in_tail = central_value_per_ton_co2 > quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) %>% 
  mutate(across(carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism, ~as.numeric(.x)))

tail_df[is.na(tail_df)] = 0
tail_df$scc_year = as.numeric(tail_df$scc_year)

feols(
  as.formula(paste("in_tail ~ discountrate + scc_year + ", paste0(names(dat)[26:36], collapse = "+")))
  , tail_df) %>% 
  tidy() %>% 
  arrange(desc(estimate))