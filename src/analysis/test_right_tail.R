# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table, janitor, magrittr, fixest, 
  broom, tidyverse, tidylog, eva
)
options("tidylog.display" = NULL)
`%notin%` = Negate(`%in%`)

# load main dataset
source("src/data_cleaining_scripts/cleaning_master.R")

# get row ids for resampling
row_samples = fread("outputs/distribution.csv")

# resample the dataset
resampled_df = dat %>% 
  slice(row_samples$row) %>% 
  select(`Central Value ($ per ton CO2)`) %>% 
  drop_na()

# test for fat tails: plot mean excess function, positive slope = fat tail
mef_df = resampled_df$`Central Value ($ per ton CO2)`

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

feols(
  me ~ x, 
  mef_df %>% 
    filter(x > quantile(resampled_df$`Central Value ($ per ton CO2)`, .5, na.rm = T)) %>% 
    filter(x < quantile(resampled_df$`Central Value ($ per ton CO2)`, .99, na.rm = T))
  )

mef_df %>% 
  filter(x > quantile(resampled_df$`Central Value ($ per ton CO2)`, .3, na.rm = T)) %>% 
  filter(x < quantile(resampled_df$`Central Value ($ per ton CO2)`, .95, na.rm = T)) %>% 
  mutate(bins = ntile(x, 100)) %>% 
  group_by(bins) %>% 
  summarise(me = mean(me, na.rm = T),
            x = mean(x, na.rm = T)) %>% 
  ggplot() +
  geom_point(aes(x = x, y = me)) +
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
    x = "Threshold Value",
    subtitle = "Mean Excess Function (MEF(x) = E[SCC - x | SCC > x])"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

