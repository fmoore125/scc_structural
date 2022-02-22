# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table, janitor, magrittr, fixest,
  broom, tidyverse, tidylog
)
options("tidylog.display" = NULL)
`%notin%` <- Negate(`%in%`)

# load main dataset
source("src/data_cleaining_scripts/cleaning_master.R")
source("src/analysis/find_distribution.R")


all.qs <- c(0,0.001,0.01, .025, .05, .1, .17, .25, .5, .75, .83, .9, .95, .975, .99,0.999, 1)
all.as.cols <- which(names(dat) == 'Min'):which(names(dat) == 'Max')

#start by generating distributions for each row
dists=list()
for (ii in 1:nrow(dat)) {
    print(ii)
    all.as <- t(dat[ii, all.as.cols])
    qs <- all.qs[!is.na(all.as)]
    as <- all.as[!is.na(all.as)]
    mu <- dat$`Central Value ($ per ton CO2)`[ii]
    if (is.na(mu) && length(qs) == 0) {
        next
    }
    
    dists[[ii]] <- generate.pdf(mu, qs, as, 1e6)
}

# select only the variables we need
dat <- dat %>%
  as_tibble() %>%
  clean_names() %>%
  rename(
    damage_tipping = tipping_points2,
    climate_tipping = tipping_points
  ) %>%
  select(
    scc_year:discountrate
  )

# for each observation, construct the difference between
# min and max values assuming a uniform distribution
dat <- dat %>%
  mutate(
    unc_minmax = (max - min) / (1 - 0),
    unc_999 = (x99_9th - x0_1th) / (.999 - .001),
    unc_990 = (x99th - x1th) / (.99 - .01),
    unc_975 = (x97_5th - x2_5th) / (.975 - .025),
    unc_950 = (x95th - x5th) / (.95 - .05),
    unc_900 = (x90th - x10th) / (.9 - .1),
    unc_830 = (x83rd - x17th) / (.83 - .17),
    unc_750 = (x75th - x25th) / (.75 - .25)
  ) %>%
  rowwise() %>%
  # this is proportional to the standard deviation (if we only had one set of percentiles per obs)
  # std dev of uniform = (max - min) / sqrt(12)
  mutate(range = mean(c_across(unc_minmax:unc_750), na.rm = T))

# variables on the RHS
vars <- names(dat)[45:58]

# ensure they are indicator variables
dat <- dat %>%
  mutate(across(tfp_growth:risk_aversion_ez_utility, ~ ifelse(is.na(.x), 0, 1)))

a <- feols(
  # uniform range
  log(range) ~ 
  # parametric uncertainty variables
  .[vars] 
  # controls that may be correlated with parametric uncertainty
  + as.numeric(scc_year) + discountrate + central_value_per_ton_co2,
  dat
) %>%
  summary()
