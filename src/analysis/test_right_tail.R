# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    data.table, janitor, magrittr, fixest, kableExtra,
    broom, tidyverse, tidylog, modelsummary
)
options("tidylog.display" = NULL)
`%notin%` = Negate(`%in%`)

# load main dataset
source("src/data_cleaining_scripts/cleaning_master.R")

dat = dat |>
    clean_names() |>
    rename(damage_tipping = tipping_points2,
           climate_tipping = tipping_points) |>
    # when doing the max we need numeric values here and not NA
    mutate(carbon_cycle = as.numeric(carbon_cycle),
           carbon_cycle = ifelse(is.na(carbon_cycle), 0, 1),
           climate_model = as.numeric(climate_model),
           climate_model = ifelse(is.na(climate_model), 0, 1)) |>
    mutate(earth_system = climate_model | carbon_cycle) |>
    select(-carbon_cycle, -climate_model) |>
    relocate(id_number:climate_tipping, earth_system)

# get row ids for resampling
row_samples = fread("outputs/distribution_v2.csv")

# resample the dataset
resampled_df = dat |>
    slice(row_samples$row) |>
    select(central_value_per_ton_co2) |>
    drop_na()

#############################
#### COMPUTE MEAN EXCESS
#############################

# test for fat tails: plot mean excess function, positive slope = fat tail
mef_df = resampled_df$central_value_per_ton_co2
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
                        me = mean(mef_df[mef_df > x] - x, na.rm = T),
                        n = length(mef_df[mef_df > x])
                    )
                })

mef_df_alt = resampled_df$central_value_per_ton_co2[resampled_df$central_value_per_ton_co2 < 50000]
mef_df_alt = map_df(unique(mef_df_alt),
                    function(x) {
                        tibble(
                            x = x,
                            me = mean(mef_df_alt[mef_df_alt > x] - x, na.rm = T),
                            n = length(mef_df_alt[mef_df_alt > x])
                        )
                    })

#############################
#### TAIL INDEX ESTIMATES
#############################

## MEF(x) = sigma / (1-xi) + xi / (1-xi) * x
## MEF =  alpha + beta * x
## 1. Regress MEF on x
## 2. Invert coefficient beta = xi / (1-xi) --> xi = beta / (1 + beta)
## 3. Tail index \equiv alpha = 1 / xi = (1 + beta) / beta
## 4. All moments > alpha are infinite


results = list()
results[["(1)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .95)),
                         weights = mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .95)])
results[["(2)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .90)),
                         weights = mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .90)])
results[["(3)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .75)),
                         weights = mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .75)])
results[["(4)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > 0),
                         weights = mef_df_alt$n[mef_df_alt$x > 0])
results[["(5)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .95)),
                         weights = 1/mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .95)])
results[["(6)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .90)),
                         weights = 1/mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .90)])
results[["(7)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > quantile(x, .75)),
                         weights = 1/mef_df_alt$n[mef_df_alt$x > quantile(mef_df_alt$x, .75)])
results[["(8)"]] = feols(me ~ x,
                         data = mef_df_alt |> filter(x > 0),
                         weights = mef_df_alt$n[mef_df_alt$x > 0])

betas = lapply(results, function(x) coefficients(x)[2]) |>
    cbind() |>
    unlist()

gpd_shape = round(betas/(betas + 1), 2)

tail_index = round((betas + 1)/betas, 2)

rows <- tribble(
    ~term, ~"(1)", ~"(2)", ~"(3)", ~"(4)", ~"(5)", ~"(6)", ~"(7)", ~"(8)",
    "Minimum Threshold Percentile", "95", "90", "75", "0", "95", "90", "75", "0",
    "Observational Weights", "Num. Obs.", "Num. Obs.", "Num. Obs.", "Num. Obs.", "1/Num. Obs.", "1/Num. Obs.", "1/Num. Obs.", "1/Num. Obs."
)

rows = rbind(rows, c("GPD Shape Parameter", gpd_shape))
rows = rbind(rows, c("Estimated Tail Index", tail_index))

msummary(results,
         coef_map = c(
             "x" = "MEF Slope"
         ),
         gof_omit = "R2|AIC|BIC|Log.|Std.|FE:",
         fmt = "%.3f",
         stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
         add_rows = rows,
         title = "Estimates of the Mean Excess Function slope, the Generalized Pareto Distribution shape parameter, and tail index. \\label{tab:tail_index}",
         output = "latex"
) |>
    kable_styling() |>
    footnote(general = "Standard errors are robust to heteroskedasticity. All estimates are from a sample that excludes Nordhaus (2019). Columns 1-4 weight observations of the mean excess by the number of SCC observations used to compute it. Columns 5-8 weight observations with the inverse.", threeparttable = TRUE,
             fixed_small_size = T) |>
    save_kable(file = "outputs/tail_index.tex")
