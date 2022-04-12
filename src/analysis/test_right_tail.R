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
           climate_tipping = tipping_points)

# get row ids for resampling
row_samples = fread("outputs/distribution.csv")

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
                        me = mean(mef_df[mef_df > x] - x, na.rm = T)
                    )
                })

mef_df_alt = resampled_df$central_value_per_ton_co2[resampled_df$central_value_per_ton_co2 < 50000]
mef_df_alt = map_df(unique(mef_df_alt),
                    function(x) {
                        tibble(
                            x = x,
                            me = mean(mef_df_alt[mef_df_alt > x] - x, na.rm = T)
                        )
                    })

#############################
#### MEAN EXCESS PLOT
#############################

# make binscattered MEF plot
mef_bins = mef_df |> 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) |> 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .99, na.rm = T)) |> 
    mutate(bins = ntile(x, 100)) |> 
    group_by(bins) |> 
    summarise(me = mean(me, na.rm = T),
              x = mean(x, na.rm = T))

mef_alt_bins = mef_df_alt |> 
    filter(x > quantile(resampled_df$central_value_per_ton_co2, .5, na.rm = T)) |> 
    filter(x < quantile(resampled_df$central_value_per_ton_co2, .99, na.rm = T)) |> 
    mutate(bins = ntile(x, 100)) |> 
    group_by(bins) |> 
    summarise(me = mean(me, na.rm = T),
              x = mean(x, na.rm = T))

ggplot() +
    geom_point(data = mef_bins, aes(x = x, y = me), size = 3, shape = 16) +
    geom_point(data = mef_alt_bins, aes(x = x, y = me), size = 3, shape = 8) +
    annotate(
        geom = "text", x = 900, y = 3700, label = "Full Sample", hjust = 0,
        size = 7
    ) +
    annotate(
        geom = "text", x = 900, y = 1000, label = "Exclude Nordhaus (2019)", hjust = 0,
        size = 7
    ) +
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

tail_prob = .95

#### Early period

early_df = dat |> 
    rename_with(~str_c("struct_", .), .cols = carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism) |> 
    rename_with(~str_c("uncert_", .), .cols = tfp_growth:risk_aversion_ez_utility) |> 
    select(-starts_with("x"), -min, -max) |> 
    select(central_value_per_ton_co2, discountrate, scc_year, struct_carbon_cycle:uncert_risk_aversion_ez_utility, id_number) |> 
    slice(row_samples$row) |>  
    filter(as.numeric(scc_year) < 2030 & as.numeric(scc_year) >= 2010) |> 
    mutate(in_tail = central_value_per_ton_co2 >= quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~as.numeric(.x))) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~replace(., is.na(.), 0))) 

early_df[is.na(early_df)] = 0

# variables to compute odds ratio
variables_vec = c(names(early_df)[grep(pattern = "uncert_", names(early_df))], names(early_df)[grep(pattern = "struct_", names(early_df))])

odds_df = early_df |> 
    pivot_longer(
        cols = all_of(variables_vec),
        names_to = "parameter",
        values_to = "on"
    ) |> 
    group_by(
        parameter, in_tail, on
    ) |> 
    summarise(n = n()) 

odds_early_df = expand.grid(
    parameter = unique(odds_df$parameter),
    in_tail = c(0, 1),
    on = c(0, 1)
) |> 
    left_join(odds_df) |> 
    mutate(across(n, ~replace(., is.na(.), 0))) |>
    pivot_wider(
        names_from = c("in_tail", "on"),
        values_from = n
    ) |> 
    # (in tail & change / not in tail & change) / (in tail & not change / not in tail & not change)
    mutate(odds_ratio_early = `1_1` / `0_1` / (`1_0` / `0_0` )) |> 
    select(parameter, odds_ratio_early) |> 
    arrange(desc(odds_ratio_early)) |> 
    mutate(across(odds_ratio_late, ~replace(., is.na(.), 0))) 

#### Mid period

mid_df = dat |> 
    rename_with(~str_c("struct_", .), .cols = carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism) |> 
    rename_with(~str_c("uncert_", .), .cols = tfp_growth:risk_aversion_ez_utility) |> 
    select(-starts_with("x"), -min, -max) |> 
    select(central_value_per_ton_co2, discountrate, scc_year, struct_carbon_cycle:uncert_risk_aversion_ez_utility, id_number) |> 
    slice(row_samples$row) |>  
    filter(as.numeric(scc_year) < 2070 & as.numeric(scc_year) >= 2030) |> 
    mutate(in_tail = central_value_per_ton_co2 >= quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~as.numeric(.x))) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~replace(., is.na(.), 0))) 

mid_df[is.na(mid_df)] = 0

odds_df = mid_df |> 
    pivot_longer(
        cols = all_of(variables_vec),
        names_to = "parameter",
        values_to = "on"
    ) |> 
    group_by(
        parameter, in_tail, on
    ) |> 
    summarise(n = n()) 

odds_mid_df = expand.grid(
    parameter = unique(odds_df$parameter),
    in_tail = c(0, 1),
    on = c(0, 1)
) |> 
    left_join(odds_df) |> 
    mutate(across(n, ~replace(., is.na(.), 0))) |>
    pivot_wider(
        names_from = c("in_tail", "on"),
        values_from = n
    ) |> 
    # (in tail & change / not in tail & change) / (in tail & not change / not in tail & not change)
    mutate(odds_ratio_mid = `1_1` / `0_1` / (`1_0` / `0_0` )) |> 
    select(parameter, odds_ratio_mid) |> 
    arrange(desc(odds_ratio_mid)) |> 
    mutate(across(odds_ratio_late, ~replace(., is.na(.), 0))) 

#### Late period

late_df = dat |> 
    rename_with(~str_c("struct_", .), .cols = carbon_cycle:alternative_ethical_approaches_not_discounted_utilitarianism) |> 
    rename_with(~str_c("uncert_", .), .cols = tfp_growth:risk_aversion_ez_utility) |> 
    select(-starts_with("x"), -min, -max) |> 
    select(central_value_per_ton_co2, discountrate, scc_year, struct_carbon_cycle:uncert_risk_aversion_ez_utility, id_number) |> 
    slice(row_samples$row) |>  
    filter(as.numeric(scc_year) <= 2100 & as.numeric(scc_year) >= 2070) |> 
    mutate(in_tail = central_value_per_ton_co2 >= quantile(central_value_per_ton_co2, tail_prob, na.rm = T)) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~as.numeric(.x))) |> 
    mutate(across(struct_carbon_cycle:uncert_risk_aversion_ez_utility, ~replace(., is.na(.), 0))) 

late_df[is.na(late_df)] = 0


odds_df = late_df |> 
    pivot_longer(
        cols = all_of(variables_vec),
        names_to = "parameter",
        values_to = "on"
    ) |> 
    group_by(
        parameter, in_tail, on
    ) |> 
    summarise(n = n()) 

odds_late_df = expand.grid(
    parameter = unique(odds_df$parameter),
    in_tail = c(0, 1),
    on = c(0, 1)
) |> 
    left_join(odds_df) |> 
    mutate(across(n, ~replace(., is.na(.), 0))) |>
    pivot_wider(
        names_from = c("in_tail", "on"),
        values_from = n
    ) |> 
    # (in tail & change / not in tail & change) / (in tail & not change / not in tail & not change)
    mutate(odds_ratio_late = `1_1` / `0_1` / (`1_0` / `0_0` )) |> 
    select(parameter, odds_ratio_late) |> 
    arrange(desc(odds_ratio_late)) |> 
    mutate(across(odds_ratio_late, ~replace(., is.na(.), 0))) 

odds_ratio_df =
    odds_early_df |> 
    left_join(odds_mid_df) |> 
    left_join(odds_late_df) |> 
    arrange(desc(odds_ratio_early))

names(odds_ratio_df) = c("Parameter", "2010-2030", "2030-2070", "2070-2100")

odds_ratio_df = odds_ratio_df |> 
    mutate(Class = 
               case_when(
                   substr(Parameter, 1, 6) == "struct" ~ "Structural Change",
                   TRUE ~ "Parametric Uncertainty"
               ),
           Parameter = substr(Parameter, 8, 50),
           Parameter = gsub("_", " ", Parameter),
           Parameter = str_to_title(Parameter)) |> 
    relocate(Class, Parameter) |> 
    arrange(Class, Parameter)

n_struct = sum(odds_ratio_df$Class == "Structural Change")
n_uncert = sum(odds_ratio_df$Class == "Parametric Uncertainty")

names(odds_ratio_df)[2] = " "

odds_table <- kbl(odds_ratio_df |> select(-Class), 
                  format = "latex",
                  caption = "Odds ratio of being in the top 5 percent of SCC values.",
                  label = "odds_table",
                  booktabs = T, 
                  digits = 1, 
                  escape = F,
                  # col.names = linebreak(c("", "Without\nMarket Adaptation", "With\nMarket Adaptation", "Without\nMarket Adaptation", "With\nMarket Adaptation"), align = "c"), 
                  align = c("l", "c", "c", "c")) %>%
    pack_rows(index = c("Parametric Uncertainty" = n_uncert, "Structural Change" = n_struct)) %>%
    add_header_above(c(" " = 1, "Odds Ratios for Being in Right Tail of Distribution" = 3)) %>%
    kable_styling() %>% 
    footnote(general = "Footnote here.", 
             threeparttable = T,
             fixed_small_size = T)

fileConn <- file("outputs/odds_ratio.tex")
writeLines(odds_table, fileConn)
close(fileConn)
