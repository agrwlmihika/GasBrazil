
# ============================================================
# 1. Summary Statistics
# ============================================================
load('Data/stations_sampleDF.RData')
load('/proj/econ/duartema/GasData/Database/Frota.RData')
library(dplyr)
library(stargazer)

# Build fleet market size
fleet_monthly <- Frota %>%
  filter(state == "DF") %>%
  group_by(month) %>%
  summarise(fleet_t = sum(AUTOMOVEL, na.rm=TRUE), .groups="drop") %>%
  mutate(market_size_fleet = fleet_t * 100)

# Build sample with fleet shares
sample_demand <- sample %>%
  filter(prod == "gas") %>%
  left_join(fleet_monthly, by="month") %>%
  group_by(month) %>%
  mutate(
    sjt       = sale.imp / market_size_fleet,
    s0t       = 1 - sum(sjt, na.rm=TRUE),
    log_ratio = log(sjt / s0t)
  ) %>%
  ungroup() %>%
  filter(is.finite(log_ratio), s0t > 0, sjt > 0)

# Summary stats table
sum_data <- sample_demand %>%
  select(
    p_sell.imp, log.nstations_group, pumps.imp,
    tank.imp, log.rent, pop.area, sjt, s0t
  ) %>%
  as.data.frame()   # stargazer requires a plain data frame, not tibble

# Print to console (useful for checking)
stargazer(sum_data, type = "text")

# Output LaTeX file
stargazer(sum_data,
          type    = "latex",
          title   = "Summary Statistics",
          label   = "tab:summary",
          digits  = 4,
          out     = "/proj/econ/duartema/GasBrazil/Final/Output/table1_summary_stats.tex",
          # Rename variables for cleaner labels in the table
          covariate.labels = c(
            "Price (sell)",
            "Log N Stations (group)",
            "Pumps",
            "Tank size",
            "Log rent",
            "Pop. density",
            "$s_{jt}$",
            "$s_{0t}$"
          )
)
