# ============================================================
# 1. Summary Statistics
# ============================================================
load('Data/stations_sampleDF.RData')
load('/proj/econ/duartema/GasData/Database/Frota.RData')
library(dplyr)

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
sum_stats <- sample_demand %>%
  select(p_sell.imp, log.nstations_group, pumps.imp, tank.imp,
         log.rent, pop.area, sjt, s0t) %>%
  summarise(across(everything(), list(
    mean   = ~round(mean(.,   na.rm=TRUE), 4),
    median = ~round(median(., na.rm=TRUE), 4),
    sd     = ~round(sd(.,     na.rm=TRUE), 4),
    min    = ~round(min(.,    na.rm=TRUE), 4),
    max    = ~round(max(.,    na.rm=TRUE), 4)
  ))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = c("variable", "stat"),
                      names_sep = "_(?=[^_]+$)") %>%
  tidyr::pivot_wider(names_from=stat, values_from=value)

cat("=== Table 1: Summary Statistics ===\n")
print(sum_stats)

# Save
sink('/proj/econ/duartema/GasBrazil/Final/Output/table1_summary_stats.txt')
print(sum_stats)
sink()
cat("Saved to Output/table1_summary_stats.txt\n")