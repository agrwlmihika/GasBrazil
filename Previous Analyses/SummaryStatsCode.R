# Summary Statistics
library(dplyr)
load('Data/stations_sampleDF.RData')

# ------------------------------------------------------------
# Variables used in demand estimation
# ------------------------------------------------------------
sum_stats <- sample_demand %>%
  filter(prod == "gas") %>%
  select(
    p_sell.imp,          # price
    log.nstations_group, # log # stations in group
    pumps.imp,           # pumps
    tank.imp,            # tank size
    log.rent,            # log rent
    pop.area,            # population density
    sjt,                 # market share
    s0t                  # outside share
  ) %>%
  summarise(across(everything(), list(
    mean   = ~round(mean(.,   na.rm=TRUE), 4),
    median = ~round(median(., na.rm=TRUE), 4),
    sd     = ~round(sd(.,     na.rm=TRUE), 4),
    min    = ~round(min(.,    na.rm=TRUE), 4),
    max    = ~round(max(.,    na.rm=TRUE), 4)
  ))) %>%
  tidyr::pivot_longer(
    everything(),
    names_to  = c("variable", "stat"),
    names_sep = "_(?=[^_]+$)"
  ) %>%
  tidyr::pivot_wider(names_from = stat, values_from = value)

print(sum_stats)