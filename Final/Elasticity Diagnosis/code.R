library(dplyr)

load('/proj/econ/duartema/GasBrazil/CodeStruct/Data/stations_sampleDF.RData')
load('/proj/econ/duartema/GasData/Database/Frota.RData')

fleet_monthly <- Frota %>%
  filter(state == "DF") %>%
  group_by(month) %>%
  summarise(fleet_t = sum(AUTOMOVEL, na.rm=TRUE), .groups="drop") %>%
  mutate(market_size_fleet = fleet_t * 100)

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

elas_summary_iv1 <- readRDS("/proj/econ/duartema/GasBrazil/Final/Elasticity Diagnosis/elas_summary_iv1.rds")

elas_outliers_iv1 <- elas_summary_iv1 %>%
  filter(own_elas < -20 | own_elas > -5)

elas_outliers_iv1_info <- elas_outliers_iv1 %>%
  left_join(sample %>% filter(prod == "gas") %>%
              select(cnpj, prod, month, sale, neigh, dist, tank, pumps, group, downtown, unbranded, n_oppo1km,
                     name, type, lat, lon, nstations_group, p_sell, std.pdev, tot.rev),
            by = c("cnpj", "month")) %>%
  left_join(sample_demand %>% select(cnpj, month, sjt),
            by = c("cnpj", "month"))

# ============================================================
# FLAGGED STATION FREQUENCY ANALYSIS
# ============================================================

# step 1: compute frequency stats from the clean outliers (no join duplicates)
flagged_frequency <- elas_outliers_iv1 %>%
  group_by(cnpj) %>%
  summarise(
    n_months_flagged  = n(),
    median_own_elas   = median(own_elas, na.rm=TRUE),
    min_own_elas      = min(own_elas,    na.rm=TRUE),
    max_own_elas      = max(own_elas,    na.rm=TRUE)
  ) %>%
  arrange(desc(n_months_flagged)) %>%
  ungroup()

# step 2: join time-varying characteristics (median across flagged months)
varying_chars <- elas_outliers_iv1_info %>%
  group_by(cnpj) %>%
  summarise(
    median_sale       = median(sale,            na.rm=TRUE),
    median_sjt        = median(sjt,             na.rm=TRUE),
    median_n_oppo1km  = median(n_oppo1km,       na.rm=TRUE),
    median_nstations  = median(nstations_group, na.rm=TRUE),
    median_p_sell     = median(p_sell,          na.rm=TRUE),
    median_std_pdev   = median(std.pdev,        na.rm=TRUE),
    median_tot_rev    = median(tot.rev,         na.rm=TRUE)
  ) %>% ungroup()

# step 3: join time-invariant characteristics (one row per station)
invariant_chars <- sample %>%
  filter(prod == "gas") %>%
  select(cnpj, neigh, group, prod, downtown, unbranded, cascol, dist, tank, pumps) %>%
  distinct(cnpj, .keep_all=TRUE)

# step 4: combine everything
flagged_frequency <- flagged_frequency %>%
  left_join(invariant_chars, by = "cnpj") %>%
  left_join(varying_chars,   by = "cnpj")

View(flagged_frequency)

# ============================================================
# SUMMARY STATS FOR FLAGGED STATIONS
# ============================================================

# 1. distribution of how often stations get flagged
cat("=== Distribution of Flagged Months ===\n")
print(summary(flagged_frequency$n_months_flagged))
print(table(flagged_frequency$n_months_flagged))

# 2. split into two groups: too inelastic (> -5) vs too elastic (< -20)
elas_outliers_iv1 <- elas_outliers_iv1 %>%
  mutate(flag_type = ifelse(own_elas > -5, "too_inelastic", "too_elastic"))

cat("\n=== Flag Type Breakdown ===\n")
print(table(elas_outliers_iv1$flag_type))

# 3. compare flagged vs non-flagged stations on key characteristics
elas_summary_iv1 <- elas_summary_iv1 %>%
  mutate(flagged = ifelse(own_elas < -20 | own_elas > -5, 1, 0))

# join characteristics to full summary (gas only)
full_info <- elas_summary_iv1 %>%
  left_join(sample %>% filter(prod == "gas") %>%
              select(cnpj, month, p_sell, n_oppo1km,
                     nstations_group, tank, pumps, unbranded, downtown,
                     cascol, dist, tot.rev, std.pdev),
            by = c("cnpj", "month")) %>%
  left_join(sample_demand %>% select(cnpj, month, sjt),
            by = c("cnpj", "month"))

# 4. compare flagged vs non-flagged
comparison <- full_info %>%
  group_by(flagged) %>%
  summarise(
    n_obs                  = n(),
    median_own_elas        = median(own_elas,        na.rm=TRUE),
    mean_own_elas        = mean(own_elas,        na.rm=TRUE),
    median_sjt             = median(sjt,             na.rm=TRUE),
    mean_sjt             = mean(sjt,             na.rm=TRUE),
    median_p_sell          = median(p_sell,          na.rm=TRUE),
    mean_p_sell            = mean(p_sell,       na.rm=TRUE),
    median_n_oppo          = median(n_oppo1km,       na.rm=TRUE),
    mean_n_oppo          = mean(n_oppo1km,       na.rm=TRUE),
    median_tank            = median(tank,            na.rm=TRUE),
    median_pumps           = median(pumps,           na.rm=TRUE),
    pct_downtown           = mean(downtown,          na.rm=TRUE),
    pct_cascol             = mean(cascol == 1,       na.rm=TRUE),
    pct_raizen             = mean(dist == "RAIZEN COMBUSTÍVEIS S.A.",          na.rm=TRUE),
    pct_petrobras          = mean(dist == "PETROBRAS DISTRIBUIDORA S.A.",      na.rm=TRUE),
    pct_ipiranga           = mean(dist == "IPIRANGA PRODUTOS DE PETRÓLEO S.A", na.rm=TRUE),
    pct_bandeira_branca    = mean(dist == "BANDEIRA BRANCA",                   na.rm=TRUE),
    median_tot_rev         = median(tot.rev,         na.rm=TRUE),
    median_std_pdev        = median(std.pdev,        na.rm=TRUE)
  )

print(comparison)

# 5. three-way split: too elastic / too inelastic / normal
full_info <- full_info %>%
  mutate(flag_type = case_when(
    own_elas < -20 ~ "too_elastic",
    own_elas > -5  ~ "too_inelastic",
    TRUE           ~ "normal"
  ))

comparison_3way <- full_info %>%
  group_by(flag_type) %>%
  summarise(
    n_obs                  = n(),
    median_own_elas        = median(own_elas,        na.rm=TRUE),
    mean_own_elas        = mean(own_elas,        na.rm=TRUE),
    median_sjt             = median(sjt,             na.rm=TRUE),
    mean_sjt             = mean(sjt,             na.rm=TRUE),
    median_p_sell          = median(p_sell,          na.rm=TRUE),
    mean_p_sell            = mean(p_sell,       na.rm=TRUE),
    median_n_oppo          = median(n_oppo1km,       na.rm=TRUE),
    mean_n_oppo          = mean(n_oppo1km,       na.rm=TRUE),
    median_tank            = median(tank,            na.rm=TRUE),
    median_pumps           = median(pumps,           na.rm=TRUE),
    pct_downtown           = mean(downtown,          na.rm=TRUE),
    pct_cascol             = mean(cascol == 1,       na.rm=TRUE),
    pct_raizen             = mean(dist == "RAIZEN COMBUSTÍVEIS S.A.",          na.rm=TRUE),
    pct_petrobras          = mean(dist == "PETROBRAS DISTRIBUIDORA S.A.",      na.rm=TRUE),
    pct_ipiranga           = mean(dist == "IPIRANGA PRODUTOS DE PETRÓLEO S.A", na.rm=TRUE),
    pct_bandeira_branca    = mean(dist == "BANDEIRA BRANCA",                   na.rm=TRUE),
    median_tot_rev         = median(tot.rev,         na.rm=TRUE),
    median_std_pdev        = median(std.pdev,        na.rm=TRUE)
  )

print(comparison_3way)

library(stargazer)

stargazer(as.data.frame(comparison),
          summary = FALSE,
          rownames = FALSE,
          out = "/proj/econ/duartema/GasBrazil/Final/Elasticity Diagnosis/comparison_table.tex",
          title = "Flagged vs Non-Flagged Stations",
          label = "tab:comparison")

stargazer(as.data.frame(comparison_3way),
          summary = FALSE,
          rownames = FALSE,
          out = "/proj/econ/duartema/GasBrazil/Final/Elasticity Diagnosis/comparison_3way_table.tex",
          title = "Too Elastic vs Too Inelastic vs Normal Stations",
          label = "tab:comparison3way")
