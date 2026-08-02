# Spatial Regression - Fleet Based Market Size
# WITHOUT month x neighborhood fixed effects

load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')
load('/proj/econ/duartema/GasData/Database/Frota.RData')

library(dplyr)
library(fixest)

# ============================================================
# STEP 1: Build fleet-based market size from Frota data
# ============================================================
fleet_monthly <- Frota %>%
  filter(state == "DF") %>%
  group_by(month) %>%
  summarise(fleet_t = sum(AUTOMOVEL, na.rm = TRUE), .groups = "drop") %>%
  mutate(market_size_fleet = fleet_t * 100) %>%
  arrange(month)

cat("Fleet monthly (head):\n");  print(head(fleet_monthly, 12))
cat("Fleet monthly (tail):\n");  print(tail(fleet_monthly, 12))
cat("Month range:", min(fleet_monthly$month), "to", max(fleet_monthly$month), "\n")
cat("Fleet range:", round(min(fleet_monthly$fleet_t)), "to",
    round(max(fleet_monthly$fleet_t)), "\n")
cat("Market size range:", round(min(fleet_monthly$market_size_fleet)), "to",
    round(max(fleet_monthly$market_size_fleet)), "\n")

sample_months  <- unique(sample$month[sample$prod == "gas"])
missing_months <- sample_months[!sample_months %in% fleet_monthly$month]
cat("\nSample months not in Frota:", length(missing_months), "\n")
if (length(missing_months) > 0) print(missing_months)

# ============================================================
# STEP 2: Define market size and compute shares
# ============================================================
sample_demand_fleet_nofe <- sample %>%
  filter(prod == "gas") %>%
  left_join(fleet_monthly, by = "month") %>%
  group_by(month) %>%
  mutate(
    total_sales_mt = sum(sale.imp, na.rm = TRUE),
    sjt            = sale.imp / market_size_fleet,
    s0t            = 1 - sum(sjt, na.rm = TRUE),
    log_ratio      = log(sjt / s0t)
  ) %>%
  ungroup() %>%
  filter(is.finite(log_ratio), s0t > 0, sjt > 0)

cat("\ns0t range:", round(min(sample_demand_fleet_nofe$s0t), 4),
    "to", round(max(sample_demand_fleet_nofe$s0t), 4), "\n")
cat("sjt range:", round(min(sample_demand_fleet_nofe$sjt), 6),
    "to", round(max(sample_demand_fleet_nofe$sjt), 6), "\n")
cat("Rows:", nrow(sample_demand_fleet_nofe), "\n")

# ============================================================
# STEP 3: OLS demand regression — NO fixed effects
# ============================================================
model_gas_fleet_nofe <- feols(
  log_ratio ~ p_sell.imp +
    log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist,
  cluster = ~cl,
  data = sample_demand_fleet_nofe
)

cat("\n=== OLS Demand (Fleet, NO FE) ===\n")
print(data.frame(
  name  = names(coef(model_gas_fleet_nofe)),
  coeff = round(coef(model_gas_fleet_nofe), 4),
  SE    = round(se(model_gas_fleet_nofe),   4),
  row.names = NULL
))

# ============================================================
# STEP 4: IV regression — NO fixed effects
# ============================================================
iv_gas_fleet_nofe <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    p_sell.imp ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_demand_fleet_nofe
)

cat("\n=== IV (Fleet, NO FE) ===\n")
cat("First stage F-stat:\n"); print(fitstat(iv_gas_fleet_nofe, "ivf"))

# ============================================================
# STEP 5: Build spatial terms
# ============================================================
gas_shares_fleet_nofe <- sample_demand_fleet_nofe %>%
  select(id, month, cnpj, sjt) %>%
  filter(!is.na(sjt))

build_Bg <- function(month_id, Dmat, shares_m) {
  gas_ids   <- shares_m$id
  D         <- Dmat[gas_ids, gas_ids]
  Sj        <- shares_m$sjt
  ratio_mat <- outer(Sj, Sj, function(sj, sk) log(1 / (1 + sk/sj)))
  diag(ratio_mat) <- 0
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)
  data.frame(
    id = gas_ids,
    B1 = rowSums(w1 * ratio_mat, na.rm = TRUE),
    B2 = rowSums(w2 * ratio_mat, na.rm = TRUE),
    B3 = rowSums(w3 * ratio_mat, na.rm = TRUE)
  )
}

months_list <- unique(gas_shares_fleet_nofe$month)
Bg_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m        <- months_list[i]
  shares_m <- gas_shares_fleet_nofe %>% filter(month == m)
  Dmat     <- Geo$Dmatrix[[i]]
  gas_ids  <- shares_m$id
  keep_ids <- gas_ids[gas_ids %in% rownames(Dmat)]
  if (length(keep_ids) < 2) next
  shares_m <- shares_m %>% filter(id %in% keep_ids)
  Bg_list[[i]] <- build_Bg(m, Dmat, shares_m)
  rm(Dmat); gc()
}

Bg_df_fleet_nofe <- bind_rows(Bg_list)
cat("Spatial terms built:", nrow(Bg_df_fleet_nofe), "rows\n")

sample_spatial_fleet_nofe <- sample_demand_fleet_nofe %>%
  left_join(Bg_df_fleet_nofe, by = "id") %>%
  filter(!is.na(B1), !is.na(B2), !is.na(B3))
cat("Sample after spatial merge:", nrow(sample_spatial_fleet_nofe), "rows\n")

# ============================================================
# STEP 6: Spatial regressions — NO fixed effects
# ============================================================

# Spatial OLS
model_spatial_ols_fleet_nofe <- feols(
  log_ratio ~ p_sell.imp +
    log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist + B1 + B2 + B3,
  cluster = ~cl,
  data = sample_spatial_fleet_nofe
)

cat("\n=== Spatial OLS (Fleet, NO FE) ===\n")
print(data.frame(
  name  = names(coef(model_spatial_ols_fleet_nofe)),
  coeff = round(coef(model_spatial_ols_fleet_nofe), 4),
  SE    = round(se(model_spatial_ols_fleet_nofe),   4),
  row.names = NULL
))

# Spatial IV1 — price instrumented, B terms exogenous, NO FE
model_spatial_iv1_fleet_nofe <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist + B1 + B2 + B3 |
    p_sell.imp ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_spatial_fleet_nofe
)

cat("\n=== Spatial IV1 (Fleet, NO FE) ===\n")
cat("First stage F-stat:\n"); print(fitstat(model_spatial_iv1_fleet_nofe, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv1_fleet_nofe)),
  coeff = round(coef(model_spatial_iv1_fleet_nofe), 4),
  SE    = round(se(model_spatial_iv1_fleet_nofe),   4),
  row.names = NULL
))

# Spatial IV2 — price AND B terms instrumented, NO FE
model_spatial_iv2_fleet_nofe <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    p_sell.imp + B1 + B2 + B3 ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_spatial_fleet_nofe
)

cat("\n=== Spatial IV2 (Fleet, NO FE) ===\n")
cat("First stage F-stats:\n"); print(fitstat(model_spatial_iv2_fleet_nofe, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv2_fleet_nofe)),
  coeff = round(coef(model_spatial_iv2_fleet_nofe), 4),
  SE    = round(se(model_spatial_iv2_fleet_nofe),   4),
  row.names = NULL
))

# ============================================================
# SAVE
# ============================================================
save(model_gas_fleet_nofe, iv_gas_fleet_nofe,
     model_spatial_ols_fleet_nofe, model_spatial_iv1_fleet_nofe,
     model_spatial_iv2_fleet_nofe,
     sample_spatial_fleet_nofe, sample_demand_fleet_nofe,
     file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_FLEET_NOFE.RData')