# Spatial Regression - Fleet Based Market Size
# Market size = fleet_t * 1000km/month * 2 months * (1/10 L/km)
# = fleet_t * 200 liters per vehicle
# ------------------------------------------------------------

load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')
load('Data/commute.RData')

library(dplyr)
library(fixest)

# ------------------------------------------------------------
# STEP 1: Build fleet-based market size
# ------------------------------------------------------------
# fleet2011 = 956,344 vehicles in Brasilia in 2011
# We apply annual growth rate to get monthly fleet estimate
# Brazil average fleet growth rate ~4% per year (DENATRAN)
fleet_base <- unique(commute$fleet2011)
cat("Base fleet (2011):", fleet_base, "\n")

# Create monthly fleet estimates with 4% annual growth
# Base year = 2011, base month = 2011-01
fleet_monthly <- data.frame(
  month = format(seq(as.Date("2010-01-01"), 
                     as.Date("2017-12-01"), by="month"), "%Y-%m")
) %>%
  mutate(
    years_from_2011 = as.numeric(as.Date(paste0(month, "-01")) - 
                                   as.Date("2011-01-01")) / 365.25,
    fleet_t         = fleet_base * (1.04 ^ years_from_2011),
    # Market size = fleet * 200 liters (2 months driving at 10km/L, 1000km/month)
    market_size_fleet = fleet_t * 200
  ) %>%
  select(month, fleet_t, market_size_fleet)

cat("Fleet market size sample:\n")
print(head(fleet_monthly))

# ------------------------------------------------------------
# STEP 2: Define market size and compute shares
# ------------------------------------------------------------
sample_demand_fleet <- sample %>%
  filter(prod == "gas") %>%
  left_join(fleet_monthly, by = "month") %>%
  group_by(month) %>%
  mutate(
    total_sales_mt    = sum(sale.imp, na.rm = TRUE),
    sjt               = sale.imp / market_size_fleet,
    s0t               = 1 - sum(sjt, na.rm = TRUE),
    log_ratio         = log(sjt / s0t)
  ) %>%
  ungroup() %>%
  filter(is.finite(log_ratio), s0t > 0, sjt > 0)

cat("s0t range:", min(sample_demand_fleet$s0t), "to", max(sample_demand_fleet$s0t), "\n")
cat("sjt range:", min(sample_demand_fleet$sjt), "to", max(sample_demand_fleet$sjt), "\n")
cat("Rows:", nrow(sample_demand_fleet), "\n")

# ------------------------------------------------------------
# STEP 3: OLS demand regression
# ------------------------------------------------------------
run_demand <- function(df) {
  feols(
    log_ratio ~ p_sell.imp +
      log.nstations_group + pumps.imp + tank.imp +
      log.rent + pop.area + dist,
    cluster = ~cl,
    data = df
  )
}

model_gas_fleet <- run_demand(sample_demand_fleet)

cat("\n=== OLS Demand (Fleet Market Size) ===\n")
print(data.frame(
  name  = names(coef(model_gas_fleet)),
  coeff = round(coef(model_gas_fleet), 4),
  SE    = round(se(model_gas_fleet),   4),
  row.names = NULL
))

# ------------------------------------------------------------
# STEP 4: IV regression
# ------------------------------------------------------------
iv_gas_fleet <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    cl |
    p_sell.imp ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_demand_fleet
)

cat("\n=== IV (Fleet Market Size) ===\n")
cat("First stage F-stat:\n"); print(fitstat(iv_gas_fleet, "ivf"))

# ------------------------------------------------------------
# STEP 5: Build spatial terms
# ------------------------------------------------------------
gas_shares_fleet <- sample_demand_fleet %>%
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

months_list <- unique(gas_shares_fleet$month)
Bg_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m        <- months_list[i]
  shares_m <- gas_shares_fleet %>% filter(month == m)
  Dmat     <- Geo$Dmatrix[[i]]
  gas_ids  <- shares_m$id
  keep_ids <- gas_ids[gas_ids %in% rownames(Dmat)]
  if (length(keep_ids) < 2) next
  shares_m <- shares_m %>% filter(id %in% keep_ids)
  Bg_list[[i]] <- build_Bg(m, Dmat, shares_m)
  rm(Dmat); gc()
}

Bg_df_fleet <- bind_rows(Bg_list)
cat("Spatial terms built:", nrow(Bg_df_fleet), "rows\n")

sample_spatial_fleet <- sample_demand_fleet %>%
  left_join(Bg_df_fleet, by = "id") %>%
  filter(!is.na(B1), !is.na(B2), !is.na(B3))
cat("Sample after merge:", nrow(sample_spatial_fleet), "rows\n")

# ------------------------------------------------------------
# STEP 6: Spatial regressions
# ------------------------------------------------------------
model_spatial_ols_fleet <- feols(
  log_ratio ~ p_sell.imp +
    log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist + B1 + B2 + B3,
  cluster = ~cl,
  data = sample_spatial_fleet
)

cat("\n=== Spatial OLS (Fleet) ===\n")
print(data.frame(
  name  = names(coef(model_spatial_ols_fleet)),
  coeff = round(coef(model_spatial_ols_fleet), 4),
  SE    = round(se(model_spatial_ols_fleet),   4),
  row.names = NULL
))

model_spatial_iv1_fleet <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist + B1 + B2 + B3 |
    cl |
    p_sell.imp ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_spatial_fleet
)

cat("\n=== Spatial IV1 (Fleet) ===\n")
cat("First stage F-stat:\n"); print(fitstat(model_spatial_iv1_fleet, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv1_fleet)),
  coeff = round(coef(model_spatial_iv1_fleet), 4),
  SE    = round(se(model_spatial_iv1_fleet),   4),
  row.names = NULL
))

model_spatial_iv2_fleet <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    cl |
    p_sell.imp + B1 + B2 + B3 ~
    refinery_price.d1:pumps.imp.mdev +
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:tank.imp.mdev +
    destillery_price.d1:tank.imp.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    Dstd.pdev.dist,
  cluster = ~cl,
  data = sample_spatial_fleet
)

cat("\n=== Spatial IV2 (Fleet) ===\n")
cat("First stage F-stats:\n"); print(fitstat(model_spatial_iv2_fleet, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv2_fleet)),
  coeff = round(coef(model_spatial_iv2_fleet), 4),
  SE    = round(se(model_spatial_iv2_fleet),   4),
  row.names = NULL
))

save(model_gas_fleet, iv_gas_fleet,
     model_spatial_ols_fleet, model_spatial_iv1_fleet, model_spatial_iv2_fleet,
     sample_spatial_fleet, sample_demand_fleet,
     file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_FLEET.RData')
cat("\nSaved to spatial_regression_FLEET.RData\n")