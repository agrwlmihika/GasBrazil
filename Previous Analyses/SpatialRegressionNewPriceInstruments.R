# Demand Model with Spatial Competition Terms
# ------------------------------------------------------------

load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')

library(dplyr)
library(fixest)

#===============Demand Model Code============================================================

# STEP 1: Define market size and compute shares
sample_demand <- sample %>%
  group_by(month, prod) %>%
  mutate(
    total_sales_mt  = sum(sale.imp, na.rm = TRUE),
    market_size_mt  = total_sales_mt * 1.5,
    sjt             = sale.imp / market_size_mt,
    s0t             = 1 - sum(sjt, na.rm = TRUE),
    log_ratio       = log(sjt / s0t)
  ) %>%
  ungroup() %>%
  filter(is.finite(log_ratio), s0t > 0)

summary(sample_demand$s0t)
summary(sample_demand$sjt)

# STEP 2: Run regression for gas
run_demand <- function(df, product) {
  df_prod <- df %>% filter(prod == product)
  feols(
    log_ratio ~ p_sell.imp +
      log.nstations_group +
      pumps.imp +
      tank.imp +
      log.rent +
      pop.area +
      dist,
    cluster = ~cl,
    data = df_prod
  )
}

model_gas <- run_demand(sample_demand, "gas")

# STEP 3: Build results table
extract_results <- function(model, product) {
  coefs <- coef(model)
  ses   <- se(model)
  data.frame(
    product = product,
    name    = names(coefs),
    coeff   = round(coefs, 4),
    SE      = round(ses,   4),
    row.names = NULL
  )
}

results_gas   <- extract_results(model_gas, "Gasoline")
results_table <- results_gas
print(results_table)

# ======================IV Regression Code===========================================================
sample_iv <- sample_demand

iv_gas <- feols(
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
  data = sample_iv %>% filter(prod == "gas")
)

iv_table <- data.frame(
  product = "Gasoline",
  name    = names(coef(iv_gas)),
  coeff   = round(coef(iv_gas), 4),
  SE      = round(se(iv_gas),   4),
  row.names = NULL
)
print(iv_table)
cat("First stage F-stat IV:\n"); print(fitstat(iv_gas, "ivf"))

#=========================Spatial Regression=====================================
gas_shares <- sample_demand %>%
  filter(prod == "gas") %>%
  select(id, month, cnpj, sjt) %>%
  filter(!is.na(sjt))

build_Bg <- function(month_id, Dmat, shares_m) {
  gas_ids  <- shares_m$id
  D        <- Dmat[gas_ids, gas_ids]
  n        <- nrow(D)
  Sj       <- shares_m$sjt
  ratio_mat <- outer(Sj, Sj, function(sj, sk) log(1 / (1 + sk/sj)))
  diag(ratio_mat) <- 0
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)
  B1 <- rowSums(w1 * ratio_mat, na.rm = TRUE)
  B2 <- rowSums(w2 * ratio_mat, na.rm = TRUE)
  B3 <- rowSums(w3 * ratio_mat, na.rm = TRUE)
  data.frame(id = gas_ids, B1 = B1, B2 = B2, B3 = B3)
}

months_list <- unique(gas_shares$month)
Bg_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m        <- months_list[i]
  shares_m <- gas_shares %>% filter(month == m)
  Dmat     <- Geo$Dmatrix[[i]]
  gas_ids  <- shares_m$id
  keep_ids <- gas_ids[gas_ids %in% rownames(Dmat)]
  if (length(keep_ids) < 2) next
  shares_m <- shares_m %>% filter(id %in% keep_ids)
  Bg_list[[i]] <- build_Bg(m, Dmat, shares_m)
  rm(Dmat); gc()
}

Bg_df <- bind_rows(Bg_list)
cat("Spatial terms built. Rows:", nrow(Bg_df), "\n")

sample_spatial <- sample_demand %>%
  filter(prod == "gas") %>%
  left_join(Bg_df, by = "id") %>%
  filter(!is.na(B1), !is.na(B2), !is.na(B3))
cat("Sample after merge:", nrow(sample_spatial), "rows\n")

# Regression 1 — Spatial OLS
model_spatial_ols <- feols(
  log_ratio ~ p_sell.imp +
    log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist +
    B1 + B2 + B3,
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Regression 1: Spatial OLS ===\n")
print(data.frame(
  name  = names(coef(model_spatial_ols)),
  coeff = round(coef(model_spatial_ols), 4),
  SE    = round(se(model_spatial_ols),   4),
  row.names = NULL
))

# Regression 2 — Spatial IV1 (price only, clean instruments)
model_spatial_iv1 <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist +
    B1 + B2 + B3 |
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
  data = sample_spatial
)

cat("\n=== Regression 2: Spatial IV1 (price only) ===\n")
cat("First stage F-stat:\n"); print(fitstat(model_spatial_iv1, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv1)),
  coeff = round(coef(model_spatial_iv1), 4),
  SE    = round(se(model_spatial_iv1),   4),
  row.names = NULL
))

# Regression 3 — Spatial IV2 (price + spatial terms)
model_spatial_iv2 <- feols(
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
  data = sample_spatial
)

cat("\n=== Regression 3: Spatial IV2 (price + spatial) ===\n")
cat("First stage F-stats:\n"); print(fitstat(model_spatial_iv2, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv2)),
  coeff = round(coef(model_spatial_iv2), 4),
  SE    = round(se(model_spatial_iv2),   4),
  row.names = NULL
))

save(model_spatial_ols, model_spatial_iv1, model_spatial_iv2,
     sample_spatial, iv_gas, model_gas,
     file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_NEW.RData')