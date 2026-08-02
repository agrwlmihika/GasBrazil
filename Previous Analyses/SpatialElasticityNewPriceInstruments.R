# Elasticity with Spatial Component
# ============================================================

load('/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_NEW.RData')
load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')

library(dplyr)
library(fixest)

# Rebuild sample_demand
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

# Rebuild iv_gas for alpha
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

# STEP 1: Extract parameters
alpha  <- coef(iv_gas)["fit_p_sell.imp"]
gamma1 <- coef(model_spatial_iv2)["fit_B1"]
gamma2 <- coef(model_spatial_iv2)["fit_B2"]
gamma3 <- coef(model_spatial_iv2)["fit_B3"]

cat("Alpha:", round(alpha,  4), "\n")
cat("Gamma1 (< 500m):       ", round(gamma1, 4), "\n")
cat("Gamma2 (500m-1.5km):   ", round(gamma2, 4), "\n")
cat("Gamma3 (1.5km-3km):    ", round(gamma3, 4), "\n")

# STEP 2: Prepare gas data
gas_data_sp <- sample_demand %>%
  filter(prod == "gas") %>%
  select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))

# STEP 3: Elasticity function (with sign fix)
compute_elasticity_month <- function(m, Dmat, gas_m, g1, g2, g3, a) {
  gas_ids  <- gas_m$id
  keep_ids <- gas_ids[gas_ids %in% rownames(Dmat)]
  if (length(keep_ids) < 2) return(NULL)
  gas_m <- gas_m %>% filter(id %in% keep_ids)
  n     <- nrow(gas_m)
  S     <- gas_m$sjt
  s0    <- unique(gas_m$s0t)[1]
  D  <- Dmat[gas_m$id, gas_m$id]
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)
  S_sum   <- outer(S, S, "+")
  gamma_w <- g1 * w1 + g2 * w2 + g3 * w3
  dG_ds   <- gamma_w / S_sum
  for (j in seq_len(n)) {
    sk_over_sum    <- S / (S[j] + S)
    sk_over_sum[j] <- 0
    spatial_term   <- -sum((g1*w1[j,] + g2*w2[j,] + g3*w3[j,]) * sk_over_sum)
    dG_ds[j, j]    <- (1 + spatial_term) / S[j]
  }
  M       <- dG_ds
  diag(M) <- diag(M) + 1/s0
  M_inv <- tryCatch(solve(M), error = function(e) NULL)
  if (is.null(M_inv)) return(NULL)
  P        <- gas_m$p_sell.imp
  elas_mat <- a * M_inv
  elas_mat <- sweep(elas_mat, 2, P, "*")
  elas_mat <- sweep(elas_mat, 1, S, "/")
  own_elas <- diag(elas_mat)
  cross_elas_b1 <- sapply(seq_len(n), function(j) {
    idx <- which(w1[j, ] & seq_len(n) != j)
    if (length(idx) == 0) return(NA)
    mean(elas_mat[j, idx], na.rm = TRUE)
  })
  cross_elas_b2 <- sapply(seq_len(n), function(j) {
    idx <- which(w2[j, ] & seq_len(n) != j)
    if (length(idx) == 0) return(NA)
    mean(elas_mat[j, idx], na.rm = TRUE)
  })
  cross_elas_b3 <- sapply(seq_len(n), function(j) {
    idx <- which(w3[j, ] & seq_len(n) != j)
    if (length(idx) == 0) return(NA)
    mean(elas_mat[j, idx], na.rm = TRUE)
  })
  data.frame(
    id            = gas_m$id,
    cnpj          = gas_m$cnpj,
    month         = m,
    own_elas      = own_elas,
    cross_elas_b1 = cross_elas_b1,
    cross_elas_b2 = cross_elas_b2,
    cross_elas_b3 = cross_elas_b3
  )
}

# STEP 4: Loop over months
months_list  <- unique(gas_data_sp$month)
elas_list_sp <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m                 <- months_list[i]
  gas_m             <- gas_data_sp %>% filter(month == m)
  Dmat              <- Geo$Dmatrix[[i]]
  elas_list_sp[[i]] <- compute_elasticity_month(m, Dmat, gas_m,
                                                gamma1, gamma2, gamma3, alpha)
  rm(Dmat); gc()
}

elas_sp <- bind_rows(elas_list_sp)

# STEP 5: Average per station across months
elas_sp_avg <- elas_sp %>%
  group_by(cnpj) %>%
  summarise(
    mean_own_elas  = mean(own_elas,      na.rm = TRUE),
    mean_cross_b1  = mean(cross_elas_b1, na.rm = TRUE),
    mean_cross_b2  = mean(cross_elas_b2, na.rm = TRUE),
    mean_cross_b3  = mean(cross_elas_b3, na.rm = TRUE)
  ) %>%
  ungroup()

cat("\n=== Spatial Elasticity Table (first 20 rows) ===\n")
print(head(elas_sp_avg, 20))

cat("\n=== Overall Medians ===\n")
cat("Median own-price elasticity:            ",
    round(median(elas_sp_avg$mean_own_elas,  na.rm=TRUE), 4), "\n")
cat("Median cross elasticity B1 (<500m):     ",
    round(median(elas_sp_avg$mean_cross_b1, na.rm=TRUE), 4), "\n")
cat("Median cross elasticity B2 (500m-1.5km):",
    round(median(elas_sp_avg$mean_cross_b2, na.rm=TRUE), 4), "\n")
cat("Median cross elasticity B3 (1.5-3km):   ",
    round(median(elas_sp_avg$mean_cross_b3, na.rm=TRUE), 4), "\n")

# Save new instrument elasticity results
save(elas_sp_avg,
     file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/all_elasticities_NEW.RData')
cat("Saved to all_elasticities_NEW.RData\n")