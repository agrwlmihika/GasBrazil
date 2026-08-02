# Demand Model with Spatial Competition Terms
# ============================================================
library(dplyr)
library(fixest)

# ------------------------------------------------------------
# STEP 1: Build spatial competition terms B^g
# For each station j and distance ring g, compute:
# B^g_j = sum over k≠j of w^g_jk * ln(1 / (1 + Sk/Sj))
# ------------------------------------------------------------

# Gas only shares needed for spatial terms
gas_shares <- sample_demand %>%
  filter(prod == "gas") %>%
  select(id, month, cnpj, sjt) %>%
  filter(!is.na(sjt))

# Function to build B^g for one month
build_Bg <- function(month_id, Dmat, shares_m) {
  
  # Filter distance matrix to gas rows only for this month
  gas_ids <- shares_m$id
  
  # Subset the distance matrix to gas-only rows and columns
  D <- Dmat[gas_ids, gas_ids]
  
  n   <- nrow(D)
  Sj  <- shares_m$sjt
  
  # Precompute ln(1/(1 + Sk/Sj)) for all j,k pairs
  # This is an n x n matrix where entry [j,k] = ln(1/(1 + Sk/Sj))
  # We compute it efficiently using outer()
  ratio_mat <- outer(Sj, Sj, function(sj, sk) log(1 / (1 + sk/sj)))
  diag(ratio_mat) <- 0  # set diagonal (own) to zero
  
  # Distance rings
  w1 <- (D > 0)    & (D < 0.5)          # < 500m
  w2 <- (D >= 0.5) & (D < 1.5)          # 500m to 1.5km
  w3 <- (D >= 1.5) & (D < 3.0)          # 1.5km to 3km
  # g=4: w=0, drops out
  
  # B^g = rowSums of w^g * ratio_mat for each ring
  B1 <- rowSums(w1 * ratio_mat, na.rm = TRUE)
  B2 <- rowSums(w2 * ratio_mat, na.rm = TRUE)
  B3 <- rowSums(w3 * ratio_mat, na.rm = TRUE)
  
  data.frame(
    id = gas_ids,
    B1 = B1,
    B2 = B2,
    B3 = B3
  )
}

# Loop over months — compute B^g month by month to manage memory
months_list <- unique(gas_shares$month)
Bg_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  
  m        <- months_list[i]
  shares_m <- gas_shares %>% filter(month == m)
  
  # Get the distance matrix for this month
  # Dmatrix is indexed by month — find the right one
  # rownames are formatted as "month|cnpj|prod"
  Dmat <- Geo$Dmatrix[[i]]
  
  # Keep only gas rows that exist in both shares and Dmatrix
  gas_ids  <- shares_m$id
  keep_ids <- gas_ids[gas_ids %in% rownames(Dmat)]
  
  if (length(keep_ids) < 2) next
  
  shares_m <- shares_m %>% filter(id %in% keep_ids)
  
  Bg_list[[i]] <- build_Bg(m, Dmat, shares_m)
  
  # Clear memory
  rm(Dmat)
  gc()
}

# Combine all months
Bg_df <- bind_rows(Bg_list)

cat("Spatial terms built. Rows:", nrow(Bg_df), "\n")

# ------------------------------------------------------------
# STEP 2: Merge spatial terms into sample_demand
# ------------------------------------------------------------
sample_spatial <- sample_demand %>%
  filter(prod == "gas") %>%
  left_join(Bg_df, by = "id") %>%
  filter(!is.na(B1), !is.na(B2), !is.na(B3))

cat("Sample after merge:", nrow(sample_spatial), "rows\n")

# ------------------------------------------------------------
# STEP 3: Regression 1 — OLS with spatial terms
# ------------------------------------------------------------
model_spatial_ols <- feols(
  log_ratio ~ p_sell.imp +
    log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist +
    B1 + B2 + B3,         # spatial competition terms
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Regression 1: OLS with Spatial Terms ===\n")
print(data.frame(
  name  = names(coef(model_spatial_ols)),
  coeff = round(coef(model_spatial_ols), 4),
  SE    = round(se(model_spatial_ols),   4),
  row.names = NULL
))

# ------------------------------------------------------------
# STEP 4: Regression 2 — IV, price instrumented only
# Z = n_oppo1km * gasC_cost
# ------------------------------------------------------------
sample_spatial <- sample_spatial %>%
  mutate(Z = n_oppo1km * gasC_cost)

#model_spatial_iv1 <- feols(
#  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
 #   log.rent + pop.area + dist +
#    B1 + B2 + B3 |
#    cl |
#    p_sell.imp ~ Z,
#  cluster = ~cl,
#  data = sample_spatial
#)

model_spatial_iv1 <- feols(
log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
  log.rent + pop.area + dist +
  B1 + B2 + B3 |
  cl |
  p_sell.imp ~ destillery_price.d1:pumps.imp.mdev+
  refinery_price.d1:close_oppo.mdev+destillery_price.d1:close_oppo.mdev+
  refinery_price.d1:n_unbranded.mdev+destillery_price.d1:n_unbranded.mdev+
  refinery_price.d1:n_oppo1km.mdev+destillery_price.d1:n_oppo1km.mdev,
cluster = ~cl,
data = sample_spatial
)



cat("\n=== Regression 2: IV (price only) with Spatial Terms ===\n")
cat("First stage F-stat: "); print(fitstat(model_spatial_iv1, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv1)),
  coeff = round(coef(model_spatial_iv1), 4),
  SE    = round(se(model_spatial_iv1),   4),
  row.names = NULL
))

# ------------------------------------------------------------
# STEP 5: Regression 3 — IV, price + spatial terms instrumented
# Y instruments for price:
#   refinery_price.d1 and destillery_price.d1 interacted with
#   competition variables + Dstd.pdev.dist interactions
# Spatial terms B1, B2, B3 treated as endogenous
# Instrumented using lagged/cost interactions
# ------------------------------------------------------------
model_spatial_iv2 <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    cl |
    p_sell.imp + B1 + B2 + B3 ~
    destillery_price.d1:pumps.imp.mdev +
    refinery_price.d1:close_oppo.mdev +
    destillery_price.d1:close_oppo.mdev +
    refinery_price.d1:n_unbranded.mdev +
    destillery_price.d1:n_unbranded.mdev +
    refinery_price.d1:n_oppo1km.mdev +
    destillery_price.d1:n_oppo1km.mdev +
    Dstd.pdev.dist +
    Dstd.pdev.dist:n_unbranded.mdev +
    Dstd.pdev.dist:n_oppo1km.mdev +
    Dstd.pdev.dist:close_oppo.mdev,
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Regression 3: IV (price + spatial) with Spatial Terms ===\n")
cat("First stage F-stats:\n"); print(fitstat(model_spatial_iv2, "ivf"))
print(data.frame(
  name  = names(coef(model_spatial_iv2)),
  coeff = round(coef(model_spatial_iv2), 4),
  SE    = round(se(model_spatial_iv2),   4),
  row.names = NULL
))