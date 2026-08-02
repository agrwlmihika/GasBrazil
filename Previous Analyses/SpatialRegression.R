# Demand Model with Spatial Competition Terms
# ------------------------------------------------------------

load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')

library(dplyr)
library(fixest)

#===============Demand Model Code============================================================

# STEP 1: Define market size and compute shares
# ------------------------------------------------------------
sample_demand <- sample %>%
  # Compute total sales per month x product
  group_by(month, prod) %>%
  mutate(
    total_sales_mt  = sum(sale.imp, na.rm = TRUE),
    market_size_mt  = total_sales_mt * 1.5,         # scale factor = 1.5
    sjt             = sale.imp / market_size_mt,     # station share
    s0t             = 1 - sum(sjt, na.rm = TRUE),    # outside share
    log_ratio       = log(sjt / s0t)                 # LHS of regression
  ) %>%
  ungroup() %>%
  # Drop any rows where LHS is undefined (log of negative/zero)
  filter(is.finite(log_ratio), s0t > 0)

#Check — s0t should be between 0 and 1
summary(sample_demand$s0t)
summary(sample_demand$sjt)

# ------------------------------------------------------------
# STEP 2: Run regression for gas
# ------------------------------------------------------------
# Model: log(sjt/s0t) = alpha*Price + B*X + error
# X = brand (dist), nstations_group, pumps, tank, log.rent, pop.area
# SEs clustered by cl (month x neighborhood)

run_demand <- function(df, product) {
  
  df_prod <- df %>% filter(prod == product)
  
  feols(
    log_ratio ~ p_sell.imp +              # price
      log.nstations_group +     # log # stations owned by retailer
      pumps.imp +               # # pumps
      tank.imp +                # tank size
      log.rent +                # log avg neighborhood rent
      pop.area +                # neighborhood population density
      dist,                     # brand dummies (factor)
    cluster = ~cl,                        # cluster SE by month x neighborhood
    data = df_prod
  )
}

model_gas <- run_demand(sample_demand, "gas")

# STEP 3: Build results table (name, coeff, SE)
# ------------------------------------------------------------
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

results_gas <- extract_results(model_gas, "Gasoline")
results_table <- results_gas
print(results_table)



# ======================IV Regression Code ===========================================================

sample_iv <- sample_demand

iv_gas <- feols(
  log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
    log.rent + pop.area + dist |
    cl |
    p_sell.imp ~ destillery_price.d1:pumps.imp.mdev+
    refinery_price.d1:close_oppo.mdev+destillery_price.d1:close_oppo.mdev+
    refinery_price.d1:n_unbranded.mdev+destillery_price.d1:n_unbranded.mdev+
    refinery_price.d1:n_oppo1km.mdev+destillery_price.d1:n_oppo1km.mdev + 
    Dstd.pdev.dist+Dstd.pdev.dist:n_unbranded.mdev + Dstd.pdev.dist:n_oppo1km.mdev+  Dstd.pdev.dist:close_oppo.mdev,
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

#=========================Spatial Regression =====================================
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
  ratio_mat <- outer(Sj, Sj, function(sj, sk) log(1 / (1 + sk/sj)))
  diag(ratio_mat) <- 0  # set diagonal (own) to zero
  
  # Distance rings
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)          
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
  
  Dmat <- Geo$Dmatrix[[i]]
  
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

save(model_spatial_iv2, sample_spatial, iv_gas, 
     file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_output.RData')