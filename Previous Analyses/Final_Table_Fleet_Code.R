# ============================================================
# FINAL TABLE - Fleet-Based Market Size
# ============================================================

library(modelsummary)
library(dplyr)

# Load fleet models and data (if not already in memory)
load('/proj/econ/duartema/GasBrazil/Mihika_Practice/spatial_regression_FLEET.RData')
load('Data/Geo.RData')

# ------------------------------------------------------------
# Updated run_elasticity: now accepts an explicit data argument
# so it works with sample_demand_fleet instead of sample_demand
# ------------------------------------------------------------
run_elasticity <- function(g1, g2, g3, alpha, model_name = "", data) {
  
  cat("Computing elasticities for", model_name, "...\n")
  
  gas_data_sp <- data %>%
    select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
    filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))
  
  months_list <- unique(gas_data_sp$month)
  n_months    <- length(months_list)
  elas_list   <- vector("list", n_months)
  
  cat("  Processing", n_months, "months...\n")
  
  for (i in seq_along(months_list)) {
    if (i %% 5 == 0) cat("    Month", i, "of", n_months, "\n")
    m     <- months_list[i]
    gas_m <- gas_data_sp %>% filter(month == m)
    Dmat  <- Geo$Dmatrix[[i]]
    elas_list[[i]] <- compute_elasticity_month(m, Dmat, gas_m, g1, g2, g3, alpha)
    rm(Dmat)
    if (i %% 10 == 0) gc(verbose = FALSE)
  }
  
  elas_df <- bind_rows(elas_list)
  
  elas_avg <- elas_df %>%
    group_by(cnpj) %>%
    summarise(
      mean_own_elas  = mean(own_elas,      na.rm = TRUE),
      mean_cross_b1  = mean(cross_elas_b1, na.rm = TRUE),
      mean_cross_b2  = mean(cross_elas_b2, na.rm = TRUE),
      mean_cross_b3  = mean(cross_elas_b3, na.rm = TRUE)
    ) %>%
    ungroup()
  
  cat("  Complete! Computed elasticities for", nrow(elas_avg), "stations.\n\n")
  return(elas_avg)
}

# ------------------------------------------------------------
# STEP 1: Check if results already computed and saved
# ------------------------------------------------------------
results_file <- '/proj/econ/duartema/GasBrazil/Mihika_Practice/all_elasticities_FLEET.RData'

if (file.exists(results_file)) {
  cat("Loading pre-computed elasticities...\n")
  load(results_file)
} else {
  
  # Extract coefficients from fleet models
  # sOLS: nothing instrumented, coefficients named plainly
  g1_sOLS <- coef(model_spatial_ols_fleet)["B1"]
  g2_sOLS <- coef(model_spatial_ols_fleet)["B2"]
  g3_sOLS <- coef(model_spatial_ols_fleet)["B3"]
  
  # sIV1: only price instrumented, B coefficients still named plainly
  g1_sIV1 <- coef(model_spatial_iv1_fleet)["B1"]
  g2_sIV1 <- coef(model_spatial_iv1_fleet)["B2"]
  g3_sIV1 <- coef(model_spatial_iv1_fleet)["B3"]
  
  # sIV2: price AND spatial terms instrumented, B coefficients have fit_ prefix
  g1_sIV2 <- coef(model_spatial_iv2_fleet)["fit_B1"]
  g2_sIV2 <- coef(model_spatial_iv2_fleet)["fit_B2"]
  g3_sIV2 <- coef(model_spatial_iv2_fleet)["fit_B3"]
  
  # Alpha: OLS from simple logit, IV from iv_gas_fleet
  alpha_ols <- coef(model_gas_fleet)["p_sell.imp"]
  alpha_iv  <- coef(iv_gas_fleet)["fit_p_sell.imp"]
  
  cat("\n=== Parameter Comparison (Fleet) ===\n")
  cat("Alpha OLS:", round(alpha_ols, 4), "  Alpha IV:", round(alpha_iv, 4), "\n")
  cat("sOLS gammas:", round(g1_sOLS, 4), round(g2_sOLS, 4), round(g3_sOLS, 4), "\n")
  cat("sIV1 gammas:", round(g1_sIV1, 4), round(g2_sIV1, 4), round(g3_sIV1, 4), "\n")
  cat("sIV2 gammas:", round(g1_sIV2, 4), round(g2_sIV2, 4), round(g3_sIV2, 4), "\n\n")
  
  start_time <- Sys.time()
  
  elas_sOLS_fleet <- run_elasticity(g1_sOLS, g2_sOLS, g3_sOLS, alpha_ols,
                                    "Spatial OLS (Fleet)", data = sample_demand_fleet)
  elas_sIV1_fleet <- run_elasticity(g1_sIV1, g2_sIV1, g3_sIV1, alpha_iv,
                                    "Spatial IV1 (Fleet)", data = sample_demand_fleet)
  elas_sIV2_fleet <- run_elasticity(g1_sIV2, g2_sIV2, g3_sIV2, alpha_iv,
                                    "Spatial IV2 (Fleet)", data = sample_demand_fleet)
  
  end_time <- Sys.time()
  cat("Total computation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n\n")
  
  save(elas_sOLS_fleet, elas_sIV1_fleet, elas_sIV2_fleet, file = results_file)
  cat("Results saved to:", results_file, "\n\n")
}

# ------------------------------------------------------------
# STEP 2: Pull elasticity stats
# ------------------------------------------------------------
med_own_sOLS  <- round(median(elas_sOLS_fleet$mean_own_elas,  na.rm = TRUE), 4)
med_b1_sOLS   <- round(median(elas_sOLS_fleet$mean_cross_b1,  na.rm = TRUE), 4)
med_b2_sOLS   <- round(median(elas_sOLS_fleet$mean_cross_b2,  na.rm = TRUE), 4)
med_b3_sOLS   <- round(median(elas_sOLS_fleet$mean_cross_b3,  na.rm = TRUE), 4)

med_own_sIV1  <- round(median(elas_sIV1_fleet$mean_own_elas,  na.rm = TRUE), 4)
med_b1_sIV1   <- round(median(elas_sIV1_fleet$mean_cross_b1,  na.rm = TRUE), 4)
med_b2_sIV1   <- round(median(elas_sIV1_fleet$mean_cross_b2,  na.rm = TRUE), 4)
med_b3_sIV1   <- round(median(elas_sIV1_fleet$mean_cross_b3,  na.rm = TRUE), 4)

med_own_sIV2  <- round(median(elas_sIV2_fleet$mean_own_elas,  na.rm = TRUE), 4)
med_b1_sIV2   <- round(median(elas_sIV2_fleet$mean_cross_b1,  na.rm = TRUE), 4)
med_b2_sIV2   <- round(median(elas_sIV2_fleet$mean_cross_b2,  na.rm = TRUE), 4)
med_b3_sIV2   <- round(median(elas_sIV2_fleet$mean_cross_b3,  na.rm = TRUE), 4)

cat("=== Elasticity Medians Comparison (Fleet) ===\n")
cat("Own elasticity:  sOLS =", med_own_sOLS, " sIV1 =", med_own_sIV1, " sIV2 =", med_own_sIV2, "\n")
cat("Cross B1:        sOLS =", med_b1_sOLS,  " sIV1 =", med_b1_sIV1,  " sIV2 =", med_b1_sIV2, "\n")
cat("Cross B2:        sOLS =", med_b2_sOLS,  " sIV1 =", med_b2_sIV1,  " sIV2 =", med_b2_sIV2, "\n")
cat("Cross B3:        sOLS =", med_b3_sOLS,  " sIV1 =", med_b3_sIV1,  " sIV2 =", med_b3_sIV2, "\n\n")

# First-stage F-stats from fleet IV models
f_iv_gas_fleet      <- round(fitstat(iv_gas_fleet,            "ivf")[[1]]$stat, 3)
f_spatial_iv1_fleet <- round(fitstat(model_spatial_iv1_fleet, "ivf")[[1]]$stat, 3)
f_spatial_iv2_fleet <- round(fitstat(model_spatial_iv2_fleet, "ivf")[[1]]$stat, 3)

# ------------------------------------------------------------
# STEP 3: Rename coefficients
# coef_map covers both plain (OLS) and fit_ (IV) prefixed names
# ------------------------------------------------------------
coef_map <- c(
  "p_sell.imp"                             = "Price",
  "fit_p_sell.imp"                         = "Price",
  "log.nstations_group"                    = "Log # stations in group",
  "pumps.imp"                              = "Pumps",
  "tank.imp"                               = "Tank size",
  "log.rent"                               = "Log rent",
  "pop.area"                               = "Population density",
  "distIPIRANGA PRODUTOS DE PETRÓLEO S.A"  = "Ipiranga",
  "distPETROBRAS DISTRIBUIDORA S.A."       = "Petrobras",
  "distRAIZEN COMBUSTÍVEIS S.A."           = "Raizen",
  "B1"                                     = "B1 (< 500m)",
  "B2"                                     = "B2 (500m-1.5km)",
  "B3"                                     = "B3 (1.5km-3km)",
  "fit_B1"                                 = "B1 (< 500m)",
  "fit_B2"                                 = "B2 (500m-1.5km)",
  "fit_B3"                                 = "B3 (1.5km-3km)"
)

# ------------------------------------------------------------
# STEP 4: Extra rows for footer
# Logit own-elasticity formula: alpha * (1 - Sj) * Pj
# Uses fleet shares and fleet alpha values
# ------------------------------------------------------------
alpha_logit_ols <- coef(model_gas_fleet)["p_sell.imp"]
alpha_logit_iv  <- coef(iv_gas_fleet)["fit_p_sell.imp"]

gas_elas_base <- sample_demand_fleet %>%
  select(cnpj, sjt, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(p_sell.imp))

med_own_logit_ols <- round(median(alpha_logit_ols * (1 - gas_elas_base$sjt) *
                                    gas_elas_base$p_sell.imp, na.rm = TRUE), 4)

med_own_logit_iv  <- round(median(alpha_logit_iv  * (1 - gas_elas_base$sjt) *
                                    gas_elas_base$p_sell.imp, na.rm = TRUE), 4)

cat("Median own elasticity - Logit OLS:", med_own_logit_ols, "\n")
cat("Median own elasticity - Logit IV: ", med_own_logit_iv,  "\n")

extra_rows <- data.frame(
  term = c(
    "Market size",
    "Median own elasticity",
    "Median cross (< 500m)",
    "Median cross (500m-1.5km)",
    "Median cross (1.5km-3km)",
    "First stage F-stat"
  ),
  OLS  = c("Fleet (200L/vehicle)", med_own_logit_ols, "—",         "—",         "—",         "—"),
  IV   = c("Fleet (200L/vehicle)", med_own_logit_iv,  "—",         "—",         "—",         f_iv_gas_fleet),
  sOLS = c("Fleet (200L/vehicle)", med_own_sOLS,      med_b1_sOLS, med_b2_sOLS, med_b3_sOLS, "—"),
  IV1  = c("Fleet (200L/vehicle)", med_own_sIV1,      med_b1_sIV1, med_b2_sIV1, med_b3_sIV1, f_spatial_iv1_fleet),
  IV2  = c("Fleet (200L/vehicle)", med_own_sIV2,      med_b1_sIV2, med_b2_sIV2, med_b3_sIV2, f_spatial_iv2_fleet)
)
attr(extra_rows, "position") <- c(25, 26, 27, 28, 29, 30)

# ------------------------------------------------------------
# STEP 5: Build table
# ------------------------------------------------------------
modelsummary(
  list("OLS"  = model_gas_fleet,
       "IV"   = iv_gas_fleet,
       "sOLS" = model_spatial_ols_fleet,
       "sIV1" = model_spatial_iv1_fleet,
       "sIV2" = model_spatial_iv2_fleet),
  coef_map   = coef_map,
  gof_map    = c("nobs", "r.squared"),
  add_rows   = extra_rows,
  stars      = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title      = "Table 4: Demand Estimates (Fleet-Based Market Size)",
  output     = "/proj/econ/duartema/GasBrazil/Final/Table3.txt"
)