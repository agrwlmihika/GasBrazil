library(modelsummary)
library(dplyr)

# ------------------------------------------------------------
# OPTIMIZED FUNCTIONS with progress tracking
# ------------------------------------------------------------

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

run_elasticity <- function(g1, g2, g3, alpha, model_name = "") {
  
  cat("Computing elasticities for", model_name, "...\n")
  
  gas_data_sp <- sample_demand %>%
    filter(prod == "gas") %>%
    select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
    filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))
  
  months_list <- unique(gas_data_sp$month)
  n_months    <- length(months_list)
  elas_list   <- vector("list", n_months)
  
  cat("  Processing", n_months, "months...\n")
  
  for (i in seq_along(months_list)) {
    if (i %% 5 == 0) cat("    Month", i, "of", n_months, "\n")  # Progress every 5 months
    
    m     <- months_list[i]
    gas_m <- gas_data_sp %>% filter(month == m)
    Dmat  <- Geo$Dmatrix[[i]]
    
    elas_list[[i]] <- compute_elasticity_month(m, Dmat, gas_m, g1, g2, g3, alpha)
    
    rm(Dmat)
    if (i %% 10 == 0) gc(verbose = FALSE)  # Garbage collect every 10 months
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
results_file <- '/proj/econ/duartema/GasBrazil/Mihika_Practice/all_elasticities.RData'

if (file.exists(results_file)) {
  cat("Loading pre-computed elasticities...\n")
  load(results_file)
} else {
  
  # Extract coefficients
  g1_sOLS <- coef(model_spatial_ols)["B1"]
  g2_sOLS <- coef(model_spatial_ols)["B2"]
  g3_sOLS <- coef(model_spatial_ols)["B3"]
  
  g1_sIV1 <- coef(model_spatial_iv1)["B1"]
  g2_sIV1 <- coef(model_spatial_iv1)["B2"]
  g3_sIV1 <- coef(model_spatial_iv1)["B3"]
  
  g1_sIV2 <- coef(model_spatial_iv2)["fit_B1"]
  g2_sIV2 <- coef(model_spatial_iv2)["fit_B2"]
  g3_sIV2 <- coef(model_spatial_iv2)["fit_B3"]
  
  alpha_ols <- coef(model_gas)["p_sell.imp"]
  alpha_iv  <- coef(iv_gas)["fit_p_sell.imp"]
  
  # Show parameter differences
  cat("\n=== Parameter Comparison ===\n")
  cat("Alpha OLS:", round(alpha_ols, 4), "  Alpha IV:", round(alpha_iv, 4), "\n")
  cat("sOLS gammas:", round(g1_sOLS, 4), round(g2_sOLS, 4), round(g3_sOLS, 4), "\n")
  cat("sIV1 gammas:", round(g1_sIV1, 4), round(g2_sIV1, 4), round(g3_sIV1, 4), "\n")
  cat("sIV2 gammas:", round(g1_sIV2, 4), round(g2_sIV2, 4), round(g3_sIV2, 4), "\n\n")
  
  # Compute elasticities
  start_time <- Sys.time()
  
  elas_sOLS <- run_elasticity(g1_sOLS, g2_sOLS, g3_sOLS, alpha_ols, "Spatial OLS")
  elas_sIV1 <- run_elasticity(g1_sIV1, g2_sIV1, g3_sIV1, alpha_iv,  "Spatial IV1")
  elas_sIV2 <- run_elasticity(g1_sIV2, g2_sIV2, g3_sIV2, alpha_iv,  "Spatial IV2")
  
  end_time <- Sys.time()
  cat("Total computation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n\n")
  
  # Save results for future use
  save(elas_sOLS, elas_sIV1, elas_sIV2, file = results_file)
  cat("Results saved to:", results_file, "\n\n")
}

# ------------------------------------------------------------
# STEP 2: Pull elasticity stats
# ------------------------------------------------------------
med_own_logit   <- round(median(elasticity_table$mean_own_elasticity,   na.rm=TRUE), 4)
med_cross_logit <- round(median(elasticity_table$mean_cross_elasticity, na.rm=TRUE), 4)

med_own_sOLS  <- round(median(elas_sOLS$mean_own_elas,  na.rm=TRUE), 4)
med_b1_sOLS   <- round(median(elas_sOLS$mean_cross_b1,  na.rm=TRUE), 4)
med_b2_sOLS   <- round(median(elas_sOLS$mean_cross_b2,  na.rm=TRUE), 4)
med_b3_sOLS   <- round(median(elas_sOLS$mean_cross_b3,  na.rm=TRUE), 4)

med_own_sIV1  <- round(median(elas_sIV1$mean_own_elas,  na.rm=TRUE), 4)
med_b1_sIV1   <- round(median(elas_sIV1$mean_cross_b1,  na.rm=TRUE), 4)
med_b2_sIV1   <- round(median(elas_sIV1$mean_cross_b2,  na.rm=TRUE), 4)
med_b3_sIV1   <- round(median(elas_sIV1$mean_cross_b3,  na.rm=TRUE), 4)

med_own_sIV2  <- round(median(elas_sIV2$mean_own_elas,  na.rm=TRUE), 4)
med_b1_sIV2   <- round(median(elas_sIV2$mean_cross_b1,  na.rm=TRUE), 4)
med_b2_sIV2   <- round(median(elas_sIV2$mean_cross_b2,  na.rm=TRUE), 4)
med_b3_sIV2   <- round(median(elas_sIV2$mean_cross_b3,  na.rm=TRUE), 4)

# Show comparison
cat("=== Elasticity Medians Comparison ===\n")
cat("Own elasticity:  sOLS =", med_own_sOLS, " sIV1 =", med_own_sIV1, " sIV2 =", med_own_sIV2, "\n")
cat("Cross B1:        sOLS =", med_b1_sOLS,  " sIV1 =", med_b1_sIV1,  " sIV2 =", med_b1_sIV2, "\n")
cat("Cross B2:        sOLS =", med_b2_sOLS,  " sIV1 =", med_b2_sIV1,  " sIV2 =", med_b2_sIV2, "\n")
cat("Cross B3:        sOLS =", med_b3_sOLS,  " sIV1 =", med_b3_sIV1,  " sIV2 =", med_b3_sIV2, "\n\n")

f_iv_gas      <- round(fitstat(iv_gas,            "ivf")[[1]]$stat, 3)
f_spatial_iv1 <- round(fitstat(model_spatial_iv1, "ivf")[[1]]$stat, 3)
f_spatial_iv2 <- round(fitstat(model_spatial_iv2, "ivf")[[1]]$stat, 3)

# ------------------------------------------------------------
# STEP 3: Rename coefficients
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
# ------------------------------------------------------------
# Compute simple logit elasticities separately for OLS and IV
alpha_logit_ols <- coef(model_gas)["p_sell.imp"]
alpha_logit_iv  <- coef(iv_gas)["fit_p_sell.imp"]

# Simple formula: alpha * (1 - Sj) * Pj
gas_elas_base <- sample_demand %>%
  filter(prod == "gas") %>%
  select(cnpj, sjt, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(p_sell.imp))

med_own_logit_ols <- round(median(alpha_logit_ols * (1 - gas_elas_base$sjt) * 
                                    gas_elas_base$p_sell.imp, na.rm=TRUE), 4)

med_own_logit_iv  <- round(median(alpha_logit_iv  * (1 - gas_elas_base$sjt) * 
                                    gas_elas_base$p_sell.imp, na.rm=TRUE), 4)

cat("Median own elasticity - Logit OLS:", med_own_logit_ols, "\n")
cat("Median own elasticity - Logit IV: ", med_own_logit_iv,  "\n")

extra_rows <- data.frame(
  term = c(
    "theta (market size)",
    "Median own elasticity",
    "Median cross (< 500m)",
    "Median cross (500m-1.5km)",
    "Median cross (1.5km-3km)",
    "First stage F-stat"
  ),
  OLS  = c("1.500", med_own_logit_ols, "—",         "—",         "—",         "—"),
  IV   = c("1.500", med_own_logit_iv,  "—",         "—",         "—",         f_iv_gas),
  sOLS = c("1.500", med_own_sOLS,      med_b1_sOLS, med_b2_sOLS, med_b3_sOLS, "—"),
  IV1  = c("1.500", med_own_sIV1,      med_b1_sIV1, med_b2_sIV1, med_b3_sIV1, f_spatial_iv1),
  IV2  = c("1.500", med_own_sIV2,      med_b1_sIV2, med_b2_sIV2, med_b3_sIV2, f_spatial_iv2)
)
attr(extra_rows, "position") <- c(25, 26, 27, 28, 29, 30)

# ------------------------------------------------------------
# STEP 5: Build table
# ------------------------------------------------------------
modelsummary(
  list("OLS" = model_gas,
       "IV"  = iv_gas,
       "sOLS" = model_spatial_ols,
       "sIV1" = model_spatial_iv1,
       "sIV2" = model_spatial_iv2),
  coef_map   = coef_map,
  gof_map    = c("nobs", "r.squared"),
  add_rows   = extra_rows,
  stars      = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  notes      = "Standard errors clustered by month x neighborhood.",
  output     = "/proj/econ/duartema/GasBrazil/Final/Final_Table.txt"
)