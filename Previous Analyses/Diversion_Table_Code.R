# Diversion x Distance Table
# Diversion ratio D(i,j) = (-dsj/dpi) / (dsi/dpi)
# ============================================================

library(dplyr)

# Pull sOLS parameters (spatial OLS - all coefficients from model_spatial_ols)
g1 <- coef(model_spatial_ols)["B1"]
g2 <- coef(model_spatial_ols)["B2"]
g3 <- coef(model_spatial_ols)["B3"]
a  <- coef(model_spatial_ols)["p_sell.imp"]

cat("Parameters for diversion calculation (Spatial OLS):\n")
cat("Alpha:", round(a, 4), "\n")
cat("Gamma1:", round(g1, 4), "Gamma2:", round(g2, 4), "Gamma3:", round(g3, 4), "\n")

# ------------------------------------------------------------
# STEP 2: Function to compute diversion ratios for one month
# ------------------------------------------------------------
compute_diversion_month <- function(m, Dmat, gas_m, g1, g2, g3, a) {
  
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
  w4 <- (D >= 3.0)             # beyond 3km
  diag(w4) <- FALSE            # exclude self
  
  # Build dG/ds matrix (same as elasticity code with sign fix)
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
  
  # dsj/dpi matrix = alpha * M_inv
  # entry [j, i] = dsj/dpi
  deriv_mat <- a * M_inv
  
  # Own price derivative for each station i = diagonal
  own_deriv <- diag(deriv_mat)  # dsi/dpi for each i
  
  # Diversion ratio D(i,j) = (-dsj/dpi) / (dsi/dpi)
  # = -deriv_mat[j,i] / own_deriv[i]
  # We want: for each station i, average diversion to stations
  # in each distance ring
  
  # Diversion matrix: entry [i,j] = D(i,j)
  # = -deriv_mat[j,i] / own_deriv[i]
  # deriv_mat[j,i] is column i of deriv_mat
  # So diversion[i,j] = -deriv_mat[j,i] / own_deriv[i]
  # = t(-deriv_mat)[i,j] / own_deriv[i]
  
  div_mat <- -t(deriv_mat)                    # [i,j] = -dsj/dpi
  div_mat <- sweep(div_mat, 1, own_deriv, "/") # divide each row i by own_deriv[i]
  diag(div_mat) <- NA                          # exclude own
  
  # Number of stations in each ring for each station i
  n_w1 <- rowSums(w1)
  n_w2 <- rowSums(w2)
  n_w3 <- rowSums(w3)
  n_w4 <- rowSums(w4)
  
  # Mean diversion to each ring for each station i
  mean_div_b1 <- sapply(seq_len(n), function(i) {
    idx <- which(w1[i, ])
    if (length(idx) == 0) return(NA)
    mean(div_mat[i, idx], na.rm = TRUE)
  })
  
  mean_div_b2 <- sapply(seq_len(n), function(i) {
    idx <- which(w2[i, ])
    if (length(idx) == 0) return(NA)
    mean(div_mat[i, idx], na.rm = TRUE)
  })
  
  mean_div_b3 <- sapply(seq_len(n), function(i) {
    idx <- which(w3[i, ])
    if (length(idx) == 0) return(NA)
    mean(div_mat[i, idx], na.rm = TRUE)
  })
  
  mean_div_b4 <- sapply(seq_len(n), function(i) {
    idx <- which(w4[i, ])
    if (length(idx) == 0) return(NA)
    mean(div_mat[i, idx], na.rm = TRUE)
  })
  
  # Sum diversion to each ring (total share diverted)
  sum_div_b1 <- sapply(seq_len(n), function(i) {
    idx <- which(w1[i, ])
    if (length(idx) == 0) return(NA)
    sum(div_mat[i, idx], na.rm = TRUE)
  })
  
  sum_div_b2 <- sapply(seq_len(n), function(i) {
    idx <- which(w2[i, ])
    if (length(idx) == 0) return(NA)
    sum(div_mat[i, idx], na.rm = TRUE)
  })
  
  sum_div_b3 <- sapply(seq_len(n), function(i) {
    idx <- which(w3[i, ])
    if (length(idx) == 0) return(NA)
    sum(div_mat[i, idx], na.rm = TRUE)
  })
  
  sum_div_b4 <- sapply(seq_len(n), function(i) {
    idx <- which(w4[i, ])
    if (length(idx) == 0) return(NA)
    sum(div_mat[i, idx], na.rm = TRUE)
  })
  
  data.frame(
    id          = gas_m$id,
    cnpj        = gas_m$cnpj,
    month       = m,
    n_b1        = n_w1,
    n_b2        = n_w2,
    n_b3        = n_w3,
    n_b4        = n_w4,
    mean_div_b1 = mean_div_b1,
    mean_div_b2 = mean_div_b2,
    mean_div_b3 = mean_div_b3,
    mean_div_b4 = mean_div_b4,
    sum_div_b1  = sum_div_b1,
    sum_div_b2  = sum_div_b2,
    sum_div_b3  = sum_div_b3,
    sum_div_b4  = sum_div_b4
  )
}

# ------------------------------------------------------------
# STEP 3: Loop over months
# ------------------------------------------------------------
gas_data_sp <- sample_demand %>%
  filter(prod == "gas") %>%
  select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))

months_list  <- unique(gas_data_sp$month)
div_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m     <- months_list[i]
  gas_m <- gas_data_sp %>% filter(month == m)
  Dmat  <- Geo$Dmatrix[[i]]
  div_list[[i]] <- compute_diversion_month(m, Dmat, gas_m, g1, g2, g3, a)
  rm(Dmat); gc()
}

div_df <- bind_rows(div_list)

# ------------------------------------------------------------
# STEP 4: Average per station across months
# ------------------------------------------------------------
div_avg <- div_df %>%
  group_by(cnpj) %>%
  summarise(
    mean_n_b1       = mean(n_b1,        na.rm = TRUE),
    mean_n_b2       = mean(n_b2,        na.rm = TRUE),
    mean_n_b3       = mean(n_b3,        na.rm = TRUE),
    mean_n_b4       = mean(n_b4,        na.rm = TRUE),
    mean_div_b1     = mean(mean_div_b1, na.rm = TRUE),
    mean_div_b2     = mean(mean_div_b2, na.rm = TRUE),
    mean_div_b3     = mean(mean_div_b3, na.rm = TRUE),
    mean_div_b4     = mean(mean_div_b4, na.rm = TRUE),
    mean_sum_div_b1 = mean(sum_div_b1,  na.rm = TRUE),
    mean_sum_div_b2 = mean(sum_div_b2,  na.rm = TRUE),
    mean_sum_div_b3 = mean(sum_div_b3,  na.rm = TRUE),
    mean_sum_div_b4 = mean(sum_div_b4,  na.rm = TRUE)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# STEP 5: Print and save summary table
# ------------------------------------------------------------
print_diversion_table <- function() {
  cat("\n=== Diversion x Distance Table ===\n")
  cat(sprintf("%-30s %10s %12s %10s %10s\n",
              "", "< 0.5km", "0.5-1.5km", "1.5-3km", "> 3km"))
  cat(strrep("-", 65), "\n")
  cat(sprintf("%-30s %10s %12s %10s %10s\n",
              "Number of stations",
              paste0(round(mean(div_avg$mean_n_b1, na.rm=T),1),
                     " (", round(sd(div_avg$mean_n_b1, na.rm=T),1), ")"),
              paste0(round(mean(div_avg$mean_n_b2, na.rm=T),1),
                     " (", round(sd(div_avg$mean_n_b2, na.rm=T),1), ")"),
              paste0(round(mean(div_avg$mean_n_b3, na.rm=T),1),
                     " (", round(sd(div_avg$mean_n_b3, na.rm=T),1), ")"),
              paste0(round(mean(div_avg$mean_n_b4, na.rm=T),1),
                     " (", round(sd(div_avg$mean_n_b4, na.rm=T),1), ")")))
  cat(sprintf("%-30s %10s %12s %10s %10s\n",
              "Median Cross-Elasticity %",
              round(median(div_avg$mean_div_b1, na.rm=T)*100, 3),
              round(median(div_avg$mean_div_b2, na.rm=T)*100, 3),
              round(median(div_avg$mean_div_b3, na.rm=T)*100, 3),
              round(median(div_avg$mean_div_b4, na.rm=T)*100, 3)))
  cat(sprintf("%-30s %10s %12s %10s %10s\n",
              "Mean Diversion %",
              paste0(round(mean(div_avg$mean_div_b1, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_div_b1, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_div_b2, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_div_b2, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_div_b3, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_div_b3, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_div_b4, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_div_b4, na.rm=T)*100, 1), ")")))
  cat(sprintf("%-30s %10s %12s %10s %10s\n",
              "Mean Diversion Sum %",
              paste0(round(mean(div_avg$mean_sum_div_b1, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_sum_div_b1, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_sum_div_b2, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_sum_div_b2, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_sum_div_b3, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_sum_div_b3, na.rm=T)*100, 1), ")"),
              paste0(round(mean(div_avg$mean_sum_div_b4, na.rm=T)*100, 1),
                     " (", round(sd(div_avg$mean_sum_div_b4, na.rm=T)*100, 1), ")")))
}

# Print to console
print_diversion_table()

# Save to file
sink('/proj/econ/duartema/GasBrazil/Final/Diversion_Table.txt')
print_diversion_table()
sink()