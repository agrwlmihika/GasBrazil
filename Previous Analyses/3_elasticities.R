# ============================================================
# 3. Spatial Elasticities — Structured by Derivation Steps
# ============================================================
load('/proj/econ/duartema/GasBrazil/Final/Data/regression_output.RData')
load('Data/Geo.RData')
library(dplyr)

# ------------------------------------------------------------
# Parameters from sIV1
# ------------------------------------------------------------
alpha  <- coef(model_sIV1)["fit_p_sell.imp"]
theta1 <- coef(model_sIV1)["B1"]
theta2 <- coef(model_sIV1)["B2"]
theta3 <- coef(model_sIV1)["B3"]

cat("Alpha:", round(alpha,4),
    " G1:", round(theta1,4),
    " G2:", round(theta2,4),
    " G3:", round(theta3,4), "\n")

# Prepare gas data
gas_data <- sample_spatial %>%
  select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))

months_list <- unique(gas_data$month)

# ============================================================
# STEP 1: dG/dS matrix for each month — save to Data
# ============================================================
dGdS_list <- vector("list", length(months_list))
names(dGdS_list) <- months_list

for (i in seq_along(months_list)) {
  m     <- months_list[i]
  gas_m <- gas_data %>% filter(month == m)
  Dmat  <- Geo$Dmatrix[[i]]
  ids   <- gas_m$id
  keep  <- ids[ids %in% rownames(Dmat)]
  if (length(keep) < 2) { rm(Dmat); gc(); next }
  gas_m <- gas_m %>% filter(id %in% keep)
  n     <- nrow(gas_m)
  S     <- gas_m$sjt
  s0    <- unique(gas_m$s0t)[1]
  D     <- Dmat[gas_m$id, gas_m$id]
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)
  
  # Build dG/dS matrix
  S_sum   <- outer(S, S, "+")
  theta_w <- theta1*w1 + theta2*w2 + theta3*w3
  dG_ds   <- -theta_w / S_sum
  
  for (j in seq_len(n)) {
    sk_sum     <- S/(S[j]+S); sk_sum[j] <- 0
    dG_ds[j,j] <- (1 + sum((theta1*w1[j,]+theta2*w2[j,]+theta3*w3[j,])*sk_sum)) / S[j]
  }
  #NOTE: DOUBLE CHECK THE SIGN OF THE ABOVE TERM
  
  rownames(dG_ds) <- gas_m$cnpj
  colnames(dG_ds) <- gas_m$cnpj
  dGdS_list[[i]]  <- dG_ds
  rm(Dmat); gc()
}

save(dGdS_list, file="/proj/econ/duartema/GasBrazil/Final/Data/dGdS_matrix.RData")
cat("Step 1 complete: dG/dS matrices saved\n")

# ============================================================
# STEP 2: dS/dP matrix for each month — save to Data
# dS/dP = alpha * [dG/dS + (1/s0)*I]^{-1}
# ============================================================
dSdP_list  <- vector("list", length(months_list))
names(dSdP_list) <- months_list

for (i in seq_along(months_list)) {
  m     <- months_list[i]
  gas_m <- gas_data %>% filter(month == m)
  dG    <- dGdS_list[[i]]
  if (is.null(dG)) next
  
  # Match stations to gas_m
  keep  <- gas_m$id[gas_m$id %in% rownames(
    Geo$Dmatrix[[i]])]
  gas_m <- gas_m %>% filter(id %in% keep)
  s0    <- unique(gas_m$s0t)[1]
  
  # Add 1/s0 to diagonal and invert
  M       <- dG
  diag(M) <- diag(M) + 1/s0
  M_inv   <- tryCatch(solve(M), error=function(e) NULL)
  if (is.null(M_inv)) next
  
  dSdP <- alpha * M_inv
  rownames(dSdP) <- rownames(dG)
  colnames(dSdP) <- colnames(dG)
  dSdP_list[[i]] <- dSdP
}

save(dSdP_list, file="/proj/econ/duartema/GasBrazil/Final/Data/dSdP_matrix.RData")
cat("Step 2 complete: dS/dP matrices saved\n")

# ============================================================
# STEP 3: Elasticity matrix for each month
# E[j,k] = dSj/dPk * Pk/Sj
# ============================================================
elas_matrix_list <- vector("list", length(months_list))
names(elas_matrix_list) <- months_list
elas_summary_list <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m     <- months_list[i]
  gas_m <- gas_data %>% filter(month == m)
  dSdP  <- dSdP_list[[i]]
  if (is.null(dSdP)) next
  
  # FIX: rownames(dSdP) are cnpj, so match on cnpj (not id)
  keep  <- gas_m$cnpj[gas_m$cnpj %in% rownames(dSdP)]
  gas_m <- gas_m %>% filter(cnpj %in% keep)
  if (nrow(gas_m) < 2) next   # FIX: removed rm(Dmat) — not yet defined here
  
  # Reorder gas_m rows to align with dSdP row order
  gas_m  <- gas_m[match(rownames(dSdP)[rownames(dSdP) %in% gas_m$cnpj],
                        gas_m$cnpj), ]
  gas_m  <- gas_m[!is.na(gas_m$cnpj), ]
  dSdP_m <- dSdP[gas_m$cnpj, gas_m$cnpj]   # subset & align dSdP
  
  S <- gas_m$sjt
  P <- gas_m$p_sell.imp
  
  elas_mat <- sweep(dSdP_m, 2, P, "*")
  elas_mat <- sweep(elas_mat, 1, S, "/")
  rownames(elas_mat) <- gas_m$cnpj
  colnames(elas_mat) <- gas_m$cnpj
  elas_matrix_list[[i]] <- elas_mat
  
  # Distance rings  (FIX: removed duplicate block)
  Dmat <- Geo$Dmatrix[[i]]
  D    <- Dmat[gas_m$id, gas_m$id]
  w1   <- (D > 0)    & (D < 0.5)
  w2   <- (D >= 0.5) & (D < 1.5)
  w3   <- (D >= 1.5) & (D < 3.0)
  n    <- nrow(gas_m)
  
  elas_summary_list[[i]] <- data.frame(
    cnpj  = gas_m$cnpj,
    month = m,
    own_elas      = diag(elas_mat),
    cross_elas_b1 = sapply(seq_len(n), function(j) {
      idx <- which(w1[j,] & seq_len(n) != j)
      if (!length(idx)) return(NA); mean(elas_mat[j, idx], na.rm = TRUE) }),
    cross_elas_b2 = sapply(seq_len(n), function(j) {
      idx <- which(w2[j,] & seq_len(n) != j)
      if (!length(idx)) return(NA); mean(elas_mat[j, idx], na.rm = TRUE) }),
    cross_elas_b3 = sapply(seq_len(n), function(j) {
      idx <- which(w3[j,] & seq_len(n) != j)
      if (!length(idx)) return(NA); mean(elas_mat[j, idx], na.rm = TRUE) })
  )
  rm(Dmat); gc()
}

elas_summary <- bind_rows(elas_summary_list)

# Average per station
elas_avg <- elas_summary %>%
  group_by(cnpj) %>%
  summarise(
    mean_own   = mean(own_elas,      na.rm=TRUE),
    mean_b1    = mean(cross_elas_b1, na.rm=TRUE),
    mean_b2    = mean(cross_elas_b2, na.rm=TRUE),
    mean_b3    = mean(cross_elas_b3, na.rm=TRUE)
  ) %>% ungroup()

cat("\n=== Elasticity Medians ===\n")
cat("Own:     ", round(median(elas_avg$mean_own, na.rm=TRUE), 4), "\n")
cat("Cross B1:", round(median(elas_avg$mean_b1,  na.rm=TRUE), 4), "\n")
cat("Cross B2:", round(median(elas_avg$mean_b2,  na.rm=TRUE), 4), "\n")
cat("Cross B3:", round(median(elas_avg$mean_b3,  na.rm=TRUE), 4), "\n")

save(elas_matrix_list, elas_summary, elas_avg,
     file="/proj/econ/duartema/GasBrazil/Final/Data/elasticity_matrix.RData")
cat("Step 3 complete: elasticity matrices saved\n")

# ============================================================
# STEP 4: Diversion table
# D(i,j) = (-dSj/dPi) / (dSi/dPi)
# ============================================================
div_list <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m     <- months_list[i]
  gas_m <- gas_data %>% filter(month == m)
  dSdP  <- dSdP_list[[i]]
  if (is.null(dSdP)) next
  
  # FIX: rownames of dSdP are cnpj (set in STEP 1), not id
  keep  <- gas_m$cnpj[gas_m$cnpj %in% rownames(dSdP)]
  gas_m <- gas_m %>% filter(cnpj %in% keep)
  if (nrow(gas_m) < 2) next
  
  # FIX: reorder gas_m rows to match dSdP row order
  # (so that D and div_mat are aligned on the same stations)
  gas_m <- gas_m[match(rownames(dSdP), gas_m$cnpj), ]
  if (any(is.na(gas_m$cnpj))) next
  
  n <- nrow(gas_m)
  
  Dmat <- Geo$Dmatrix[[i]]
  D    <- Dmat[gas_m$id, gas_m$id]
  w1   <- (D > 0)    & (D < 0.5)
  w2   <- (D >= 0.5) & (D < 1.5)
  w3   <- (D >= 1.5) & (D < 3.0)
  w4   <- (D >= 3.0); diag(w4) <- FALSE
  
  own_deriv     <- diag(dSdP)
  div_mat       <- -t(dSdP)
  div_mat       <- sweep(div_mat, 1, own_deriv, "/")
  diag(div_mat) <- NA
  
  ring_stats <- function(w) {
    list(
      n    = rowSums(w),
      mean = sapply(seq_len(n), function(i) {
        idx <- which(w[i,]); if (!length(idx)) return(NA)
        mean(div_mat[i,idx], na.rm=TRUE) }),
      sum  = sapply(seq_len(n), function(i) {
        idx <- which(w[i,]); if (!length(idx)) return(NA)
        sum(div_mat[i,idx], na.rm=TRUE) })
    )
  }
  
  r1 <- ring_stats(w1); r2 <- ring_stats(w2)
  r3 <- ring_stats(w3); r4 <- ring_stats(w4)
  
  div_list[[i]] <- data.frame(
    cnpj=gas_m$cnpj, month=m,
    n_b1=r1$n, n_b2=r2$n, n_b3=r3$n, n_b4=r4$n,
    mean_div_b1=r1$mean, mean_div_b2=r2$mean,
    mean_div_b3=r3$mean, mean_div_b4=r4$mean,
    sum_div_b1=r1$sum,   sum_div_b2=r2$sum,
    sum_div_b3=r3$sum,   sum_div_b4=r4$sum
  )
  rm(Dmat); gc()
}

div_avg <- bind_rows(div_list) %>%
  group_by(cnpj) %>%
  summarise(across(starts_with("mean_") | starts_with("sum_") |
                     starts_with("n_"), ~mean(.x, na.rm=TRUE))) %>%
  ungroup()

print_div_table <- function(d) {
  cat("\n=== Diversion x Distance Table ===\n")
  cat(sprintf("%-30s %10s %12s %10s %10s\n","","<0.5km","0.5-1.5km","1.5-3km",">3km"))
  cat(strrep("-", 65), "\n")
  fmt <- function(x) paste0(round(mean(x,na.rm=T),1)," (",round(sd(x,na.rm=T),1),")")
  cat(sprintf("%-30s %10s %12s %10s %10s\n","Number of stations",
              fmt(d$n_b1),fmt(d$n_b2),fmt(d$n_b3),fmt(d$n_b4)))
  cat(sprintf("%-30s %10s %12s %10s %10s\n","Median Cross-Elasticity %",
              round(median(d$mean_div_b1,na.rm=T)*100,3),
              round(median(d$mean_div_b2,na.rm=T)*100,3),
              round(median(d$mean_div_b3,na.rm=T)*100,3),
              round(median(d$mean_div_b4,na.rm=T)*100,3)))
  cat(sprintf("%-30s %10s %12s %10s %10s\n","Mean Diversion %",
              fmt(d$mean_div_b1*100),fmt(d$mean_div_b2*100),
              fmt(d$mean_div_b3*100),fmt(d$mean_div_b4*100)))
  cat(sprintf("%-30s %10s %12s %10s %10s\n","Mean Diversion Sum %",
              fmt(d$sum_div_b1*100),fmt(d$sum_div_b2*100),
              fmt(d$sum_div_b3*100),fmt(d$sum_div_b4*100)))
}

print_div_table(div_avg)
sink('/proj/econ/duartema/GasBrazil/Final/Output/table3_diversion.txt')
print_div_table(div_avg)
sink()

save(div_avg, file="/proj/econ/duartema/GasBrazil/Final/Data/diversion_table.RData")
cat("Step 4 complete: diversion table saved\n")

# ============================================================
# STEP 5: Positive Semi-Definite Test on dG/dS matrices
# For consumer theory consistency, dG/dS must be PSD
# (stations must be substitutes, not complements)
# A matrix is PSD if all eigenvalues >= 0
# ============================================================
psd_results <- data.frame(
  month          = months_list,
  n_stations     = NA_integer_,
  min_eigenvalue = NA_real_,
  is_PSD         = NA,
  pct_positive_offdiag = NA_real_
)

for (i in seq_along(months_list)) {
  dG <- dGdS_list[[i]]
  if (is.null(dG)) next
  
  # Eigenvalues — all must be >= 0 for PSD
  eigs <- eigen(dG, symmetric=TRUE, only.values=TRUE)$values
  min_eig <- min(eigs)
  
  # Check off-diagonal signs — positive = substitutes
  off_diag <- dG[lower.tri(dG)]
  pct_pos  <- mean(off_diag > 0, na.rm=TRUE) * 100
  
  psd_results$n_stations[i]          <- nrow(dG)
  psd_results$min_eigenvalue[i]      <- round(min_eig, 6)
  psd_results$is_PSD[i]              <- min_eig >= -1e-8  # tolerance for numerical error
  psd_results$pct_positive_offdiag[i] <- round(pct_pos, 1)
}

cat("\n=== Positive Semi-Definite Test: dG/dS Matrix ===\n")
cat("Total months tested:", sum(!is.na(psd_results$min_eigenvalue)), "\n")
cat("Months passing PSD:", sum(psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Months failing PSD:", sum(!psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Min eigenvalue overall:", round(min(psd_results$min_eigenvalue, na.rm=TRUE), 6), "\n")
cat("% off-diagonal positive (substitutes) — median across months:",
    round(median(psd_results$pct_positive_offdiag, na.rm=TRUE), 1), "%\n")

cat("\n=== PSD Summary Table (by month) ===\n")
print(psd_results, row.names=FALSE)

# Save PSD results
sink('/proj/econ/duartema/GasBrazil/Final/Output/table4_psd_test.txt')
cat("=== Positive Semi-Definite Test: dG/dS Matrix ===\n")
cat("Months passing PSD:", sum(psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Months failing PSD:", sum(!psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Min eigenvalue overall:", round(min(psd_results$min_eigenvalue, na.rm=TRUE), 6), "\n")
cat("% off-diagonal positive (substitutes) — median:",
    round(median(psd_results$pct_positive_offdiag, na.rm=TRUE), 1), "%\n\n")
print(psd_results, row.names=FALSE)
sink()

save(psd_results, file="/proj/econ/duartema/GasBrazil/Final/Data/psd_test.RData")
cat("Step 5 complete: PSD test saved to Output/table4_psd_test.txt\n")