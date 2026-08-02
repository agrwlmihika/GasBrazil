# ============================================================
# 3. Spatial Elasticities — Structured by Derivation Steps
# ============================================================
load('/proj/econ/duartema/GasBrazil/Final/Data/regression_output.RData')
load('Data/Geo.RData')
library(dplyr)

# ------------------------------------------------------------
# Parameters from sOLS and sIV1
# ------------------------------------------------------------
alpha_ols  <- coef(model_sOLS)["p_sell.imp"]
theta1_ols <- coef(model_sOLS)["B1"]
theta2_ols <- coef(model_sOLS)["B2"]
theta3_ols <- coef(model_sOLS)["B3"]

alpha_iv1  <- coef(model_sIV1)["fit_p_sell.imp"]
theta1_iv1 <- coef(model_sIV1)["B1"]
theta2_iv1 <- coef(model_sIV1)["B2"]
theta3_iv1 <- coef(model_sIV1)["B3"]

cat("sOLS — Alpha:", round(alpha_ols,4),
    " G1:", round(theta1_ols,4),
    " G2:", round(theta2_ols,4),
    " G3:", round(theta3_ols,4), "\n")
cat("sIV1 — Alpha:", round(alpha_iv1,4),
    " G1:", round(theta1_iv1,4),
    " G2:", round(theta2_iv1,4),
    " G3:", round(theta3_iv1,4), "\n")

# ------------------------------------------------------------
# Prepare gas data
# ------------------------------------------------------------
gas_data <- sample_spatial %>%
  select(id, month, cnpj, sjt, s0t, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(s0t), !is.na(p_sell.imp))

months_list <- unique(gas_data$month)

# ============================================================
# FUNCTIONS
# ============================================================

# STEP 1: dG/dS matrix
compute_dGdS <- function(gas_data, months_list, t1, t2, t3) {
  out <- vector("list", length(months_list))
  names(out) <- months_list
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
    D     <- Dmat[gas_m$id, gas_m$id]
    w1 <- (D > 0)    & (D < 0.5)
    w2 <- (D >= 0.5) & (D < 1.5)
    w3 <- (D >= 1.5) & (D < 3.0)
    S_sum   <- outer(S, S, "+")
    theta_w <- t1*w1 + t2*w2 + t3*w3
    dG_ds   <- -theta_w / S_sum
    for (j in seq_len(n)) {
      sk_sum     <- S/(S[j]+S); sk_sum[j] <- 0
      dG_ds[j,j] <- (1 + sum((t1*w1[j,]+t2*w2[j,]+t3*w3[j,])*sk_sum)) / S[j]
    }
    rownames(dG_ds) <- gas_m$cnpj
    colnames(dG_ds) <- gas_m$cnpj
    out[[i]] <- dG_ds
    rm(Dmat); gc()
  }
  out
}

# STEP 2: dS/dP matrix
compute_dSdP <- function(gas_data, months_list, dGdS_list, alpha) {
  out <- vector("list", length(months_list))
  names(out) <- months_list
  for (i in seq_along(months_list)) {
    m     <- months_list[i]
    gas_m <- gas_data %>% filter(month == m)
    dG    <- dGdS_list[[i]]
    if (is.null(dG)) next
    keep  <- gas_m$id[gas_m$id %in% rownames(Geo$Dmatrix[[i]])]
    gas_m <- gas_m %>% filter(id %in% keep)
    s0    <- unique(gas_m$s0t)[1]
    M     <- dG; diag(M) <- diag(M) + 1/s0
    M_inv <- tryCatch(solve(M), error=function(e) NULL)
    if (is.null(M_inv)) next
    dSdP  <- alpha * M_inv
    rownames(dSdP) <- rownames(dG)
    colnames(dSdP) <- colnames(dG)
    out[[i]] <- dSdP
  }
  out
}

# STEP 3: Elasticity matrix + summary
compute_elas <- function(gas_data, months_list, dSdP_list) {
  elas_matrix_list  <- vector("list", length(months_list))
  elas_summary_list <- vector("list", length(months_list))
  names(elas_matrix_list) <- months_list
  for (i in seq_along(months_list)) {
    m     <- months_list[i]
    gas_m <- gas_data %>% filter(month == m)
    dSdP  <- dSdP_list[[i]]
    if (is.null(dSdP)) next
    keep  <- gas_m$cnpj[gas_m$cnpj %in% rownames(dSdP)]
    gas_m <- gas_m %>% filter(cnpj %in% keep)
    if (nrow(gas_m) < 2) next
    gas_m  <- gas_m[match(rownames(dSdP)[rownames(dSdP) %in% gas_m$cnpj],
                          gas_m$cnpj), ]
    gas_m  <- gas_m[!is.na(gas_m$cnpj), ]
    dSdP_m <- dSdP[gas_m$cnpj, gas_m$cnpj]
    S <- gas_m$sjt; P <- gas_m$p_sell.imp
    elas_mat <- sweep(dSdP_m, 2, P, "*")
    elas_mat <- sweep(elas_mat, 1, S, "/")
    rownames(elas_mat) <- gas_m$cnpj
    colnames(elas_mat) <- gas_m$cnpj
    elas_matrix_list[[i]] <- elas_mat
    Dmat <- Geo$Dmatrix[[i]]
    D    <- Dmat[gas_m$id, gas_m$id]
    w1 <- (D > 0) & (D < 0.5)
    w2 <- (D >= 0.5) & (D < 1.5)
    w3 <- (D >= 1.5) & (D < 3.0)
    n  <- nrow(gas_m)
    elas_summary_list[[i]] <- data.frame(
      cnpj=gas_m$cnpj, month=m,
      own_elas      = diag(elas_mat),
      cross_elas_b1 = sapply(seq_len(n), function(j) {
        idx <- which(w1[j,] & seq_len(n)!=j)
        if (!length(idx)) return(NA); mean(elas_mat[j,idx], na.rm=TRUE) }),
      cross_elas_b2 = sapply(seq_len(n), function(j) {
        idx <- which(w2[j,] & seq_len(n)!=j)
        if (!length(idx)) return(NA); mean(elas_mat[j,idx], na.rm=TRUE) }),
      cross_elas_b3 = sapply(seq_len(n), function(j) {
        idx <- which(w3[j,] & seq_len(n)!=j)
        if (!length(idx)) return(NA); mean(elas_mat[j,idx], na.rm=TRUE) })
    )
    rm(Dmat); gc()
  }
  list(matrices=elas_matrix_list,
       summary=bind_rows(elas_summary_list))
}

# STEP 4: Diversion table
compute_diversion <- function(gas_data, months_list, dSdP_list) {
  div_list <- vector("list", length(months_list))
  for (i in seq_along(months_list)) {
    m     <- months_list[i]
    gas_m <- gas_data %>% filter(month == m)
    dSdP  <- dSdP_list[[i]]
    if (is.null(dSdP)) next
    keep  <- gas_m$cnpj[gas_m$cnpj %in% rownames(dSdP)]
    gas_m <- gas_m %>% filter(cnpj %in% keep)
    if (nrow(gas_m) < 2) next
    gas_m <- gas_m[match(rownames(dSdP), gas_m$cnpj), ]
    if (any(is.na(gas_m$cnpj))) next
    n    <- nrow(gas_m)
    Dmat <- Geo$Dmatrix[[i]]
    D    <- Dmat[gas_m$id, gas_m$id]
    w1 <- (D > 0) & (D < 0.5)
    w2 <- (D >= 0.5) & (D < 1.5)
    w3 <- (D >= 1.5) & (D < 3.0)
    w4 <- (D >= 3.0); diag(w4) <- FALSE
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
    r1<-ring_stats(w1); r2<-ring_stats(w2)
    r3<-ring_stats(w3); r4<-ring_stats(w4)
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
  bind_rows(div_list) %>%
    group_by(cnpj) %>%
    summarise(across(starts_with("mean_") | starts_with("sum_") |
                       starts_with("n_"), ~mean(.x, na.rm=TRUE))) %>%
    ungroup()
}

# ============================================================
# RUN FOR sOLS
# ============================================================
cat("\n--- Running sOLS ---\n")
dGdS_ols <- compute_dGdS(gas_data, months_list, theta1_ols, theta2_ols, theta3_ols)
dSdP_ols <- compute_dSdP(gas_data, months_list, dGdS_ols, alpha_ols)
elas_ols  <- compute_elas(gas_data, months_list, dSdP_ols)
div_ols   <- compute_diversion(gas_data, months_list, dSdP_ols)
elas_summary_ols <- elas_ols$summary

elas_avg_ols <- elas_ols$summary %>%
  group_by(cnpj) %>%
  summarise(median_own=median(own_elas,na.rm=TRUE),
            median_b1=median(cross_elas_b1,na.rm=TRUE),
            median_b2=median(cross_elas_b2,na.rm=TRUE),
            median_b3=median(cross_elas_b3,na.rm=TRUE)) %>% ungroup()
cat("\n=== sOLS Elasticity Medians ===\n")
cat("Own:     ", round(median(elas_avg_ols$median_own, na.rm=TRUE), 4), "\n")
cat("Cross B1:", round(median(elas_avg_ols$median_b1,  na.rm=TRUE), 4), "\n")
cat("Cross B2:", round(median(elas_avg_ols$median_b2,  na.rm=TRUE), 4), "\n")
cat("Cross B3:", round(median(elas_avg_ols$median_b3,  na.rm=TRUE), 4), "\n")

# ============================================================
# RUN FOR sIV1
# ============================================================
cat("\n--- Running sIV1 ---\n")
dGdS_iv1 <- compute_dGdS(gas_data, months_list, theta1_iv1, theta2_iv1, theta3_iv1)
dSdP_iv1 <- compute_dSdP(gas_data, months_list, dGdS_iv1, alpha_iv1)
elas_iv1  <- compute_elas(gas_data, months_list, dSdP_iv1)
div_iv1   <- compute_diversion(gas_data, months_list, dSdP_iv1)
elas_summary_iv1 <- elas_iv1$summary

elas_avg_iv1 <- elas_iv1$summary %>%
  group_by(cnpj) %>%
  summarise(median_own=median(own_elas,na.rm=TRUE),
            median_b1=median(cross_elas_b1,na.rm=TRUE),
            median_b2=median(cross_elas_b2,na.rm=TRUE),
            median_b3=median(cross_elas_b3,na.rm=TRUE)) %>% ungroup()
cat("\n=== sIV1 Elasticity Medians ===\n")
cat("Own:     ", round(median(elas_avg_iv1$median_own, na.rm=TRUE), 4), "\n")
cat("Cross B1:", round(median(elas_avg_iv1$median_b1,  na.rm=TRUE), 4), "\n")
cat("Cross B2:", round(median(elas_avg_iv1$median_b2,  na.rm=TRUE), 4), "\n")
cat("Cross B3:", round(median(elas_avg_iv1$median_b3,  na.rm=TRUE), 4), "\n")
# ============================================================
# STEP 5: PSD Test — run on sIV1 dG/dS (preferred spec)
# ============================================================
cat("\n--- PSD Test (sIV1) ---\n")
psd_results <- data.frame(
  month                = months_list,
  n_stations           = NA_integer_,
  min_eigenvalue       = NA_real_,
  is_PSD               = NA,
  pct_positive_offdiag = NA_real_
)

for (i in seq_along(months_list)) {
  dG <- dGdS_iv1[[i]]
  if (is.null(dG)) next
  eigs    <- eigen(dG, symmetric=TRUE, only.values=TRUE)$values
  min_eig <- min(eigs)
  off_diag <- dG[lower.tri(dG)]
  pct_pos  <- mean(off_diag > 0, na.rm=TRUE) * 100
  psd_results$n_stations[i]           <- nrow(dG)
  psd_results$min_eigenvalue[i]       <- round(min_eig, 6)
  psd_results$is_PSD[i]               <- min_eig >= -1e-8
  psd_results$pct_positive_offdiag[i] <- round(pct_pos, 1)
}

cat("Months passing PSD:", sum(psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Months failing PSD:", sum(!psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Min eigenvalue overall:", round(min(psd_results$min_eigenvalue, na.rm=TRUE), 6), "\n")
cat("% off-diagonal positive — median:",
    round(median(psd_results$pct_positive_offdiag, na.rm=TRUE), 1), "%\n")
print(psd_results, row.names=FALSE)

# ============================================================
# PRINT DIVERSION TABLE
# ============================================================
print_div_table <- function(d, label) {
  cat("\n=== Diversion x Distance Table —", label, "===\n")
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

print_div_table(div_iv1,  "sIV1")

# ============================================================
# SAVE ALL
# ============================================================
save(dGdS_ols, dSdP_ols, elas_ols, elas_avg_ols, div_ols,
     dGdS_iv1, dSdP_iv1, elas_iv1, elas_avg_iv1, div_iv1,
     psd_results,
     file="/proj/econ/duartema/GasBrazil/Final/Data/elasticity_matrix.RData")

sink('/proj/econ/duartema/GasBrazil/Final/Output/table3_diversion.txt')
print_div_table(div_iv1, "sIV1")
sink()

sink('/proj/econ/duartema/GasBrazil/Final/Output/table4_psd_test.txt')
cat("=== PSD Test (sIV1) ===\n")
cat("Months passing:", sum(psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Months failing:", sum(!psd_results$is_PSD, na.rm=TRUE), "\n")
cat("Min eigenvalue:", round(min(psd_results$min_eigenvalue, na.rm=TRUE), 6), "\n")
print(psd_results, row.names=FALSE)
sink()

saveRDS(elas_summary_iv1, "/proj/econ/duartema/GasBrazil/Final/Data/elas_summary_iv1.rds")
saveRDS(elas_summary_ols, "/proj/econ/duartema/GasBrazil/Final/Data/elas_summary_ols.rds")

cat("\nAll results saved.\n")