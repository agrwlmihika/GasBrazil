# Elasticity Calculation - Gasoline Only
# Own-price elasticity:  dSj/Pj * Pj/Sj = alpha * (1 - Sj) * Pj
# Cross-price elasticity: dSj/Pk * Pk/Sj = -alpha * Sk * Pk
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# STEP 1: Extract alpha from OLS gasoline model
# ------------------------------------------------------------
alpha <- coef(iv_gas)["fit_p_sell.imp"]
cat("Alpha (price coefficient):", round(alpha, 4), "\n")

# ------------------------------------------------------------
# STEP 2: Prepare gasoline data with shares and prices
# ------------------------------------------------------------
gas_data <- sample_demand %>%
  filter(prod == "gas") %>%
  select(month, cnpj, sjt, p_sell.imp) %>%
  filter(!is.na(sjt), !is.na(p_sell.imp))

# ------------------------------------------------------------
# STEP 3: Own-price elasticity
# dSj/Pj * Pj/Sj = alpha * (1 - Sj) * Pj
# This is one number per station per month — no big matrix needed
# ------------------------------------------------------------
own_elas <- gas_data %>%
  mutate(own_elasticity = alpha * (1 - sjt) * p_sell.imp)

# Average own-price elasticity per station across months
own_elas_avg <- own_elas %>%
  group_by(cnpj) %>%
  summarise(mean_own_elasticity = mean(own_elasticity, na.rm = TRUE)) %>%
  ungroup()

cat("\n=== Own-Price Elasticity (avg across months) ===\n")
summary(own_elas_avg$mean_own_elasticity)

# ------------------------------------------------------------
# STEP 4: Cross-price elasticity — handle memory carefully
# dSj/Pk * Pk/Sj = -alpha * Sk * Pk
# This produces a J x J matrix PER month (300+ x 300+ = 90,000+ cells)
# We compute month by month and immediately summarise to avoid storing
# the full matrix in memory
# ------------------------------------------------------------

# Get list of months
months <- unique(gas_data$month)

# Storage for average cross elasticity per station pair
cross_elas_list <- vector("list", length(months))

for (i in seq_along(months)) {
  
  m <- months[i]
  df_m <- gas_data %>% filter(month == m)
  
  n <- nrow(df_m)
  if (n < 2) next  # skip months with only 1 station
  
  # Cross elasticity matrix: entry [j,k] = -alpha * Sk * Pk
  # This is the elasticity of station j's share with respect to station k's price
  # Note: does not depend on j at all, only on k
  # So each column k is constant: -alpha * Sk * Pk
  
  Sk <- df_m$sjt          # vector of shares
  Pk <- df_m$p_sell.imp   # vector of prices
  
  # Compute -alpha * Sk * Pk for each k (one value per competitor)
  cross_k <- -alpha * Sk * Pk  # length n vector
  
  # Store as a summary — mean cross elasticity that station j faces
  # (average over all k ≠ j competitors in that month)
  cross_elas_list[[i]] <- data.frame(
    month = m,
    cnpj  = df_m$cnpj,
    mean_cross_elasticity = sapply(seq_len(n), function(j) {
      mean(cross_k[-j], na.rm = TRUE)  # exclude own j
    })
  )
  
  # Clear intermediates to protect memory
  rm(Sk, Pk, cross_k)
  gc()
}

# Combine all months
cross_elas <- bind_rows(cross_elas_list)

# Average cross elasticity per station across all months
cross_elas_avg <- cross_elas %>%
  group_by(cnpj) %>%
  summarise(mean_cross_elasticity = mean(mean_cross_elasticity, na.rm = TRUE)) %>%
  ungroup()

cat("\n=== Cross-Price Elasticity (avg across months) ===\n")
summary(cross_elas_avg$mean_cross_elasticity)

# ------------------------------------------------------------
# STEP 5: Final summary table
# ------------------------------------------------------------
elasticity_table <- own_elas_avg %>%
  left_join(cross_elas_avg, by = "cnpj") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

cat("\n=== Elasticity Table (first 20 rows) ===\n")
print(head(elasticity_table, 20))

cat("\n=== Overall Averages ===\n")
cat("Mean own-price elasticity: ",  round(mean(elasticity_table$mean_own_elasticity,  na.rm=TRUE), 4), "\n")
cat("Mean cross-price elasticity:", round(mean(elasticity_table$mean_cross_elasticity, na.rm=TRUE), 4), "\n")

save(elasticity_table, file = '/proj/econ/duartema/GasBrazil/Mihika_Practice/elasticity_output.RData')