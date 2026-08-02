# Demand Model
library(dplyr)
library(fixest)

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