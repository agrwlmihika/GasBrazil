# Construct instrument
sample_iv <- sample_demand %>%
  mutate(Z = n_oppo1km * gasC_cost)

# IV regression for gasoline
# iv_gas <- feols(
#   log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
#     log.rent + pop.area + dist |
#     cl |
#     p_sell.imp ~ destillery_price.d1:pumps.imp.mdev+
#                   refinery_price.d1:close_oppo.mdev+destillery_price.d1:close_oppo.mdev+
#                   refinery_price.d1:n_unbranded.mdev+destillery_price.d1:n_unbranded.mdev+
#                   refinery_price.d1:n_oppo1km.mdev+destillery_price.d1:n_oppo1km.mdev,
#   cluster = ~cl,
#   data = sample_iv %>% filter(prod == "gas")
# )
# 
# 
# iv_gas <- feols(
#   log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
#     log.rent + pop.area + dist |
#     cl |
#     p_sell.imp ~ Dstd.pdev.dist+Dstd.pdev.dist:n_unbranded.mdev + Dstd.pdev.dist:n_oppo1km.mdev+  Dstd.pdev.dist:close_oppo.mdev,
#   cluster = ~cl,
#   data = sample_iv %>% filter(prod == "gas")
# )


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

# # IV regression for ethanol
# iv_eth <- feols(
#   log_ratio ~ log.nstations_group + pumps.imp + tank.imp +
#     log.rent + pop.area + dist |
#     cl |
#     p_sell.imp ~ Z,
#   cluster = ~cl,
#   data = sample_iv %>% filter(prod == "eth")
# )

# IV results table only
# iv_table <- bind_rows(
#   data.frame(product="Gasoline", name=names(coef(iv_gas)),
#              coeff=round(coef(iv_gas),4), SE=round(se(iv_gas),4), row.names=NULL),
#   data.frame(product="Ethanol",  name=names(coef(iv_eth)),
#              coeff=round(coef(iv_eth),4), SE=round(se(iv_eth),4), row.names=NULL)
# )

iv_table <- data.frame(
  product = "Gasoline",
  name    = names(coef(iv_gas)),
  coeff   = round(coef(iv_gas), 4),
  SE      = round(se(iv_gas),   4),
  row.names = NULL
)
print(iv_table)