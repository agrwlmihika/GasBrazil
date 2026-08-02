# ============================================================
# 2. Spatial Regressions — Spec 3 (sOLS) and Spec 4 (sIV1)
# Fleet market size, new price instruments, station + time FE
# ============================================================
load('Data/stations_sampleDF.RData')
load('Data/Geo.RData')
load('/proj/econ/duartema/GasData/Database/Frota.RData')
library(dplyr)
library(fixest)
library(modelsummary)

# ------------------------------------------------------------
# STEP 1: Fleet market size + shares
# ------------------------------------------------------------
fleet_monthly <- Frota %>%
  filter(state == "DF") %>%
  group_by(month) %>%
  summarise(fleet_t = sum(AUTOMOVEL, na.rm=TRUE), .groups="drop") %>%
  mutate(market_size_fleet = fleet_t * 100)

sample_demand <- sample %>%
  filter(prod == "gas") %>%
  left_join(fleet_monthly, by="month") %>%
  group_by(month) %>%
  mutate(
    sjt       = sale.imp / market_size_fleet,
    s0t       = 1 - sum(sjt, na.rm=TRUE),
    log_ratio = log(sjt / s0t)
  ) %>%
  ungroup() %>%
  filter(is.finite(log_ratio), s0t > 0, sjt > 0)
sample_demand <- sample_demand %>%
  mutate(station_id = paste(lat, lon, sep = "_"))

cat("s0t range:", round(min(sample_demand$s0t),4),
    "to", round(max(sample_demand$s0t),4), "\n")
cat("Rows:", nrow(sample_demand), "\n")

# ------------------------------------------------------------
# STEP 2: Build spatial terms B1 B2 B3
# ------------------------------------------------------------
gas_shares <- sample_demand %>%
  select(id, month, station_id, sjt) %>%
  filter(!is.na(sjt))

build_Bg <- function(Dmat, shares_m) {
  ids       <- shares_m$id
  D         <- Dmat[ids, ids]
  Sj        <- shares_m$sjt
  ratio_mat <- outer(Sj, Sj, function(sj,sk) log(1+sk/sj))
  diag(ratio_mat) <- 0
  w1 <- (D > 0)    & (D < 0.5)
  w2 <- (D >= 0.5) & (D < 1.5)
  w3 <- (D >= 1.5) & (D < 3.0)
  data.frame(id=ids,
             B1=rowSums(w1*ratio_mat, na.rm=TRUE),
             B2=rowSums(w2*ratio_mat, na.rm=TRUE),
             B3=rowSums(w3*ratio_mat, na.rm=TRUE))
}

months_list <- unique(gas_shares$month)
Bg_list     <- vector("list", length(months_list))

for (i in seq_along(months_list)) {
  m        <- months_list[i]
  shares_m <- gas_shares %>% filter(month == m)
  Dmat     <- Geo$Dmatrix[[i]]
  keep     <- shares_m$id[shares_m$id %in% rownames(Dmat)]
  if (length(keep) < 2) { rm(Dmat); gc(); next }
  Bg_list[[i]] <- build_Bg(Dmat, shares_m %>% filter(id %in% keep))
  rm(Dmat); gc()
}

sample_spatial <- sample_demand %>%
  left_join(bind_rows(Bg_list), by="id") %>%
  filter(!is.na(B1), !is.na(B2), !is.na(B3))
cat("Spatial sample rows:", nrow(sample_spatial), "\n")

# ------------------------------------------------------------
# STEP 3: Instruments
# ------------------------------------------------------------
instruments <- ~ refinery_price.d1:pumps.imp.mdev +
  destillery_price.d1:pumps.imp.mdev +
  refinery_price.d1:tank.imp.mdev +
  destillery_price.d1:tank.imp.mdev +
  refinery_price.d1:n_unbranded.mdev +
  destillery_price.d1:n_unbranded.mdev +
  Dstd.pdev.dist

# ------------------------------------------------------------
# STEP 4: Spec 3 — Spatial OLS with station + time FE
# ------------------------------------------------------------
model_sOLS <- feols(
  log_ratio ~ p_sell.imp + pumps.imp +
    tank.imp + B1 + B2 + B3 |
    station_id + month,
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Spec 3: Spatial OLS ===\n")
print(data.frame(name=names(coef(model_sOLS)),
                 coeff=round(coef(model_sOLS),4),
                 SE=round(se(model_sOLS),4), row.names=NULL))

# ------------------------------------------------------------
# STEP 5: Spec 4 — Spatial IV1 with station + time FE
# ------------------------------------------------------------
model_sIV1 <- feols(
  as.formula(paste(
    "log_ratio ~ pumps.imp + tank.imp +",
    "B1 + B2 + B3 | station_id + month |",
    "p_sell.imp ~",
    paste(attr(terms(instruments), "term.labels"), collapse=" + ")
  )),
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Spec 4: Spatial IV1 ===\n")
cat("First stage F-stat:\n"); print(fitstat(model_sIV1, "ivf"))
print(data.frame(name=names(coef(model_sIV1)),
                 coeff=round(coef(model_sIV1),4),
                 SE=round(se(model_sIV1),4), row.names=NULL))

# ------------------------------------------------------------
# STEP 5b: Spec 5 — Spatial IV2 with station + time FE
# (price AND B1, B2, B3 all endogenous)
# ------------------------------------------------------------
model_sIV2 <- feols(
  as.formula(paste(
    "log_ratio ~ pumps.imp + tank.imp | station_id + month |",
    "p_sell.imp + B1 + B2 + B3 ~",
    paste(attr(terms(instruments), "term.labels"), collapse=" + ")
  )),
  cluster = ~cl,
  data = sample_spatial
)

cat("\n=== Spec 5: Spatial IV2 (price + spatial endogenous) ===\n")
print(data.frame(name=names(coef(model_sIV2)),
                 coeff=round(coef(model_sIV2),4),
                 SE=round(se(model_sIV2),4), row.names=NULL))


# ------------------------------------------------------------
# STEP 6: Output table
# ------------------------------------------------------------
coef_map <- c(
  "p_sell.imp"                            = "Price",
  "fit_p_sell.imp"                        = "Price",
  "log.nstations_group"                   = "Log # stations in group",
  "pumps.imp"                             = "Pumps",
  "tank.imp"                              = "Tank size",
  "log.rent"                              = "Log rent",
  "pop.area"                              = "Population density",
  "distIPIRANGA PRODUTOS DE PETRÓLEO S.A" = "Ipiranga",
  "distPETROBRAS DISTRIBUIDORA S.A."      = "Petrobras",
  "distRAIZEN COMBUSTÍVEIS S.A."          = "Raizen",
  "B1"                                    = "B1 (< 500m)",
  "B2"                                    = "B2 (500m-1.5km)",
  "B3"                                    = "B3 (1.5km-3km)",
  "fit_B1"                                = "B1 (< 500m)",
  "fit_B2"                                = "B2 (500m-1.5km)",
  "fit_B3"                                = "B3 (1.5km-3km)"
)

modelsummary(
  list("sOLS (3)" = model_sOLS, "sIV1 (4)" = model_sIV1, "sIV2 (5)" = model_sIV2),
  coef_map   = coef_map,
  gof_map    = c("nobs", "r.squared"),
  stars      = c("*"=0.1, "**"=0.05, "***"=0.01),
  notes      = "SE clustered by station and month. Station + time FE included.",
  output     = "/proj/econ/duartema/GasBrazil/Final/Output/table2_regressions.txt"
)

# ------------------------------------------------------------
# STEP 7: Save
# ------------------------------------------------------------
save(model_sOLS, model_sIV1, model_sIV2, sample_spatial, sample_demand,
     file = "/proj/econ/duartema/GasBrazil/Final/Data/regression_output.RData")
cat("\nSaved to Data/regression_output.RData\n")