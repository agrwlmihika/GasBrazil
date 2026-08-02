library(dplyr)
library(fixest)
library(geosphere)

# Load the saved sample
load('Data/stations_sampleDF.RData')  # loads object called 'sample'

# Step 1: Compute n_competitors_3km ------------------------------------------

stations_geo <- sample %>%
  filter(prod == 'gas') %>%
  select(cnpj, lat, lon) %>%
  distinct()

coords <- as.matrix(stations_geo[, c('lon', 'lat')])

dist_matrix <- distm(coords, fun = distHaversine)
rownames(dist_matrix) <- stations_geo$cnpj
colnames(dist_matrix) <- stations_geo$cnpj

n_competitors_3km <- apply(dist_matrix, 1, function(x) sum(x > 0 & x <= 3000, na.rm = TRUE))

n_competitors_3km_df <- data.frame(
  cnpj = names(n_competitors_3km),
  n_competitors_3km = as.integer(n_competitors_3km)
)

# Step 2: Merge into sample --------------------------------------------------

reg_sample <- sample %>%
  filter(prod == 'gas') %>%
  left_join(n_competitors_3km_df, by = 'cnpj')

# Sanity check
summary(reg_sample$n_competitors_3km)
hist(reg_sample$n_competitors_3km, main = "Competitors within 3km", xlab = "N competitors")

# Step 3: Run regression -----------------------------------------------------

reg <- feols(
  p_sell ~ n_competitors_3km | month,
  data    = reg_sample,
  cluster = ~cl
)

summary(reg)