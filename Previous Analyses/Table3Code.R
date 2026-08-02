library(dplyr)
library(tidyr)

# 1️⃣ Load data
obj_name <- load("/proj/econ/duartema/GasData/Database/PriceData.RData")
PriceData <- get(obj_name)

# 2️⃣ Clean and prepare data
PriceData_clean <- PriceData %>%
  filter(!is.na(p_buy_gas)) %>%
  mutate(
    week = as.Date(week),
    branded = ifelse(dist == "BRANCA", 0, 1)
  )

# 3️⃣ Define periods
periods <- list(
  '2007-2010' = c(as.Date("2007-01-01"), as.Date("2010-12-31")),
  '2011-2015' = c(as.Date("2011-01-01"), as.Date("2015-12-31")),
  '2016-2020' = c(as.Date("2016-01-01"), as.Date("2020-12-31"))
)

# Initialize results
row1 <- c()
row2 <- c()
row3 <- c()

# 4️⃣ Loop through periods
for (period_name in names(periods)) {
  
  start_date <- periods[[period_name]][1]
  end_date <- periods[[period_name]][2]
  
  data_period <- PriceData_clean %>%
    filter(week >= start_date & week <= end_date)
  
  # --- Row 1: Average price (weekly avg → mean) ---
  avg_price <- data_period %>%
    group_by(week) %>%
    summarise(weekly_avg = mean(p_buy_gas, na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_price = mean(weekly_avg, na.rm = TRUE)) %>%
    pull(mean_price)
  
  row1[period_name] <- avg_price * 100
  
  # --- Row 2: Weekly standard deviation ---
  avg_sd <- data_period %>%
    group_by(week) %>%
    summarise(weekly_sd = sd(p_buy_gas, na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_sd = mean(weekly_sd, na.rm = TRUE)) %>%
    pull(mean_sd)
  
  row2[period_name] <- avg_sd * 100
  
  # --- Row 3: Branded - Unbranded difference ---
  diff_val <- data_period %>%
    group_by(week, branded) %>%
    summarise(avg_price = mean(p_buy_gas, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = branded, values_from = avg_price) %>%
    mutate(diff = `0` - `1`) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE)) %>%
    pull(mean_diff)
  
  row3[period_name] <- diff_val * 100
}

# 5️⃣ Combine into final Table 3
table3 <- rbind(
  Average_Wholesale_Price = row1,
  Weekly_SD = row2,
  Branded_minus_Unbranded = row3
)

print(table3)
