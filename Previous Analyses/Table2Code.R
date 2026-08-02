# Step 0: Load data
load("/proj/econ/duartema/GasData/Database/SellDist.RData")
load("/proj/econ/duartema/GasData/Database/SellTotState.RData")

# Step 1: Extract year
SellDist$year <- as.numeric(substr(SellDist$month, 1, 4))
SellTotState$year <- as.numeric(substr(SellTotState$month, 1, 4))

# Step 2: Restrict to 2011–2015
SellDist <- SellDist[SellDist$year >= 2011 & SellDist$year <= 2015, ]
SellTotState <- SellTotState[SellTotState$year >= 2011 & SellTotState$year <= 2015, ]

# Step 3: Restrict to Federal District
SellDist_DF <- SellDist[SellDist$state == "DF", ]
SellTotState_DF <- SellTotState[SellTotState$state == "DF", ]
total_DF_sales <- sum(SellTotState_DF$vol_sell)

# Step 4: Normalize distributor names
normalize_name <- function(name) {
  name <- toupper(name)                    # uppercase
  name <- gsub("[[:punct:]]", "", name)   # remove punctuation
  name <- gsub("\\b(SA|S A|S/A|LTDA|LTD)\\b", "", name) # remove common suffixes
  name <- trimws(name)                     # remove leading/trailing spaces
  return(name)
}

SellDist_DF$dist_norm <- sapply(SellDist_DF$dist, normalize_name)

# Step 5: Aggregate by normalized distributor name
dist_totals_DF <- aggregate(vol_sell ~ dist_norm, data = SellDist_DF, sum)

# Step 6: Identify top 3 distributors automatically
dist_totals_DF <- dist_totals_DF[order(-dist_totals_DF$vol_sell), ]
top3_DF <- head(dist_totals_DF, 3)

# Step 7: Compute Gas Sale (%) for each top distributor
top3_DF$GasSale_pct <- round((top3_DF$vol_sell / total_DF_sales) * 100, 5)

# Step 8: Compute total share of top 3
top3_total <- round(sum(top3_DF$vol_sell) / total_DF_sales * 100, 5)

# Step 9: Final output
result <- rbind(
  top3_DF[, c("dist_norm", "GasSale_pct")],
  data.frame(dist_norm = "Total", GasSale_pct = top3_total)
)

print(result)