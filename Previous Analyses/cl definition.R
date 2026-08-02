# Check what cl looks like
sample_spatial %>%
  select(cl, cnpj, month) %>%
  distinct() %>%
  head(20)

# How many unique values?
cat("Unique cl values:    ", n_distinct(sample_spatial$cl),   "\n")
cat("Unique cnpj values:  ", n_distinct(sample_spatial$cnpj), "\n")
cat("Unique months:       ", n_distinct(sample_spatial$month),"\n")
cat("cnpj x month combos: ", n_distinct(paste(sample_spatial$cnpj,
                                              sample_spatial$month)), "\n")