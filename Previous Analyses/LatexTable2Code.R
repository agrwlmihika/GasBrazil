options("modelsummary_format_numeric_latex" = "plain")
modelsummary(
  list("OLS" = model_gas,
       "IV"  = iv_gas,
       "sOLS" = model_spatial_ols,
       "sIV1" = model_spatial_iv1,
       "sIV2" = model_spatial_iv2),
  coef_map   = coef_map,
  gof_map    = c("nobs", "r.squared"),
  add_rows   = extra_rows,
  stars      = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  fmt        = 3,
  escape     = FALSE,
  output     = "latex_tabular"  # This forces booktabs format
) %>%
  # Add title, label, and notes manually
  {
    paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      "\\caption{Demand Estimates}\n",
      "\\label{tab:demand}\n",
      "\\begin{threeparttable}\n",
      "\\small\n",
      .,
      "\\begin{tablenotes}\n",
      "\\small\n",
      "\\item \\textit{Notes:} Dependent variable is $\\ln(s_{jt}/s_{0t})$, where $s_{jt}$ is station $j$'s market share in month $t$ and $s_{0t}$ is the outside share (constant at 1/3). ",
      "Standard errors (in parentheses) clustered by month $\\times$ neighborhood. ",
      "Spatial terms $B^g_{jt}$ capture competition from nearby stations in three distance bands. ",
      "IV regressions instrument price using interactions of cost shifters with station characteristics. ",
      "IV2 additionally instruments spatial terms using distance-based price deviations. ",
      "* $p < 0.1$, ** $p < 0.05$, *** $p < 0.01$.\n",
      "\\end{tablenotes}\n",
      "\\end{threeparttable}\n",
      "\\end{table}\n"
    )
  } %>%
  cat()
