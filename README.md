# Aggregate utility

This repository contains the source code for the [R Shiny app](
https://janewliang.shinyapps.io/agg_utility/) associated with the paper 

> **Liang, J. W., Christensen, K. D., Green, R. C., and Kraft, P. (2023). A framework evaluating the utility of multi-gene, multi-disease population-based panel testing that accounts for uncertainty in penetrance estimates.**

The app uses the `rhandsontable` package to incorporate interactive data tables, as well as the `plotly` and `ggplot2` packages for visualization. Please see the `ui.R` and `server.R` files for more information. 

Those interested in using the aggregate utility directly should peruse the code and documentation in `functions.R`. Some example usage follows: 

```r
# Net utility for a gene with prevalence 0.005, disease with penetrance 0.6, 
# G+D- (d_01) and G-D+ (d_10) disutilities of 1, and test utility (K) of 0
Delta_ij(prev = 0.005, pen = 0.6, d_01 = 1, d_10 = 1, K = 0)
#>[1] 0.001
```

```r
# 5% quantile, where the uncertainty distribution has precision 100
Delta_ij_quantile(prev = 0.005, pen = 0.6, d_01 = 1, d_10 = 1, K = 0, 
	prob = 0.05, precision = 100)
#> [1] 0.0001863664
```

```r
# Probability of a positive net utility
Delta_ij_prob_pos(prev = 0.005, pen = 0.6, d_01 = 1, d_10 = 1, K = 0, 
	precision = 100)
#> [1] 0.9780624
```