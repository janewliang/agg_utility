#' Net utility for a given gene and disease
#' 
#' @param prev Prevalence of carrying a pathogenic variant of the gene
#' @param pen Cumulative lifetime penetrance for the disease, conditional on 
#' carrying a pathogenic variant of the gene
#' @param d_01 Disutility associated with a G+D- result for the gene, in 
#' relation to the disease
#' @param d_10 Disutility associated with a G-D+ result for the gene, in 
#' relation to the disease
#' @param K Disutility (financial and psychological) associated with 
#' conducting the test, regardless of results
#' 
#' @return Net utility as a single output value
Delta_ij = function(prev, pen, d_01, d_10, K) {
    prev * (- d_01 + (d_01 + d_10) * pen) + K
}


#' Net utility quantile for a given gene and disease
#' 
#' @param prev Prevalence of carrying a pathogenic variant of the gene
#' @param pen Cumulative lifetime penetrance for the disease, conditional on 
#' carrying a pathogenic variant of the gene
#' @param d_01 Disutility associated with a G+D- result for the gene, in 
#' relation to the disease
#' @param d_10 Disutility associated with a G-D+ result for the gene, in 
#' relation to the disease
#' @param K Disutility (financial and psychological) associated with 
#' conducting the test, regardless of results
#' @param prob Quantile probability desired
#' @param precision Hypothetical trial size used to motivate the choice of 
#' parameters in the Beta uncertainty distribution
#' 
#' @details Uncertainty distribution used to obtain quantile is a Beta 
#' distribution with parameters chosen to be the expected number of carriers 
#' who develop and do not develop the disease in a study with trial size equal 
#' to `precision` argument. 
#' @return Net utility quantile as a single output value
Delta_ij_quantile = function(prev, pen, d_01, d_10, K, prob, precision) {
    p = qbeta(prob, 
              shape1 = pen * precision, 
              shape2 = precision - pen * precision)
    Delta_ij(prev, p, d_01, d_10, K)
}


#' Probability that the net utility for a given gene and disease is positive
#' 
#' @param prev Prevalence of carrying a pathogenic variant of the gene
#' @param pen Cumulative lifetime penetrance for the disease, conditional on 
#' carrying a pathogenic variant of the gene
#' @param d_01 Disutility associated with a G+D- result for the gene, in 
#' relation to the disease
#' @param d_10 Disutility associated with a G-D+ result for the gene, in 
#' relation to the disease
#' @param K Disutility (financial and psychological) associated with 
#' conducting the test, regardless of results
#' @param precision Hypothetical trial size used to motivate the choice of 
#' parameters in the Beta uncertainty distribution
#' 
#' @details Uncertainty distribution used to calculate probability is a Beta 
#' distribution with parameters chosen to be the expected number of carriers 
#' who develop and do not develop the disease in a study with trial size equal 
#' to `precision` argument. 
#' @return Probability that the net utility is positive
Delta_ij_prob_pos = function(prev, pen, d_01, d_10, K, precision) {
    1 - pbeta(d_01 / (d_01 + d_10), 
              shape1 = pen * precision, 
              shape2 = precision - pen * precision)
}
