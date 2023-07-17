#' Generalized Variance Empirical Distribution
#'
#' This function calculates the empirical distribution of the pivotal random
#' variable that can be used to perform inferential procedures
#' for the Generalized Variance of the released Single Synthetic dataset
#' generated under Plug-in Sampling, assuming that the original distribution
#' is normally distributed.
#'
#' We define
#' \deqn{T_1^\star = (n-1)\frac{|\boldsymbol{S}^*|}{|\boldsymbol{\Sigma}|},}
#' where \eqn{\boldsymbol{S}^\star = \sum_{i=1}^n (v_i - \bar{v})(v_i - \bar{v})^{\top}}, \eqn{\boldsymbol{\Sigma}} is the population covariance matrix
#' and \eqn{v_i} is the \eqn{i}th observation of the synthetic dataset.
#' Its distribution is stochastic equivalent to
#' \deqn{ \prod_{i=1}^n \chi_{n-i}^2 \prod_{i=1}^p \chi_{n-i}^2}
#' where \eqn{\chi_{n-i}^2} are all independent chi-square random variables.
#' The \eqn{(1-\alpha)} level confidence interval for \eqn{|\boldsymbol{\Sigma}|}
#' is given by
#' \deqn{\left(\frac{(n-1)^p|\tilde{\boldsymbol{S}}^\star|}{t^\star_{1,1-\alpha/2}},
#' \frac{(n-1)^p|\tilde{\boldsymbol{S}}^\star|}{t^\star_{1,\alpha/2}} \right)}
#' where \eqn{\tilde{\boldsymbol{S}}^\star} is the observed value of
#' \eqn{\boldsymbol{S}^\star} and \eqn{t^\star_{1,\gamma}}
#' is the \eqn{\gamma}th percentile of \eqn{T_1}.
#'
#' @param nsample Sample size.
#' @param pvariates Number of variables.
#' @param iterations Number of iterations for simulating values from the distribution and finding the quantiles. Default is \code{10000}.
#'
#' @references
#' Klein, M., Moura, R. and Sinha, B. (2021). Multivariate Normal Inference based on Singly Imputed Synthetic Data under Plug-in Sampling. Sankhya B 83, 273–287.
#'
#' @examples
#'
#'# Original data creation
#'library(MASS)
#'mu <- c(1,2,3,4)
#'Sigma <- matrix(c(1, 0.5, 0.5, 0.5,
#'                   0.5, 1, 0.5, 0.5,
#'                   0.5, 0.5, 1, 0.5,
#'                   0.5, 0.5, 0.5, 1), nrow = 4, ncol = 4, byrow = TRUE)
#' seed = 1
#' n_sample = 100
#' # Create original simulated dataset
#' df = mvrnorm(n_sample, mu = mu, Sigma = Sigma)
#'
#'# Synthetic data created
#'
#'df_s = simSynthData(df)
#'
#'
#'# Gather the 0.025 and 0.975 quantiles and construct confident interval for sigma^2
#'# Check that sigma^2 is inside in both cases
#'p = dim(df_s)[2]
#'
#'T <- GVdist(100, p, 10000)
#'q975 <- quantile(T, 0.975)
#'q025 <- quantile(T, 0.025)
#'
#'left <- (n_sample-1)^p * det(cov(df_s)*(n_sample-1))/q975
#'right <- (n_sample-1)^p * det(cov(df_s)*(n_sample-1))/q025
#'
#'cat(left,right,'\n')
#'print(det(Sigma))
#'# The synthetic value is inside the confidence interval of GV
#'
#' @export


GVdist <- function(nsample, pvariates, iterations = 10000) {
  q <- chiprod(pvariates, nsample - 1,iterations)
  p <- chiprod(pvariates, nsample - 1,iterations)
  return(p * q)
}

