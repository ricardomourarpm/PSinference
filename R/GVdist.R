#' GVdist  is
#'
#' This function calculates the empirical distribution of the pivotal random variable that can be used to perform inferential procedures
#' for the Generalized Variance based on the released Single Synthetic data generated under Plug-in Sampling.
#' \deqn{T_1^\star = (n-1)\frac{|S^*|}{|\Sigma|}}
#' where \eqn{S^\star = \sum_{i=1}^n (v_i - \bar{v})(v_i - \bar{v})^{\top}} and \eqn{v_i} is the \eqn{i}th observation of the synthetic dataset.
#' Its distribution is
#' \deqn{T_1^* \sim  \prod_{i=1}^n \chi_{n-i}^2 \prod_{i=1}^p \chi_{n-i}^2}
#' where \eqn{\chi_{n-i}^2} are all independent chi-square random variables.
#' The \eqn{(1-\alpha)} level confidence interval for \eqn{|\Sigma|} is given by
#'
#' \deqn{\left(\frac{(n-1)^p|\tilde{S}^\star|}{t^\star_{1,1-\alpha/2}},\left(\frac{(n-1)^p|\tilde{S}^\star|}{t^\star_{1,\alpha/2}} \right)}
#' where \eqn{|\tilde{S}^\star|} is the observed value of \eqn{|\{S}^\star|} and \eqn{t^\star_{1,\gamma}} is the \eqn{\gamma}th percentile of \eqn{T_1}.
#'
#' @param nsample Sample size.
#' @param pvariates Number of variables.
#' @param iterations Number of iterations for simulating values from the distribution and finding the quantiles. Default is \code{10000}.
#'
#' @references
#'  ref
#'
#' @examples
#' library(MASS)
#'
#'# Original data created
#'# two different population variances (determinants)
#'Sigma1 <- matrix(c(1, 0.5, 0.5, 0.5,
#'                   0.5, 1, 0.5, 0.5,
#'                   0.5, 0.5, 1, 0.5,
#'                   0.5, 0.5, 0.5, 1), nrow = 4, ncol = 4, byrow = TRUE)
#'
#'Sigma2 <- matrix(c(1, 0.5, 0, 0,
#'                   0.5, 2, 0, 0,
#'                   0, 0, 3, 0.2,
#'                   0, 0, 0.2, 4), nrow = 4, ncol = 4, byrow = TRUE)
#' n_sample = 100
#' # Create original simulated dataset
#' df_o1 = mvrnorm(n_sample, mu = c(1,0,5,0), Sigma = Sigma1)
#' df_o2 = mvrnorm(n_sample, mu = c(1,0,5,0), Sigma = Sigma2)
#'
#'# Sim data created
#'
#'df_s1 = simSynthData(df_o1)
#'df_s2 = simSynthData(df_o2)
#'
#'
#'# Gather the 0.025 and 0.975 quantiles and construct confident interval for sigma^2
#'# Check that sigma^2 is inside in both cases
#'p1 = dim(df_s1)[2]
#'p2 = dim(df_s2)[2]
#'
#'T1 <- GVdist(100, p1, 10000)
#'q975 <- quantile(T1, 0.975)
#'q025 <- quantile(T1, 0.025)
#'
#'left <- (100-1)^p1*det(cov(df_s1)*99)/q975
#'right <- (100-1)^p1*det(cov(df_s1)*99)/q025
#'
#'cat(left,right,'\n')
#'print(det(Sigma1))
#'# The synthetic value is inside the confidence interval of GV
#'
#' @export


GVdist <- function(nsample, pvariates, iterations = 10000) {
  q <- chiprod(pvariates, nsample - 1,iterations)
  p <- chiprod(pvariates, nsample - 1,iterations)
  return(p * q)
}

