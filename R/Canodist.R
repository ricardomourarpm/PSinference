#' Canonical Empirical Distribution

#'
#' This function calculates the empirical distribution of the pivotal random variable that can be used to perform inferential procedures
#' for the test for the regression of one Set of variables on the other based on the released Single Synthetic data generated under Plug-in Sampling.
#' \deqn{T_4^\star = \frac{|S^{\star}_{12}(S^{\star}_{22})^{-1}-\Delta)S^{\star}_{22}(S^{\star}_{12})(S^{\star}_22)^{-1}-\Delta)'|}{|S^{\star}_{11.2}|}}
#' where \eqn{S^\star = \sum_{i=1}^n (v_i - \bar{v})(v_i - \bar{v})^{\top}} and \eqn{v_i} is the \eqn{i}th observation of the synthetic dataset. It is partinoed as 
#' \deqn{S^{\star}=\left[\begin{array}{lll} S^{\star}_{11}& S^{\star}_{12}\\ \S^{\star}_{21} & S^{\star}_{22} \end{array}\right]}
#' 
#' For \eqn{\Delta = \Sigma_{12}\Sigma_{22}^{-1}}, its distribution is
#' \deqn{T_4^* \sim  \frac{|\Omega_{12}\Omega_{22}^{-1}\Omega_{21}|}{|\Omega_{11}-\Omega_{12}\Omega_{22}^{-1}\Omega_{21}|}
#' where \eqn{\Omega} are all independent Wishart random variables, \eqn{\Omega \sim W_p(n-1, \frac{W}{n-1})} where \eqn{W \sim W_P(n-1, I_p)} and \eqn{\omega} partioned in the same way as \eqn{S^{\star}}.
#' The \eqn{(1-\alpha)} level confidence interval for \eqn{|\Sigma|} is given by
#' \deqn{\left(\frac{(n-1)^p|\tilde{S}^\star|}{t^\star_{1,1-\alpha/2}},\frac{(n-1)^p|\tilde{S}^\star|}{t^\star_{1,\alpha/2}} \right)}
#' where \eqn{\tilde{S}^\star} is the observed value of \eqn{S^\star} and \eqn{t^\star_{1,\gamma}} is the \eqn{\gamma}th percentile of \eqn{T_1}.
#'
#' @param part Number of partitions.
#' @param nsample Sample size.
#' @param pvariates Number of variables.
#' @param iterations Number of iterations for simulating values from the distribution and finding the quantiles. Default is \code{10000}.
#'
#'
#' @references
#'  Klein, M., Moura, R. and Sinha, B. (2021). Multivariate Normal Inference based on Singly Imputed Synthetic Data under Plug-in Sampling. Sankhya B 83, 273–287.
#'
#' @importFrom stats rWishart
#' @examples
#' canodist(part = 2, nsample = 100, pvariates = 4, iterations = 2)
#'
#' @export

canodist <- function(part, nsample, pvariates, iterations) {
  T <- rep(NA, iterations)
  W1 <- stats::rWishart(iterations, nsample - 1, diag(pvariates))
  for (i in 1:iterations) {
    W2 <- stats::rWishart(1, nsample - 1, W1[, , i] / (nsample - 1))
    A <- partition(W2[, , 1], part, part)
    W2_11 <- A[[1]]
    W2_12 <- A[[2]]
    W2_21 <- A[[3]]
    W2_22 <- A[[4]]
    Q <- W2_12 %*% solve(W2_22) %*% W2_21
    T[i] <- det(Q) / det(W2_11 - Q)
  }
  return(T)
}
