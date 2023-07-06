#' Independence Empirical Distribution
#'
#' This function calculates
#'
#' @param part .
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#'
#' @references
#' Klein, M., Moura, R. and Sinha, B. (2021). Multivariate Normal Inference based on Singly Imputed Synthetic Data under Plug-in Sampling. Sankhya B 83, 273–287.
#'
#' @importFrom stats rWishart
#' @examples
#' Inddist(part = 2, nsample = 100, pvariates = 4, iterations = 2)
#'
#' @export

Inddist <- function(part, nsample, pvariates, iterations) {
  T <- rep(NA, iterations)
  W1 <- stats::rWishart(iterations, nsample - 1, diag(pvariates))
  if (part == 1) {
    for (i in 1:iterations) {
      W2 <- stats::rWishart(1, nsample - 1, W1[, , i] / (nsample - 1))
      A <- partition(W2[, , 1], part, part)
      W2_11 <- A[[1]]
      W2_22 <- A[[4]]
      T[i] <- det(W2[, , 1]) / (W2_11 * det(W2_22))
    }
  } else {
    for (i in 1:iterations) {
      W2 <- stats::rWishart(1, nsample - 1, W1[, , i] / (nsample - 1))
      A <- partition(W2[, , 1], part, part)
      W2_11 <- A[[1]]
      W2_22 <- A[[4]]
      T[i] <- det(W2[, , 1]) / (det(W2_11) * det(W2_22))
    }
  }
  return(T)
}
