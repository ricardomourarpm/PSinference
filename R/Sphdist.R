#' sphdist
#'
#' This function calculates
#'
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#'
#' @references
#'  ref
#'
#' @examples
#' Sphdist(nsample = 100, pvariates = 4, iterations = 2)
#'
#' @export

Sphdist <- function(nsample, pvariates, iterations) {
  T <- rep(NA, iterations)
  W1 <- rWishart(iterations, nsample - 1, diag(pvariates) / (nsample - 1))
  W2 <- rWishart(iterations, nsample - 1, diag(pvariates))
  for (i in 1:iterations) {
    inner_prod <- crossprod(W1[,,i], W2[,,i])
    T[i] <- det(inner_prod)^(1/pvariates) / sum(inner_prod)  }
  return(T)
}
