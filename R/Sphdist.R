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
#' sphdist(nsample = 100, pvariates = 4, iterations = 4)
#'
#' @export

sphdist <- function(nsample=100, pvariates=4, iterations=2) {
  T <- rep(NA, iterations)
  W1 <- rWishart(iterations, nsample - 1, diag(pvariates) / (nsample - 1))
  W2 <- rWishart(iterations, nsample - 1, diag(pvariates))
  for (i in 1:iterations) {
    T[i] <- det(W1[,,i] %*% W2[,,i])^(1/pvariates) / sum(W1[,,i] %*% W2[,,i])
  }
  return(T)
}
