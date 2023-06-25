#' Canodist
#'
#' This function calculates the positive-rule Stein estimator. This estimator is
#'
#' part, nsample, pvariates, iterations
#' @param part .
#' @param nsample .
#' @param pvariates .
#' @param iterations .
#' @references
#'  ref
#'
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
