#' GVdist
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
#' GVdist(nsample = 10, pvariates = 4, iterations = 5)
#'
#' @export


GVdist <- function(nsample, pvariates, iterations) {
  T <- rep(NA, iterations)
  for (i in 1:iterations) {
    q <- chiprod(pvariates, nsample - 1)
    p <- chiprod(pvariates, nsample - 1)
    T[i] <- p * q
  }
  return(T)
}
