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
#' @importFrom stats rchisq
#' @examples
#' GVdist(nsample = 10, pvariates = 4, iterations = 5)
#'
#' @export


GVdist <- function(nsample, pvariates, iterations) {
  q <- chiprod(pvariates, nsample - 1,iterations)
  p <- chiprod(pvariates, nsample - 1,iterations)
  return(p * q)
}
