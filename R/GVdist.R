GVdist <- function(nsample, pvariates, iterations) {
  T <- c()
  for (i in 1:iterations) {
    q <- chiprod(pvariates, nsample - 1)
    p <- chiprod(pvariates, nsample - 1)
    T <- c(T, p * q)
  }
  return(T)

}
