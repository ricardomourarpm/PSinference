Sphdist <- function(nsample, pvariates, iterations) {
  T <- c()
  for (i in 1:iterations) {
    W1 <- rWishart(1, nsample - 1, diag(pvariates) / (nsample - 1))
    W2 <- rWishart(1, nsample - 1, diag(pvariates))
    T <- c(T, det(W1 %% W2)^(1/pvariates) / sum(W1 %% W2))
  }
  return(T)
}
