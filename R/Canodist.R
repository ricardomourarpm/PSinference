library(stats)

Canodist <- function(part, nsample, pvariates, iterations) {
  T <- c()
  for (i in 1:iterations) {
    W1 <- rWishart(1, nsample - 1, diag(pvariates))
    W2 <- rWishart(1, nsample - 1, W1 / (nsample - 1))
    W2_11 <- W2[part, part]
    W2_12 <- W2[part, -part]
    W2_21 <- W2[-part, part]
    W2_22 <- W2[-part, -part]
    Q <- W2_12 %% solve(W2_22) %% W2_21
    T <- c(T, det(Q) / det(W2_11 - Q))
  }
  return(T)

}
