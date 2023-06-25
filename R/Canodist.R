library(stats)
part <- 2
nsample <- 100
pvariates <- 4
iterations <- 2
Canodist <- function(part, nsample, pvariates, iterations) {
  T <- c()
  for (i in 1:iterations) {
    W1 <- rWishart(1, nsample - 1, diag(pvariates))
    W2 <- rWishart(1, nsample - 1, W1[,,1] / (nsample - 1))
    A <- partition(W2[,,1],part,part)
    W2_11 <- A[[1]]
    W2_12 <- A[[2]]
    W2_21 <- A[[3]]
    W2_22 <- A[[4]]
    Q <- W2_12 %*% solve(W2_22) %*% W2_21
    T[i] <-  det(Q) / det(W2_11 - Q)
  }
  return(T)

}
Canodist(2,100,4,2)
