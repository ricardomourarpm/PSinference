library(Matrix)
library(mvtnorm)
library(ggplot2)

set.seed(4321)
sim <- 100

# Sample size and partition size
n <- 100
p1 <- 1
p2 <- 2

# Some population mean
mu <- c(1, 2, 3, 4)
# Number of covariates
p <- length(mu)

# Three different population variances and respective partitions
Sigma1 = diag(p)
S1<- partition(Sigma1,p1,p1)

Sigma1_11 <- S1[[1]]
Sigma1_12 <- S1[[2]]
Sigma1_21 <- S1[[3]]
Sigma1_22 <- S1[[4]]

Sigma1_112 <- Sigma1_11 - Sigma1_12 %*% (solve(Sigma1_22) %*% Sigma1_21)

Sigma2 <- matrix(c(1, 0.5, 0, 0, 0.5, 2, 0, 0, 0, 0, 3, 0.2, 0, 0, 0.2, 4), nrow = p, ncol = p, byrow = TRUE)
S2 <- partition(Sigma2,p2,p2)
Sigma2_11 <- S2[[1]]
Sigma2_12 <- S2[[2]]
Sigma2_21 <- S2[[3]]
Sigma2_22 <- S2[[4]]

Sigma2_112 <- Sigma2_11 - Sigma2_12 %*% (solve(Sigma2_22) %*% Sigma2_21)

T1 <- Inddist(p1, n, p, sim)
T2 <- Inddist(p2, n, p, sim)


q051 <- quantile(T1, 0.05)
q052 <- quantile(T2, 0.05)

T1_1 <- c()
T1_2 <- c()

for (i in 1:sim) {

  # Generate original data sample from normal with mu and Sigma
  X1 <- MASS::mvrnorm(n, mu, Sigma1)
  X2 <- MASS::mvrnorm(n, mu, Sigma2)
  mean1 <- colMeans(X1)
  mean2 <- colMeans(X2)
  S1 <- (t(X1) - mean1)%*%t(t(X1) - mean1)
  S2 <- (t(X2) - mean2)%*%t(t(X2) - mean2)

  # Generate PLS synthetic single data
  V1 <- MASS::mvrnorm(n, mean1, S1/(n-1))
  V2 <- MASS::mvrnorm(n, mean2, S2/(n-1))

  # PLS estimates of mu and Sigma
  meanV1 <- colMeans(V1)
  meanV2 <- colMeans(V2)
  S_star1 <- (t(V1) - meanV1)%*%t(t(V1) - meanV1)
  S_star2 <- (t(V2) - meanV2)%*%t(t(V2) - meanV2)

  s3 <- partition(S_star1,p1,p1)
  S_star1_11 <- as.matrix(s3[[1]])
  S_star1_12 <- s3[[2]]
  S_star1_21 <- s3[[3]]
  S_star1_22 <- s3[[4]]

  S_star1_112 <- S_star1_11 - S_star1_12 %*% (solve(S_star1_22) %*% S_star1_21)

  s4 <- partition(S_star2,p2,p2)
  S_star2_11 <- s4[[1]]
  S_star2_12 <- s4[[2]]
  S_star2_21 <- s4[[3]]
  S_star2_22 <- s4[[4]]

  S_star2_112 <- S_star2_11 - S_star2_12 %*% (solve(S_star2_22) %*% S_star2_21)

  T1temp <- det(S_star1) / (det(S_star1_11) * det(S_star1_22))
  T2temp <- det(S_star2) / (det(S_star2_11) * det(S_star2_22))

  T1_1[i] <- T1temp
  T1_2[i] <- T2temp
}

print(c( min(T1_1), quantile(T1_1, 0.10), quantile(T1_1, 0.50), quantile(T1_1, 0.90), max(T1_1)))
print(c( min(T1_2), quantile(T1_2, 0.10), quantile(T1_2, 0.50), quantile(T1_2, 0.90), max(T1_2)))
print(c( min(T1), quantile(T1, 0.10), quantile(T1, 0.50), quantile(T1, 0.90), max(T1)))
print(c( min(T2), quantile(T2, 0.10), quantile(T2, 0.50), quantile(T2, 0.90), max(T2)))

windows(10)
par(mfrow=c(2,2))
plot(density(T1) ,col=1)
plot(density(T2) ,col=1)
plot(density(T1_1),col=2)
plot(density(T1_2),col=3)

par(mfrow=c(1,1))
plot(density(T1),col=1)
lines(density(T1_1),col=2)

plot(density(T2),col=1)
lines(density(T1_2),col=2)


cov1 = mean(T1_1>q051)
cov2 = mean(T1_2>q052)

print(c(cov1,cov2))

