library(MASS)
set.seed(1234)
sim <- 100
#Sample size and partition size
n <- 100
#Some population mean
mu <- c(1, 2, 3, 4)
#Number of covariates
p <- length(mu)
#Three different population variances and respective partitions
Sigma1 <- diag(p)
Sigma2 <- 5 * diag(p)
T <- Sphdist(n, p, sim)
q05 <- quantile(T, 0.05)
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

  T1temp <- det(S_star1)^(1/p) / sum(diag(S_star1))
  T2temp <- det(S_star2)^(1/p) / sum(diag(S_star2))
  T1_1[i] <- T1temp
  T1_2[i] <- T2temp
}

print(c(min(T1_1), quantile(T1_1, 0.1), quantile(T1_1, 0.5), quantile(T1_1, 0.9), max(T1_1)))
print(c(min(T1_2), quantile(T1_2, 0.1), quantile(T1_2, 0.5), quantile(T1_2, 0.9), max(T1_2)))
print(c(min(T), quantile(T, 0.1), quantile(T, 0.5), quantile(T, 0.9), max(T)))
# Plotting KDE
windows()
par(mfrow=c(2,2))
plot(density(T ^ (1 / p)),col=1)
plot(density(T1_1 ^ (1 / p)),col=2)
plot(density(T1_2 ^ (1 / p)),col=3)
par(mfrow=c(1,1))

plot(density(T ),col=1)
lines(density(T1_1),col=2)
lines(density(T1_2) ,col=3)



cov1 <- mean(T1_1 > q05)
cov2 <- mean(T1_2 > q05)

print(c(cov1, cov2))
