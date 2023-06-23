library(Matrix)
library(mvtnorm)
library(ggplot2)

set.seed(4321)
sim <- 100000

# Sample size and partition size
n <- 100
p1 <- 1
p2 <- 2

# Some population mean
mu <- c(1, 2, 3, 4)
# Number of covariates
p <- length(mu)

# Three different population variances and respective partitions

Sigma1 <- diag(p)

Sigma1_11 <- Sigma1[1:p1, 1:p1]
Sigma1_12 <- Sigma1[1:p1, (p1+1):p]
Sigma1_21 <- Sigma1[(p1+1):p, 1:p1]
Sigma1_22 <- Sigma1[(p1+1):p, (p1+1):p]

Sigma1_112 <- Sigma1_11 - Sigma1_12 %*% solve(Sigma1_22) %*% Sigma1_21

Sigma2 <- matrix(c(1, 0.5, 0, 0, 0.5, 2, 0, 0, 0, 0, 3, 0.2, 0, 0, 0.2, 4), nrow = p, ncol = p, byrow = TRUE)

Sigma2_11 <- Sigma2[1:p2, 1:p2]
Sigma2_12 <- Sigma2[1:p2, (p2+1):p]
Sigma2_21 <- Sigma2[(p2+1):p, 1:p2]
Sigma2_22 <- Sigma2[(p2+1):p, (p2+1):p]

Sigma2_112 <- Sigma2_11 - Sigma2_12 %*% solve(Sigma2_22) %*% Sigma2_21

start_time <- Sys.time()
T1 <- Inddist(p1, n, p, sim)
T2 <- Inddist(p2, n, p, sim)
end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time1:", dtime))

q051 <- quantile(T1, 0.05)
q052 <- quantile(T2, 0.05)

T1_1 <- c()
T1_2 <- c()

start_time <- Sys.time()
for (i in 1:sim) {

  # Generate original data sample from multivariate normal with mu and Sigma

  X1 <- rmvnorm(n, mu, Sigma1)
  X2 <- rmvnorm(n, mu, Sigma2)

  mean1 <- colMeans(X1)
  mean2 <- colMeans(X2)

  S1 <- t(X1 - mean1) %*% (X1 - mean1)
  S2 <- t(X2 - mean2) %*% (X2 - mean2)

  # Generate PLS synthetic single data

  V1 <- rmvnorm(n, mean1, S1 / (n - 1))
  V2 <- rmvnorm(n, mean2, S2 / (n - 1))

  # PLS estimates of mu and Sigma

  meanV1 <- colMeans(V1)
  meanV2 <- colMeans(V2)

  S_star1 <- t(V1 - meanV1) %*% (V1 - meanV1)
  S_star2 <- t(V2 - meanV2) %*% (V2 - meanV2)

  S_star1_11 <- S_star1[1:p1, 1:p1]
  S_star1_12 <- S_star1[1:p1, (p1+1):p]
  S_star1_21 <- S_star1[(p1+1):p, 1:p1]
  S_star1_22 <- S_star1[(p1+1):p, (p1+1):p]

  S_star1_112 <- S_star1_11 - S_star1_12 %*% solve(S_star1_22) %*% S_star1_21

  S_star2_11 <- S_star2[1:p2, 1:p2]
  S_star2_12 <- S_star2[1:p2, (p2+1):p]
  S_star2_21 <- S_star2[(p2+1):p, 1:p2]
  S_star2_22 <- S_star2[(p2+1):p, (p2+1):p]

  S_star2_112 <- S_star2_11 - S_star2_12 %*% solve(S_star2_22) %*% S_star2_21

  T1temp <- det(S_star1) / (det(S_star1_11) * det(S_star1_22))
  T2temp <- det(S_star2) / (det(S_star2_11) * det(S_star2_22))

  T1_1 <- c(T1_1, T1temp)
  T1_2 <- c(T1_2, T2temp)
}

end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time:", dtime))

print(paste("T1_1:", min(T1_1), quantile(T1_1, 0.10), quantile(T1_1, 0.50), quantile(T1_1, 0.90), max(T1_1)))
print(paste("T1_2:", min(T1_2), quantile(T1_2, 0.10), quantile(T1_2, 0.50), quantile(T1_2, 0.90), max(T1_2)))
print(paste("T1:", min(T1), quantile(T1, 0.10), quantile(T1, 0.50), quantile(T1, 0.90), max(T1)))
print(paste("T2:", min(T2), quantile(T2, 0.10), quantile(T2, 0.50), quantile(T2, 0.90), max(T2)))

plt1 <- ggplot(data.frame(T1 = T1), aes(x = T1)) +
  geom_density(color = "blue") +
  labs(title = "KDE of T1")

plt2 <- ggplot(data.frame(T2 = T2), aes(x = T2)) +
  geom_density(color = "blue") +
  labs(title = "KDE of T2")

plt3 <- ggplot(data.frame(T1_1 = T1_1), aes(x = T1_1)) +
  geom_density(color = "green") +
  labs(title = "KDE of T1_1")

plt4 <- ggplot(data.frame(T1_2 = T1_2), aes(x = T1_2)) +
  geom_density(color = "yellow") +
  labs(title = "KDE of T1_2")

plt <- ggplot() +
  geom_density(data = data.frame(T = T1), aes
