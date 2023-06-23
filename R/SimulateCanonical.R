library(MASS)
library(matrixcalc)

set.seed(1234)
sim <- 100000

# Sample size and partition size
n <- 100
p1 <- 2
p2 <- 1

# Some population mean
mu <- c(1, 2, 3, 4)
# Number of covariates
p <- length(mu)

# Three different population variances and respective partitions

Sigma1 <- matrix(c(1, 0.5, 0.5, 0.5,
                   0.5, 1, 0.5, 0.5,
                   0.5, 0.5, 1, 0.5,
                   0.5, 0.5, 0.5, 1), nrow = p, ncol = p, byrow = TRUE)

Sigma1_11 <- Sigma1[1:p1, 1:p1]
Sigma1_12 <- Sigma1[1:p1, (p1 + 1):p]
Sigma1_21 <- Sigma1[(p1 + 1):p, 1:p1]
Sigma1_22 <- Sigma1[(p1 + 1):p, (p1 + 1):p]

Sigma1_112 <- Sigma1_11 - Sigma1_12 %*% solve(Sigma1_22) %*% Sigma1_21
Delta1 <- Sigma1_12 %*% solve(Sigma1_22)

Sigma2 <- matrix(c(1, 0.5, 0, 0,
                   0.5, 2, 0, 0,
                   0, 0, 3, 0.2,
                   0, 0, 0.2, 4), nrow = p, ncol = p, byrow = TRUE)

Sigma2_11 <- Sigma2[1:p2, 1:p2]
Sigma2_12 <- Sigma2[1:p2, (p2 + 1):p]
Sigma2_21 <- Sigma2[(p2 + 1):p, 1:p2]
Sigma2_22 <- Sigma2[(p2 + 1):p, (p2 + 1):p]

Sigma2_112 <- Sigma2_11 - Sigma2_12 %*% solve(Sigma2_22) %*% Sigma2_21
Delta2 <- Sigma2_12 %*% solve(Sigma2_22)

T1 <- Canodist(p1, n, p, sim)
T2 <- Canodist(p2, n, p, sim)

start_time <- Sys.time()

q951 <- quantile(T1, 0.95)
q952 <- quantile(T2, 0.95)

T1_1 <- NULL
T1_2 <- NULL

for (i in 1:sim) {
  # Generate original data sample from normal with mu and Sigma
  X1 <- mvrnorm(n, mu, Sigma1)
  X2 <- mvrnorm(n, mu, Sigma2)

  mean1 <- colMeans(X1)
  mean2 <- colMeans(X2)

  S1 <- t(X1 - mean1) %*% (X1 - mean1)
  S2 <- t(X2 - mean2) %*% (X2 - mean2)

  # Generate PLS synthetic single data
  V1 <- mvrnorm(n, mean1, S1 / (n - 1))
  V2 <- mvrnorm(n, mean2, S2 / (n - 1))

  # PLS estimates of mu and Sigma
  meanV1 <- colMeans(V1)
  meanV2 <- colMeans(V2)

  S_star1 <- t(V1 - meanV1) %*% (V1 - meanV1)
  S_star2 <- t(V2 - meanV2) %*% (V2 - meanV2)

  S_star1_11 <- S_star1[1:p1, 1:p1]
  S_star1_12 <- S_star1[1:p1, (p1 + 1):p]
  S_star1_21 <- S_star1[(p1 + 1):p, 1:p1]
  S_star1_22 <- S_star1[(p1 + 1):p, (p1 + 1):p]

  S_star1_112 <- S_star1_11 - S_star1_12 %*% solve(S_star1_22)
  Delta_star1 <- S_star1_12 %*% solve(S_star1_22)

  S_star2_11 <- S_star2[1:p2, 1:p2]
  S_star2_12 <- S_star2[1:p2, (p2 + 1):p]
  S_star2_21 <- S_star2[(p2 + 1):p, 1:p2]
  S_star2_22 <- S_star2[(p2 + 1):p, (p2 + 1):p]

  S_star2_112 <- S_star2_11 - S_star2_12 %*% solve(S_star2_22)
  Delta_star2 <- S_star2_12 %*% solve(S_star2_22)

  T1temp <- det((Delta_star1 - Delta1) %*% (S_star1_22 %*% t(Delta_star1 - Delta1))) / det(S_star1_112)
  T2temp <- det((Delta_star2 - Delta2) %*% (S_star2_22 %*% t(Delta_star2 - Delta2))) / det(S_star2_112)

  T1_1 <- c(T1_1, T1temp)
  T1_2 <- c(T1_2, T2temp)
}

end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time1:", dtime))

start_time <- Sys.time()

print(paste(min(T1_1), quantile(T1_1, 0.10), quantile(T1_1, 0.50), quantile(T1_1, 0.90), max(T1_1)))
print(paste(min(T1_2), quantile(T1_2, 0.10), quantile(T1_2, 0.50), quantile(T1_2, 0.90), max(T1_2)))
print(paste(min(T1), quantile(T1, 0.10), quantile(T1, 0.50), quantile(T1, 0.90), max(T1)))
print(paste(min(T2), quantile(T2, 0.10), quantile(T2, 0.50), quantile(T2, 0.90), max(T2)))

plot_data <- data.frame(T1 = T1, T2 = T2, T1_1 = T1_1, T1_2 = T1_2)

par(mfrow = c(3, 2))
hist(T1, main = "T1", col = "blue")
hist(T2, main = "T2", col = "blue")
hist(T1_1, main = "T1_1", col = "green")
hist(T1_2, main = "T1_2", col = "yellow")
legend("topright", legend = c("T1", "T2", "T1_1", "T1_2"), fill = c("blue", "blue", "green", "yellow"))

par(mfrow = c(2, 2))
plot(density(T1^(1/p)), main = "T1^(1/p)", col = "blue", lwd = 2)
lines(density(T1_1^(1/p)), col = "green", lwd = 2)
plot(density(T2^(1/p)), main = "T2^(1/p)", col = "blue", lwd = 2)
lines(density(T1_2^(1/p)), col = "yellow", lwd = 2)

cov1 <- mean(T1_1 < q951)
cov2 <- mean(T1_2 < q952)

print(c(cov1, cov2))
