library(MASS)
library(matrixStats)
library(ggplot2)

set.seed(1234)
sim <- 100000

# Sample size and partition size
n <- 100

# Some population mean
mu <- c(1, 2, 3, 4)
# Number of covariates
p <- length(mu)

# Three different population variances and respective partitions
Sigma1 <- matrix(c(1, 0.5, 0.5, 0.5,
                   0.5, 1, 0.5, 0.5,
                   0.5, 0.5, 1, 0.5,
                   0.5, 0.5, 0.5, 1), nrow = 4, ncol = 4, byrow = TRUE)

Sigma2 <- matrix(c(1, 0.5, 0, 0,
                   0.5, 2, 0, 0,
                   0, 0, 3, 0.2,
                   0, 0, 0.2, 4), nrow = 4, ncol = 4, byrow = TRUE)

start_time <- Sys.time()
T <- GVdist(n, p, sim)
q975 <- quantile(T, 0.975)
q025 <- quantile(T, 0.025)
T1_1 <- c()
T1_2 <- c()
end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time1:", dtime))

start_time <- Sys.time()
for (i in 1:sim) {
  # Generate original data sample from normal with mu and Sigma
  X1 <- mvrnorm(n, mu, Sigma1)
  X2 <- mvrnorm(n, mu, Sigma2)

  mean1 <- colMeans(X1)
  mean2 <- colMeans(X2)

  S1 <- crossprod(X1 - mean1) / (n - 1)
  S2 <- crossprod(X2 - mean2) / (n - 1)

  # Generate PLS synthetic single data
  V1 <- mvrnorm(n, mean1, S1)
  V2 <- mvrnorm(n, mean2, S2)

  # PLS estimates of mu and Sigma
  meanV1 <- colMeans(V1)
  meanV2 <- colMeans(V2)

  S_star1 <- crossprod(V1 - meanV1) / (n - 1)
  S_star2 <- crossprod(V2 - meanV2) / (n - 1)

  T1temp <- (n - 1) ^ p * det(S_star1) / det(Sigma1)
  T2temp <- (n - 1) ^ p * det(S_star2) / det(Sigma2)

  T1_1 <- c(T1_1, T1temp)
  T1_2 <- c(T1_2, T2temp)
}

end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time2:", dtime))

print(paste(min(T1_1), quantile(T1_1, 0.1), quantile(T1_1, 0.5), quantile(T1_1, 0.9), max(T1_1)))
print(paste(min(T1_2), quantile(T1_2, 0.1), quantile(T1_2, 0.5), quantile(T1_2, 0.9), max(T1_2)))
print(paste(min(T), quantile(T, 0.1), quantile(T, 0.5), quantile(T, 0.9), max(T)))

# Plotting
ggplot() +
  geom_density(aes(T ^ (1 / p), color = "T"), cut = 0) +
  geom_density(aes(T1_1 ^ (1 / p), color = "T1_1"), cut = 0) +
  geom_density(aes(T1_2 ^ (1 / p), color = "T1_2"), cut = 0) +
  labs(x = "T ^ (1 / p)", y = "Density") +
  scale_color_manual(values = c("T" = "blue", "T1_1" = "green", "T1_2" = "yellow")) +
  facet_wrap(~ 1, ncol = 2) +
  theme_minimal()

ggplot() +
  geom_density(aes(T, color = "T"), cut = 0) +
  geom_density(aes(T1_1, color = "T1_1"), cut = 0) +
  geom_density(aes(T1_2, color = "T1_2"), cut = 0) +
  labs(x = "T", y = "Density") +
  scale_color_manual(values = c("T" = "blue", "T1_1" = "green", "T1_2" = "yellow")) +
  theme_minimal()

cov1 <- mean(T1_1 < q975 & T1_1 > q025)
cov2 <- mean(T1_2 < q975 & T1_2 > q025)

print(c(cov1, cov2))
