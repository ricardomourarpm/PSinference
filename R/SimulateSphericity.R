library(Matrix)
library(mvtnorm)
library(ggplot2)
set.seed(1234)
sim <- 100000
#Sample size and partition size
n <- 100
#Some population mean
mu <- c(1, 2, 3, 4)
#Number of covariates
p <- length(mu)
#Three different population variances and respective partitions
Sigma1 <- diag(p)
Sigma2 <- 5 * diag(p)
start_time <- Sys.time()
T <- Sphdist(n, p, sim)
end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time1:", dtime))
q05 <- quantile(T, 0.05)
T1_1 <- c()
T1_2 <- c()
start_time <- Sys.time()
for (i in 1:sim) {
# Generate original data sample from multivariate normal with mu and Sigma
  X1 <- rmvnorm(n, mu, Sigma1)
  X2 <- rmvnorm(n, mu, Sigma2)

  mean1 <- colMeans(X1)
  mean2 <- colMeans(X2)

  S1 <- t(X1 - mean1) %% (X1 - mean1)
  S2 <- t(X2 - mean2) %% (X2 - mean2)
#  Generate PLS synthetic single data

  V1 <- rmvnorm(n, mean1, S1 / (n - 1))
  V2 <- rmvnorm(n, mean2, S2 / (n - 1))
 # PLS estimates of mu and Sigma

  meanV1 <- colMeans(V1)
  meanV2 <- colMeans(V2)

  S_star1 <- t(V1 - meanV1) %% (V1 - meanV1)
  S_star2 <- t(V2 - meanV2) %% (V2 - meanV2)

  T1temp <- det(S_star1)^(1/p) / sum(diag(S_star1))
  T2temp <- det(S_star2)^(1/p) / sum(diag(S_star2))

  T1_1 <- c(T1_1, T1temp)
  T1_2 <- c(T1_2, T2temp)
}

end_time <- Sys.time()
dtime <- end_time - start_time
print(paste("time2:", dtime))

print(paste(min(T1_1), quantile(T1_1, 0.1), quantile(T1_1, 0.5), quantile(T1_1, 0.9), max(T1_1)))
print(paste(min(T1_2), quantile(T1_2, 0.1), quantile(T1_2, 0.5), quantile(T1_2, 0.9), max(T1_2)))
print(paste(min(T), quantile(T, 0.1), quantile(T, 0.5), quantile(T, 0.9), max(T)))
Plotting KDE

df <- data.frame(T = T, T1_1 = T1_1, T1_2 = T1_2)

plt1 <- ggplot(df, aes(x = T)) +
  geom_density(color = "blue") +
  labs(title = "KDE of T")

plt2 <- ggplot(df, aes(x = T1_1)) +
  geom_density(color = "green") +
  labs(title = "KDE of T1_1")

plt3 <- ggplot(df, aes(x = T1_2)) +
  geom_density(color = "yellow") +
  labs(title = "KDE of T1_2")

plt <- ggplot(df) +
  geom_density(aes(x = T, color = "T"), cut = 0) +
  geom_density(aes(x = T1_1, color = "T1_1"), cut = 0) +
  geom_density(aes(x = T1_2, color = "T1_2"), cut = 0) +
  labs(title = "KDE Comparison", x = "T") +
  scale_color_manual(values = c("T" = "blue", "T1_1" = "green", "T1_2" = "yellow")) +
  theme(legend.title = element_blank())

print(plt)

cov1 <- mean(T1_1 > q05)
cov2 <- mean(T1_2 > q05)

print(c(cov1, cov2))
