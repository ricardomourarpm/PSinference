## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 8,
  fig.height = 6.5,
  warning    = FALSE,
  message    = FALSE
)
library(PSinference)
library(MASS)
set.seed(2026)

## ----load_data----------------------------------------------------------------
data(brittany_soil_ps)
cat("Dimensions:", nrow(brittany_soil_ps), "x", ncol(brittany_soil_ps), "\n")
cat("Variables :", colnames(brittany_soil_ps), "\n\n")

round(head(brittany_soil_ps, 4), 3)

## ----corr---------------------------------------------------------------------
round(cor(brittany_soil_ps), 2)

## ----mvn, fig.cap = "Multivariate normality diagnostic panel for `brittany_soil_ps`. Histograms (blue = Shapiro-Wilk test not rejected) with fitted normal curves, and chi-square Q-Q plot.", fig.height = 6, fig.width = 6.9----
mvn_result <- mvn_test(brittany_soil_ps, hz_nsim = 1000)

## ----mvn_summary--------------------------------------------------------------
print(mvn_result)

## ----gen_synth----------------------------------------------------------------
set.seed(2026)
X <- brittany_soil_ps
n <- nrow(X) # 37
p <- ncol(X) # 6

# M = 1: single release
V1 <- simSynthData(X)
cat("Single release: ", nrow(V1), "x", ncol(V1), "\n")

# M = 3: three stacked releases
V3 <- simSynthData(X, M = 3)
cat("Three releases: ", nrow(V3), "x", ncol(V3), "\n")

## ----gv, fig.cap = "GV null distribution (log10 scale) for M = 1 (left) and M = 3 (right). The null distribution narrows with M; the observed statistic falls within the non-rejection region in both cases.", fig.height = 4, fig.width = 6.9----
op <- par(mfrow = c(1, 2))
gv1 <- gv_test(V1, M = 1, Sigma = cov(X), iterations = 5000)
plot(gv1, main = "Gen. Variance  (M = 1)")
gv3 <- gv_test(V3, M = 3, Sigma = cov(X), iterations = 5000)
plot(gv3, main = "Gen. Variance  (M = 3)")
par(op)

## ----gv_print-----------------------------------------------------------------
cat("M = 1: p =", gv1$p.value, "|", gv1$decision, "\n")
cat("M = 3: p =", gv3$p.value, "|", gv3$decision, "\n")

## ----sph, fig.cap = "Sphericity null distribution for M = 1 (left) and M = 3 (right). The observed statistic is far into the left rejection region in both cases; larger M concentrates the null distribution.", fig.height = 4, fig.width = 6.9----
op <- par(mfrow = c(1, 2))
sph1 <- sphericity_test(V1, M = 1, iterations = 5000)
plot(sph1, main = "Sphericity  (M = 1)")
sph3 <- sphericity_test(V3, M = 3, iterations = 5000)
plot(sph3, main = "Sphericity  (M = 3)")
par(op)

## ----sph_print----------------------------------------------------------------
cat(
  "M = 1: stat =", round(sph1$statistic, 4),
  "| p =", sph1$p.value, "|", sph1$decision, "\n"
)
cat(
  "M = 3: stat =", round(sph3$statistic, 4),
  "| p =", sph3$p.value, "|", sph3$decision, "\n"
)

## ----ind, fig.cap = "Independence null distribution for M = 1 (left) and M = 3 (right). Rejection is clear in both cases; the critical value moves right as M increases.", fig.height = 4, fig.width = 6.7----
op <- par(mfrow = c(1, 2))
ind1 <- independence_test(V1,
  M = 1,
  group_a = c("pH_water", "pH_KCl", "log_CEC_Metson"),
  group_b = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
  iterations = 5000
)
plot(ind1, main = "Independence  (M = 1)")
ind3 <- independence_test(V3,
  M = 3,
  group_a = c("pH_water", "pH_KCl", "log_CEC_Metson"),
  group_b = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
  iterations = 5000
)
plot(ind3, main = "Independence  (M = 3)")
par(op)

## ----ind_print----------------------------------------------------------------
cat(
  "M = 1: stat =", round(ind1$statistic, 4),
  "| p =", ind1$p.value, "|", ind1$decision, "\n"
)
cat(
  "M = 3: stat =", round(ind3$statistic, 4),
  "| p =", ind3$p.value, "|", ind3$decision, "\n"
)

## ----reg, fig.cap = "Regression null distribution (log10 scale) for M = 1 (left) and M = 3 (right). Both reject the zero-regression null.", fig.height = 4, fig.width = 6.9----
op <- par(mfrow = c(1, 2))
reg1 <- regression_test(V1,
  M = 1,
  response = c("pH_water", "pH_KCl", "log_CEC_Metson"),
  predictors = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
  Delta0 = matrix(0, 3, 3),
  iterations = 5000
)
plot(reg1, main = "Regression  (M = 1)")
reg3 <- regression_test(V3,
  M = 3,
  response = c("pH_water", "pH_KCl", "log_CEC_Metson"),
  predictors = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
  Delta0 = matrix(0, 3, 3),
  iterations = 5000
)
plot(reg3, main = "Regression  (M = 3)")
par(op)

## ----reg_print----------------------------------------------------------------
cat(
  "M = 1: stat =", round(reg1$statistic, 5),
  "| p =", reg1$p.value, "|", reg1$decision, "\n"
)
cat(
  "M = 3: stat =", round(reg3$statistic, 5),
  "| p =", reg3$p.value, "|", reg3$decision, "\n"
)

## ----M_effect, fig.cap = "Sphericity null distribution for M = 1, 2, 3, 5. The distribution concentrates as M grows, lowering the critical value and strengthening evidence against the false null.", fig.height = 5, fig.width = 6.9----
op <- par(mfrow = c(2, 2))
for (m in c(1L, 2L, 3L, 5L)) {
  Vm <- simSynthData(X, M = m)
  res <- sphericity_test(Vm, M = m, iterations = 2000)
  plot(res, main = sprintf("Sphericity: M = %d  (N = %d)", m, m * n))
}
par(op)

## ----M_table------------------------------------------------------------------
# Numerical summary across M
results <- lapply(c(1L, 2L, 3L, 5L), function(m) {
  Vm <- simSynthData(X, M = m)
  s <- sphericity_test(Vm, M = m, iterations = 3000)
  i <- independence_test(Vm,
    M = m,
    group_a = c("pH_water", "pH_KCl", "log_CEC_Metson"),
    group_b = c("log_Organic_C", "log_Total_N", "log_P_Olsen"),
    iterations = 3000
  )
  data.frame(
    M = m,
    N = m * n,
    sph_stat = round(s$statistic, 5),
    sph_pval = round(s$p.value, 5),
    ind_stat = round(i$statistic, 5),
    ind_pval = round(i$p.value, 5)
  )
})
do.call(rbind, results)

