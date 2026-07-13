## ============================================================================
## R code: Classical tests on the original Brittany soil data
## Produces the "Original" column in Table (tab:realdata) of Section 6.3
##
## Tests performed on X = brittany_soil_ps  (n=37, p=6, standardised)
## Partition: block 1 = {pH_water, pH_KCl}           (p1 = 2)
##            block 2 = {log_CEC_Metson, log_Organic_C,
##                       log_Total_N, log_P_Olsen}   (p2 = 4)
## ============================================================================

library(PSinference)

# ── 0. Load and inspect data ─────────────────────────────────────────────────
X <- brittany_soil_ps
n <- nrow(X)   # 37
p <- ncol(X)   # 6
p1 <- 2
p2 <- p - p1   # 4

cat("Dimensions: n =", n, "| p =", p, "\n\n")

# ── 1. Sample covariance matrix (plug-in parameter) ──────────────────────────
S_hat <- cov(X)                  # p x p; equals correlation matrix (standardised)
cat("Sample covariance matrix (Sigma_hat):\n")
print(round(S_hat, 4))

# Block partition
blk  <- partition(S_hat, part1 = c("pH_water", "pH_KCl"))
S11  <- blk[["A"]]               # 2 x 2
S12  <- blk[["B"]]               # 2 x 4
S21  <- blk[["C"]]               # 4 x 2
S22  <- blk[["D"]]               # 4 x 4

# ── 2. Null-hypothesis reference values ──────────────────────────────────────

## GV: sigma_0 = |Sigma_hat|
sigma0_gv  <- det(S_hat)
cat("\n--- Null-hypothesis reference values ---\n")
cat("GV    sigma_0 = |Sigma_hat| =", formatC(sigma0_gv, format = "e", digits = 4), "\n")

## Sphericity: sigma^2 = tr(S_hat)/p  (= 1 for standardised data)
sigma2_sph <- sum(diag(S_hat)) / p
cat("Sph   sigma^2 = tr(S_hat)/p =", round(sigma2_sph, 6), "\n")

## Regression: Delta_0 = S12 %*% solve(S22)
Delta0 <- S12 %*% solve(S22)
cat("Reg   Delta_0 = S12 * S22^{-1}:\n")
print(round(Delta0, 4))

# ── 3. Classical test on original data ───────────────────────────────────────

## ── 3a. Generalized variance ─────────────────────────────────────────────────
## H0: |Sigma| = |Sigma_hat|  => stat = |S_hat| itself; p-value = 1 by construction.
## (The plug-in null is the MLE, so the LR statistic equals 1.)
gv_stat  <- det(S_hat)
gv_pval  <- 1.0
cat("\n--- Test 1: Generalized Variance ---\n")
cat("Stat = |S_hat| =", formatC(gv_stat, format = "e", digits = 4), "\n")
cat("p-value        =", gv_pval, "(H0 is the plug-in MLE => LR = 1)\n")

## ── 3b. Sphericity ───────────────────────────────────────────────────────────
## Bartlett-Box test: H0: Sigma = sigma^2 * I_p
## W  = |S_hat| / (tr(S_hat)/p)^p
## chi2 = -[(n-1) - (2p^2+p+2)/(6p)] * log(W)
## df  = p*(p+1)/2 - 1

W_sph    <- det(S_hat) / (sum(diag(S_hat)) / p)^p
c_sph    <- (n - 1) - (2 * p^2 + p + 2) / (6 * p)
chi2_sph <- -c_sph * log(W_sph)
df_sph   <- p * (p + 1) / 2 - 1
pval_sph <- pchisq(chi2_sph, df = df_sph, lower.tail = FALSE)

cat("\n--- Test 2: Sphericity (Bartlett-Box) ---\n")
cat("W   =", round(W_sph,    6), "\n")
cat("chi2=", round(chi2_sph, 3), "| df =", df_sph,
    "| p-value =", formatC(pval_sph, format = "e", digits = 3), "\n")

## ── 3c. Independence ─────────────────────────────────────────────────────────
## Bartlett factored-likelihood test: H0: Sigma_12 = 0
## Lambda = |S_hat| / (|S11| * |S22|)
## chi2   = -[(n-1) - (p+3)/2] * log(Lambda)
## df     = p1 * p2

Lambda_ind  <- det(S_hat) / (det(S11) * det(S22))
c_ind       <- (n - 1) - (p + 3) / 2
chi2_ind    <- -c_ind * log(Lambda_ind)
df_ind      <- p1 * p2
pval_ind    <- pchisq(chi2_ind, df = df_ind, lower.tail = FALSE)

cat("\n--- Test 3: Independence (Bartlett factored-likelihood) ---\n")
cat("Lambda =", round(Lambda_ind, 6), "\n")
cat("chi2   =", round(chi2_ind,   3), "| df =", df_ind,
    "| p-value =", formatC(pval_ind, format = "e", digits = 3), "\n")

## ── 3d. Regression ───────────────────────────────────────────────────────────
## H0: Delta = Delta_0 = S12 * S22^{-1}
## Wilks Lambda = |S11.2| / |S11|   where S11.2 = S11 - S12*S22^{-1}*S21
##                                  and numerator has Delta = Delta_0
##
## Since Delta_0 is the plug-in MLE, the numerator of the Wilks statistic
## equals |S11.2| exactly, and the residual Schur complement is also S11.2.
## => The test statistic T*_4 = 0, Wilks Lambda = 1, F = 0 when Delta=Delta_0.
##
## For a meaningful comparison we test H0: Delta = 0 (zero regression),
## which is the natural classical counterpart assessed on the original data.

S11_2    <- S11 - S12 %*% solve(S22) %*% S21   # Schur complement

## Wilks Lambda for H0: Delta = 0
## Numerator: |S11.2|  (residual SS after regressing block1 on block2)
## Denominator: |S11|  (total SS of block1)
Lambda_reg <- det(S11_2) / det(S11)

## F-approximation (Rao 1951 / Anderson 1984 §8.4):
## s = sqrt((p1^2*p2^2-4)/(p1^2+p2^2-5))   [= p2 when p1=1 or p2=1]
## df1 = p1*p2
## df2 = s*((n-1) - (p1+p2+1)/2) - (p1*p2-2)/2 - 1
## F   = ((1 - Lambda^{1/s}) / Lambda^{1/s}) * (df2 / df1)

if (p1^2 + p2^2 - 5 > 0) {
  s_reg <- sqrt((p1^2 * p2^2 - 4) / (p1^2 + p2^2 - 5))
} else {
  s_reg <- 1
}
df1_reg  <- p1 * p2
df2_reg  <- s_reg * ((n - 1) - (p1 + p2 + 1) / 2) - (p1 * p2 - 2) / 2 - 1
Lam_s    <- Lambda_reg^(1 / s_reg)
F_reg    <- ((1 - Lam_s) / Lam_s) * (df2_reg / df1_reg)
pval_reg <- pf(F_reg, df1 = df1_reg, df2 = df2_reg, lower.tail = FALSE)

cat("\n--- Test 4: Regression H0: Delta = 0 (Wilks F-approximation) ---\n")
cat("(Note: H0: Delta=Delta_0 with Delta_0=MLE gives stat=0, p=1 by construction)\n")
cat("(Reported here: H0: Delta=0, the natural classical comparison)\n")
cat("Lambda =", round(Lambda_reg, 6), "\n")
cat("F      =", round(F_reg,      3),
    "| df1 =", round(df1_reg, 1),
    "| df2 =", round(df2_reg, 1),
    "| p-value =", round(pval_reg, 4), "\n")

## ── 4. Summary table ─────────────────────────────────────────────────────────
cat("\n")
cat("========================================================\n")
cat(" SUMMARY: Classical tests on original data (n=37, p=6)\n")
cat("========================================================\n")
cat(sprintf("%-14s  %-20s  %s\n", "Test", "Statistic", "p-value"))
cat("--------------------------------------------------------\n")
cat(sprintf("%-14s  %-20s  %s\n",
            "GV",
            paste0("|S_hat| = ", formatC(gv_stat, format="e", digits=3)),
            "1.000 (plug-in null)"))
cat(sprintf("%-14s  %-20s  %s\n",
            "Sphericity",
            paste0("chi2 = ", round(chi2_sph, 1), " (df=", df_sph, ")"),
            formatC(pval_sph, format="e", digits=3)))
cat(sprintf("%-14s  %-20s  %s\n",
            "Independence",
            paste0("chi2 = ", round(chi2_ind, 1), " (df=", df_ind, ")"),
            formatC(pval_ind, format="e", digits=3)))
cat(sprintf("%-14s  %-20s  %s\n",
            "Regression",
            paste0("F = ", round(F_reg, 3),
                   " (df1=", df1_reg, ", df2=", round(df2_reg,1), ")"),
            round(pval_reg, 4)))
cat("========================================================\n")
cat("Note: Regression tests H0: Delta=0 on original data.\n")
cat("      PS regression tests H0: Delta=Delta_0 on synthetic data.\n")
