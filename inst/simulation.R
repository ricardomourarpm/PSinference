# =============================================================================
#  PSinference -- Monte Carlo Simulation Study
#
#  Exact PS procedures and comparison with approximate Reiter-type rule
#
#  Fixes included:
#   1. Correct use of M inside .coverage_cell()
#   2. Correct regression partition using Sigma12 %*% solve(Sigma22)
#   3. Correct upper-tail chi-square p-values for Reiter-type tests
#   4. Reiter synthetic datasets generated with simSynthData(X, M = 1L)
#   5. Exact null distributions precomputed once per cell where possible
# =============================================================================

library(PSinference)
library(MASS)

set.seed(20260101)

## ---------------------------------------------------------------------------
## Simulation parameters
## ---------------------------------------------------------------------------

N_inner <- 1e4L      # increase for better calibration, e.g. 2000 or 5000
N_outer <- 1e4L      # use 1e5 for publication-quality tables
ALPHA   <- 0.05
P       <- 4L
MU      <- 1:P

## Covariance matrices
Sigma1 <- diag(P)

Sigma2 <- 5 * diag(P)

Sigma3 <- matrix(0.5, P, P)
diag(Sigma3) <- 1

Sigma4 <- matrix(c(1, .5,  0,   0,
                   .5, 2,  0,   0,
                   0,  0,  3,  .2,
                   0,  0, .2,   4), P, P)

n_vals <- c(10L, 20L, 50L, 100L, 500L)
M_vals <- c(1L,  2L,  5L,  10L)

## ---------------------------------------------------------------------------
## Helper: compute population regression Delta0
## ---------------------------------------------------------------------------

.compute_Delta0 <- function(Sigma, part) {
  p <- ncol(Sigma)

  if (part <= 0L || part >= p) {
    stop("'part' must satisfy 1 <= part < p.", call. = FALSE)
  }

  idx1 <- seq_len(part)
  idx2 <- seq.int(part + 1L, p)

  Sigma12 <- Sigma[idx1, idx2, drop = FALSE]
  Sigma22 <- Sigma[idx2, idx2, drop = FALSE]

  Sigma12 %*% solve(Sigma22)
}

## ---------------------------------------------------------------------------
## Helper: simulate coverage for one exact PS cell
## ---------------------------------------------------------------------------

.coverage_cell <- function(n, M, Sigma, test,
                           part = 2L,
                           N_sim = N_outer,
                           n_inner = N_inner) {

  test <- match.arg(test, c("gv", "sph", "ind", "reg"))

  ## Pre-compute critical value once
  cv <- switch(
    test,

    gv = {
      nd <- GVdist(
        nsample    = n,
        pvariates  = P,
        M          = M,
        iterations = n_inner
      )

      c(
        quantile(nd, ALPHA / 2),
        quantile(nd, 1 - ALPHA / 2)
      )
    },

    sph = {
      quantile(
        Sphdist(
          nsample    = n,
          pvariates  = P,
          M          = M,
          iterations = n_inner
        ),
        ALPHA
      )
    },

    ind = {
      quantile(
        Inddist(
          part       = part,
          nsample    = n,
          pvariates  = P,
          M          = M,
          iterations = n_inner
        ),
        ALPHA
      )
    },

    reg = {
      quantile(
        canodist(
          part       = part,
          nsample    = n,
          pvariates  = P,
          M          = M,
          iterations = n_inner
        ),
        1 - ALPHA
      )
    }
  )

  ## Population Delta0 for regression null
  if (test == "reg") {
    D0 <- .compute_Delta0(Sigma, part)
  }

  in_ci <- logical(N_sim)

  for (r in seq_len(N_sim)) {

    X <- MASS::mvrnorm(n, MU, Sigma)

    V <- simSynthData(X, M = M)

    S_star <- crossprod(sweep(V, 2, colMeans(V)))

    T_obs <- switch(
      test,

      gv = {
        ## Pivot for generalised variance
        (n - 1L)^P * det(S_star) / det(Sigma)
      },

      sph = {
        ## Sphericity statistic
        det(S_star)^(1 / P) / (sum(diag(S_star)) / P)
      },

      ind = {
        ## Independence statistic
        idx1 <- seq_len(part)
        idx2 <- seq.int(part + 1L, P)

        det(S_star) /
          (
            det(S_star[idx1, idx1, drop = FALSE]) *
              det(S_star[idx2, idx2, drop = FALSE])
          )
      },

      reg = {
        ## Regression statistic
        idx1 <- seq_len(part)
        idx2 <- seq.int(part + 1L, P)

        S11  <- S_star[idx1, idx1, drop = FALSE]
        S12  <- S_star[idx1, idx2, drop = FALSE]
        S22  <- S_star[idx2, idx2, drop = FALSE]
        S22i <- solve(S22)

        dif <- S12 %*% S22i - D0

        det(dif %*% S22 %*% t(dif)) /
          det(S11 - S12 %*% S22i %*% t(S12))
      }
    )

    in_ci[r] <- switch(
      test,
      gv  = T_obs >= cv[1] && T_obs <= cv[2],
      sph = T_obs >= cv,
      ind = T_obs >= cv,
      reg = T_obs <= cv
    )
  }

  round(mean(in_ci), 4)
}

## ---------------------------------------------------------------------------
## Study 1A: Generalised Variance and Sphericity
## ---------------------------------------------------------------------------

cat("\n=== Study 1A: Generalised Variance and Sphericity ===\n")
cat(sprintf("Outer replications: %d   Inner MC: %d\n\n",
            N_outer, N_inner))

## GV: Sigma3 and Sigma4
gv_results <- matrix(
  NA_real_,
  nrow = length(n_vals),
  ncol = length(M_vals) * 2,
  dimnames = list(
    paste0("n=", n_vals),
    c(
      paste0("GV_Sigma3_M", M_vals),
      paste0("GV_Sigma4_M", M_vals)
    )
  )
)

for (ni in seq_along(n_vals)) {
  n <- n_vals[ni]

  for (mi in seq_along(M_vals)) {
    m <- M_vals[mi]

    if (n * m <= P + 1L) next

    cat(sprintf("  GV Sigma3 n=%d M=%d ... ", n, m))
    gv_results[ni, mi] <- .coverage_cell(n, m, Sigma3, "gv")
    cat(sprintf("coverage = %.4f\n", gv_results[ni, mi]))

    cat(sprintf("  GV Sigma4 n=%d M=%d ... ", n, m))
    gv_results[ni, mi + length(M_vals)] <- .coverage_cell(
      n,
      m,
      Sigma4,
      "gv"
    )
    cat(sprintf("coverage = %.4f\n",
                gv_results[ni, mi + length(M_vals)]))
  }
}

## Sphericity: Sigma1 and Sigma2
sph_results <- matrix(
  NA_real_,
  nrow = length(n_vals),
  ncol = length(M_vals) * 2,
  dimnames = list(
    paste0("n=", n_vals),
    c(
      paste0("Sph_Sigma1_M", M_vals),
      paste0("Sph_Sigma2_M", M_vals)
    )
  )
)

for (ni in seq_along(n_vals)) {
  n <- n_vals[ni]

  for (mi in seq_along(M_vals)) {
    m <- M_vals[mi]

    if (n * m <= P + 1L) next

    cat(sprintf("  Sph Sigma1 n=%d M=%d ... ", n, m))
    sph_results[ni, mi] <- .coverage_cell(n, m, Sigma1, "sph")
    cat(sprintf("coverage = %.4f\n", sph_results[ni, mi]))

    cat(sprintf("  Sph Sigma2 n=%d M=%d ... ", n, m))
    sph_results[ni, mi + length(M_vals)] <- .coverage_cell(
      n,
      m,
      Sigma2,
      "sph"
    )
    cat(sprintf("coverage = %.4f\n",
                sph_results[ni, mi + length(M_vals)]))
  }
}

## ---------------------------------------------------------------------------
## Study 1B: Independence and Regression
## ---------------------------------------------------------------------------

cat("\n=== Study 1B: Independence and Regression ===\n")

## Independence: Sigma1 with p1 = 1, Sigma4 with p1 = 2
ind_results <- matrix(
  NA_real_,
  nrow = length(n_vals),
  ncol = length(M_vals) * 2,
  dimnames = list(
    paste0("n=", n_vals),
    c(
      paste0("Ind_Sigma1p1_M", M_vals),
      paste0("Ind_Sigma4p2_M", M_vals)
    )
  )
)

for (ni in seq_along(n_vals)) {
  n <- n_vals[ni]

  for (mi in seq_along(M_vals)) {
    m <- M_vals[mi]

    if (n * m <= P + 1L) next

    cat(sprintf("  Ind Sigma1(p1=1) n=%d M=%d ... ", n, m))
    ind_results[ni, mi] <- .coverage_cell(
      n,
      m,
      Sigma1,
      "ind",
      part = 1L
    )
    cat(sprintf("coverage = %.4f\n", ind_results[ni, mi]))

    cat(sprintf("  Ind Sigma4(p1=2) n=%d M=%d ... ", n, m))
    ind_results[ni, mi + length(M_vals)] <- .coverage_cell(
      n,
      m,
      Sigma4,
      "ind",
      part = 2L
    )
    cat(sprintf("coverage = %.4f\n",
                ind_results[ni, mi + length(M_vals)]))
  }
}

## Regression: Sigma3 with p1 = 2, Sigma4 with p1 = 1
reg_results <- matrix(
  NA_real_,
  nrow = length(n_vals),
  ncol = length(M_vals) * 2,
  dimnames = list(
    paste0("n=", n_vals),
    c(
      paste0("Reg_Sigma3p2_M", M_vals),
      paste0("Reg_Sigma4p1_M", M_vals)
    )
  )
)

for (ni in seq_along(n_vals)) {
  n <- n_vals[ni]

  for (mi in seq_along(M_vals)) {
    m <- M_vals[mi]

    if (n * m <= P + 1L) next

    cat(sprintf("  Reg Sigma3(p1=2) n=%d M=%d ... ", n, m))
    reg_results[ni, mi] <- .coverage_cell(
      n,
      m,
      Sigma3,
      "reg",
      part = 2L
    )
    cat(sprintf("coverage = %.4f\n", reg_results[ni, mi]))

    cat(sprintf("  Reg Sigma4(p1=1) n=%d M=%d ... ", n, m))
    reg_results[ni, mi + length(M_vals)] <- .coverage_cell(
      n,
      m,
      Sigma4,
      "reg",
      part = 1L
    )
    cat(sprintf("coverage = %.4f\n",
                reg_results[ni, mi + length(M_vals)]))
  }
}

## ---------------------------------------------------------------------------
## Study 2: Exact PS vs approximate Reiter-type rule
## ---------------------------------------------------------------------------

cat("\n=== Study 2: Exact PS vs approximate Reiter-type rule ===\n")
cat("(Sphericity test, M in {2,5,10}, n in {10,20,30,50})\n\n")

.reiter_sph_pval <- function(X, M) {

  n <- nrow(X)
  p <- ncol(X)

  ## Generate M separate synthetic datasets using the same PS generator
  V_list <- lapply(seq_len(M), function(j) {
    simSynthData(X, M = 1L)
  })

  S_list <- lapply(V_list, cov)
  S_bar  <- Reduce("+", S_list) / M

  lambda <- det(S_bar) / (sum(diag(S_bar)) / p)^p

  ## Numerical protection
  lambda <- min(max(lambda, .Machine$double.eps), 1)

  ## Sphericity degrees of freedom
  df <- p * (p + 1) / 2 - 1

  ## Approximate chi-square statistic
  chi_sq <- -(
    M * (n - 1) -
      (2 * p^2 + p + 2) / (6 * p)
  ) * log(lambda)

  ## Correct upper-tail p-value
  pchisq(
    chi_sq,
    df = df,
    lower.tail = FALSE
  )
}

n_comp <- c(10L, 20L, 30L, 50L)
M_comp <- c(2L, 5L, 10L)
N_comp <- 2000L

comp_results <- data.frame()

for (n in n_comp) {
  for (m in M_comp) {

    if (n * m <= P + 1L) next

    cat(sprintf("  n=%d M=%d ... ", n, m))

    ## Pre-compute exact null distribution once
    nd_exact <- Sphdist(
      nsample    = n,
      pvariates  = P,
      M          = m,
      iterations = N_inner
    )

    cv_exact <- quantile(nd_exact, ALPHA)

    cvg_exact  <- rep(NA_real_, N_comp)
    cvg_reiter <- rep(NA_real_, N_comp)

    pv_exact   <- rep(NA_real_, N_comp)
    pv_reiter  <- rep(NA_real_, N_comp)

    for (r in seq_len(N_comp)) {

      X <- MASS::mvrnorm(n, MU, Sigma1)

      ## Exact PS stacked synthetic data
      V <- simSynthData(X, M = m)

      S_star <- crossprod(sweep(V, 2, colMeans(V)))

      T_obs <- det(S_star)^(1 / P) /
        (sum(diag(S_star)) / P)

      ## Exact PS p-value
      pv_e <- mean(nd_exact <= T_obs)

      pv_exact[r]  <- pv_e
      cvg_exact[r] <- T_obs >= cv_exact

      ## Approximate Reiter-type p-value
      pv_r <- tryCatch(
        .reiter_sph_pval(X, m),
        error = function(e) NA_real_
      )

      pv_reiter[r]  <- pv_r
      cvg_reiter[r] <- !is.na(pv_r) && pv_r > ALPHA
    }

    row <- data.frame(
      n               = n,
      M               = m,
      N_eff           = n * m,
      exact_coverage  = round(mean(cvg_exact, na.rm = TRUE), 4),
      reiter_coverage = round(mean(cvg_reiter, na.rm = TRUE), 4),
      exact_meanpval  = round(mean(pv_exact, na.rm = TRUE), 4),
      reiter_meanpval = round(mean(pv_reiter, na.rm = TRUE), 4)
    )

    comp_results <- rbind(comp_results, row)

    cat(sprintf(
      "exact=%.4f  reiter=%.4f\n",
      row$exact_coverage,
      row$reiter_coverage
    ))
  }
}

## ---------------------------------------------------------------------------
## Print summary tables
## ---------------------------------------------------------------------------

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("FINAL SUMMARY TABLES\n")
cat(strrep("=", 70), "\n")

MC_SE <- round(sqrt(0.95 * 0.05 / N_outer), 4)

cat(sprintf(
  "Nominal level = %.2f  |  MC SE ≈ %.4f  |  Outer reps = %d\n\n",
  ALPHA,
  MC_SE,
  N_outer
))

cat("Table 1: Generalised Variance coverage  (target 0.95)\n")
print(gv_results)

cat("\nTable 2: Sphericity coverage  (target 0.95)\n")
print(sph_results)

cat("\nTable 3: Independence coverage  (target 0.95)\n")
print(ind_results)

cat("\nTable 4: Regression coverage  (target 0.95)\n")
print(reg_results)

cat("\nTable 5: Exact PS vs approximate Reiter-type rule, sphericity\n")
cat("Coverage probability under H0: Sigma = I4, nominal = 0.95\n\n")
print(comp_results, row.names = FALSE)

## ---------------------------------------------------------------------------
## Save results
## ---------------------------------------------------------------------------

sim_results <- list(
  gv           = gv_results,
  sphericity   = sph_results,
  independence = ind_results,
  regression   = reg_results,
  comparison   = comp_results,
  params       = list(
    N_outer = N_outer,
    N_inner = N_inner,
    N_comp  = N_comp,
    alpha   = ALPHA,
    p       = P,
    mu      = MU,
    M_vals  = M_vals,
    n_vals  = n_vals
  )
)

saveRDS(sim_results, "simulation_results.rds")

cat("\nResults saved to simulation_results.rds\n")
