#' @title Multivariate Normality Assessment
#'
#' @description
#' Assesses multivariate normality of a data set using five complementary
#' approaches: (1) univariate Shapiro-Wilk tests on each variable,
#' (2) Mardia's multivariate skewness test, (3) Mardia's multivariate
#' kurtosis test, (4) the Henze-Zirkler omnibus test, and (5) Royston's
#' multivariate extension of the Shapiro-Wilk test. A visual diagnostic
#' panel shows one histogram with a fitted normal curve for each variable
#' and a chi-square Q-Q plot of squared Mahalanobis distances.
#'
#' @param X A numeric matrix or data frame with dimension \eqn{n \times p}.
#' @param alpha Significance level for all tests. The default is \code{0.05}.
#' @param plot Logical. If \code{TRUE}, the default, a diagnostic plot
#' panel is produced.
#' @param hz_nsim Integer. Number of Monte Carlo draws used to
#' calibrate the Henze-Zirkler null distribution. The default is \code{2000L}.
#' Increasing this value gives more accurate \eqn{p}-values at the cost
#' of additional computation time.
#' @param verbose Logical. If \code{TRUE}, the default, a formatted summary of
#' all test results is printed.
#'
#' @return
#' A list of class \code{mvn_test}, returned invisibly, with components:
#' \describe{
#'   \item{shapiro}{Data frame of per-variable Shapiro-Wilk statistics
#'     and \eqn{p}-values.}
#'   \item{mardia_skewness}{Named list with components \code{statistic},
#'        \code{df}, \code{p.value}, and \code{decision}.}
#'   \item{mardia_kurtosis}{Named list with components \code{statistic},
#'     \code{p.value}, \code{decision}.}
#'   \item{henze_zirkler}{Named list Named list with components
#'       \code{statistic}, \code{p.value}, \code{decision}.}
#'   \item{royston}{Named list Named list with components \code{statistic}, the
#'        Royston \eqn{H} statistic; \code{df}, the effective degrees of freedom;
#'        \code{p.value}, and \code{decision}.}
#'   \item{mahal_distances}{Numeric vector of squared Mahalanobis
#'     distances.}
#'   \item{overall}{Character string giving the overall conclusion based on
#'     all tests.}
#' }
#'
#' @details
#' \strong{Mardia's skewness test} evaluates the null hypothesis of zero
#' multivariate skewness:
#' \deqn{
#'   \kappa = \frac{n}{6} b_{1,p}
#'   \sim \chi^2\!\left(\frac{p(p+1)(p+2)}{6}\right),
#'   \qquad
#'   b_{1,p}
#'   =
#'   \frac{1}{n^2}
#'   \sum_{a=1}^n
#'   \sum_{b=1}^n
#'   d_{ab}^3.
#' }
#' Here
#' \deqn{
#'   d_{ab}
#'   =
#'   (\boldsymbol{x}_a-\bar{\boldsymbol{x}})'
#'   S^{-1}
#'   (\boldsymbol{x}_b-\bar{\boldsymbol{x}})
#' }
#' is a Mahalanobis inner product between observations \eqn{a} and
#' \eqn{b}. The indices \eqn{a} and \eqn{b} run over observations,
#' not variables.
#'
#' \strong{Mardia's kurtosis test} evaluates whether the multivariate
#' kurtosis equals \eqn{p(p+2)}:
#' \deqn{
#'   z =
#'   \frac{b_{2,p} - p(p+2)}
#'   {\sqrt{8p(p+2)/n}}
#'   \sim N(0,1),
#'   \qquad
#'   b_{2,p}
#'   =
#'   \frac{1}{n}
#'   \sum_{a=1}^n
#'   d_{aa}^2.
#' }
#' The quantity \eqn{d_{aa}} is the squared Mahalanobis distance of
#' observation \eqn{a} from the sample mean.
#'
#' \strong{The Henze-Zirkler omnibus test} is based on a weighted
#' \eqn{L^2} distance between the empirical and theoretical
#' multivariate normal characteristic functions. The statistic is
#' \deqn{
#'   \mathrm{HZ} = \frac{1}{n}\sum_{i=1}^n\sum_{j=1}^n
#'     e^{-\frac{\beta^2}{2}\|\bm{x}_i-\bm{x}_j\|^2_S}
#'   - 2(1+\beta^2)^{-p/2}\frac{1}{n}\sum_{i=1}^n
#'     e^{-\frac{\beta^2}{2(1+\beta^2)}d_i^2}
#'   + (1+2\beta^2)^{-p/2},
#' }
#' where
#' \deqn{
#'    \beta =
#'      \frac{1}{\sqrt{2}}
#'        \left(\frac{2p+1}{4}\right)^{1/(p+4)}
#'        n^{1/(p+4)}.
#' }
#' The null distribution of \eqn{\mathrm{HZ}} is
#' approximated by a log-normal distribution whose parameters are
#' estimated by Monte Carlo simulation of size \code{hz_nsim} from
#' \eqn{\mathcal{N}_p(\bm{0}, \bm{I}_p)}. This test is particularly
#' powerful useful against heavy-tailed and skewed alternatives.
#'
#' \strong{Royston's H test} extends the univariate Shapiro-Wilk statistic
#' to the multivariate setting. For each variable, the Shapiro-Wilk
#' \eqn{p}-value \eqn{p_j} is transformed to \eqn{Z_j = \Phi^{-1}(1 - p_j)}.
#' The test statistic is
#' \deqn{
#'   H = e^{-1} \sum_{j=1}^p Z_j^2 \sim \chi^2_e,
#' }
#' where \eqn{e = p / \bigl[1 + (p-1)\hat\rho_z\bigr]} is an effective
#' degree of freedom parameter that accounts for correlation among the
#' \eqn{Z_j} values. The quantity \eqn{\hat\rho_z} estimated from the average
#' squared pairwise correlation \eqn{\hat\rho_z} of the original variables.
#'
#' \strong{Diagnostic panel:} The diagnostic panel contains \eqn{p + 1} plots
#' arranged in a grid. The first \eqn{p} panels show histograms with fitted
#' \eqn{\mathcal{N}(\bar{x}_j, s_j^2)} density curves. The bar color is
#' steel-blue when the Shapiro-Wilk test fails to reject normality and tomato-red
#' when it rejects. The final panel shows the chi-square Q-Q plot of squared
#' Mahalanobis distances.
#'
#' @references
#' Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis
#' with applications. \emph{Biometrika}, \strong{57}, 519--530.
#'
#' Henze, N. and Zirkler, B. (1990). A class of invariant consistent tests for
#' multivariate normality. \emph{Communications in Statistics: Theory
#' and Methods}, \strong{19}, 3595--3617.
#'
#' Royston, J. P. (1992). Approximating the Shapiro-Wilk W test for
#' non-normality. \emph{Statistics and Computing}, \strong{2}, 117--119.
#'
#' @seealso
#' \code{\link{simSynthData}},
#' \code{\link{ps_attitude}},
#' \code{\link{ps_mtcars}}
#'
#' @export
#'
#' @examples
#' data(attitude)
#' mvn_test(attitude)
#'
#' data(mtcars)
#' mvn_test(mtcars)
mvn_test <- function(X, alpha = 0.05,
                     plot = TRUE,
                     hz_nsim = 2000L,
                     verbose = TRUE) {
  X <- .validate_X(X)
  n <- nrow(X)
  p <- ncol(X)
  vnames <- if (!is.null(colnames(X))) {
    colnames(X)
  } else {
    paste0("V", seq_len(p))
  }

  if (n < 8L) {
    stop("At least 8 observations are needed for MVN testing.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------------
  # 1. Univariate Shapiro-Wilk
  # ------------------------------------------------------------------
  sw_res <- lapply(seq_len(p), function(j) {
    tst <- stats::shapiro.test(X[, j])
    data.frame(
      variable = vnames[j],
      statistic = round(tst$statistic, 4),
      p.value = round(tst$p.value, 4),
      decision = if (tst$p.value < alpha) {
        "Reject H0"
      } else {
        "Fail to Reject H0"
      },
      stringsAsFactors = FALSE
    )
  })
  sw_df <- do.call(rbind, sw_res)
  rownames(sw_df) <- NULL

  # ------------------------------------------------------------------
  # 2 & 3. Mardia's skewness and kurtosis
  # ------------------------------------------------------------------
  Xc <- sweep(X, 2L, colMeans(X), "-")
  # Use MLE covariance (divide by n) for Mardia's statistics
  S_mle <- crossprod(Xc) / n
  Si <- tryCatch(
    solve(S_mle),
    error = function(e) {
      stop("Covariance matrix is singular; cannot run MVN tests.",
        call. = FALSE
      )
    }
  )

  D <- Xc %*% Si %*% t(Xc) # n x n pairwise Mahalanobis products
  d2 <- diag(D) # squared Mahalanobis distances (MLE)

  # Mardia skewness
  b1p <- sum(D^3) / n^2
  kappa <- n * b1p / 6
  df_sk <- p * (p + 1L) * (p + 2L) / 6L
  p_sk <- stats::pchisq(kappa, df = df_sk, lower.tail = FALSE)

  mardia_sk <- list(
    b1p       = round(b1p, 4),
    statistic = round(kappa, 4),
    df        = df_sk,
    p.value   = round(p_sk, 4),
    decision  = if (p_sk < alpha) "Reject H0" else "Fail to Reject H0"
  )

  # Mardia kurtosis
  b2p <- sum(d2^2) / n
  z_ku <- (b2p - p * (p + 2L)) / sqrt(8L * p * (p + 2L) / n)
  p_ku <- 2 * stats::pnorm(abs(z_ku), lower.tail = FALSE)

  mardia_ku <- list(
    b2p       = round(b2p, 4),
    statistic = round(z_ku, 4),
    p.value   = round(p_ku, 4),
    decision  = if (p_ku < alpha) "Reject H0" else "Fail to Reject H0"
  )

  # ------------------------------------------------------------------
  # 4. Henze-Zirkler omnibus test
  #    Uses unbiased S (divide by n-1) as in the original paper.
  # ------------------------------------------------------------------
  S_hz <- crossprod(Xc) / (n - 1L) # unbiased covariance
  Si_hz <- tryCatch(solve(S_hz), error = function(e) NULL)

  hz_result <- if (is.null(Si_hz)) {
    list(
      statistic = NA_real_, p.value = NA_real_,
      decision = "Cannot compute (singular covariance)"
    )
  } else {
    beta <- (1 / sqrt(2)) *
      ((2 * p + 1) / 4)^(1 / (p + 4)) *
      n^(1 / (p + 4))
    b2 <- beta^2

    # Squared Mahalanobis distances using unbiased S
    d2_hz <- stats::mahalanobis(X, colMeans(X), S_hz)

    # Pairwise Mahalanobis distances: ||xi - xj||^2_{S^{-1}}
    # Use Cholesky decomposition for numerical stability
    L <- tryCatch(chol(Si_hz), error = function(e) NULL)
    if (is.null(L)) {
      list(
        statistic = NA_real_, p.value = NA_real_,
        decision = "Cannot compute (Cholesky failed)"
      )
    } else {
      XL <- X %*% t(L) # n x p, xi transformed to ||.||^2 space
      pw_dist <- as.matrix(stats::dist(XL))^2 # pairwise squared distances

      # HZ statistic
      T1 <- mean(exp(-b2 / 2 * pw_dist))
      T2 <- 2 * (1 + b2)^(-p / 2) *
        mean(exp(-b2 / (2 * (1 + b2)) * d2_hz))
      T3 <- (1 + 2 * b2)^(-p / 2)
      HZ <- T1 - T2 + T3

      # Calibrate null distribution by Monte Carlo simulation
      set.seed(42L)
      hz_null <- replicate(as.integer(hz_nsim), {
        Znull <- MASS::mvrnorm(
          n = n,
          mu = rep(0, p),
          Sigma = diag(p)
        )
        Sz <- stats::var(Znull)
        Szi <- solve(Sz)
        d2z <- stats::mahalanobis(Znull, rep(0, p), Sz)
        Lz <- chol(Szi)
        ZnL <- Znull %*% t(Lz)
        pw_z <- as.matrix(stats::dist(ZnL))^2
        t1z <- mean(exp(-b2 / 2 * pw_z))
        t2z <- 2 * (1 + b2)^(-p / 2) *
          mean(exp(-b2 / (2 * (1 + b2)) * d2z))
        t1z - t2z + T3
      })

      mu_sim <- mean(hz_null)
      var_sim <- stats::var(hz_null)
      log_mu <- log(mu_sim^2 / sqrt(var_sim + mu_sim^2))
      log_sig <- sqrt(log(1 + var_sim / mu_sim^2))
      p_hz <- stats::plnorm(HZ, log_mu, log_sig, lower.tail = FALSE)

      list(
        statistic = round(HZ, 5),
        p.value   = round(p_hz, 4),
        beta      = round(beta, 4),
        decision  = if (p_hz < alpha) "Reject H0" else "Fail to Reject H0"
      )
    }
  }

  # ------------------------------------------------------------------
  # 5. Royston's H test (multivariate Shapiro-Wilk extension)
  #    Royston (1992) Statistics and Computing 2:117-119
  # ------------------------------------------------------------------
  sw_pvals <- sw_df$p.value

  # Protect against numerical boundary values
  sw_pvals_safe <- pmin(pmax(sw_pvals, 1e-6), 1 - 1e-6)

  # Transform to z-scores: Z_j = Phi^{-1}(1 - p_j)
  z_vec <- stats::qnorm(1 - sw_pvals_safe)

  # Effective degrees of freedom: e = p / (1 + (p-1)*rho_z)
  # where rho_z = average squared pairwise correlation of original vars
  if (p > 1L) {
    R2_sum <- sum(stats::cor(X)^2) - p # off-diagonal sum of r^2
    rho_z <- R2_sum / (p * (p - 1L))
    e_df <- p / (1 + (p - 1L) * rho_z)
  } else {
    rho_z <- 0
    e_df <- 1
  }

  # H statistic and chi-square p-value
  H_stat <- e_df * mean(z_vec^2)
  p_H <- stats::pchisq(H_stat, df = e_df, lower.tail = FALSE)

  royston <- list(
    statistic = round(H_stat, 4),
    df        = round(e_df, 4),
    p.value   = round(p_H, 4),
    decision  = if (p_H < alpha) "Reject H0" else "Fail to Reject H0"
  )

  # ------------------------------------------------------------------
  # 6. Overall conclusion
  # ------------------------------------------------------------------
  test_pvals <- c(
    min(sw_df$p.value),
    mardia_sk$p.value,
    mardia_ku$p.value,
    if (!is.na(hz_result$p.value)) hz_result$p.value else 1,
    royston$p.value
  )

  n_reject <- sum(test_pvals < alpha)
  overall <- if (n_reject == 0L) {
    "Multivariate normality not rejected"
  } else if (n_reject >= 3L) {
    "Multivariate normality rejected (multiple tests)"
  } else {
    "Multivariate normality marginal (some tests reject)"
  }

  # ------------------------------------------------------------------
  # 7. Diagnostic plots
  # ------------------------------------------------------------------
  if (plot) {
    n_panels <- p + 1L
    nc <- if (p <= 3L) p + 1L else ceiling(sqrt(n_panels))
    nr <- ceiling(n_panels / nc)

    op <- graphics::par(
      mfrow = c(nr, nc),
      mar   = c(3.8, 3.8, 3.8, 1.0),
      mgp   = c(2.3, 0.6, 0),
      oma   = c(0, 0, 2.8, 0)
    )
    on.exit(graphics::par(op), add = TRUE)

    pass_col <- grDevices::adjustcolor("steelblue", 0.55)
    fail_col <- grDevices::adjustcolor("tomato", 0.55)
    crv_col <- "firebrick"

    # (a) Histogram + normal curve per variable
    for (j in seq_len(p)) {
      xj <- X[, j]
      sw_pval_j <- sw_df$p.value[j]
      sw_pass_j <- sw_pval_j > alpha
      bar_col <- if (sw_pass_j) pass_col else fail_col

      h <- graphics::hist(xj, plot = FALSE, breaks = "Sturges")
      y_max <- max(
        h$density,
        stats::dnorm(
          stats::median(xj),
          mean(xj), stats::sd(xj)
        )
      ) * 1.15

      graphics::hist(xj,
        freq   = FALSE,
        breaks = "Sturges",
        col    = bar_col,
        border = "white",
        main   = "",
        xlab   = "",
        ylab   = "Density",
        ylim   = c(0, y_max),
        las    = 1L
      )

      xseq <- seq(min(xj) - diff(range(xj)) * 0.1,
        max(xj) + diff(range(xj)) * 0.1,
        length.out = 200L
      )
      graphics::lines(xseq,
        stats::dnorm(xseq,
          mean = mean(xj),
          sd   = stats::sd(xj)
        ),
        col = crv_col, lwd = 2.2
      )

      graphics::mtext(vnames[j],
        side = 3L, line = 1.75,
        cex = 0.90, col = "black", font = 2L
      )
      graphics::mtext(sprintf("SW  p = %s", .fmt_p(sw_pval_j, 3)),
        side = 3L, line = 0.45,
        cex = 0.78,
        col = if (sw_pass_j) "steelblue4" else "firebrick",
        font = 2L
      )
    }

    # (b) Chi-square Q-Q plot (using MLE Mahalanobis distances)
    q_theo <- stats::qchisq(stats::ppoints(n), df = p)
    q_obs <- sort(d2)
    qq_col <- grDevices::adjustcolor("steelblue", 0.75)

    graphics::plot(q_theo, q_obs,
      main = "",
      xlab = expression(chi[p]^2 ~ "quantiles"),
      ylab = expression("Mahalanobis" ~ d[i]^2),
      pch = 16L, col = qq_col, las = 1L, cex = 0.85
    )
    graphics::abline(0, 1, col = "firebrick", lwd = 2L, lty = 2L)

    graphics::mtext(expression("Chi-square Q-Q  (" * d[i]^2 * ")"),
      side = 3L, line = 1.75,
      cex = 0.90, col = "black", font = 2L
    )

    # Summarize Mardia results in the Q-Q subtitle
    graphics::mtext(
      sprintf("Mardia: skew p=%s  kurt p=%s", .fmt_p(p_sk, 3), .fmt_p(p_ku, 3)),
      side = 3L, line = 0.45, cex = 0.75,
      col = if (p_sk > alpha && p_ku > alpha) {
        "steelblue4"
      } else {
        "firebrick"
      },
      font = 2L
    )

    # Overall title
    graphics::mtext(
      sprintf(
        "MVN Assessment: %s  (n=%d, p=%d)",
        overall, n, p
      ),
      outer = TRUE, side = 3L, line = 0.8,
      cex = 0.88, font = 2L,
      col = if (n_reject == 0L) "darkgreen" else "firebrick"
    )
  }

  # ------------------------------------------------------------------
  # 8. Print summary
  # ------------------------------------------------------------------
  if (verbose) {
    bar <- strrep("-", 60)
    cat("\nMultivariate Normality Assessment\n")
    cat(bar, "\n")
    cat(sprintf("  n = %d observations | p = %d variables\n", n, p))
    cat(bar, "\n\n")

    cat("1. Univariate Shapiro-Wilk Tests\n")
    sw_df_print <- sw_df
    sw_df_print$p.value <- vapply(sw_df$p.value, .fmt_p, character(1))
    print(sw_df_print, row.names = FALSE)

    cat("\n2. Mardia Skewness Test\n")
    cat(sprintf(
      "   b1p = %.4f | chi-sq(df=%d) = %.4f | p = %s | %s\n",
      mardia_sk$b1p, mardia_sk$df, mardia_sk$statistic,
      .fmt_p(mardia_sk$p.value), mardia_sk$decision
    ))

    cat("\n3. Mardia Kurtosis Test\n")
    cat(sprintf(
      "   b2p = %.4f | z = %.4f | p = %.4f | %s\n",
      mardia_ku$b2p, mardia_ku$statistic,
      mardia_ku$p.value, mardia_ku$decision
    ))

    cat("\n4. Henze-Zirkler Omnibus Test\n")
    if (!is.na(hz_result$statistic)) {
      cat(sprintf(
        "   HZ = %.5f | beta = %.4f | p = %.4f | %s\n",
        hz_result$statistic, hz_result$beta,
        hz_result$p.value, hz_result$decision
      ))
    } else {
      cat(sprintf("   %s\n", hz_result$decision))
    }

    cat("\n5. Royston H Test\n")
    cat(sprintf(
      "   H = %.4f | eff. df = %.2f | p = %.4f | %s\n",
      royston$statistic, royston$df,
      royston$p.value, royston$decision
    ))

    cat(sprintf("\n%s\n", bar))
    cat(sprintf("Overall: %s\n\n", overall))
  }

  # ------------------------------------------------------------------
  # 9. Return
  # ------------------------------------------------------------------
  result <- list(
    shapiro         = sw_df,
    mardia_skewness = mardia_sk,
    mardia_kurtosis = mardia_ku,
    henze_zirkler   = hz_result,
    royston         = royston,
    mahal_distances = d2,
    n               = n,
    p               = p,
    alpha           = alpha,
    overall         = overall
  )
  class(result) <- "mvn_test"
  invisible(result)
}

#' @title Print Method for \code{mvn_test} Objects
#' @exportS3Method print mvn_test
print.mvn_test <- function(x, ...) {
  cat(sprintf(
    "MVN Assessment: %s (n = %d, p = %d)\n",
    x$overall, x$n, x$p
  ))
  invisible(x)
}

#' @title Plot Method for \code{mvn_test} Objects
#'
#' @description
#' Re-draws the chi-square Q-Q diagnostic. For the full histogram panel,
#' call \code{mvn_test(X, plot = TRUE)} on the original data directly.
#'
#' @param x An object of class \code{mvn_test}.
#' @param ... Further arguments (currently ignored).
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @exportS3Method plot mvn_test
plot.mvn_test <- function(x, ...) {
  message(
    "For the full panel call mvn_test(X, plot = TRUE). ",
    "Showing chi-square Q-Q plot only."
  )

  n <- x$n
  p <- x$p
  d2 <- x$mahal_distances
  p_sk <- x$mardia_skewness$p.value
  p_ku <- x$mardia_kurtosis$p.value
  mvn_pass <- p_sk > x$alpha && p_ku > x$alpha

  q_theo <- stats::qchisq(stats::ppoints(n), df = p)
  q_obs <- sort(d2)
  qq_col <- grDevices::adjustcolor("steelblue", 0.75)

  graphics::plot(q_theo, q_obs,
    main = expression("Chi-square Q-Q  (" * d[i]^2 * ")"),
    xlab = expression(chi[p]^2 ~ "quantiles"),
    ylab = expression("Mahalanobis" ~ d[i]^2),
    pch = 16L, col = qq_col, las = 1L, cex = 0.85
  )
  graphics::abline(0, 1, col = "firebrick", lwd = 2L, lty = 2L)
  graphics::mtext(
    sprintf("Mardia: skew p=%.3f  kurt p=%.3f", p_sk, p_ku),
    side = 3L, line = 0.15, cex = 0.75,
    col = if (mvn_pass) "steelblue4" else "firebrick", font = 2L
  )

  invisible(x)
}
