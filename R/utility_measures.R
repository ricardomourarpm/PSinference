#' @title Utility Measures for Plug-in Sampling Synthetic Data
#'
#' @description
#' Computes five complementary utility measures quantifying how well plug-in
#' sampling (PS) synthetic data preserve the statistical properties of the
#' original confidential data. The reported headline measures are Frobenius
#' distance, mean standardized mean difference, variance ratio range,
#' propensity score MSE (pMSE) ratio, and mean confidence interval overlap.
#' For \eqn{M > 1} releases, per-release statistics are computed and averaged
#' where appropriate.
#'
#' @param X A numeric matrix or data frame with dimension \eqn{n \times p},
#'   containing the original confidential data.
#' @param V A numeric matrix or data frame with dimension \eqn{Mn \times p},
#'   containing the stacked PS synthetic data returned by
#'   \code{\link{simSynthData}}.
#' @param M Positive integer giving the number of synthetic releases. The
#'   default is \code{1L}.
#' @param alpha Significance level used for confidence interval overlap. The
#'   default is \code{0.05}.
#' @param verbose Logical. If \code{TRUE}, the default, a formatted summary
#'   is printed.
#'
#' @return
#' A list of class \code{ps_utility}, returned invisibly, with components:
#' \describe{
#'   \item{frobenius}{Numeric. Frobenius distance
#'     \eqn{\|\hat{\boldsymbol{\Sigma}}_X -
#'     \bar{\hat{\boldsymbol{\Sigma}}}_V\|_F}.}
#'   \item{smd_mean}{Numeric. Mean standardized mean difference,
#'     \eqn{p^{-1}\sum_j |\bar{v}_j - \bar{x}_j|/s_{x_j}}.}
#'   \item{var_ratio_range}{Named numeric vector with the minimum and maximum
#'     per-variable variance ratios.}
#'   \item{pmse_ratio}{Numeric. Ratio \code{pmse / pmse_null}.}
#'   \item{ci_overlap_mean}{Numeric. Mean confidence interval overlap across
#'     all variables.}
#'   \item{M}{Integer. Number of synthetic releases.}
#'   \item{n}{Integer. Original sample size.}
#'   \item{p}{Integer. Number of variables.}
#'   \item{alpha}{Numeric. Significance level used.}
#' }
#'
#' @details
#' \subsection{Frobenius distance}{
#' The Frobenius distance measures covariance matrix preservation:
#' \deqn{
#'   d_F =
#'   \|\hat{\boldsymbol{\Sigma}}_X -
#'   \bar{\hat{\boldsymbol{\Sigma}}}_V\|_F.
#' }
#' For \eqn{M > 1}, \eqn{\bar{\hat{\boldsymbol{\Sigma}}}_V} is the
#' average per-release sample covariance matrix.
#' }
#'
#' \subsection{Mean standardized mean difference}{
#' The mean standardized mean difference is
#' \deqn{
#'   \frac{1}{p}
#'   \sum_{j=1}^p
#'   \frac{|\bar{v}_j - \bar{x}_j|}{s_{x_j}}.
#' }
#' Values below 0.10 indicate negligible marginal mean differences.
#' }
#'
#' \subsection{Variance ratio range}{
#' The variance ratio range is
#' \deqn{
#'   \left[
#'   \min_j \frac{s^2_{v_j}}{s^2_{x_j}},
#'   \max_j \frac{s^2_{v_j}}{s^2_{x_j}}
#'   \right].
#' }
#' Values close to 1 indicate good variance preservation.
#' }
#'
#' \subsection{Propensity score MSE ratio}{
#' The original data, labeled 0, and the synthetic data, labeled 1, are
#' combined. A logistic classifier is fitted to distinguish original from
#' synthetic records. The pMSE is
#' \deqn{
#'   \mathrm{pMSE}
#'   =
#'   \frac{1}{n + Mn}
#'   \sum_{i=1}^{n+Mn}
#'   (\hat{p}_i - c)^2,
#'   \qquad
#'   c = \frac{Mn}{n + Mn}.
#' }
#' The expected value under a correctly specified synthesis model is
#' \deqn{
#'   \mathrm{pMSE}_{\mathrm{null}}
#'   =
#'   \frac{c(1-c)(p+1)}{n+Mn}.
#' }
#' The pMSE ratio is
#' \deqn{
#'   \mathrm{pMSE}/\mathrm{pMSE}_{\mathrm{null}}.
#' }
#' Values near 1 indicate good utility, whereas values well below 1 indicate
#' that the synthetic data blend in well with the original data.
#' }
#'
#' \subsection{Mean confidence interval overlap}{
#' For each variable \eqn{j}, a \eqn{(1-\alpha)} confidence interval (CI) is
#' computed from the original data and from the synthetic data. The overlap
#' coefficient is
#' \deqn{
#'   J_j =
#'   \frac{
#'   \max\left[
#'   0,\;
#'   \min(u^X_j, u^V_j) - \max(l^X_j, l^V_j)
#'   \right]
#'   }{
#'   \max(u^X_j - l^X_j,\; u^V_j - l^V_j)
#'   }.
#' }
#' where \eqn{[l^X_j, u^X_j]} and \eqn{[l^V_j, u^V_j]} are the CIs
#' from original and synthetic data. The reported measure is the average
#' of \eqn{J_j} across all variables.
#' }
#'
#' @references
#' Karr, A. F., Kohnen, C. N., Oganian, A., Reiter, J. P., and
#' Sanil, A. P. (2006). A framework for evaluating the utility of data
#' altered to protect confidentiality. \emph{The American Statistician},
#' \strong{60}, 224--232.
#'
#' Snoke, J., Raab, G. M., Nowok, B., Dibben, C., and Slavkovic, A.
#' (2018). General and specific utility measures for synthetic data.
#' \emph{Journal of the Royal Statistical Society: Series A},
#' \strong{181}, 663--688.
#'
#' Woo, M.-J., Reiter, J. P., Oganian, A., and Karr, A. F. (2009).
#' Global measures of data utility for microdata masked for disclosure
#' limitation. \emph{Journal of Privacy and Confidentiality},
#' \strong{1}, 111--124.
#'
#' @seealso
#' \code{\link{simSynthData}},
#' \code{\link{ps_test}}
#'
#' @export
#'
#' @examples
#' data(brittany_soil_ps)
#'
#' set.seed(1)
#' V3 <- simSynthData(brittany_soil_ps, M = 3)
#'
#' utility_measures(brittany_soil_ps, V3, M = 3)
utility_measures <- function(X, V, M = 1L,
                             alpha = 0.05,
                             verbose = TRUE) {

  ## ------------------------------------------------------------------
  ## Input validation
  ## ------------------------------------------------------------------

  X <- .validate_X(X)
  V <- .validate_X(V)
  M <- .validate_M(M)
  alpha <- .validate_alpha_utility(alpha)

  n <- nrow(X)
  p <- ncol(X)
  N <- nrow(V)

  if (N != M * n) {
    stop(
      sprintf(
        "nrow(V) = %d is not equal to M * nrow(X) = %d * %d = %d.",
        N, M, n, M * n
      ),
      call. = FALSE
    )
  }

  if (ncol(V) != p) {
    stop(
      sprintf(
        "ncol(V) = %d must equal ncol(X) = %d.",
        ncol(V), p
      ),
      call. = FALSE
    )
  }

  if (!is.null(colnames(X)) && !is.null(colnames(V))) {
    if (!identical(colnames(X), colnames(V))) {
      stop(
        "Column names of 'X' and 'V' must match and be in the same order.",
        call. = FALSE
      )
    }
  }

  ## ------------------------------------------------------------------
  ## Per-release statistics
  ## ------------------------------------------------------------------

  rel_idx <- lapply(seq_len(M), function(m) {
    (m - 1L) * n + seq_len(n)
  })

  vnames <- if (!is.null(colnames(X))) {
    colnames(X)
  } else {
    paste0("V", seq_len(p))
  }

  S_list <- lapply(rel_idx, function(idx) {
    stats::cov(V[idx, , drop = FALSE])
  })

  mu_list <- lapply(rel_idx, function(idx) {
    colMeans(V[idx, , drop = FALSE])
  })

  S_avg <- Reduce("+", S_list) / M
  mu_avg <- Reduce("+", mu_list) / M

  ## ------------------------------------------------------------------
  ## Original statistics
  ## ------------------------------------------------------------------

  S_orig <- stats::cov(X)
  mu_orig <- colMeans(X)
  sd_orig <- apply(X, 2L, stats::sd)

  if (any(sd_orig <= .Machine$double.eps)) {
    stop(
      "All variables in 'X' must have positive variance.",
      call. = FALSE
    )
  }

  ## ------------------------------------------------------------------
  ## 1. Frobenius distance
  ## ------------------------------------------------------------------

  frob <- as.numeric(norm(S_orig - S_avg, "F"))


  ## ------------------------------------------------------------------
  ## 2. Mean standardized mean difference
  ## ------------------------------------------------------------------

  smd <- abs(mu_avg - mu_orig) / sd_orig
  smd_mean <- mean(smd)


  ## ------------------------------------------------------------------
  ## 3. Variance ratio range
  ## ------------------------------------------------------------------

  var_syn <- diag(S_avg)
  var_orig <- diag(S_orig)

  var_ratio <- var_syn / var_orig

  var_ratio_range <- range(var_ratio, na.rm = TRUE)
  names(var_ratio_range) <- c("min", "max")


  ## ------------------------------------------------------------------
  ## 4. Propensity score MSE ratio
  ## ------------------------------------------------------------------

  c_prop <- M * n / (n + M * n)

  df_pool <- data.frame(
    rbind(X, V),
    check.names = TRUE
  )

  label_name <- ".ps_label"

  while (label_name %in% names(df_pool)) {
    label_name <- paste0(label_name, "_")
  }

  df_pool[[label_name]] <- c(rep(0L, n), rep(1L, M * n))

  form_pmse <- stats::as.formula(paste(label_name, "~ ."))

  fit_pmse <- tryCatch(
    suppressWarnings(
      stats::glm(
        form_pmse,
        data = df_pool,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_pmse) && isTRUE(fit_pmse$converged)) {
    phat <- stats::fitted(fit_pmse)

    if (all(is.finite(phat))) {
      pmse <- mean((phat - c_prop)^2)
      pmse_null <- c_prop * (1 - c_prop) * (p + 1L) / (n + M * n)
      pmse_ratio <- pmse / pmse_null
    } else {
      pmse_ratio <- NA_real_
    }
  } else {
    pmse_ratio <- NA_real_
  }

  ## ------------------------------------------------------------------
  ## 5. Mean confidence interval overlap
  ## ------------------------------------------------------------------

  t_crit_x <- stats::qt(1 - alpha / 2, df = n - 1L)
  t_crit_v <- stats::qt(1 - alpha / 2, df = M * n - 1L)

  mu_syn_pooled <- colMeans(V)
  sd_syn_pooled <- apply(V, 2L, stats::sd)

  ci_overlap <- vapply(seq_len(p), function(j) {
    mx <- mu_orig[j]
    sx <- sd_orig[j]

    mv <- mu_syn_pooled[j]
    sv <- sd_syn_pooled[j]

    lx <- mx - t_crit_x * sx / sqrt(n)
    ux <- mx + t_crit_x * sx / sqrt(n)

    lv <- mv - t_crit_v * sv / sqrt(M * n)
    uv <- mv + t_crit_v * sv / sqrt(M * n)

    overlap <- max(0, min(ux, uv) - max(lx, lv))
    width <- max(ux - lx, uv - lv)

    if (width < .Machine$double.eps) {
      1
    } else {
      overlap / width
    }
  }, numeric(1L))

  ci_overlap_mean <- mean(ci_overlap)

  ## ------------------------------------------------------------------
  ## Return object
  ## ------------------------------------------------------------------

  result <- list(
    frobenius       = frob,
    smd_mean        = smd_mean,
    var_ratio_range = var_ratio_range,
    pmse_ratio      = pmse_ratio,
    ci_overlap_mean = ci_overlap_mean,
    M               = M,
    n               = n,
    p               = p,
    alpha           = alpha
  )

  class(result) <- "ps_utility"

  if (verbose) {
    print(result)
  }

  invisible(result)
}

#' @title Print Method for \code{ps_utility} Objects
#'
#' @description
#' Prints a formatted summary of the five headline utility measures returned
#' by \code{\link{utility_measures}}.
#'
#' @param x An object of class \code{ps_utility}.
#' @param ... Further arguments, currently ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @exportS3Method print ps_utility
print.ps_utility <- function(x, ...) {

  bar <- strrep("-", 68)

  cat("\nPS Synthetic Data Utility Assessment\n")
  cat(bar, "\n")
  cat(sprintf(
    "  n = %d | p = %d | M = %d releases (N = %d) | alpha = %.3f\n",
    x$n, x$p, x$M, x$M * x$n, x$alpha
  ))
  cat(bar, "\n\n")

  summary_df <- data.frame(
    Measure = c(
      "Frobenius distance",
      "Mean SMD",
      "Variance ratio range",
      "pMSE ratio",
      "Mean CI overlap"
    ),
    Value = c(
      sprintf("%.4f", x$frobenius),
      sprintf("%.4f", x$smd_mean),
      sprintf(
        "[%.4f, %.4f]",
        x$var_ratio_range["min"],
        x$var_ratio_range["max"]
      ),
      if (is.na(x$pmse_ratio)) {
        "NA"
      } else {
        sprintf("%.4f", x$pmse_ratio)
      },
      sprintf("%.4f", x$ci_overlap_mean)
    ),
    Flag = c(
      "",
      .ps_utility_flag_smd(x$smd_mean),
      .ps_utility_flag_vr(x$var_ratio_range),
      .ps_utility_flag_pmse(x$pmse_ratio),
      .ps_utility_flag_ci(x$ci_overlap_mean)
    ),
    stringsAsFactors = FALSE
  )

  print(summary_df, row.names = FALSE)

  cat("\n")
  cat("Flag rules:\n")
  cat("  Mean SMD:              * > 0.10, ** > 0.20\n")
  cat("  Variance ratio range:  * outside [0.90, 1.10], ** outside [0.80, 1.20]\n")
  cat("  pMSE ratio:            * > 1.50, ** > 2.00\n")
  cat("  Mean CI overlap:       * < 0.90, ** < 0.70\n")
  cat(bar, "\n")

  invisible(x)
}


