#' @title S3 Class for PS Inference Test Results
#'
#' @description
#' The \code{ps_test} class is the unified output object returned by the
#' inferential functions in \pkg{PSinference}. It stores the test result,
#' the simulated null distribution, and relevant metadata, and provides
#' \code{print}, \code{summary}, and \code{plot} methods for convenient
#' inspection and reporting.
#'
#' @section Slots:
#' \describe{
#'   \item{statistic}{Numeric. Observed value of the test statistic.}
#'   \item{p.value}{Numeric. Monte Carlo p-value.}
#'   \item{alpha}{Numeric. Significance level used.}
#'   \item{decision}{Character. \code{"Reject H0"} or
#'     \code{"Fail to reject H0"}.}
#'   \item{null.dist}{Numeric vector. Simulated null distribution.}
#'   \item{test}{Character. One of \code{"gv"}, \code{"sphericity"},
#'     \code{"independence"}, \code{"regression"}.}
#'   \item{n}{Integer. Original sample size.}
#'   \item{M}{Integer. Number of synthetic releases.}
#'   \item{N}{Integer. Effective sample size \eqn{N = Mn}.}
#'   \item{p}{Integer. Number of variables.}
#'   \item{conf.int}{Numeric vector of length 2 or \code{NULL}.
#'     Confidence interval (generalized variance only).}
#'   \item{sigma2.hat}{Numeric or \code{NULL}. Plug-in estimator of
#'     \eqn{\sigma^2} (sphericity only).}
#'   \item{Delta.hat}{Matrix or \code{NULL}. Plug-in estimator of
#'     \eqn{\Delta}, used for the regression test.}
#'   \item{lbl1}{Character or \code{NULL}. Label for the first variable
#'     block, used by block-based tests.}
#'   \item{lbl2}{Character or \code{NULL}. Label for the second variable
#'     block, used by block-based tests.}
#'   \item{iterations}{Integer. Number of Monte Carlo iterations used
#'     to calibrate the null distribution.}
#' }
#'
#' @name ps_test-class
NULL

#' Construct a \code{ps_test} object.
#'
#' @noRd
new_ps_test <- function(statistic,
                        p.value,
                        alpha,
                        decision,
                        null.dist,
                        test,
                        n,
                        M,
                        p,
                        conf.int = NULL,
                        sigma2.hat = NULL,
                        Delta.hat = NULL,
                        lbl1 = NULL,
                        lbl2 = NULL) {
  stopifnot(
    length(statistic) == 1L,
    is.numeric(statistic) || is.na(statistic),
    length(p.value) == 1L,
    is.numeric(p.value) || is.na(p.value),
    length(alpha) == 1L,
    is.numeric(alpha),
    length(decision) == 1L,
    is.character(decision),
    is.numeric(null.dist),
    length(test) == 1L,
    is.character(test),
    length(n) == 1L,
    is.numeric(n),
    length(M) == 1L,
    is.numeric(M),
    length(p) == 1L,
    is.numeric(p)
  )

  structure(
    list(
      statistic  = statistic,
      p.value    = p.value,
      alpha      = alpha,
      decision   = decision,
      null.dist  = null.dist,
      test       = test,
      n          = as.integer(n),
      M          = as.integer(M),
      N          = as.integer(M * n),
      p          = as.integer(p),
      conf.int   = conf.int,
      sigma2.hat = sigma2.hat,
      Delta.hat  = Delta.hat,
      lbl1       = lbl1,
      lbl2       = lbl2,
      iterations = length(null.dist)
    ),
    class = "ps_test"
  )
}

#' @title Print a \code{ps_test} Object
#'
#' @description
#' Prints a concise, human-readable summary of the test result stored in
#' a \code{ps_test} object.
#'
#' @param x An object of class \code{ps_test}.
#' @param ... Further arguments, currently ignored.
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @exportS3Method print ps_test
#'
#' @examples
#' data(ps_attitude)
#'
#' set.seed(1)
#' V <- simSynthData(ps_attitude, M = 3)
#'
#' \donttest{
#' res <- sphericity_test(V, M = 3, iterations = 1000L)
#' print(res)
#' }
print.ps_test <- function(x, ...) {
  test_label <- .test_label(x$test)
  bar <- strrep("-", 56)

  cat("\n")
  cat("PSinference:", test_label, "\n")
  cat(bar, "\n")
  cat(sprintf("  Original sample size  n = %d\n", x$n))
  cat(sprintf("  Number of variables   p = %d\n", x$p))
  cat(sprintf("  Number of releases    M = %d\n", x$M))
  cat(sprintf("  Effective sample size N = Mn = %d\n", x$N))
  cat(bar, "\n")
  cat(sprintf("  Test statistic  : %.6g\n", x$statistic))
  cat(sprintf(
    "  p-value         : %s\n",
    .ps_fmt_pvalue(x$p.value, x$iterations)
  ))
  cat(sprintf("  alpha           : %.2f\n", x$alpha))
  cat(sprintf("  Decision        : %s\n", x$decision))

  ## Test-specific output
  if (x$test == "gv") {
    cat("  H0 : |Sigma| = |Sigma_0|  (two-sided test)\n")

    if (!is.null(x$conf.int)) {
      cat(sprintf(
        "  %.0f%% confidence interval for |Sigma|: (%.4e, %.4e)\n",
        100 * (1 - x$alpha),
        x$conf.int[1],
        x$conf.int[2]
      ))
    }
  }

  if (x$test == "sphericity" && !is.null(x$sigma2.hat)) {
    cat("  H0 : Sigma = sigma^2 * I_p\n")
    cat(sprintf("  sigma2_hat = %.4f  (under H0)\n", x$sigma2.hat))
  }

  if (x$test == "independence") {
    if (!is.null(x$lbl1) && !is.null(x$lbl2)) {
      cat(sprintf(
        "  H0 : {%s} independent of {%s}\n",
        x$lbl1,
        x$lbl2
      ))
    } else {
      cat("  H0 : Sigma_12 = 0\n")
    }
  }

  if (x$test == "regression") {
    if (!is.null(x$lbl1) && !is.null(x$lbl2)) {
      cat(sprintf(
        "  H0 : Delta = Delta_0  (regress {%s} on {%s})\n",
        x$lbl1,
        x$lbl2
      ))
    } else {
      cat("  H0 : Delta = Delta_0\n")
    }

    if (!is.null(x$Delta.hat)) {
      cat("  Delta_hat (plug-in slope):\n")
      print(round(x$Delta.hat, 4))
    }
  }

  cat(bar, "\n")
  cat(sprintf("  Monte Carlo iterations: %d\n", x$iterations))
  cat("\n")

  invisible(x)
}

#' @title Summarize a \code{ps_test} Object
#'
#' @description
#' Prints a detailed summary including the null-distribution quantiles
#' and a comparison with the observed statistic.
#'
#' @param object An object of class \code{ps_test}.
#' @param ... Further arguments, currently ignored.
#'
#' @return Invisibly returns \code{object}.
#'
#' @exportS3Method summary ps_test
#'
#' @examples
#' data(attitude)
#' V <- simSynthData(attitude, M = 3)
#' res <- sphericity_test(V, M = 3)
#' summary(res)
summary.ps_test <- function(object, ...) {
  print(object)

  nd <- object$null.dist

  cat("Null distribution summary (Monte Carlo):\n")

  q <- stats::quantile(
    nd,
    probs = c(
      0.01, 0.025, 0.05,
      0.25, 0.50, 0.75,
      0.95, 0.975, 0.99
    )
  )

  print(round(q, 6))

  cat(sprintf("\nObserved statistic : %.6g\n", object$statistic))

  pct <- round(100 * mean(nd <= object$statistic), 1)
  cat(sprintf("Percentile in null : %.1f%%\n\n", pct))

  invisible(object)
}


#' @title Plot a \code{ps_test} Object
#'
#' @description
#' Produces a density plot of the simulated null distribution with the
#' observed test statistic and critical value(s) marked. The rejection
#' region is shaded.
#'
#' The x-axis always includes both the null distribution and the observed
#' statistic. For the generalized variance and regression tests, a
#' log10 scale is used automatically because these statistics may span
#' several orders of magnitude. Key information is placed below the plot
#' as text so it does not overlap the density curve. The function is
#' multi-panel aware: inside \code{par(mfrow = ...)}, it uses compact
#' in-plot annotations and does not modify the outer margins.
#'
#' @param x An object of class \code{ps_test}.
#' @param main Optional title string. If \code{NULL}, a title is generated
#'   automatically.
#' @param shade_col Color for the rejection-region shading.
#' @param dist_col Color for the null-distribution density fill.
#' @param stat_col Color for the observed-statistic line.
#' @param crit_col Color for the critical-value line(s).
#' @param ... Further arguments passed to \code{plot()}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @exportS3Method plot ps_test
#'
#' @examples
#' data(ps_attitude)
#' V <- simSynthData(ps_attitude, M = 3)
#' plot(sphericity_test(V, M = 3))
plot.ps_test <- function(x,
                         main = NULL,
                         shade_col = grDevices::adjustcolor("tomato", 0.45),
                         dist_col = grDevices::adjustcolor("steelblue", 0.22),
                         stat_col = "firebrick",
                         crit_col = "steelblue4",
                         ...) {
  nd <- x$null.dist
  obs <- x$statistic

  ## Log10 scale for generalized variance and regression tests
  use_log <- x$test %in% c("gv", "regression")

  if (use_log) {
    nd_pos <- nd[is.finite(nd) & nd > 0]

    if (length(nd_pos) < 2L) {
      stop(
        "The null distribution must contain at least two positive finite values for a log-scale plot.",
        call. = FALSE
      )
    }

    nd_plot <- log10(nd_pos)
    obs_plot <- if (!is.na(obs) && is.finite(obs) && obs > 0) {
      log10(obs)
    } else {
      NA_real_
    }

    xlab_str <- expression(log[10](italic(T)^"*"))
  } else {
    nd_plot <- nd[is.finite(nd)]
    obs_plot <- obs
    xlab_str <- expression(italic(T)^"*")
  }

  ## Tail type
  two_sided <- x$test == "gv"
  left_tail <- x$test %in% c("sphericity", "independence")

  ## Critical values on the plot scale
  .cv <- function(prob) {
    v <- as.numeric(stats::quantile(nd, prob, na.rm = TRUE))

    if (use_log) {
      log10(max(v, .Machine$double.eps))
    } else {
      v
    }
  }

  .cv_raw <- function(prob) {
    as.numeric(stats::quantile(nd, prob, na.rm = TRUE))
  }

  if (two_sided) {
    cl_plot <- .cv(x$alpha / 2)
    ch_plot <- .cv(1 - x$alpha / 2)

    cl_raw <- .cv_raw(x$alpha / 2)
    ch_raw <- .cv_raw(1 - x$alpha / 2)
  } else if (left_tail) {
    cl_plot <- .cv(x$alpha)
    ch_plot <- NA_real_

    cl_raw <- .cv_raw(x$alpha)
    ch_raw <- NA_real_
  } else {
    cl_plot <- NA_real_
    ch_plot <- .cv(1 - x$alpha)

    cl_raw <- NA_real_
    ch_raw <- .cv_raw(1 - x$alpha)
  }

  ## Density on the plot scale
  d <- stats::density(nd_plot, bw = "SJ")

  ## x-axis: include null distribution, observed statistic, and critical values
  finite_x <- c(d$x, obs_plot, cl_plot, ch_plot)
  finite_x <- finite_x[is.finite(finite_x)]

  x_range <- range(finite_x)
  pad <- diff(x_range) * 0.08

  if (!is.finite(pad) || pad == 0) {
    pad <- max(abs(x_range), 1) * 0.08
  }

  x_lo <- x_range[1L] - pad
  x_hi <- x_range[2L] + pad

  ylim <- c(0, max(d$y) * 1.08)

  ## Title and annotation strings
  main_ <- if (!is.null(main)) {
    main
  } else {
    sprintf("%s  (M = %d, N = %d)", .test_label(x$test), x$M, x$N)
  }

  .fmt <- function(v_plot) {
    if (is.na(v_plot) || !is.finite(v_plot)) {
      return("NA")
    }

    if (use_log) {
      .ps_fmt_sci(10^v_plot)
    } else {
      sprintf("%.4g", v_plot)
    }
  }

  ann_obs <- sprintf("Observed T* = %s", .fmt(obs_plot))
  ann_pval <- sprintf(
    "p-value = %s",
    .ps_fmt_pvalue(x$p.value, x$iterations)
  )
  ann_dec <- x$decision

  if (two_sided) {
    ann_crit <- sprintf(
      "Reject if T* < %s or T* > %s",
      .fmt(cl_plot),
      .fmt(ch_plot)
    )
  } else if (left_tail) {
    ann_crit <- sprintf("Reject if T* < %s", .fmt(cl_plot))
  } else {
    ann_crit <- sprintf("Reject if T* > %s", .fmt(ch_plot))
  }

  ## Multi-panel awareness
  multi_panel <- prod(graphics::par("mfrow")) > 1L

  if (multi_panel) {
    op <- graphics::par(
      mar = c(5.4, 4.0, 2.5, 0.8),
      mgp = c(2.4, 0.6, 0)
    )

    cex_main <- 0.85
    cex_lab <- 0.80
    cex_axis <- 0.70
    cex_ann <- 0.50

  } else {
    op <- graphics::par(
      oma = c(6.0, 0, 0, 0),
      mar = c(4.0, 4.5, 2.8, 1.2),
      mgp = c(2.8, 0.7, 0)
    )

    cex_main <- 0.90
    cex_lab <- 0.85
    cex_axis <- 0.75
    cex_ann <- 0.75
  }

  on.exit(graphics::par(op), add = TRUE)

  ## Base plot
  plot(
    d,
    main = main_,
    xlab = xlab_str,
    ylab = "Density",
    xlim = c(x_lo, x_hi),
    ylim = ylim,
    lwd = 2.5,
    col = "steelblue",
    axes = FALSE,
    zero.line = FALSE,
    cex.main = cex_main,
    cex.lab = cex_lab,
    ...
  )

  graphics::axis(1L, las = 1L, cex.axis = cex_axis)
  graphics::axis(2L, las = 1L, cex.axis = cex_axis)
  graphics::box()

  ## Null distribution fill
  graphics::polygon(
    c(d$x, rev(d$x)),
    c(d$y, rep(0, length(d$y))),
    col = dist_col,
    border = NA
  )

  graphics::lines(d, col = "steelblue", lwd = 2.5)


  ## Shade rejection region(s)
  .shade <- function(thresh, left) {
    if (is.na(thresh) || !is.finite(thresh)) {
      return(invisible(NULL))
    }

    idx <- if (left) {
      d$x <= thresh
    } else {
      d$x >= thresh
    }

    if (sum(idx) < 2L) {
      return(invisible(NULL))
    }

    graphics::polygon(
      c(d$x[idx], rev(d$x[idx])),
      c(d$y[idx], rep(0, sum(idx))),
      col = shade_col,
      border = NA
    )

    invisible(NULL)
  }

  if (two_sided) {
    .shade(cl_plot, TRUE)
    .shade(ch_plot, FALSE)
  } else if (left_tail) {
    .shade(cl_plot, TRUE)
  } else {
    .shade(ch_plot, FALSE)
  }

  graphics::lines(d, col = "steelblue", lwd = 2.5)

  ## Draw vertical lines inside the plot frame only
  .vline_inside <- function(v, col, lwd = 2, lty = 1L) {
    if (is.na(v) || !is.finite(v)) {
      return(invisible(NULL))
    }

    graphics::segments(
      x0 = v,
      y0 = 0,
      x1 = v,
      y1 = ylim[2L],
      col = col,
      lwd = lwd,
      lty = lty,
      xpd = FALSE
    )

    invisible(NULL)
  }

  .vline_inside(obs_plot, stat_col, lwd = 2.2, lty = 2L)
  .vline_inside(cl_plot, crit_col, lwd = 1.6, lty = 3L)
  .vline_inside(ch_plot, crit_col, lwd = 1.6, lty = 3L)

  ## Annotation below the plot
  dec_col <- if (grepl("^Reject", x$decision)) {
    "firebrick"
  } else {
    "darkgreen"
  }

  if (multi_panel) {
    sub1 <- sprintf("%s   |   %s", ann_obs, ann_pval)
    sub2 <- ann_dec
    sub3 <- ann_crit

    graphics::mtext(
      sub1,
      side = 1L,
      line = 3.0,
      cex = cex_ann,
      col = "black",
      font = 1L
    )

    graphics::mtext(
      sub2,
      side = 1L,
      line = 3.8,
      cex = cex_ann,
      col = dec_col,
      font = 1L
    )
  } else {
    ann <- list(
      list(ann_obs, stat_col, 2L),
      list(ann_pval, "black", 1L),
      list(ann_dec, dec_col, 2L),
      list(ann_crit, crit_col, 1L)
    )

    old_xpd <- graphics::par("xpd")
    graphics::par(xpd = NA)
    on.exit(graphics::par(xpd = old_xpd), add = TRUE)

    for (k in seq_along(ann)) {
      graphics::mtext(
        ann[[k]][[1L]],
        side = 1L,
        outer = TRUE,
        line = (k - 1L) * 1.15 + 0.8,
        cex = cex_ann,
        adj = 0.5,
        font = ann[[k]][[3L]],
        col = ann[[k]][[2L]]
      )
    }
  }

  invisible(x)
}

#' @title Test Whether an Object Has Class \code{ps_test}
#'
#' @description
#' Checks whether an object inherits from class \code{ps_test}.
#'
#' @param x Any R object.
#'
#' @return
#' A logical value: \code{TRUE} if \code{x} inherits from class
#' \code{ps_test}, and \code{FALSE} otherwise.
#'
#' @export
#'
#' @examples
#' data(ps_attitude)
#'
#' set.seed(1)
#' V <- simSynthData(ps_attitude, M = 3)
#'
#' \donttest{
#' res <- sphericity_test(V, M = 3, iterations = 1000L)
#' is.ps_test(res)
#' }
is.ps_test <- function(x) {
  inherits(x, "ps_test")
}


#' Convert an internal test code to a display label.
#'
#' @param test Character string identifying the test.
#'
#' @return A character label.
#'
#' @noRd
.test_label <- function(test) {
  switch(test,
         gv           = "Generalized Variance Test",
         sphericity   = "Sphericity Test",
         independence = "Independence Test",
         regression   = "Regression Test",
         test
  )
}

#' Format a number using compact or scientific notation.
#'
#' @param v Numeric value to format.
#'
#' @return A character string.
#'
#' @noRd
.ps_fmt_sci <- function(v) {
  if (is.na(v) || !is.finite(v)) {
    return("NA")
  }

  if (abs(v) >= 0.001 && abs(v) < 1e5) {
    return(sprintf("%.4g", v))
  }

  sprintf("%.3e", v)
}

#' Format a Monte Carlo p-value for display.
#'
#' This helper never returns \code{"0.0000"}. If the Monte Carlo count is
#' zero, it reports a conservative scientific-notation bound based on the
#' simulation resolution. For example, with \code{B = 2000}, it reports
#' \code{"< 1e-3"}.
#'
#' @param pval Numeric \eqn{p}-value in \eqn{[0, 1]}.
#' @param B Number of Monte Carlo iterations used.
#'
#' @return A character string.
#'
#' @noRd
.ps_fmt_pvalue <- function(pval, B = NULL) {
  if (is.na(pval) || !is.finite(pval)) {
    return("NA")
  }

  if (is.null(B) || is.na(B) || !is.finite(B) || B < 1L) {
    B <- 10000L
  }

  B <- as.integer(B)

  min_detectable <- 1 / B

  if (pval <= 0 || pval < min_detectable) {
    exponent <- ceiling(log10(min_detectable))
    threshold <- 10^exponent
    return(sprintf("< 1e%d", exponent))
  }

  if (pval >= 0.001) {
    return(sprintf("%.4f", pval))
  }

  sprintf("%.2e", pval)
}
