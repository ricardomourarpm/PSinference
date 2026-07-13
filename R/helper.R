## ------------------------------------------------------------------
## Internal validators
## ------------------------------------------------------------------

#' Validate a numeric input matrix.
#'
#' @noRd
.validate_X <- function(X, name = "X") {
  if (is.data.frame(X)) {
    non_numeric <- !vapply(X, is.numeric, logical(1L))

    if (any(non_numeric)) {
      stop(
        sprintf(
          "'%s' must contain only numeric columns. Non-numeric columns: %s.",
          name,
          paste(names(X)[non_numeric], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop(
      sprintf("'%s' must be a numeric matrix or data frame.", name),
      call. = FALSE
    )
  }

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(n) || is.null(p)) {
    stop(
      sprintf("'%s' must be a two-dimensional numeric object.", name),
      call. = FALSE
    )
  }

  if (n <= 1L) {
    stop(
      sprintf("'%s' must have at least 2 rows.", name),
      call. = FALSE
    )
  }

  if (p < 1L) {
    stop(
      sprintf("'%s' must have at least 1 column.", name),
      call. = FALSE
    )
  }

  if (n <= p) {
    stop(
      sprintf(
        "Sample size n = %d must exceed the number of variables p = %d.",
        n, p
      ),
      call. = FALSE
    )
  }

  if (any(!is.finite(X))) {
    stop(
      sprintf("'%s' contains non-finite values: NA, NaN, Inf, or -Inf.", name),
      call. = FALSE
    )
  }

  X
}


#' Validate the number of synthetic releases.
#'
#' @noRd
.validate_M <- function(M) {
  if (length(M) != 1L || !is.numeric(M) || !is.finite(M)) {
    stop("'M' must be a single positive integer.", call. = FALSE)
  }

  if (M < 1L || M != floor(M)) {
    stop("'M' must be a single positive integer.", call. = FALSE)
  }

  as.integer(M)
}


#' Validate an integer block partition.
#'
#' @noRd
.validate_part <- function(part, p) {
  if (length(part) != 1L || !is.numeric(part) || !is.finite(part)) {
    stop("'part' must be a single integer.", call. = FALSE)
  }

  if (part != floor(part)) {
    stop("'part' must be a single integer.", call. = FALSE)
  }

  p1 <- as.integer(part)

  if (p1 < 1L || p1 >= p) {
    stop(
      sprintf(
        "'part' must be an integer in {1, ..., p - 1}. Got %s with p = %d.",
        part, p
      ),
      call. = FALSE
    )
  }

  p1
}


#' Check whether a square matrix is positive definite.
#'
#' @noRd
.check_pd <- function(A, name = "matrix") {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop(
      sprintf("The %s must be a square matrix.", name),
      call. = FALSE
    )
  }

  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  tol <- .Machine$double.eps * max(abs(ev)) * nrow(A)

  if (any(ev <= tol)) {
    stop(
      sprintf("The %s is not positive definite.", name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Check the effective sample size.
#'
#' @noRd
.check_N <- function(N, p) {
  if (N <= p + 1L) {
    stop(
      sprintf(
        paste(
          "Effective sample size N = Mn = %d must exceed p + 1 = %d.",
          "Increase n or M."
        ),
        N, p + 1L
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Compute the centered cross-product matrix.
#'
#' @noRd
.compute_S_star <- function(V) {
  V <- as.matrix(V)
  Vc <- sweep(V, 2L, colMeans(V), "-")

  crossprod(Vc)
}


## ------------------------------------------------------------------
## Internal block-resolution helper
## ------------------------------------------------------------------

#' Resolve variable blocks for independence and regression tests.
#'
#' @noRd
.resolve_blocks <- function(V, p, part, group_a, group_b, fun_name) {
  cn <- colnames(V)

  ## Integer interface: first part columns versus remaining columns
  if (!is.null(part) && is.null(group_a) && is.null(group_b)) {
    p1 <- .validate_part(part, p)

    idx1 <- seq_len(p1)
    idx2 <- seq_len(p - p1) + p1

    lbl1 <- if (!is.null(cn)) {
      paste(cn[idx1], collapse = ", ")
    } else {
      paste0("cols 1-", p1)
    }

    lbl2 <- if (!is.null(cn)) {
      paste(cn[idx2], collapse = ", ")
    } else {
      paste0("cols ", p1 + 1L, "-", p)
    }

    return(list(
      V = V,
      p1 = p1,
      idx1 = idx1,
      idx2 = idx2,
      lbl1 = lbl1,
      lbl2 = lbl2
    ))
  }

  ## Named or indexed interface
  if (!is.null(group_a) && !is.null(group_b)) {
    .to_idx <- function(g, which_arg) {
      if (length(g) < 1L) {
        stop(
          sprintf("'%s' must contain at least one variable.", which_arg),
          call. = FALSE
        )
      }

      if (is.character(g)) {
        if (is.null(cn)) {
          stop(
            sprintf(
              "In '%s': column names are required when '%s' uses names.",
              fun_name, which_arg
            ),
            call. = FALSE
          )
        }

        bad <- setdiff(g, cn)

        if (length(bad) > 0L) {
          stop(
            sprintf(
              "'%s' contains names not found in colnames(V): %s.",
              which_arg,
              paste(bad, collapse = ", ")
            ),
            call. = FALSE
          )
        }

        match(g, cn)
      } else if (is.numeric(g)) {
        if (any(!is.finite(g)) || any(g != floor(g))) {
          stop(
            sprintf("'%s' must contain integer column indices.", which_arg),
            call. = FALSE
          )
        }

        g <- as.integer(g)

        if (any(g < 1L | g > p)) {
          stop(
            sprintf(
              "'%s' must contain column indices between 1 and p = %d.",
              which_arg, p
            ),
            call. = FALSE
          )
        }

        g
      } else {
        stop(
          sprintf(
            "'%s' must be either a character vector of column names or an integer vector of column indices.",
            which_arg
          ),
          call. = FALSE
        )
      }
    }

    idx1 <- .to_idx(group_a, "group_a")
    idx2 <- .to_idx(group_b, "group_b")

    if (anyDuplicated(idx1)) {
      stop("'group_a' must not contain duplicate variables.", call. = FALSE)
    }

    if (anyDuplicated(idx2)) {
      stop("'group_b' must not contain duplicate variables.", call. = FALSE)
    }

    if (anyDuplicated(c(idx1, idx2))) {
      stop(
        sprintf(
          "In '%s': 'group_a' and 'group_b' must not overlap.",
          fun_name
        ),
        call. = FALSE
      )
    }

    if (!setequal(c(idx1, idx2), seq_len(p))) {
      stop(
        sprintf(
          "In '%s': 'group_a' and 'group_b' must cover all %d columns.",
          fun_name, p
        ),
        call. = FALSE
      )
    }

    p1 <- length(idx1)

    V <- V[, c(idx1, idx2), drop = FALSE]

    lbl1 <- if (!is.null(cn)) {
      paste(cn[idx1], collapse = ", ")
    } else {
      paste(idx1, collapse = ", ")
    }

    lbl2 <- if (!is.null(cn)) {
      paste(cn[idx2], collapse = ", ")
    } else {
      paste(idx2, collapse = ", ")
    }

    return(list(
      V = V,
      p1 = p1,
      idx1 = seq_len(p1),
      idx2 = seq_len(p - p1) + p1,
      lbl1 = lbl1,
      lbl2 = lbl2
    ))
  }

  stop(
    sprintf(
      "In '%s': supply either 'part' or both 'group_a' and 'group_b'.",
      fun_name
    ),
    call. = FALSE
  )
}

#' Validate the significance level.
#'
#' @noRd
.validate_alpha <- function(alpha) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha)) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }

  alpha
}

#' Validate the Monte Carlo sample size.
#'
#' @noRd
.validate_iterations <- function(iterations) {
  if (
    length(iterations) != 1L ||
      !is.numeric(iterations) ||
      !is.finite(iterations)
  ) {
    stop("'iterations' must be a single positive integer.", call. = FALSE)
  }

  if (iterations < 1L || iterations != floor(iterations)) {
    stop("'iterations' must be a single positive integer.", call. = FALSE)
  }

  as.integer(iterations)
}


#' Validate a user-supplied null distribution.
#'
#' @noRd
.validate_null_dist <- function(null_dist) {
  if (!is.numeric(null_dist) || length(null_dist) < 1L) {
    stop("'null_dist' must be a non-empty numeric vector.", call. = FALSE)
  }

  if (any(!is.finite(null_dist))) {
    stop(
      "'null_dist' contains non-finite values: NA, NaN, Inf, or -Inf.",
      call. = FALSE
    )
  }

  null_dist
}


#' Resolve stacked PS dimensions.
#'
#' @noRd
.resolve_ps_dimensions <- function(V, M) {
  N <- nrow(V)
  p <- ncol(V)

  if (N %% M != 0L) {
    stop(
      sprintf(
        paste(
          "The number of rows in 'V' must be divisible by M.",
          "Got nrow(V) = %d and M = %d."
        ),
        N, M
      ),
      call. = FALSE
    )
  }

  n <- N %/% M

  if (n <= p) {
    stop(
      sprintf(
        paste(
          "The implied original sample size n = nrow(V) / M = %d",
          "must exceed the number of variables p = %d."
        ),
        n, p
      ),
      call. = FALSE
    )
  }

  .check_N(N, p)

  list(
    N = N,
    n = n,
    p = p
  )
}

#' Validate a positive integer scalar.
#'
#' @noRd
.validate_positive_integer <- function(x, name, min_value = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    stop(
      sprintf("'%s' must be a single positive integer.", name),
      call. = FALSE
    )
  }

  if (x != floor(x) || x < min_value) {
    stop(
      sprintf("'%s' must be a single positive integer.", name),
      call. = FALSE
    )
  }

  as.integer(x)
}


#' Validate common arguments for null-distribution simulations.
#'
#' @noRd
.validate_distribution_args <- function(nsample, pvariates, iterations, M) {
  n <- .validate_positive_integer(nsample, "nsample", min_value = 2L)
  p <- .validate_positive_integer(pvariates, "pvariates", min_value = 1L)
  it <- .validate_positive_integer(iterations, "iterations", min_value = 1L)
  M <- .validate_positive_integer(M, "M", min_value = 1L)

  if (n <= p) {
    stop("'nsample' must exceed 'pvariates'.", call. = FALSE)
  }

  list(
    n = n,
    p = p,
    iterations = it,
    M = M,
    N = n * M
  )
}


#' Validate the block size argument.
#'
#' @noRd
.validate_distribution_part <- function(part, p, require_p1_le_p2 = FALSE) {
  p1 <- .validate_positive_integer(part, "part", min_value = 1L)

  if (p1 >= p) {
    stop("'part' must satisfy 1 <= part < pvariates.", call. = FALSE)
  }

  p2 <- p - p1

  if (require_p1_le_p2 && p1 > p2) {
    stop("Regression test requires p1 <= p2.", call. = FALSE)
  }

  list(
    p1 = p1,
    p2 = p2
  )
}


.wilks_f_params <- function(p1, p2, nu) {
  # nu = n - 1 (original-data Wishart df)
  # s  = sqrt((p1^2 * p2^2 - 4) / (p1^2 + p2^2 - 5))  [= p2 if p1=1, p1 if p2=1]
  denom2 <- p1^2 + p2^2 - 5
  s <- if (denom2 > 0) sqrt((p1^2 * p2^2 - 4) / denom2) else 1.0
  df1 <- p1 * p2
  df2 <- s * (nu - (p1 + p2 + 1) / 2) - (p1 * p2 - 2) / 2 - 1
  list(s = s, df1 = df1, df2 = df2)
}

#' Validate a significance level for utility measures.
#'
#' @param alpha Numeric significance level.
#'
#' @return A numeric scalar.
#'
#' @noRd
.validate_alpha_utility <- function(alpha) {
  if (
    !is.numeric(alpha) ||
    length(alpha) != 1L ||
    is.na(alpha) ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    stop(
      "'alpha' must be a single numeric value in the interval (0, 1).",
      call. = FALSE
    )
  }

  alpha
}


#' Flag mean standardized mean difference.
#'
#' @param x Numeric mean SMD.
#'
#' @return A character string.
#'
#' @noRd
.ps_utility_flag_smd <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("")
  }

  if (x > 0.20) {
    return("**")
  }

  if (x > 0.10) {
    return("*")
  }

  ""
}

#' Flag variance ratio range.
#'
#' @param x Numeric vector of length two with minimum and maximum variance
#'   ratios.
#'
#' @return A character string.
#'
#' @noRd
.ps_utility_flag_vr <- function(x) {
  if (length(x) != 2L || any(is.na(x)) || any(!is.finite(x))) {
    return("")
  }

  if (x["min"] < 0.80 || x["max"] > 1.20) {
    return("**")
  }

  if (x["min"] < 0.90 || x["max"] > 1.10) {
    return("*")
  }

  ""
}


#' Flag pMSE ratio.
#'
#' @param x Numeric pMSE ratio.
#'
#' @return A character string.
#'
#' @noRd
.ps_utility_flag_pmse <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("")
  }

  if (x > 2.00) {
    return("**")
  }

  if (x > 1.50) {
    return("*")
  }

  ""
}


#' Flag mean confidence interval overlap.
#'
#' @param x Numeric mean confidence interval overlap.
#'
#' @return A character string.
#'
#' @noRd
.ps_utility_flag_ci <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("")
  }

  if (x < 0.70) {
    return("**")
  }

  if (x < 0.90) {
    return("*")
  }

  ""
}


#' correct p_values
#'
#' @param p p-value
#'
#' @return numeric
#'
#' @noRd
.fmt_p <- function(p, digits = 4) {
  if (is.na(p)) return("NA")
  thresh <- 10^(-digits)
  if (p < thresh) {
    sprintf("< %s", format(thresh, scientific = FALSE))
  } else {
    sprintf(paste0("%.", digits, "f"), p)
  }
}
