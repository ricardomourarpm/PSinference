#' @title Unified Wrapper for PS Inference Tests
#'
#' @description
#' Dispatches to the appropriate exact inferential procedure based on the
#' \code{test} argument and can optionally produce the diagnostic plot
#' immediately. This is the main entry point for users who prefer a single
#' function instead of calling the four individual test functions directly.
#'
#' @param V Stacked synthetic data set, given as an \eqn{Mn \times p}
#'   numeric matrix, as returned by \code{\link{simSynthData}}.
#' @param M Positive integer giving the number of synthetic releases. The
#'   default is \code{1L}. Setting \code{M = 1} recovers the
#'   single-release procedures of Klein et al. (2021).
#' @param test Character string specifying the test. One of
#'   \code{"gv"} for generalized variance, \code{"sphericity"},
#'   \code{"independence"}, or \code{"regression"}.
#' @param plot Logical. If \code{TRUE}, \code{plot()} is called on the
#'   result before it is returned, so the null-distribution diagnostic is
#'   displayed automatically. The default is \code{FALSE}.
#' @param ... Additional arguments passed to the corresponding test
#'   function or, when \code{plot = TRUE}, to \code{plot.ps_test()}.
#'   Arguments intended for the test function, such as \code{part},
#'   \code{Sigma}, \code{Delta0}, \code{alpha}, or \code{iterations}, and
#'   graphical arguments, such as \code{main}, \code{shade_col},
#'   \code{dist_col}, \code{stat_col}, and \code{crit_col}, are separated
#'   automatically.
#'
#' @return
#' An object of class \code{\link{ps_test-class}}, invisibly when
#' \code{plot = TRUE}.
#'
#' @seealso
#' \code{\link{gv_test}},
#' \code{\link{sphericity_test}},
#' \code{\link{independence_test}},
#' \code{\link{regression_test}}
#'
#' @export
#'
#' @examples
#' data(attitude)
#'
#' set.seed(1)
#' V <- simSynthData(attitude, M = 3)
#'
#' \donttest{
#' # Run and print only
#' ps_test(V, M = 3, test = "sphericity", iterations = 1000L)
#'
#' # Run and plot in one call
#' ps_test(V,
#'   M = 3, test = "sphericity",
#'   iterations = 1000L, plot = TRUE
#' )
#'
#' # Independence with named blocks
#' ps_test(
#'   V,
#'   M = 3,
#'   test = "independence",
#'   group_a = c("rating", "complaints"),
#'   group_b = c("privileges", "learning"),
#'   iterations = 1000L,
#'   plot = TRUE
#' )
#'
#' # Generalized variance with a reference covariance matrix
#' ps_test(
#'   V,
#'   M = 3,
#'   test = "gv",
#'   Sigma = cov(attitude),
#'   iterations = 1000L,
#'   plot = TRUE
#' )
#' }
ps_test <- function(V,
                    M = 1L,
                    test = c(
                      "gv", "sphericity",
                      "independence", "regression"
                    ),
                    plot = FALSE,
                    ...) {
  test <- match.arg(test)

  ## Separate test arguments from plot arguments
  plot_arg_names <- c(
    "main",
    "shade_col",
    "dist_col",
    "stat_col",
    "crit_col"
  )

  dots <- list(...)

  plot_args <- dots[intersect(names(dots), plot_arg_names)]
  test_args <- dots[setdiff(names(dots), plot_arg_names)]

  ## Run the selected test
  res <- switch(test,
    gv = do.call(
      gv_test,
      c(list(V = V, M = M), test_args)
    ),
    sphericity = do.call(
      sphericity_test,
      c(list(V = V, M = M), test_args)
    ),
    independence = do.call(
      independence_test,
      c(list(V = V, M = M), test_args)
    ),
    regression = do.call(
      regression_test,
      c(list(V = V, M = M), test_args)
    )
  )

  ## Optional plot
  if (plot) {
    do.call(graphics::plot, c(list(x = res), plot_args))
    return(invisible(res))
  }

  res
}

#' @title Partition a Matrix into Four Blocks
#'
#' @description
#' Splits a numeric matrix \eqn{\mathbf{M}} into four submatrices according
#' to a two-group partition of its rows and columns:
#' \deqn{
#'   \mathbf{M}
#'   =
#'   \begin{bmatrix}
#'   \mathbf{A} & \mathbf{B} \\
#'   \mathbf{C} & \mathbf{D}
#'   \end{bmatrix}.
#' }
#'
#' The blocks are:
#' \describe{
#'   \item{\eqn{\mathbf{A}}}{Rows in the first group and columns in the
#'     first group.}
#'   \item{\eqn{\mathbf{B}}}{Rows in the first group and columns in the
#'     second group.}
#'   \item{\eqn{\mathbf{C}}}{Rows in the second group and columns in the
#'     first group.}
#'   \item{\eqn{\mathbf{D}}}{Rows in the second group and columns in the
#'     second group.}
#' }
#'
#' Two interfaces are available:
#' \describe{
#'   \item{Integer interface}{
#'     Supply \code{nrows} and \code{ncols}. The first \code{nrows} rows
#'     and the first \code{ncols} columns form the first block.
#'   }
#'   \item{Name interface}{
#'     Supply \code{part1} as a character vector of row and column names.
#'     This interface is intended for square matrices, such as covariance
#'     or correlation matrices. The matrix is reordered so that the
#'     \code{part1} rows and columns appear first. The optional
#'     \code{part2} argument is checked against the complement of
#'     \code{part1}.
#'   }
#' }
#'
#' @param Matrix A numeric matrix. When using the name interface, it must
#'   be a square matrix with matching row and column names.
#' @param nrows Integer giving the number of rows in the first row block.
#'   Ignored when \code{part1} is supplied.
#' @param ncols Integer giving the number of columns in the first column
#'   block. Ignored when \code{part1} is supplied.
#' @param part1 Character vector of names forming the first group. Used
#'   for both row and column reordering in the name interface.
#' @param part2 Optional character vector naming the second group. If
#'   supplied, it is checked against the complement of \code{part1}.
#'
#' @return
#' A named list of class \code{"ps_partition"} with elements \code{A},
#' \code{B}, \code{C}, and \code{D}. Numeric indexing with
#' \code{[[1]]} through \code{[[4]]} also works.
#'
#' @export
#'
#' @examples
#' M <- matrix(
#'   1:16,
#'   4,
#'   4,
#'   dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
#' )
#'
#' # Integer interface
#' b <- partition(M, nrows = 2, ncols = 2)
#' b$A
#'
#' # Name interface
#' b2 <- partition(M, part1 = c("A", "B"))
#' b2$A
#' b2$D
#'
#' # Covariance matrix example
#' data(attitude)
#' b3 <- partition(
#'   cov(attitude),
#'   part1 = c("rating", "complaints", "raises"),
#'   part2 = c("privileges", "learning", "advance", "critical")
#' )
#' b3$A
#' b3$D
partition <- function(Matrix,
                      nrows = NULL,
                      ncols = NULL,
                      part1 = NULL,
                      part2 = NULL) {
  ## Input checks
  if (!is.matrix(Matrix) || !is.numeric(Matrix)) {
    stop("'Matrix' must be a numeric matrix.", call. = FALSE)
  }

  if (any(!is.finite(Matrix))) {
    stop(
      "'Matrix' contains non-finite values: NA, NaN, Inf, or -Inf.",
      call. = FALSE
    )
  }

  r <- nrow(Matrix)
  cc <- ncol(Matrix)
  rn <- rownames(Matrix)
  cn <- colnames(Matrix)

  ## Name interface
  if (!is.null(part1)) {
    if (!is.character(part1) || length(part1) < 1L) {
      stop(
        "'part1' must be a non-empty character vector of names.",
        call. = FALSE
      )
    }

    if (anyDuplicated(part1)) {
      stop("'part1' must not contain duplicate names.", call. = FALSE)
    }

    if (!is.null(part2) && anyDuplicated(part2)) {
      stop("'part2' must not contain duplicate names.", call. = FALSE)
    }

    if (r != cc) {
      stop(
        "The name interface requires a square matrix.",
        call. = FALSE
      )
    }

    if (is.null(rn) || is.null(cn)) {
      stop(
        "'Matrix' must have row and column names when 'part1' is supplied.",
        call. = FALSE
      )
    }

    if (!identical(rn, cn)) {
      stop(
        "The name interface requires identical row and column names.",
        call. = FALSE
      )
    }

    bad <- setdiff(part1, rn)

    if (length(bad) > 0L) {
      stop(
        "'part1' contains names not found in 'Matrix': ",
        paste(bad, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    part2_expected <- setdiff(rn, part1)

    if (length(part2_expected) == 0L) {
      stop(
        "'part1' covers all rows and columns; the second block would be empty.",
        call. = FALSE
      )
    }

    if (!is.null(part2)) {
      if (!is.character(part2)) {
        stop("'part2' must be a character vector of names.", call. = FALSE)
      }

      if (!setequal(part2, part2_expected)) {
        stop(
          "'part2' must be the complement of 'part1'. Got: ",
          paste(sort(part2), collapse = ", "),
          ". Expected: ",
          paste(sort(part2_expected), collapse = ", "),
          ".",
          call. = FALSE
        )
      }

      part2_expected <- part2
    }

    Matrix <- Matrix[
      c(part1, part2_expected),
      c(part1, part2_expected),
      drop = FALSE
    ]

    nrows <- length(part1)
    ncols <- length(part1)
    r <- nrow(Matrix)
    cc <- ncol(Matrix)
  } else {
    ## Integer interface
    if (is.null(nrows) || is.null(ncols)) {
      stop(
        "Supply either 'part1' or both 'nrows' and 'ncols'.",
        call. = FALSE
      )
    }

    if (
      !is.numeric(nrows) ||
        length(nrows) != 1L ||
        is.na(nrows) ||
        nrows < 1L ||
        nrows != floor(nrows) ||
        nrows >= r
    ) {
      stop(
        sprintf(
          "'nrows' must be a single integer in {1, ..., %d}. Got: %s.",
          r - 1L,
          as.character(nrows)
        ),
        call. = FALSE
      )
    }

    if (
      !is.numeric(ncols) ||
        length(ncols) != 1L ||
        is.na(ncols) ||
        ncols < 1L ||
        ncols != floor(ncols) ||
        ncols >= cc
    ) {
      stop(
        sprintf(
          "'ncols' must be a single integer in {1, ..., %d}. Got: %s.",
          cc - 1L,
          as.character(ncols)
        ),
        call. = FALSE
      )
    }

    nrows <- as.integer(nrows)
    ncols <- as.integer(ncols)
  }

  ## Split
  idx_r1 <- seq_len(nrows)
  idx_r2 <- seq.int(nrows + 1L, r)

  idx_c1 <- seq_len(ncols)
  idx_c2 <- seq.int(ncols + 1L, cc)

  A <- Matrix[idx_r1, idx_c1, drop = FALSE]
  B <- Matrix[idx_r1, idx_c2, drop = FALSE]
  C <- Matrix[idx_r2, idx_c1, drop = FALSE]
  D <- Matrix[idx_r2, idx_c2, drop = FALSE]

  structure(
    list(A = A, B = B, C = C, D = D),
    class = "ps_partition"
  )
}

#' Extract a block from a \code{ps_partition} object.
#'
#' @noRd
#' @export
`[[.ps_partition` <- function(x, i) {
  unclass(x)[[i]]
}

#' Print a \code{ps_partition} object.
#'
#' @noRd
#' @exportS3Method print ps_partition
print.ps_partition <- function(x, ...) {
  cat(
    "ps_partition:",
    " A (", nrow(x$A), "x", ncol(x$A), ")",
    " B (", nrow(x$B), "x", ncol(x$B), ")",
    " C (", nrow(x$C), "x", ncol(x$C), ")",
    " D (", nrow(x$D), "x", ncol(x$D), ")\n",
    sep = ""
  )

  invisible(x)
}
