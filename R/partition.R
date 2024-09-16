#' Split a matrix into blocks
#'
#' This function split a matrix into a list of blocks (either by rows and columns).
#'
#' @param Matrix a matrix to split .
#' @param nrows positive integer indicating the number of rows blocks.
#' @param ncols positive integer indicating the number of columns blocks.
#'
#' @return a list of partitioned submatrices
#' @examples
#' Mat = matrix(c(1,0.5,0,0,
#'                   0.5,2,0,0,
#'                   0,0,3,0.2,
#'                   0, 0, 0.2,4), nrow = 4, ncol = 4, byrow = TRUE)
#' partition(Matrix = Mat, nrows = 2, ncols = 2)
#' @export

partition <- function(Matrix, nrows, ncols) {
  # split the matrix into submatrices
  r <- dim(Matrix)[1]
  c <- dim(Matrix)[2]
  if (r < nrows & c < ncols) {
    stop("ncols or nrows is larger than array .")
  }
  A <- Matrix[1:nrows, 1:ncols]
  B <- Matrix[1:nrows, (ncols + 1):c]
  C <- Matrix[(nrows + 1):r, 1:ncols]
  D <- Matrix[(nrows + 1):r, (ncols + 1):c]
  return(list(A, B, C, D))
}
