partition <- function(mat, nrows, ncols) {
  # split matrix into submatrices
  r <- dim(mat)[1]
  c <- dim(mat)[2]
  if (r < nrows & c < ncols) {
    stop("ncols or nrows is larger than array .")
  }
  A <- mat[1:nrows, 1:ncols]
  B <- mat[1:nrows, (ncols + 1):c]
  C <- mat[(nrows + 1):r, 1:ncols]
  D <- mat[(nrows + 1):r, (ncols + 1):c]
  return(list(A, B, C, D))
}


chiprod <- function(dimension, degrees) {
  product <- prod(rchisq(dimension, (degrees - dimension + 1):degrees))
  return(product)
}
