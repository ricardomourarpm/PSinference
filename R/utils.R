partition <- function(array, nrows, ncols) {
  #split matrix into submatrices
  r <- dim(array)[1]
  c <- dim(array)[2]
  if(r>nrows & c>ncols){
    stop("ncols or nrows is larger than array .")
  }
  A <- array[1:nrows, 1:ncols]
  B <- array[1:nrows, (ncols + 1):c]
  C <- array[(nrows + 1):r, 1:ncols]
  D <- array[(nrows + 1):r, (ncols + 1):c]
  return(list(A, B, C, D))
}

