chiprod <- function(dimension, degrees) {
  product <- 1
  for (i in 1:dimension) {
    product <- product * rchisq(1, degrees - i + 1)
  }
  return(product)
}
