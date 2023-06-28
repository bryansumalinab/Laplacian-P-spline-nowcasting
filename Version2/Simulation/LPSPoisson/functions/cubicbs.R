# Function to construct cubic B-spline matrix
cubicbs <- function (x, lower, upper, K) 
{
  if (!is.vector(x, mode = "numeric")) 
    stop("x must be a numeric vector")
  if (anyNA(x)) 
    stop("x cannot contain NA or NaN values")
  if (!is.vector(lower, mode = "numeric") || !is.vector(upper, 
                                                        mode = "numeric")) 
    stop("Lower bound and/or upper bound is not a numeric vector")
  if (length(lower) > 1 || length(upper) > 1) 
    stop("Lower bound and/or upper bound must be of length 1")
  if (is.infinite(lower) || is.infinite(upper)) 
    stop("Lower bound and/or upper bound must be finite")
  if (lower >= upper) 
    stop("Lower bound must be smaller than upper bound")
  if (any(x < lower) || any(x > upper)) 
    stop("values in x must be between lower and upper")
  if (!is.vector(K, mode = "numeric") || length(K) > 1) 
    stop("K must be a numeric vector of length 1")
  if (is.na(K)) 
    stop("K cannot be NA or NaN")
  if (floor(K) <= 3 || is.infinite(K)) 
    stop("K must be a finite integer larger than 3")
  nx <- length(x)
  B <- matrix(0, nrow = nx, ncol = K)
  dimB <- dim(B)
  ndx <- K - 3
  dx <- (upper - lower)/ndx
  nknots <- ndx + 2 * 3 + 1
  knots <- seq(lower - 3 * dx, upper + 3 * dx, by = dx)
  for (i in 1:nx) {
    for (j in 1:(nknots - 4)) {
      temp <- 0
      cub <- x[i] - knots[j]
      if (cub > 0) {
        temp <- temp + cub^3
        cub <- x[i] - knots[j + 1]
        if (cub > 0) {
          temp <- temp - 4 * cub^3
          cub <- x[i] - knots[j + 2]
          if (cub > 0) {
            temp <- temp + 6 * cub^3
            cub <- x[i] - knots[j + 3]
            if (cub > 0) {
              temp <- temp - 4 * cub^3
              cub <- x[i] - knots[j + 4]
              if (cub > 0) {
                temp <- temp + cub^3
              }
            }
          }
        }
      }
      B[i, j] <- temp/(6 * dx^3)
      if (abs(B[i, j]) < 1e-10) 
        B[i, j] <- 0
    }
  }
  listout <- list(x = x, lower = lower, upper = upper, K = K, 
                  knots = knots, nknots = nknots, dimbasis = dimB, Bmatrix = B)
  attr(listout, "class") <- "cubicbs"
  listout
}