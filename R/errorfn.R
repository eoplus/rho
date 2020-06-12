
#' Error functions for minimization
#'
#' Those functions calculate the (weighted) difference between two vectors.
#'
#' @param x Reference vector.
#' @param y Test vector.
#' @param w Weights vector.
#'
#' @details Those functions are commonly used in the literature as error 
#' functions for minimization in inversion problems. \code{mlrt} computes the
#' mean log ratio, \code{mard} computes the mean absolute relative difference, 
#' \code{rmsd} computes the root mean squared difference, and \code{rrse} 
#' computes the relative root squared error between two vectors. All function 
#' accept a vector of weights, as in \code{weighted.mean}.
#'
#' @return A scalar with the difference between the two vectors.
#'
#' @seealso \code{\link{weighted.mean}}
#'
#' @export

mlrt <- function(x, y, w = rep(1, length(x))) {
  log(y / x) %>%
  abs() %>%
  weighted.mean(w, na.rm = T)
}

#' @rdname mlrt
#'
#' @export

mard <- function(x, y, w = rep(1, length(x))) {
  abs((y - x) / x) %>%
  weighted.mean(w, na.rm = T)
}

#' @rdname mlrt
#'
#' @export

rmsd <- function(x, y, w = rep(1, length(x))) {
  (y - x)^2 %>%
  weighted.mean(w, na.rm = T) %>%
  sqrt()
}

#' @rdname mlrt
#'
#' @export

rrsd <- function(x, y, w = rep(1, length(x))) {
  ((y - x)^2 * w) %>%
  sum() %>%
  sqrt() %>%
  `/`(sum(x * w))
}

