
#' Convert between radians and degrees
#'
#' This functions convert between angles represented in radians and degrees.
#'
#' @param x Angle (radians or degrees)
#'
#' @return A numeric vector
#' 
#' @seealso \code{\link{snell}}, \code{\link{snell_decomp}}
#'
#' @export

rad <- function(x) {
  x * pi / 180
}

#' @rdname rad
#'
#' @export

deg <- function(x) {
  x * 180 / pi
}
