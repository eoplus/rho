
#' Calculate retrieval statistics 
#'
#' The function calculates useful descriptive statistics to evaluate the 
#' estimation and create a vector of strings for plotting.
#'
#' @param x     Reference data.
#' @param y     Estimated data.
#' @param units Units of the data.
#'
#' @details Those functions are wrappers for descriptive statistics and 
#' convenient text strings for addition to a plot. They provide the statistics 
#' of the agreement of the estimated data against the reference data; i.e.
#' evaluates estimations against the 1:1 relation. Statistics include the mean 
#' absolute percentage error (MAPE), the mean percentage error (bias), the root 
#' mean squared error (RMSE), the coefficient of determination (R^2), the number
#' of data points, and the range of the reference data.
#'
#' The parameter units is also available for convenience of use with the
#' \code{lstats} function, to add those statistics to a plot. It should be a 
#' character string that can be used in expressions. See examples.
#'
#' @return A named list with component statistics and units for \code{rstat} and 
#' a vector of expressions for plot legend.
#'
#' @examples
#' x <- 1:10
#' y <- x + rnorm(10)
#' (rs <- rstat(x, y, units = bquote(m^-1)))
#' plot(x, y)
#' legend("bottomright", lstat(rs), bty = "n")
#'
#' @export

rstat <- function(x, y, yp, units = '') {
  ymx  <- y - x
  mape <- mean(abs(ymx / x), na.rm = TRUE) * 100
  bias <- mean((ymx) / x, na.rm = TRUE) * 100
  rmse <- sqrt(mean((ymx)^2, na.rm = TRUE))
  r2   <- 1 - (sum(ymx^2, na.rm = T)/(sum((x - mean(x, na.rm = T))^2, na.rm = T)))
  n    <- nrow(na.omit(cbind(x, y)))
  rang <- range(na.omit(cbind(x, y))[, 1], na.rm = T) 
  res  <- list(mape = mape, bias = bias, rmse = rmse, r2 = r2, n = n, 
               rang = rang, units = units)
  return(res)
}

#' @rdname rstat
#'
#' @param stats  object created by \code{rstat}.
#' @param digits number of digits to be passed to \code{round}.
#'
#' @export

lstat <- function(stats, digits = 3) {

  r2   <- substitute(expression(R^2 == r2),
            env = list(r2 = round(stats$r2, digits))) %>% 
          eval()

  rmse <- substitute(expression(RMSE == rmse~units),
            list(rmse = round(stats$rmse, digits), 
            units = stats$units)) %>%
          eval()

  mape <- substitute(expression(MAPE == mape~"%"),
            list(mape = round(stats$mape, digits))) %>%
          eval()

  bias <- substitute(expression(Bias == bias~"%"),
            list(bias = round(stats$bias, digits))) %>%
          eval()

  n    <- substitute(expression(N == n), list(n = stats$n)) %>%
          eval()

  rang <- substitute(expression(Range == a - b~units),
            list(a = round(stats$rang[1], digits), 
            b = round(stats$rang[2], digits), units = stats$units)) %>%
          eval()

  ltext <- c(r2, rmse, mape, bias, n, rang)

  return(ltext)

}


