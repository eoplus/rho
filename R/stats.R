
#' Calculate retrieval statistics 
#'
#' The function calculates useful descriptive statistics to evaluate the 
#' estimation and create a vector of strings for plotting.
#'
#' @param x     Reference data.
#' @param y     Estimated data.
#' @param log   Logical. Should log based performance statistics be computed?
#' @param ux    Logical. Calculate unbiased differences?
#' @param units Units of the data.
#'
#' @details Those functions are wrappers for descriptive statistics and 
#' convenient text strings for addition to a plot. They provide the statistics 
#' of the agreement of the estimated data against the reference data; i.e.
#' evaluates estimations against the 1:1 relation. Statistics include the mean 
#' absolute percentage difference (MAPD), the mean percentage error (bias), the 
#' root mean squared difference (RMSD), the coefficient of determination (R^2), 
#' the number of data points, and the range of the reference data.
#'
#' The parameter 'units' is also available for convenience of use with the
#' \code{lstats} function, to add those statistics to a plot. It should be a 
#' character string that can be used in expressions. See examples.
#'
#' When x does not present smaller uncertainty than y, ux can be set to TRUE to 
#' calculate unbiased statistics in linear space. In that case, the normalization 
#' of MAPD and bias is not on x by on the average of x and y. This approach 
#' might also be used if some values are very close to zero.
#'
#' For heterocedastic data or data spanning multiple orders of magnitude, log 
#' can be set to TRUE to calculate performance statistics in log space. Note 
#' in this case the units are set to '' (unitless). The MAPD and bias are 
#' retrieved in % with:
#'
#' MAPD = 100 * (10^((1 / n) * sum(|log10(x) - log10(y)|)) - 1)
#' Bias = 100 * (10^((1 / n) * sum(log10(x) - log10(y))) - 1)
#'
#' @return A named list with component statistics and units for \code{rstat} and 
#' a vector of expressions for plot legend.
#'
#' @references
#' IOCCG. 2019. Uncertainties in Ocean Colour Remote Sensing. Mélin F. (ed.), 
#'   IOCCG Report Series, No. 18, International Ocean Colour Coordinating Group, 
#'   Dartmouth, Canada. http://dx.doi.org/10.25607/OBP-696
#' Seegers, B. N.; Stumpf, R. P.; Schaeffer, B. A.; Loftin, K. A.; Werdell, P. 
#'   J. 2018. Performance metrics for the assessment of satellite data products: 
#'   an ocean color case study. Optics Express 26, 6, 7404. 
#'   DOI: 10.1364/OE.26.007404
#'
#' @examples
#' x <- 1:10
#' y <- x + rnorm(10)
#' (rs <- rstat(x, y, units = bquote(m^-1)))
#' plot(x, y)
#' legend("bottomright", lstat(rs), bty = "n")
#'
#' @export

rstat <- function(x, y, log = FALSE, ux = FALSE, units = '') {

  if(log) {

    if(ux)
      warning("Performance statistics in log space do not include normalization")

    x  <- log10(x)
    y  <- log10(y)
    xy <- na.omit(cbind(x, y))
    d  <- xy[, 2] - xy[, 1]

    mapd <- 100 * (10^mean(abs(d), na.rm = TRUE) - 1)
    bias <- 100 * (10^mean(d, na.rm = TRUE) - 1)

  } else {

    xy <- na.omit(cbind(x, y))
    d  <- xy[, 2] - xy[, 1]
    if(ux) {
      norm <- apply(xy, 1, mean)
    } else {
      norm <- xy[, 1]
    }

    mapd <- mean(abs(d / norm), na.rm = TRUE) * 100
    bias <- mean(d / norm, na.rm = TRUE) * 100
  }

  rmsd <- sqrt(mean(d^2, na.rm = TRUE))
  rss  <- sum(d^2)
  tss  <- sum((xy[, 1] - mean(xy[, 1]))^2)
  r2   <- 1 - rss/tss
  n    <- nrow(xy)
  rang <- range(xy[, 1], na.rm = T) 
  res  <- list(mapd = mapd, bias = bias, rmsd = rmsd, r2 = r2, n = n, 
               rang = rang, units = units)
  res

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

  rmsd <- substitute(expression(RMSD == rmsd~units),
            list(rmsd = round(stats$rmsd, digits), 
            units = stats$units)) %>%
          eval()

  mapd <- substitute(expression(MAPD == mapd~"%"),
            list(mapd = round(stats$mapd, digits))) %>%
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

  ltext <- c(r2, rmsd, mapd, bias, n, rang)

  return(ltext)

}


