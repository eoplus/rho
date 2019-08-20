
#' Gaussian deconvolution
#'
#' This function deconvolves a measured signal into a sum of Gaussian functions.
#'
#' @param spec   A vector with the spectra to be deconvolved.
#' @param sr     A vector with the spectral reference of the spectra.
#' @param start  A matrix with first guess parameters. See details.
#' @param lower  A matrix with the lower limit of the parameters to be fitted.
#'               See details.
#' @param upper  A matrix with the upper limit of the parameters to be fitted.
#'               See details.
#' @param ef     Error function to be minimized.
#' @param maxit  Maximum number of iteractions.
#' @param reltol Relative tolerance for convergence criteria.
#'
#' @details This function is a wrapper for the \code{optim} function, with the 
#' Nelder-Mead method, optimized for deconvolution of a spectrum as sum of 
#' Gaussian curves. It is appropriate for deconvolution of absorption peaks if 
#' the peaks are present in or just beyond the limits of the measured spectrum. 
#' When applied for deconvolution of pigment absorption, it is important that 
#' other absorption signals have been removed form the spectrum for a meaningful
#' result. The spectral reference scale may be wavelength or wavenumber, 
#' but the fitting space, \code{start} values and \code{sr} must be specified in 
#' the correct scale.
#'
#' The arguments \code{start}, \code{lower}, and \code{upper}, must be three 
#' column matrices, with center, sigma and scale, respectively and peaks per
#' row. The number of Gaussian curves used depends on the number of rows and the 
#' fitting space is completely defined by the lower and upper boundaries (may be
#' infinite). The boundaries are enforced by returning NA on the error function 
#' evaluation when optim test parameter values outside the fitting space.
#'
#' The spectra to be deconvolved may contain NA, but not \code{sr}, 
#' \code{start}, \code{lower}, and \code{upper}. \code{maxit} and \code{reltol}
#' are passed to optim. A set of error functions is provided: \code{mard}, 
#' \code{rmsd}, \code{mlrt}, and \code{rrsd}.
#'
#' @return The return of optim, with fitted coefficients converted to matrix
#' format for convenience.
#'
#' @references
#' Hoepffner, N.; Sathyendranath, S. 1991. Effect of pigment composition on 
#' absorption properties of phytoplankton. Marine Ecology Progress Series 73, 
#' 11-23. DOI: 10.3354/meps073011
#' 
#' @seealso \code{\link{gpeaks}}, \code{\link{mard}}, \code{\link{rmsd}}, 
#' \code{\link{mlrt}}, and \code{\link{rrsd}}
#' 
#' @export

gdeconv <- function(spec, sr, start, lower, upper, ef, maxit = 5E4, 
                    reltol = 1E-12) {

  if(any(lower > start))
    stop("Lower limits higher than start values")
  if(any(upper < start))
    stop("Upper limits lower than start values")

  # Function to be minimized:
  fn <- function(x, center, sigma, sr, spec, ef, l, u, pre = FALSE) {

    if(pre) {
      if(any(x < l[, 3]) || any(x > u[, 3])) return(NA)
      est <- gpeaks(center = center, sigma = sigma, scale = x, sr = sr, 
                   peaks = FALSE)
    } else {
      if(any(x < as.vector(l)) || any(x > as.vector(u))) return(NA)
      tp  <- 
      est <- gpeaks(center = matrix(x, ncol = 3), sr = sr, peaks = FALSE)
    }

    dif <- ef(spec, est)
    return(dif)
  }

  # Prefit:
  start[, 3] <- optim(par = start[, 3], fn = fn, center = start[, 1], 
                      sigma = start[, 2], sr = sr, spec = spec, ef = ef, 
                      l = lower, u = upper, pre = TRUE, method = "Nelder-Mead", 
                      control = list(maxit = maxit, reltol = reltol))$par

  # Fit:
  opt        <- optim(par = as.vector(start), fn = fn, sr = sr, spec = spec, 
                      ef = ef, l = lower, u = upper, pre = FALSE, method = 
                      "Nelder-Mead", control = list(maxit = maxit, 
                      reltol = reltol))

  opt$par   <- matrix(opt$par, ncol = 3)
  
  return(opt)
}

#' Calculate Gaussian absorption spectra
#'
#' This function evaluates Gaussian function with the input parameters at the 
#' desired spectral reference points.
#'
#' @param center Center spectral reference of each Gaussian curve. 
#'               Alternatively, a matrix in which columns are center, sigma and
#'               scale. See details.
#' @param sigma  Standard deviation in the appropriate spectral reference units.
#' @param scale  Scale height at the center of each Gaussian curve.
#' @param sr     Spectral reference at which Gaussians should be evaluated.
#' @param peaks  Logical. Should the Gaussian curves for each absorption band be 
#'               returned?
#' 
#' @details Units for center and sigma depend on the spectral scale used during
#' parameter fit and may be any convenient scale were absorption bands can be 
#' described by a Gaussian curve (e.g., wavelength in nm or wavenumber in 1/cm).
#'
#' If \code{center} is matrix, it must have three columns (center, sigma, 
#' scale), with absorption peaks by row. Note that if sigma and scale parameters
#' take precedence and will be used instead of the matrix columns if specified.
#'
#' @return If peaks is FALSE, a numeric vector with the summed peaks. If peaks 
#' is true, a matrix is returned, with peaks per row.
#'
#' @seealso \code{\link{gdeconv}}
#'
#' @export

gpeaks <- function(center, sigma, scale, sr, peaks = FALSE) {

  if(is.matrix(center)) {
    if(dim(center)[2] == 3) {
      if(missing(sigma)) sigma <- center[, 2]
      if(missing(scale)) scale <- center[, 3]
      center <- center[, 1]
    } else if(missing(sigma) || missing(scale)) {
      stop("sigma and scale must be specified if center is not a 3-column matrix")
    }
  } else if(missing(sigma) || missing(scale)) {
      stop("sigma and scale must be specified if center is not a 3-column matrix")
  }

  gauss <- matrix(NA, ncol = length(sr), nrow = length(center))
  for(i in 1:length(center)) {
    gauss[i, ]  <- (scale[i]  * exp(-0.5 * ((sr - center[i]) / sigma[i])^2))
  }

  if(peaks) return(gauss)
  else return(apply(gauss, 2, sum))

}

