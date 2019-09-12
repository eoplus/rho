
#' Absorption band deconvolution
#'
#' This function deconvolves a measured signal into a sum of model functions.
#'
#' @param spec   A vector with the spectra to be deconvolved.
#' @param sr     A vector with the spectral reference of the spectra.
#' @param start  A matrix with first guess parameters. See details.
#' @param lower  A matrix with the lower limit of the parameters to be fitted.
#'               See details.
#' @param upper  A matrix with the upper limit of the parameters to be fitted.
#'               See details.
#' @param ef     Error function to be minimized.
#' @param model  One of 'Gauss', 'Lorentz' or 'Voigt'. See details.
#' @param maxit  Maximum number of iterations.
#' @param reltol Relative tolerance for convergence criteria.
#'
#' @details This function is a wrapper for the \code{optim} function, with the 
#' Nelder-Mead method, for deconvolution of a spectrum as sum of Gaussian, 
#' Lorentzian or Voigt curves. It is appropriate for deconvolution of absorption 
#' peaks if the peaks are present in or just beyond the limits of the 
#' measured spectrum. When applied for deconvolution of pigment absorption, it 
#' is important that other absorption signals have been removed form the 
#' spectrum for a meaningful result. The spectral reference scale may be any 
#' desired (e.g., energy, wavenumber, wavelength), but the fitting space, 
#' \code{start} values and \code{sr} must be specified in the correct scale.
#'
#' The arguments \code{start}, \code{lower}, and \code{upper}, must be three 
#' column matrices for models 'Gauss' and 'Lorentz', with center, sigma 
#' ('Gauss') or hwhm ('Lorentz') and scale, respectively and peaks per
#' row. For the 'Voigt' model, it should have four columns, with center, sigma, 
#' hwhm and scale (in this order). The number of absorption bands used depends 
#' on the number of rows and the fitting space is completely defined by the 
#' lower and upper boundaries (may be infinite). The boundaries are enforced by 
#' returning NA on the error function evaluation when optim test parameter 
#' values outside the fitting space.
#'
#' The spectra to be deconvolved may contain NA, but not \code{sr}, 
#' \code{start}, \code{lower}, and \code{upper}. \code{maxit} and \code{reltol}
#' are passed to optim. A set of error functions is provided: \code{mard}, 
#' \code{rmsd}, \code{mlrt}, and \code{rrsd}. For more details on the Gaussian, 
#' Lorentzian and Voigt parameters, see \code{apeaks}.
#'
#' @return The return of optim, with fitted coefficients converted to matrix
#' format for convenience.
#'
#' @references
#' Hoepffner, N.; Sathyendranath, S. 1991. Effect of pigment composition on 
#' absorption properties of phytoplankton. Marine Ecology Progress Series 73, 
#' 11-23. DOI: 10.3354/meps073011
#' 
#' @seealso \code{\link{apeaks}}, \code{\link{mard}}, \code{\link{rmsd}}, 
#' \code{\link{mlrt}}, and \code{\link{rrsd}}
#' 
#' @export

adeconv <- function(spec, sr, start, lower, upper, ef, model = c('Gauss', 
  'Lorentz', 'Voigt'), maxit = 5E4, reltol = 1E-12) {

  if(any(lower > start))
    stop("Lower limits higher than start values")
  if(any(upper < start))
    stop("Upper limits lower than start values")

  if(model == 'Voigt' && (ncol(start) != 4 || ncol(lower) != 4 || ncol(upper) != 4))
    stop("start, upper and lower must have four columns for Voigt model")
  if(model != 'Voigt' && (ncol(start) != 3 || ncol(lower) != 3 || ncol(upper) != 3))
    stop("start, upper and lower must have three columns for Gauss and Lorentz models")

  if(model == 'Voigt') {
    start <- start[, c(1, 2, 4, 3)]
    lower <- lower[, c(1, 2, 4, 3)]
    upper <- upper[, c(1, 2, 4, 3)]
  }

  # Function to be minimized:
  fn <- function(x, spec, ef, l, u, pre = FALSE, center, sigma, hwhm, sr, model) {

    if(pre) {
      if(any(x < l[, 3]) || any(x > u[, 3])) return(NA)
      est <- apeaks(center = center, sigma = sigma, hwhm = hwhm, scale = x, 
               sr = sr, model = model)

    } else {
      tp  <- matrix(x, ncol = ifelse(model == 'Voigt', 4, 3))
      if(any(tp < l) || any(tp > u)) return(NA)
      est <- apeaks(center = tp[, 1], sigma = tp[, 2], hwhm = 
               tp[, ifelse(model == 'Voigt', 4, 2)], scale = tp[, 3], sr = sr, 
               model = model)
    }

    dif <- ef(spec, est)
    return(dif)
  }

  # Prefit:
  start[, 3] <- optim(par = start[, 3], fn = fn, center = start[, 1], 
                      sigma = start[, 2], sr = sr, spec = spec, ef = ef, 
                      l = lower, u = upper, model = model, pre = TRUE, 
                      method = "Nelder-Mead", control = list(maxit = maxit, 
                      reltol = reltol), hwhm = start[, ifelse(model == 'Voigt', 
                      4, 2)])$par

  # Fit:
  opt        <- optim(par = as.vector(start), fn = fn, sr = sr, spec = spec, 
                      ef = ef, l = lower, u = upper, model = model, pre = FALSE, 
                      method = "Nelder-Mead", control = list(maxit = maxit, 
                      reltol = reltol))

  if(model == 'Voigt') {
    opt$par   <- matrix(opt$par, ncol = 4)[, c(1, 2, 4, 3)]
  } else {
    opt$par   <- matrix(opt$par, ncol = 3)
  }
  return(opt)
}

#' Calculate absorption spectra
#'
#' This function evaluates an absorption band model at the desired reference 
#' positions.
#'
#' @param center Center position reference of each absorption band. 
#' @param sigma  Standard deviation for the Gaussian model.
#' @param hwhm   Half width at half maximum for the Lorentzian model.
#' @param scale  Scale height at the center of each absorption band.
#' @param sr     Position reference at which the absorption models should be 
#'               evaluated in any suitable scale (energy, frequency, wavelength).
#' @param model  One of 'Gauss', 'Lorentz' or 'Voigt'. See details.
#' @param peaks  Logical. Should the peak curves for each absorption band be 
#'               returned?
#' 
#' @details Observed absorption bands are not sharp due to several line 
#' broadening process grouped into homogeneous processes represented by a 
#' Lorentzian profile and inhomogeneous processes represented by a Gaussian 
#' profile. The Voigt profile combine both processes. 
#'
#' The line width parameters for each profile are requested in their most common
#' formats, standard deviation for the Gaussian model and half width at half 
#' maximum (HWHM) for the Lorentzian model. The HWHM is related to the standard 
#' deviation for the Gaussian function by:
#' HWHM = sqrt(2*log(2)) * sigma
#'
#' The Voigt profile is calculated with the analytical approximation of Liu 
#' et al. (2001) and the empirical approximation of Olivero & Longbothum (1977) 
#' for the HWHM of the Voigt profile.
#'
#' Units for the parameters depend on the scale used to represent the absorption 
#' process (e.g., energy, frequency, wavenumber, wavelength).
#'
#' @return If peaks is FALSE, a numeric vector with the summed peaks. If peaks 
#' is true, a matrix is returned, with peaks per row.
#'
#' @references
#' Liu, Y.; Lin, J.; Huang, G.; Guo, Y.; Duan, C. 2001. Simple empirical 
#' analytical approximation to the Voigt profile. Journal of the Optical Society 
#' of America B 18, 5, 666-672. DOI: 10.1364/JOSAB.18.000666
#'
#' Olivero, J. J.; Longbothum, R. L. 1977. Empirical fits to the Voigt line 
#' width: A brief review. Journal of Quantitative Spectroscopy and Radiative 
#' Transfer 17, 2, 233-236. DOI: 10.1016/0022-4073(77)90161-3
#'
#' @seealso \code{\link{adeconv}}
#'
#' @examples
#' sr <- seq(-10, 10, 0.01)
#' pg <- apeaks(0, 1, 0, 1, sr = sr, 'Gauss')
#' pl <- apeaks(0, 0, sqrt(2*log(2)), 1, sr = sr, 'Lorentz')
#' pv <- apeaks(0, 0.7638, 0.4 * sqrt(2*log(2)), 1, sr = sr, 'Voigt')
#' plot(sr, pg, type = 'l', col = 2, ylab = 'Height')
#' lines(sr, pl, col = 3)
#' lines(sr, pv, col = 4)
#' legend('topright', c('Gauss', 'Lorentz', 'Voigt'), col = 2:4, lty = 1, bty = 'n')
#' legend('topright', expression(HWHM == sqrt(2*log(2))), bty = 'n')
#
#' @export

apeaks <- function(center, sigma, hwhm, scale, sr, model = c('Gauss', 'Lorentz', 
  'Voigt'), peaks = FALSE) {

  if((model == 'Gauss' || model == 'Voigt') && missing(sigma))
    stop('sigma must be specified for the Gauss and Voigt models')
  if((model == 'Lorentz' || model == 'Voigt') && missing(hwhm))
    stop('hwhm must be specified for the Lorentz and Voigt models')

  if(model == 'Gauss') {
    x <- cbind(center, sigma, scale)
    apeaks <- t(apply(x, 1, .gpeak, sr = sr))
  } else if(model == 'Lorentz') {
    x <- cbind(center, hwhm, scale)
    apeaks <- t(apply(x, 1, .lpeak, sr = sr))
  } else if(model == 'Voigt') {
    x <- cbind(center, sigma, hwhm, scale)
    apeaks <- t(apply(x, 1, .vpeak, sr = sr))
  }

  if(peaks) return(apeaks)
  else return(apply(apeaks, 2, sum))

}

.gpeak <- function(x, sr) {
  # x[1] = center, x[2] = sigma, x[3] = scale
  x[3] * exp(-0.5 * ((sr - x[1]) / x[2])^2)
}

.lpeak <- function(x, sr) {
  # x[1] = center, x[2] = hwhm, x[3] = scale
  x[3] * x[2]^2 / ((x[1] - sr)^2 + x[2]^2) 
}

.vpeak <- function(x, sr) {
  # x[1] = center, x[2] = sigma, x[3] = hwhm, x[4] = scale
  hwhm_g <- x[2] * sqrt(2 * log(2))
  hwhm_l <- x[3]
  hwhm_v <- 0.5 * (1.0692 * hwhm_l + sqrt(0.86638 * hwhm_l^2 + 4 * hwhm_g^2))
  
  d   <- (hwhm_l - hwhm_g) / (hwhm_l + hwhm_g)
  cg  <- 0.32460 - 0.61825 * d + 0.17681 * d^2 - 0.12109 * d^3
  cl  <- 0.68188 + 0.61293 * d - 0.18384 * d^2 - 0.11568 * d^3

  gpk <- .gpeak(c(x[1], hwhm_v / sqrt(2 * log(2)), x[4]), sr)
  lpk <- .lpeak(c(x[1], hwhm_v, x[4]), sr)
  
  vpk <- cg * gpk + cl * lpk
  as.numeric(vpk * x[4] / max(vpk))
}

