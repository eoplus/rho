
#' Spectral convolution
#'
#' This function calculates the sum approximation for the convolution of a 
#' hyperspectral signal with a given sensor relative spectral response function.
#' 
#' @param spec A matrix containing spectral data, where the first column is the 
#'             wavelength reference of each measurement. 
#' @param rsrf A matrix containing the Relative Spectral Response Function of 
#'             the sensor, where the first column is the wavelength reference.
#'
#' @details The function calculates the sum approximation of the convolution. 
#' The RSRF is linearly interpolated to the spectral reference of the spec data, 
#' therefore spectral data beyond the interval of spec is implicitly taken as 
#' zero. 
#'
#' @return A matrix with number of rows equal to the number of bands and number 
#' of columns equal to number of spectral data in spec.
#'
#' @seealso \code{rsrf}
#'
#' @export

specconv <- function(spec, rsrf) {
  if(is.matrix(spec) & is.matrix(rsrf)) {
    if(dim(spec)[2] < 2 | dim(rsrf)[2] < 2) {
      stop('spec and srf must have at least two columns')
    }
  } else {
    stop('spec (and rsrf) must be a matrix')
  }

  specs_c <- matrix(NA, nrow = ncol(rsrf) - 1, ncol = ncol(spec))
  colnames(specs_c) <- colnames(spec)
  for(i in 2:ncol(rsrf)) {
    specs_c[i-1, 1] <- sum(rsrf[, 1] * rsrf[, i], na.rm = T) / 
                       sum(rsrf[, i], na.rm = T)
   }

  rsrf <- .approxt(rsrf, spec[, 1])
  for(j in 2:ncol(spec)) {
    for(i in 2:ncol(rsrf)) {
      if(sum(is.na(spec[, j] * rsrf[, i])) == nrow(spec)) {
        specs_c[i-1, j] <- NA
       } else {
        specs_c[i-1, j] <- sum(spec[, j] * rsrf[, i], na.rm = T) / 
                           sum(rsrf[, i], na.rm = T)
       }
     }
   }
   
  return(specs_c)
 }

.approxt <- function(xys, xout) {
  nc <- ncol(xys)
  xys_i <- matrix(NA, ncol = ncol(xys), nrow = length(xout))
  xys_i[, 1] <- xout
  for(i in 2:nc)
    xys_i[, i] <- stats::approx(x = xys[, 1], y = xys[, i], xout)$y
  colnames(xys_i) <- colnames(xys)
  return(xys_i)
 }

