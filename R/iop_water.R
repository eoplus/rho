
#' Absorption coefficient of pure water
#'
#' This function calculates the temperature and salinity dependent absorption of
#' pure (saline) water.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param Tc     Temperature [0,30] (ºC).
#' @param S      Salinity [0,40] (parts per thousand).
#'
#' @details The harmonized pure water absorption data between 300 and 4000 nm 
#' from Röttgers et al. (2011) and the absorption derivatives with respect to 
#' temperature and salinity (Sullivan et al., 2006) are used to calculate the 
#' absorption at the desired conditions. For more information see ?a_water_wopp.
#'
#' Coefficients are linear interpolated from tabulated values (2 nm steps) to 
#' the requested wavelength.
#'
#' @return A numeric vector with the absorption coefficient (1/m) of pure 
#' (saline) water.
#'
#' @references
#' Röttgers, R.; Doerffer, R.; McKee, D.; Schönfeld, W. 2011. The Water Optical 
#' Properties Processor (WOPP). Pure water spectral absorption, scattering, and 
#' real part of refractive index model. Algorithm Theoretical Basis Document.
#' 
#' Sullivan, J. M.; Twardowski, M. S.; Zaneveld, J. R. V.; Moore, C. M.; 
#' Barnard, A. H.; Donaghay, P. L.; Rhoades. B. 2006. The hyperspectral 
#' temperature and salt dependencies of absorption by water and heavy water in 
#' the 400–750 nm spectral range. Applied Optics 45, 21, 5294-5309. DOI:
#' 10.1364/AO.45.005294
#'
#' @seealso \code{\link{n_water}}, \code{\link{b_water}}, \code{\link{pf_water}}
#'
#' @examples
#' # Retrieve the absorption of average seawater in the visible range:
#' a_water(lambda = 400:700, S = 34.72, Tc = 3.5)
#'
#' @export

a_water <- function(lambda, S = 0, Tc = 20) {

  if(.check_vect_args(arg_ls = list(lambda, S, Tc)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  if(any(lambda < 300) || any(lambda > 4000))
    stop("Requested wavelength outside data domain: [300,4000] (nm)")

  a    <- approx(a_water_wopp[, 1], a_water_wopp[, 2], xout = lambda)$y
  dads <- approx(a_water_wopp[, 1], a_water_wopp[, 3], xout = lambda)$y
  dadt <- approx(a_water_wopp[, 1], a_water_wopp[, 4], xout = lambda)$y

  a <- a + (Tc - 20) * dadt + S * dads

  return(a)
}
