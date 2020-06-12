
#' Absorption coefficient of pure water
#'
#' This function calculates the temperature and salinity dependent absorption of
#' pure (saline) water.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param S      Salinity [0,40] (parts per thousand).
#' @param Tc     Temperature [0,30] (ºC).
#'
#' @details The harmonized pure water absorption data between 300 and 4000 nm 
#' from Röttgers et al. (2011) and the absorption derivatives with respect to 
#' temperature and salinity (Röttgers, 2010) are used to calculate the 
#' absorption at the desired conditions, following Röttgers et al. (2011). For 
#' more information see ?a_water_wopp.
#'
#' Coefficients are linear interpolated from tabulated values (2 nm steps) to 
#' the requested wavelength.
#'
#' @return A numeric vector with the absorption coefficient (1/m) of pure 
#' (saline) water.
#'
#' @references
#' Röttgers, R. 2010. Measurements of inherent optical properties of pure water. 
#' Technical note, STSE-WaterRadiance D5, ESA.
#'
#' Röttgers, R.; Doerffer, R.; McKee, D.; Schönfeld, W. 2011. The Water Optical 
#' Properties Processor (WOPP). Pure water spectral absorption, scattering, and 
#' real part of refractive index model. Algorithm Theoretical Basis Document.
#'
#' @seealso \code{\link{n_water}}, \code{\link{vsf_water}}, \code{\link{b_water}}
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

#' Scattering of pure (saline) water
#'
#' This functions calculate the angular integral (scattering coefficient) and 
#' the angular distribution (volume scattering function) of  pure (saline) water 
#' scattering.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param S      Salinity [0,40] (parts per thousand).
#' @param Tc     Temperature [0,30] (ºC).
#' @param psi    Polar scattering angle (radians).
#' @param P      Pressure (bar in excess of 1 ATM).
#'
#' @details The formulation of scattering of pure (saline) water is based on the 
#' Einstein-Smoluchowski theory of scattering in a particle free media due to 
#' fluctuations of the dieletric constant (square of the refractive index) 
#' caused by random motion of molecules. This implementation follows Zhang & Hu
#' (2009A) that revisited the equations formulated with respect to the density 
#' derivative of the refractive index, after the new determinations of Proutiere 
#' et al. (1992). The excess scattering caused by diluted salts is based on salt 
#' concentration fluctuations by incorporating the salinity derivative of the 
#' refractive index (Zhang & Hu 2009B).
#' 
#' The density fluctuations are given by thermodynamic statistics and since both 
#' the density, the refractive index and the isothermal compressibility of water
#' are a function of pressure, the pressure is allowed to vary. Note however 
#' that the model is only validated for surface pressure. The default pressure 
#' therefore is 0 bar, with P defined as excess pressure to 1 ATM (in bar).
#'
#' The volume scattering function is provided per polar scattering angle. In 
#' natural waters, the angular distribution of scattering is only a function of 
#' the polar angle of scattering, i.e, it is constant for all azimuths at a 
#' given polar scattering angle (Jonaz & Fournier, 2007).
#'
#' Finally, this function will provide values slightly higher than the ones 
#' provided by the Water Optical Properties Processor (WOPP; Röttgers et al., 
#' 2011) due to different equations for water density. See \code{d_water} for 
#' details. The required physical constants are taken from the 2018 edition of
#' CODATA, available from: \url{http://physics.nist.gov/cuu/Constants/}
#'
#' @return A numeric vector with the volume scattering function (m^-1 sr^-1) of 
#' pure (saline) water for function \code{vsf_water} and a numeric vector with 
#' the scattering coefficient.
#'
#' @seealso \code{\link{a_water}}, \code{\link{n_water}}
#'
#' @examples
#' # Retrieve the volume scattering function of average seawater at 550 nm:
#' psi <- seq(0, pi, length.out = 180)
#' vsf_water(lambda = 550, S = 34.72, Tc = 3.5, psi = psi)
#'
#' # Retrieve the scattering coefficient of average seawater in the visible 
#' range:
#' b_water(lambda = 400:700, S = 34.72, Tc = 3.5)
#'
#' @export

vsf_water <- function(lambda, S = 0, Tc = 20, psi = pi/2, P = 0) {

  # Depolarization ratio from Farinato & Rowel (1976):
  delta   <- 0.039

  beta_90 <- .beta_water_90(lambda = lambda, S = S, Tc = Tc, P = P)
  beta_w  <- beta_90 * (1 + cos(psi)^2 * (1 - delta) / (1 + delta))

  return(beta_w)

}

#' @rdname vsf_water
#'
#' @export

b_water <- function(lambda, S = 0, Tc = 20, P = 0) {

  # Depolarization ratio from Farinato & Rowel (1976):
  delta   <- 0.039 

  beta_90 <- .beta_water_90(lambda = lambda, S = S, Tc = Tc, P = P)
  b       <- 8 * pi * beta_90 * (2 + delta) / (1 + delta) / 3

  return(b)

}

#' Volume scattering at 90 degrees (m^-1 sr^-1)

.beta_water_90 <- function(lambda, S = 0, Tc = 20, P = 0) {

  K      <- 1.380649E-23   # Boltzmann constant (joule/Kelvin)
  Na     <- 6.02214076E23  # Avogrado's constant (1/mol)
  M0_H2O <- 0.01801528     # Molecular weight of water (kg/mol)

  # Depolarization ratio from Farinato & Rowel (1976):
  delta   <- 0.039  
  f_delta <- (6 + 6 * delta) / (6 - 7 * delta)  # Cabbanes factor of water

  n       <- n_water(S = S, Tc = Tc, lambda = lambda)
  d       <- d_water(S = S, Tc = Tc, P = P)   
  dfri    <- liquid_DFRI(n = n, d = d * 1E-3)  # kg m^-3 to kg dm^-3
  dnds    <- .dnds(lambda = lambda, Tc = Tc)
  IsoComp <- 1E-5 / .K_w(S = S, Tc = Tc, P = P)
  dlnawds <- .dlna0dS(Tc = Tc, S = S)

  lambda <- lambda * 1E-9  # Wavelength in m
  L4     <- lambda^4
  Tk     <- Tc + 273.15    # Temperature in Kelvin

  beta_d    <- (pi^2 / (2 * L4)) * dfri^2 * K * Tk * IsoComp * f_delta
  flu_con   <- S * M0_H2O * dnds^2 / d / -dlnawds / Na
  beta_c    <- 2 * pi^2 * n^2 * flu_con * f_delta / L4
  beta_w_90 <- beta_d + beta_c

  return(beta_w_90)

}

#' Partial derivative of seawater refractive index with respect to salinity, 
#' dn/dS. Based on Quan & Fry (1995).

.dnds <- function(lambda, Tc) {

  Tc2  <- Tc^2

  k1   <-  1.779E-4
  k2   <- -1.05E-6
  k3   <-  1.6E-8
  k6   <-  0.01155

  dnds <- (k1 + k2 * Tc + k3 * Tc2) + (k6 * n_air(lambda) / lambda) 

  return(dnds)

}

#' Partial derivative of the natural logarithm of solvent activity with respect 
#' to salinity, dln(a0)/dS is from Millero & Leung (1976). Has an estimated 
#' error of 0.04%.

.dlna0dS <- function(S, Tc) {

  Tc2 <- Tc^2
  Tc3 <- Tc^3

  a0  <- -5.58651E-4
  a1  <-  2.40452E-7
  a2  <- -3.12165E-9
  a3  <-  2.40808E-11
  a4  <-  1.79613E-5
  a5  <- -9.9422E-8
  a6  <-  2.08919E-9
  a7  <- -1.39872E-11
  a8  <- -2.31065E-6
  a9  <- -1.37674E-9
  a10 <- -1.93316E-11

  dlna0dS <- (a0 + a1 * Tc + a2 * Tc2 + a3 * Tc3) + 
             1.5 * S^0.5 * (a4 + a5 * Tc + a6 * Tc2 + a7 * Tc3)  + 
             2 * S * (a8 + a9 * Tc + a10 * Tc2)

  return(dlna0dS)

}


