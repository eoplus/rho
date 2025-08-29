
#' Dispersion formula for the refractive index of air
#'
#' This function calculates the real part of the refractive index of standard 
#' dry air relative to vaccum. See details.
#'
#' @param lambda Wavelength in vacuum (nm).
#'
#' @details The dispersion model implemented is that by Ciddor (1996), for the 
#' refractive index of standard dry air. It is reference for dry air (0 
#' moisture), 15 ºC, 101.325 Pa and 450 ppm CO2. This equation is recomended by 
#' Zhang et al. (2009) for the convertion of (saline) water refractive index 
#' relative to air to that relative to vacuum.
#'
#' Ciddor (1996) also provides equations to account for different levels of CO2 
#' and water vapor. Those are not implemented in the current version.
#'
#' The valid range is 200 nm to 1100 nm. If the refractive index is requested 
#' outside this range, the function will return a warning.
#'
#' @return A numeric vector with the real part of refractive index of standard 
#' dry air (unitless).
#'
#' @references
#' Ciddor, P. E. 1996. Refractive index of air: new equations for the visible 
#' and near infrared. Applied Optics 35, 9, 1566-73. DOI: 10.1364/AO.35.001566
#'
#' Zhang, X.; Hu, L.; He, M.-X. 2009. Scattering by pure seawater: effect of 
#' salinity. Optics Express 17, 7, 5698-710. DOI: 10.1364/OE.17.005698
#'
#' @seealso \code{\link{n_water}}, \code{\link{n_cellulose}}, 
#' \code{\link{n_calcite}}, \code{\link{n_quartz}}, \code{\link{snell}}, 
#' \code{\link{fresnel}}
#'
#' @examples
#' # Get the refractive index for teh visible range: 
#' n_air(400:700) 
#'
#' @export

n_air <- function(lambda) {

  if(any(lambda < 200) || any(lambda > 1100))
    warning("Requested wavelength outside model domain: [200,1100]")

  l2  <- (lambda * 1E-3)^-2  # Convert to wavlength in microns and square

  k0  <- 5792105
  k1  <- 238.0185
  k2  <- 167917
  k3  <- 57.362

  n <- 1 + ((k0 / (k1 - l2) + k2 / (k3 - l2)) * 1E-8)

  return(n)

}

#' Dispersion formula for the refractive index of liquid water
#'
#' This function calculates the real part of the refractive index of pure water
#' relative to vaccum. The complex refractive index can also be returned. See 
#' details.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param Tc     Temperature (ºC).
#' @param S      Salinity (parts per thousand).
#' @param P      Pressure (bar in excess of 1 ATM).
#' @param type   One of: "real" or "complex". See details.
#'
#' @details The temperature and salinity dependent spectral refractive index of 
#' (saline) water is based on the empirical dispersion equation parametrized by 
#' Quan & Fry (1995) for the visible range. Huibers (1997) analyzed the 
#' available data and models and attested that the model could be extrapolated 
#' to the UV and NIR (220 to 1100 nm).
#' 
#' The Quan & Fry (1995) formulation is for standard atmosphere pressure. It can 
#' be converted to any desired pressure with the isothermic density derivative 
#' of the refractive index of liquids parametrized by Proutiere et al. (1992). 
#' That pressure, in bar, should be in excess of surface pressure, i.e., to 
#' retrieve surface values, P = 0.
#'
#' The imaginary part of the refractive index of pure water is calculated from
#' the temperature and salinity dependent absorption of pure (saline) water (
#' function \code{a_water}).
#'
#' The original parametrization of Quan & Fry (1995) is based on refractive 
#' index relative to air measured by Austin & Halikas (1976), and the function 
#' will return the values relative to vacuum by using the refractive index of 
#' standard dry air (\code{n_air}), as recommended by Zhang et al. (2009).
#' 
#' Another parametrization with comparable to that of Quan & Fry (1995) in the 
#' visible range Huibers, 1997) is that of Schiebener et al. (1990). This model
#' has an extended range (200 to 2500 nm) and could be added as an option in a
#' future version.
#'
#' @return A numeric (or complex) vector of the refractive index (unitless) of 
#' pure (saline) water relative to vacuum .
#'
#' @references
#' Austin, R. W.; Halikas, G. 1976. The index of refraction of seawater. SIO 
#' Ref. 76-1, Scripps Institution of Oceanography, La Jolla, California.
#'
#' Huibers, P. D. T. 1997. Models for the wavelength dependence of the index of 
#' refraction of water. Applied Optics 36, 16, 3785-3787. DOI: 
#' 10.1364/ao.36.003785
#'
#' Jonasz, M.; Fournier, G. R. 2007. Light Scattering by Particles in Water - 
#' Theoretical and Experimental Foundations. Academic Press, 715 pp.
#'
#' Proutiere, A.; Megnassan, E.; Hucteau, H. 1992. Refractive index and density 
#' variations in pure liquids: A new theoretical relation. The Journal of 
#' Physical Chemistry 96, 3485-3489. DOI: 10.1021/j100187a058
#'
#' Quan, X.; Fry, E. S. 1995. Empirical equation for the index of refraction of 
#' seawater. Applied Optics 34, 3477–3480. DOI: 10.1364/AO.34.003477
#'
#' Schiebener, P.; Straub, J.; Levelt Sengers, J. M. H.; Gallagher, J. S. 1990. 
#' Refractive index of water and steam as a function of wavelength, temperature 
#' and density. Journal of Physical and Chemical Reference Data 19, 677–717. 
#' DOI: 10.1063/1.555859
#'
#' Zhang, X. Hu, L. 2009. Estimating scattering of pure water from density 
#' fluctuation of the refractive index. Optics Express 17, 3, 1671-1678. DOI: 
#' 10.1364/OE.17.001671
#'
#' @seealso \code{\link{n_air}}, \code{\link{n_cellulose}}, 
#' \code{\link{n_calcite}}, \code{\link{n_quartz}}, \code{\link{snell}}, 
#' \code{\link{fresnel}}
#'
#' @examples
#' # Real part of the refractive index of average seawater at 1 atmospheric 
#' # pressure and in the visible range:
#' n_water(lambda = 400:700, S = 34.72, Tc = 3.5, P = 0)
#'
#' @export

n_water <- function(lambda, Tc = 20, S = 0, P = 0, type = c('real', 'complex')) {

  if(any(lambda < 200) || any(lambda > 1100))
    warning("Requested wavelength outside model domain: [200,1100] (nm)")

  if(any(Tc < 0) || any(Tc > 30))
    warning("Requested temperature outside model domain: [0,30] (degree Celsius)")

  if(any(S < 0) || any(S > 40))
    warning("Requested salinity outside model domain: [0,40] (parts per thousand)")

  if(any(P < 0) || any(P > 1028))
    warning("Requested pressure outside model domain: [0,1028] (bar)")

  type <- match.arg(type)
 
  Tc2 <- Tc^2

  k0  <- 1.31405
  k1  <- 1.779e-4
  k2  <- -1.05e-6
  k3  <- 1.6e-8
  k4  <- -2.02e-6
  k5  <- 15.868
  k6  <- 0.01155
  k7  <- -0.00423
  k8  <- -4382
  k9  <- 1.1455e6

  n <- k0 + (k1 + k2 * Tc + k3 * Tc2) * S + k4 * Tc2 + (k5 + k6 * S + k7 * 
         Tc) / lambda + k8 / lambda^2 + k9 / lambda^3

  # Convert n_water relative to n_air to relative to vaccum:
  n <- n_air(lambda = lambda) * n

  if(any(P > 0)) {
    d0   <- d_water(Tc = Tc, S = S, P = 0) * 1E-3 # kg m^-3 to kg dm^-3
    d    <- d_water(Tc = Tc, S = S, P = P) * 1E-3 # kg m^-3 to kg dm^-3
    dfri <- liquid_DFRI(n, d0)
    n  <- sqrt(n^2 + dfri * (d - d0))
  }

  if(type == "complex") {
    a  <- a_water(lambda = lambda, Tc = Tc, S = S)
    ni <- lambda * 1E-9 * a / 4 / pi
    n  <- complex(real = n, imaginary = ni)
  }

  return(n)

}

#' Calculate the density fluctuation of the refractive index (DFRI) for liquids
#'
#' This function calculates the density fluctuation of the refractive index 
#' (DFRI) for liquids. This is the partial derivative of the squared refractive 
#' index with respect to density when keeping temperature (isothermic) or 
#' pressure (isobaric) constant.
#'
#' @param n Refractive index of the medium (unitless).
#' @param d Density of the medium (kg dm^-3).
#'
#' @details Between the different formulations for the DFRI, the formulation of 
#' Proutiere et al. (1992) show the best experimental comparison with observed 
#' data of scattering, when deriving scattering by pure water from density 
#' fluctuations (Zhang et al., 2009). Strictly, the formulation was derived for 
#' isothermal conditions, but theory supports that it should be the same as 
#' isobaric (at least approximately).
#'
#' @return A numeric vector with the partial derivative of the squared 
#' refractive index with respect to density (dn^2 / dd, dm^3 kg^-1), with 
#' density given in kg dm^-3.
#'
#' @references
#' Proutiere, A.; Megnassan, E.; Hucteau, H. 1992. Refractive index and density 
#' variations in pure liquids: A new theoretical relation. The Journal of 
#' Physical Chemistry 96, 3485-3489. DOI: 10.1021/j100187a058
#'
#' Zhang, X. Hu, L. 2009. Estimating scattering of pure water from density 
#' fluctuation of the refractive index. Optics Express 17, 3, 1671-1678. DOI: 
#' 10.1364/OE.17.001671
#'
#' @seealso \code{\link{n_water}}
#'
#' @examples
#' # Water entry on Table II of Proutiere et al. (1992):
#' liquid_DFRI(1.3330, 0.9982)
#'
#' @export

liquid_DFRI <- function(n, d) {

  if(.check_vect_args(arg_ls = list(n, d)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  n2     <- n^2

  dn2_dd <- (1 + ((n2 - 1) / 3 / n)^2 * 2 * (n2 + 2) / 3) * (n2 - 1) / d

  return(dn2_dd)

}

#' Dispersion formula for the refractive index of cellulose
#'
#' This function calculates the real part of the refractive index of cellulose 
#' relative to vaccum.
#'
#' @param lambda Wavelength in vacuum (nm).
#'
#' @return A numeric vector with the refractive index (unitless) of cellulose 
#' relative to vaccum.
#'
#' @references
#' Fournier, G. R.; Ardouin, J.-P.; Levesque, M. 2018. Modeling Sea Bottom 
#' Hyperspectral Reflectance. Applied Sciences 8, 12, 23 p. DOI: 
#' 10.3390/app8122680
#'
#' Sultanova, N.; Kasarova, S.; Nikolov, I. 2009. Dispersion properties of 
#' optical polymers. Acta Phys. Pol. A 116, 585–587. DOI: 
#'
#' @seealso \code{\link{n_air}}, \code{\link{n_water}}, \code{\link{n_calcite}}, 
#' \code{\link{n_quartz}}, \code{\link{snell}}, \code{\link{fresnel}}
#'
#' @examples
#' # Retrieve the real part of the refractive index of cellulose in the visible
#' # range:
#' n_cellulose(lambda = 400:700)
#'
#' @export

n_cellulose <- function(lambda) {

  if(any(lambda < 400) || any(lambda > 1060))
    warning("Requested wavelength outside model domain: [400,1060] (nm)")

  lambda <- lambda * 1E-3  # Conversion from nm to microns

  n <- sqrt(1 + 1.124 * lambda^2 / (lambda^2 - 0.011087))

  return(n)

}

#' Dispersion formula for the refractive index of calcite
#'
#' This function calculates the ordinary and extraordinary refractive index of 
#' calcite relative to vaccum. See details.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param comp   Logical. Should the ordinary and extraordinary indexes be 
#'               returned?
#'
#' @details Calcite is birefringent with two orientation dependent indexes of 
#' refraction, the extraordinary index of refraction (propagation along the 
#' optical axis) and the ordinary index of refraction (propagation orthogonal to 
#' optical axis). Values are provided at "room temperature".
#' 
#' If "comp" is FALSE (default), the average refractive index is returned. 
#' Otherwise a matrix is returned.
#'
#' @return A numeric vector with the refractive index (unitless) of calcite 
#' relative to vaccum.
#'
#' @references
#' Ballard, S. S.; Browder, J. S.; Ebersole, J. F. 1982. Index of Refraction. 
#' In: Gray D. E. Ed., American Institute of Physics Handbook, McGraw-Hill, New 
#' York.
#'
#' Fournier, G. R.; Ardouin, J.-P.; Levesque, M. 2018. Modeling Sea Bottom 
#' Hyperspectral Reflectance. Applied Sciences 8, 12, 23 p. DOI: 
#' 10.3390/app8122680
#'
#' Ghosh, G. 1999. Dispersion-equation coefficients for the refractive index and 
#' birefringence of calcite and quartz crystals. Opt. Commun. 163, 95–102.
#'
#' @seealso \code{\link{n_air}}, \code{\link{n_water}}, \code{\link{n_quartz}}, 
#' \code{\link{n_cellulose}}, \code{\link{snell}}, \code{\link{fresnel}}
#'
#' @examples
#' # Retrieve the average real part of the refractive index of calcite in the 
#' # visible range:
#' n_calcite(lambda = 400:700)
#'
#' @export

n_calcite <- function(lambda, comp = F) {

  if(any(lambda < 200) || any(lambda > 2170))
    warning("Requested wavelength outside model domain: [200,2170] (nm)")

  lambda <- lambda * 1E-3 # Conversion from nm to mum

  t1_o <- 1.73358749
  t2_o <- 0.96464345 * lambda^2 / (lambda^2 - 0.0194325203)
  t3_o <- 1.8283145 * lambda^2 / (lambda^2 - 120)
  n_o  <- sqrt(t1_o + t2_o + t3_o)

  t1_e <- 1.35859695
  t2_e <- 0.82427830 * lambda^2 / (lambda^2 - 0.0106689543)
  t3_e <- 0.14429128 * lambda^2 / (lambda^2 - 120)
  n_e  <- sqrt(t1_e + t2_e + t3_e)

  n <- cbind(n_o, n_e)

  if(comp) {
    return(n)
  } else { 
    return(apply(n, 1, mean))
  }

}

#' Dispersion formula for the refractive index of quartz
#'
#' This function calculates the ordinary and extraordinary refractive index of 
#' quartz relative to vaccum. See details.
#'
#' @param lambda Wavelength in vacuum (nm).
#' @param comp   Logical. Should the ordinary and extraordinary indexes be 
#'               returned?
#'
#' @details Quatrz is birefringent with two orientation dependent indexes of 
#' refraction, the extraordinary index of refraction (propagation along the 
#' optical axis) and the ordinary index of refraction (propagation orthogonal to 
#' optical axis). Values are provided at "room temperature".
#' 
#' If "comp" is FALSE (default), the average refractive index is returned. 
#' Otherwise a matrix is returned.
#'
#' @return A numeric vector with the refractive index (unitless) of quartz 
#' relative to vaccum.
#'
#' @references
#' Ballard, S. S.; Browder, J. S.; Ebersole, J. F. 1982. Index of Refraction. 
#' In: Gray D. E. Ed., American Institute of Physics Handbook, McGraw-Hill, New 
#' York.
#'
#' Fournier, G. R.; Ardouin, J.-P.; Levesque, M. 2018. Modeling Sea Bottom 
#' Hyperspectral Reflectance. Applied Sciences 8, 12, 23 p. DOI: 
#' 10.3390/app8122680
#'
#' Ghosh, G. 1999. Dispersion-equation coefficients for the refractive index and 
#' birefringence of calcite and quartz crystals. Opt. Commun. 163, 95–102.
#'
#' @seealso \code{\link{n_air}}, \code{\link{n_water}}, \code{\link{n_quartz}}, 
#' \code{\link{n_cellulose}}, \code{\link{snell}}, \code{\link{fresnel}}
#'
#' @examples
#' # Retrieve the average real part of the refractive index of quartz in the 
#' # visible range:
#' n_quartz(lambda = 400:700)
#'
#' @export

n_quartz <- function(lambda, comp = F) {

  if(any(lambda < 200) || any(lambda > 2050))
    warning("Requested wavelength outside model domain: [200,2050] (nm)")

  lambda <- lambda * 1E-3  # Conversion from nm to mum

  t1_o <- 1.28604141
  t2_o <- 1.07044083 * lambda^2 / (lambda^2 - 0.0100585997)
  t3_o <- 1.10202242 * lambda^2 / (lambda^2 - 100)
  n_o  <- sqrt(t1_o + t2_o + t3_o)

  t1_e <- 1.28851804
  t2_e <- 1.09509924 * lambda^2 / (lambda^2 - 0.0102101864)
  t3_e <- 1.15662475 * lambda^2 / (lambda^2 - 100)
  n_e  <- sqrt(t1_e + t2_e + t3_e)

  n <- cbind(n_o, n_e)
 
  if(comp) {
    return(n)
  } else {
    return(apply(n, 1, mean))
  }

}
