
#' Density of water
#'
#' This function calculates the water density as a function of pressure, 
#' temprature and salinity.
#'
#' @param Tc     Temperature (ºC).
#' @param S      Salinity (parts per thousand).
#' @param P      Pressure (bar in excess of 1 ATM).
#'
#' @details The density of pure (saline) water at diferent pressures is based on 
#' the UNESCO (1981) report, as recomended by Zhang et al. (2009) and Libes 
#' (2009). The UNESCO equation has a relative error of 0.0004% (Zhang et al.,
#' 2009). The actual function was implemented with corrections by from 
#' Forfonoff (1985), as recomended by Libes (2009).
#'
#' @references
#' Fofonoff, N. P. 1985. Physical properties of seawater: A new salinity scale 
#' and equation of state for seawater. Journal of Geophysical Research: Oceans 
#' 90, C2, 3332-3342. DOI: 10.1029/JC090iC02p03332
#'
#' Libes, S. 2009. Introduction to Marine Biogeochemistry. 2nd Edition. Academic 
#' Press, 928 pp.
#'
#' Zhang, X. Hu, L. 2009. Estimating scattering of pure water from density 
#' fluctuation of the refractive index. Optics Express 17, 3, 1671-1678. DOI: 
#' 10.1364/OE.17.001671
#'
#' @examples
#' # Average seawater density at surface under 1 ATM:
#' d_water(Tc = 3.5, S = 34.72, P = 0)
#' 
#' @export

d_water <- function(Tc = 20, S = 0, P = 0) {

  if(.check_vect_args(arg_ls = list(Tc, S, P)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  T2 <- Tc^2
  T3 <- Tc^3
  T4 <- Tc^4
  T5 <- Tc^5

  b0 <-  999.842594
  b1 <-  6.793952e-2
  b2 <- -9.09529e-3
  b3 <-  1.001685e-4
  b4 <- -1.120083e-6
  b5 <-  6.536332e-9
  c0 <-  8.24493e-1
  c1 <- -4.0899e-3
  c2 <-  7.6438e-5
  c3 <- -8.2467e-7
  c4 <-  5.3875e-9
  c5 <- -5.72466e-3
  c6 <-  1.0227e-4
  c7 <- -1.6546e-6
  c8 <-  4.8314e-4

  d <- (b0+(b1*Tc)+(b2*T2)+(b3*T3)+(b4*T4)+(b5*T5))+
       (S*(c0+(c1*Tc)+(c2*T2)+(c3*T3)+(c4*T4)))+((c5+(c6*Tc)+(c7*T2))*(S^1.5))+
       (c8*(S^2))
  d <- d / (1 - P / .K_w(Tc = Tc, S = S, P = P))

  return(d)

}

#' Secant bulk modus of pure seawater:

.K_w <- function(Tc, S, P) {

  T2 <- Tc^2
  T3 <- Tc^3
  T4 <- Tc^4

  e0 <-  19652.21
  e1 <-  148.4206
  e2 <- -2.327105
  e3 <-  1.360477e-2
  e4 <- -5.155288e-5
  f0 <-  54.6746
  f1 <- -0.603459
  f2 <-  1.09987e-2
  f3 <- -6.1670e-5
  g0 <-  7.944e-2
  g1 <-  1.6483e-2
  g2 <- -5.3009e-4
  h0 <-  3.239
  h1 <-  1.43716e-3
  h2 <-  1.16092e-4 
  h3 <- -5.77905e-7
  i0 <-  2.2838e-3 
  i1 <- -1.0981e-5
  i2 <- -1.6078e-6 
  j0 <-  1.91075e-4 
  m0 <-  8.50935e-5
  m1 <- -6.12293e-6 
  m2 <-  5.2787e-8
  n0 <- -9.9348e-7
  n1 <-  2.0816e-8
  n2 <-  9.1697e-10

  kw <- (e0+(e1*Tc)+(e2*T2)+(e3*T3)+(e4*T4))+((f0+(f1*Tc)+(f2*T2)+
        (f3*T3))*S)+((g0+(g1*Tc)+(g2*T2))*(S^1.5))
  A   <- ((h0+(h1*Tc)+(h2*T2)+(h3*T3))+((i0+(i1*Tc)+(i2*T2))*S)+
         (j0*(S^1.5)))*P
  B   <- ((m0+(m1*Tc)+(m2*T2))+(n0+(n1*Tc)+(n2*T2)))*(P^2)

  kw <- kw + A + B

  return(kw)

}


