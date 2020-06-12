#' Calculate the refraction angle 
#'
#' The function calculates the (complex) angle of refraction for eletromagnetic 
#' radiation transmitted to a medium.
#'
#' @param thetai Incidence angle relative to the normal of the interface between
#'               the two mediums (radians).
#' @param ni     Refractive index of the medium of incidence (unitless).
#' @param nr     Refractive index of the medium of refraction (unitless).
#'
#' @details The angle of incidence is limited to be between 0 and pi/2. Although
#' the sine function is periodical and the modulus will be correct for multiples
#' of a given angle, this serves as a protection in case incidence angles > pi/2
#' are inadvertently provided. \code{NA}s are allowed in all arguments.
#'
#' The indexes of refraction may be complex numbers. If so, the function will 
#' return the complex refraction angle. This complex refractive angle should be 
#' passed to the \code{fresnel} function to calculate the reflectance at the 
#' interface with an absorptive medium. Complex refractive angles can be 
#' converted to angles of constant phase or constant amplitude with the 
#' appropriate function (\code{snell_decomp}). Note that the refractive indexes 
#' are normalized to each other internally, so if passing refractive indexes 
#' that are already relative, one of the refractive indexes should be equal 
#' unity. See examples. 
#'
#' If one of the mediums have a negative refractive index, the refracted ray 
#' will occur on the same side of the normal as the incident ray, i.e., the 
#' azimuth is rotated 180&deg; to that of the incident ray. Since this is a 
#' special condition not occurring for hydrological optics, the functions are 
#' simplified by requiring that the refractive indexes be positive.
#'
#' For real indexes of refraction, in case of incidence angle is higher than the 
#' critical angle for total internal reflection (thetai > asin(nr / ni)), the 
#' function will return a real refraction angle of pi/2, i.e., absence of 
#' refraction. 
#'
#' @return A numeric or complex vector.
#'
#' @references
#' Bohren, C. F.; Huffman, D. R. 1983. Absorption and Scattering of Light by 
#' Small Particles. Wiley, New York.
#'
#' @seealso \code{\link{snell_decomp}}, \code{\link{fresnel}}
#'
#' @examples
#' # Refraction from standard dry air to average seawater at 550 nm:
#' snell(c(0, pi/4, pi/2), ni = 1.000278, nr = 1.342033)
#'
#' # Same as before but seawater index already relative to air:
#' snell(c(0, pi/4, pi/2), ni = 1, nr = 1.342033 / 1.000278)
#'
#' # Same as before, but reversed direction:
#' snell(c(0, pi/4, pi/2), ni = 1.342033 / 1.000278, nr = 1)
#'
#' # Same as before, but with the complex index of refraction of water:
#' snell(c(0, pi/4, pi/2), ni = complex(real = 1.342033, imaginary = 2.442E-09) 
#'       / 1.000278, nr = 1)
#' @export

snell <- function(thetai, ni, nr) {

  if(.check_vect_args(arg_ls = list(thetai, ni, nr)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  if(any(thetai < 0) || any(thetai > pi / 2))
    stop("Incidence angle must be between 0 and pi/2, included")
	
  if(any(Re(ni) < 0) || any(Re(nr) < 0))
    stop("Handling of negative refractive indexes is not implemented")

#  When with just one angle and several refractive indexes, it can break because
#  the initial angle will not be replicated and a vectors of NA will result...
#  if(!is.complex(ni) && !is.complex(nr)) {
#    ca <- suppressWarnings(asin(nr / ni))
#    thetai[thetai > ca] <- NaN
#  }

  thetar <- suppressWarnings(asin(sin(thetai) * ni / nr))
  thetar[is.nan(thetar)] <- pi / 2
  return(thetar)

}

#' Decompose a complex refracted angle 
#'
#' The function calculates the angle of constant phase and the angle of constant
#' amplitude for complex refracted angles of eletromagnetic radiation 
#' transmitted to a medium.
#'
#' @param thetar Complex refracted angle (radians).
#' @param nr     Refractive index of the medium of refraction (unitless).
#'
#' @details Refractive index of the medium of refraction may be real or complex.
#' The refracted angle is the angle at constant phase.
#'
#' @return A matrix with the constant-phase (theta_kp) and constant-amplitude
#' (theta_ka) angles.
#'
#' @references
#' Liu, Y.; Qian, J.; Tian, Y. 2003. Succinct formulas for decomposition of 
#' complex refraction angle. IEEE Antennas and Propagation Society International 
#' Symposium. Digest. Held in conjunction with: USNC/CNC/URSI North American 
#' Radio Sci. Meeting (Cat. No.03CH37450), Columbus, OH, pp. 487-490, vol. 3. 
#' DOI: 10.1109/APS.2003.1219892
#'
#' @seealso \code{\link{snell}}, \code{\link{fresnel}}
#'
#' @examples
#' # Refraction from average seawater to standard dry air at 550 nm:
#' refr <- snell(c(0, pi/4, pi/2), ni = complex(real = 1.342033, 
#'               imaginary = 2.442E-09), nr = 1.000278)
#' snell_decomp(thetar = refr, nr = 1)
#' 
#' @export

snell_decomp <- function(thetar, nr) {

  if(!is.complex(thetar))
    stop("Refracted angle must be complex")

  if(.check_vect_args(arg_ls = list(thetar, nr)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  theta_kp <- atan(Re(nr * sin(thetar)) / Re(nr * cos(thetar)))
  theta_ka <- atan(Im(nr * sin(thetar)) / Im(nr * cos(thetar)))

  theta_ka[theta_ka < 1E-15] <- 0  # Account for rounding errors...
  theta_ka[is.nan(theta_ka)] <- 0

  return(cbind(theta_kp, theta_ka))

}

#' Calculate the reflectance at a flat interface
#'
#' The function calculates the reflectance, i.e., the squared modulus of the 
#' (complex) amplitude reflection coefficient, of (polarized) electromagnetic 
#' radiation incident at a flat interface.
#'
#' @param thetai Incidence angle relative to the normal of the interface between
#'               the two mediums (radians).
#' @param thetar Refracted angle relative to the normal of the interface between
#'               the two mediums (radians). If missing, function \code{snell} 
#'               will be called internally.
#' @param ni     Refractive index of the medium of incidence (unitless).
#' @param nr     Refractive index of the medium of refraction (unitless).
#' @param fp     Fraction of parallel polarization (unitless).
#' @param ave    Logical. Should the reflections for S and P polarizations be 
#'               returned separately?
#'
#' @details Uses Fresnel's formulations to calculate the reflectance based on 
#' the incident and (complex) refracted angles and the (complex) indexes of 
#' refraction. \code{NA}s allowed in all arguments.
#'
#' If \code{ave} is set to TRUE, the weighted average reflectance will be 
#' calculated, with the weights of each polarization given by \code{fp}.
#'
#' @return A numeric vector if average reflectance is desired or a matrix with
#' the reflectance for parallel (R_p) and perpendicular polarizations (R_s).
#'
#' @references
#' Bohren, C. F.; Huffman, D. R. 1983. Absorption and Scattering of Light by 
#' Small Particles. Wiley, New York.
#'
#' @seealso \code{\link{snell}}, \code{\link{snell_decomp}}
#'
#' @examples
#' # Reflection from standard dry air to average seawater at 550 nm:
#' fresnel(thetai = c(0, pi/4, pi/2), ni = 1.000278, nr = 
#' complex(real = 1.342033, imaginary = 2.442E-09), ave = TRUE)
#' 
#' @export

fresnel <- function(thetai, thetar, ni, nr, fp = 0.5, ave = FALSE) {

  if(missing(thetar))
    thetar <- snell(thetai = thetai, ni = ni, nr = nr)

  if(.check_vect_args(arg_ls = list(thetai, thetar, ni, nr)))
    stop("For vector applications, the length of all vectors must be integer ", 
         "multiples of the longer vector")

  if(any(thetai < 0) || any(thetai > pi / 2))
    stop("Incidence angle must be between 0 and pi/2, included")

  refr <- snell_decomp(thetar = as.complex(thetar), nr = nr)[, 1]
  if(any(refr < 0) || any(refr > pi / 2))
    stop("Refracted angle must be between 0 and pi/2, included")

  if(any(fp < 0) || any(fp > 1))
    stop("Fraction of parallel polarization radiation must be between 0 and 1")

  R_p <- .fresnel_p(thetai = thetai, thetar = thetar, ni = ni, nr = nr)
  R_s <- .fresnel_p(thetai = thetai, thetar = thetar, ni = nr, nr = ni)

  if(ave) { 
    return(R_p * fp + R_s * (1 - fp))
  } else {
    return(cbind(R_p, R_s)) 
  }
}

.fresnel_p <- function(thetai, thetar, ni, nr) {
  r <- (nr * cos(thetar) - ni * cos(thetai)) / 
       (nr * cos(thetar) + ni * cos(thetai))
  r[is.nan(r)] <- 1
  return(abs(r)^2)
}


