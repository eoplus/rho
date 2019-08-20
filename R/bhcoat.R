

#' Calculate the efficiency factors for coated spherical particles
#'
#' This function calculates the efficiency factors of extinction, absorption, 
#' scattering and backscattering for particle consisting of two concentric 
#' spheres.
#'
#' @param ri     Radius of the inner sphere.
#' @param ro     Radius of the outer sphere.
#' @param ni     Refractive index of the inner sphere relative to surrounding 
#'               medium.
#' @param no     Refractive index of the outer sphere relative to surrounding 
#'               medium.
#' @param lambda Wavelength in free space, in the same units as the radius.
#'
#' @detail This function is a translation of the BHCOAT FORTRAN subroutine 
#' published by Bohren & Huffman (1983). The function should not be used for 
#' large, highly absorbing spheres and will issue a warning if input parameters
#' describe this condition.
#' 
#' Note that all parameters must be scalars, while \code{ni} and \code{no} can 
#' be complex numbers. The refractive indexes should be relative to the 
#' surrounding medium (water) and \code{lambda} be in the same units as the 
#' radius of the concentric spheres.
#'
#' @return A numerical vector with the efficiency factors for extinction, 
#' absorption, scattering and backscattering.
#'
#' @references
#' Bohren, C. F.; Huffman, D. R. 1983. Absorption and Scattering of Light by 
#' Small Particles. John Wiley & Sons, New York, 533 pp.
#'
#' @examples
#' ni  <-  complex(real = 1.02, imaginary = 0.006)
#' no  <-  complex(real = 1.10, imaginary = 0.0)
#' bhcoat(ri = 5, ro = 5.2, ni = ni, no = no, lambda = 0.675)
#'
#' 

bhcoat <- function(ri, ro, ni, no, lambda) {

  if(any(lengths(list(ri, ro, ni, no, lambda)) > 1))
    stop("All arguments must be scalars")

  tol <- 1E-8 # Inner sphere convergency criterion

  x <- 2 * pi * ri * Re(nm) / lambda
  y <- 2 * pi * ro * Re(nm) / lambda

  x1 <- ni * x
  x2 <- no * x
  y2 <- no * y

  if(any(Im(c(x1, x2, y2)) > 30))
    warning("This function should not be used for large, highly absorbing",  
            " spheres")

  maxn <- y + 4 * y^(1/3) + 2 # Maximum number of terms
  noi  <- no / ni

  d0x1   <-  cos(x1) / sin(x1)
  d0x2   <-  cos(x2) / sin(x2)
  d0y2   <-  cos(y2) / sin(y2)
  psi0y  <-  cos(y)
  psi1y  <-  sin(y)
  chi0y  <- -sin(y)
  chi1y  <-  cos(y)
  xi0y   <-  complex(real = psi0y, imaginary = -chi0y)
  xi1y   <-  complex(real = psi1y, imaginary = -chi1y)
  chi0y2 <- -sin(y2)
  chi1y2 <-  cos(y2)
  chi0x2 <- -sin(x2)
  chi1x2 <-  cos(x2)

  QSCA   <- 0.0
  QEXT   <- 0.0
  XBACK  <- complex(real = 0.0, imaginary = 0.0)
  iflag  <- 0

  for(n in 1:maxn) {
    psiy <- (2 * n - 1) * psi1y / y - psi0y
    chiy <- (2 * n - 1) * chi1y / y - chi0y
    xiy  <- complex(real = psiy, imaginary = -chiy)
    d1y2 <- 1 / (n / y2 - d0y2) - n / y2

    if(iflag == 0) {
      d1x1   <- 1 / (n / x1 - d0x1) - n / x1
      d1x2   <- 1 / (n / x2 - d0x2) - n / x2
      chix2  <- (2 * n - 1) * chi1x2 / x2 - chi0x2
      chiy2  <- (2 * n - 1) * chi1y2 / y2 - chi0y2
      chipx2 <- chi1x2 - n * chix2 / x2
      chipy2 <- chi1y2 - n * chiy2 / y2
      ancap  <- noi * d1x1 - d1x2
      ancap  <- ancap / (noi * d1x1 * chix2 - chipx2)
      ancap  <- ancap / (chix2 * d1x2 - chipx2)
      brack  <- ancap * (chiy2 * d1y2 - chipy2)
      bncap  <- noi * d1x2 - d1x1
      bncap  <- bncap / (noi * chipx2 - d1x1 * chix2)
      bncap  <- bncap / (chix2 * d1x2 - chipx2)
      crack  <- bncap * (chiy2 * d1y2 - chipy2)

      amess1 <- brack * chipy2
      amess2 <- brack * chiy2
      amess3 <- crack * chipy2
      amess4 <- crack * chiy2
    }

    if(abs(amess1) < tol*abs(d1y2) && abs(amess2) < tol && abs(amess3) < tol*abs(d1y2) && abs(amess4) < tol) {
      brack <- complex(real = 0., imaginary = 0.)
      crack <- complex(real = 0., imaginary = 0.)
      iflag <- 1
   } else {
      iflag <- 0
   }

   dnbar <- d1y2 - brack * chipy2
   dnbar <- dnbar / (1 - brack * chiy2)
   gnbar <- d1y2 - crack * chipy2
   gnbar <- gnbar / (1 - crack * chiy2)

   a_n <- (dnbar / no + n / y) * psiy - psi1y
   a_n <- a_n / ((dnbar / no + n / y) * xiy - xi1y)
   b_n <- (no * gnbar + n / y) * psiy - psi1y
   b_n <- b_n / ((no * gnbar + n / y) * xiy - xi1y)

   # calculate sums for qsca, qext, xback
   QSCA  <- QSCA  + (2 * n + 1) * (abs(a_n)^2 + abs(b_n)^2)
   QEXT  <- QEXT  + (2 * n + 1) * (Re(a_n) + Re(b_n))
   XBACK <- XBACK + (2 * n + 1) * (-1)^n * (a_n - b_n)

   psi0y  <- psi1y
   psi1y  <- psiy
   chi0y  <- chi1y
   chi1y  <- chiy
   xi1y   <- complex(real = psi1y,  imaginary = -chi1y)
   chi0x2 <- chi1x2
   chi1x2 <- chix2
   chi0y2 <- chi1y2
   chi1y2 <- chiy2
   d0x1   <- d1x1
   d0x2   <- d1x2
   d0y2   <- d1y2
  }

  QSCA  <- (2 / y^2) * QSCA
  QEXT  <- (2 / y^2) * QEXT
  QBACK <- (1/ y^2) * XBACK * Conj(XBACK)
  QABS  <- QEXT - QSCA

  Q <- c(QEXT, QABS, QSCA, QBACK)
  names(Q) <- c("Qext", "Qabs", "Qsca", "Qback")
  return(Q)
}

