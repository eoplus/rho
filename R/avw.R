
#' Calculate the Apparent Visible Wavelength (AVW) 
#'
#' The function calculates the Apparent Visible Wavelength (AVW) for hyperspectral 
#' and selected multispectral sensors.  
#'
#' @param wave   Reference centre wavelength (nm).
#' @param refl   Reflectance data (1/sr or unitless).
#' @param sensor One of: 'hyperspectral' or 'aeronet-oc'. See details.
#'
#' @details The Apparent Visible Wavelength (AVW) is the weighted harmonic mean
#' wavelength, with the weights given by the reflectance for all wavelengths in
#' the visible range (400 nm to 700 nm).
#'
#' For multispectral sensors, the calculated AVW can be concerted to the hyperspectral
#' equivalent AVW by using a regression function for a given multispectral sensor waveband 
#' set using simulated hyperspectral data. Currently the only multispectral sensor supported
#' is the AERONET-OC. The base set of wavebands should be used as input (412, 443, 488, 547, 667).
#'
#' @references
#' Vandermeulen, R.; Mannino, A.; Craig, S.; Werdell, P.G. 2020. 150 shades of green: Using the full 
#' spectrum of remote sensing reflectance to elucidate color shifts in the ocean. Remote Sensing of 
#' Environment 247, 111900. DOI: https://doi.org/10.1016/j.rse.2020.111900
#'
#' @export

avw  <- function(wave, refl, sensor) {

    if(sensor == "aeronet-oc") {
   
        # minimal set for AERONET-OC
        wave_aoc <- c(412, 443, 488, 547, 667)
        if(length(wave) > 5) stop("Data should be provided at the minimum AERONET-OC waveband set: 412, 443, 488, 547, 667")
        if(any(abs(wave - wave_oc) > 5)) stop("Data should be provided at wavelengths closest to the nominal 412, 443, 488, 547, 667")
        
        if(any(refl < 0)) warning("Negative reflectances in the selected range")
        avw_m <- sum(refl, na.rm = T) / sum(refl/wave, na.rm = T)
        
        # Multispectral coefficients calculated as in Vandermeulen et al. 2020:
        coef     <- c(
            "c0" = -2.409993e-09,
            "c1" = 5.802722e-06,
            "c2" = -5.575237e-03,
            "c3" = 2.669700e+00,
            "c4" = -6.353280e+02,
            "c5" = 6.031505e+04
        )
        avw <- coef["c0"] * avw_m^5 +
               coef["c1"] * avw_m^4 + 
               coef["c2"] * avw_m^3 + 
               coef["c3"] * avw_m^2 + 
               coef["c4"] * avw_m^1 + 
               coef["c5"]
        return(avw)
        
    } else if(sensor == "hyperspectral") {
        id_vis <- which(wave > 400 & wave < 700)
        if(any(refl[id_vis] < 0)) warning("Negative reflectances in the selected range")
        avw    <- sum(refl[id_vis]) / sum(refl[id_vis]/wave[id_vis])
        return(avw)
    } else {
        stop("sensor must be 'aeronet-oc' or 'hyperspectral'")
    }

}

