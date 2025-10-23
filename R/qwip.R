
#' Calculate the QWIP score 
#'
#' The function calculates the Quality Water Index Polynomial (QWIP) for hyperspectral 
#' and selected multispectral sensors.
#'
#' @param wave   Reference centre wavelength (nm).
#' @param refl   Reflectance data (1/sr or unitless).
#' @param sensor One of: 'hyperspectral' or 'aeronet-oc'. See details.
#'
#' @details The QWIP score the difference between the observed reflectance Normalized Difference 
#' Index at 492 nm and 667 nm and the predicted NDI based on an empirical relation against the 
#' Apparent Visible Wavelength (AVW).
#'
#' For multispectral sensors, the calculated AVW can be concerted to the hyperspectral
#' equivalent AVW by using a regression function for a given multispectral sensor waveband 
#' set using simulated hyperspectral data. Currently the only multispectral sensor supported
#' is the AERONET-OC. The base set of wavebands should be used as input (412, 443, 488, 547, 667).
#'
#' @references
#' Dierssen, H.; Vandermeulen, R.; Barnes, B.; Castagna, A.; Knaeps, E.; Vanhellemont, Q. 2022. 
#' QWIP: A Quantitative Metric for Quality Control of Aquatic Reflectance Spectral Shape Using the 
#' Apparent Visible Wavelength. Frontiers in Remote Sensing 3. DOI: 
#' https://doi.org/10.3389/frsen.2022.869611
#'
#' @export

qwip_score  <- function(wave, refl, sensor) {

    avw     <- avw(wave, refl)
    ndi     <- .qwip_ndi(wave, refl)
    qwip    <- -8.399885E-9 * avw^4 + 1.715532E-5 * avw^3 - 1.301670E-2 * avw^2 + 4.357838 * avw - 5.449532E2
    score   <- ndi - qwip
    return(score)

}

#' Calculate the QWIP quality control flag 
#'
#' The function performs the quality control of reflectance spectra using the QWIP approach.
#'
#' @param wave      Reference centre wavelength (nm).
#' @param refl      Reflectance data (1/sr or unitless).
#' @param sensor    One of: 'hyperspectral' or 'aeronet-oc'. See details.
#' @param score_lim A magnitude of maximum QWIP score to accept in the quality control.
#'
#' @details The QWIP quality control calculates the QWIP score and evaluate if it is within the user
#' defined limits. It also evaluates if the AVW occurs is the expected range for the classified 
#' water type using a simple three types classification.
#'
#' For multispectral sensors, the calculated AVW can be concerted to the hyperspectral
#' equivalent AVW by using a regression function for a given multispectral sensor waveband 
#' set using simulated hyperspectral data. Currently the only multispectral sensor supported
#' is the AERONET-OC. The base set of wavebands should be used as input (412, 443, 488, 547, 667).
#'
#' @references
#' Dierssen, H.; Vandermeulen, R.; Barnes, B.; Castagna, A.; Knaeps, E.; Vanhellemont, Q. 2022. 
#' QWIP: A Quantitative Metric for Quality Control of Aquatic Reflectance Spectral Shape Using the 
#' Apparent Visible Wavelength. Frontiers in Remote Sensing 3. DOI: 
#' https://doi.org/10.3389/frsen.2022.869611
#'
#' @export

qwip_qc <- function(wave, refl, sensor, score_lim = 0.2) {

    avw  <- avw(wave, refl, sensor)
    owt  <- .qwip_owt(wave, refl)
    qwip <- qwip_score(wave, refl, sensor)
    
    if(is.na(owt)) return(NA)
    
    if(owt == "owt_1") {
        avw_owt_range <- avw >= 440 & avw <= 530
    } else if(owt == "owt_2") {
        avw_owt_range <- avw >= 510 & avw <= 580
    } else {
        avw_owt_range <- avw >= 550 & avw <= 590
    }
    
    qwip_score_range <- abs(qwip) <= score_lim
    qwip_qc <- avw_owt_range & qwip_score_range
    return(qwip_qc)
    
}

#' Dont export
.qwip_ndi  <- function(wave, refl) {

    id_blue <- which.min(abs(wave - 492))
    id_red  <- which.min(abs(wave - 665))
    ndi     <- (refl[id_red] - refl[id_blue]) / (refl[id_red] + refl[id_blue])
    return(ndi)
    
}

#' Dont export
.qwip_owt <- function(wave, refl) {

    id_blue  <- which.min(abs(wave - 492))
    id_green <- which.min(abs(wave - 560))
    id_red   <- which.min(abs(wave - 665))

    if(any(is.na(refl[c(id_blue, id_green, id_red)]))) return(NA)

    if((refl[id_red] > refl[id_green]) | (refl[id_red] > 0.025)) {
        owt <- "owt_3"
    } else if (refl[id_green] < refl[id_blue]) {
        owt <- "owt_1"
    } else {
        owt <- "owt_2"
    }
    return(owt)

}

