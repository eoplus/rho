#' Absorption coefficient of liquid water 
#'
#' A harmonized dataset containing multi-source data of absorption coefficient
#' of liquid water at 0 salinity and 20 degreesC, from 300 nm to 4000 nm in 
#' steps of 2 nm. This dataset combines the measurements of several sources 
#' listed below, harmonized to a standard temperature by using the absorption 
#' derivative with respect to temperature (Röttgers, 2010). The derivative with 
#' respect to temperature and salinity (Röttgers, 2010) are also provided. The 
#' dataset is taken from the absorption coefficient Version 3 of the Water 
#' Optical Properties Processor (WOPP) and fully documented in Röttgers et al. 
#' (2011).
#'
#' This is the base data used by functions \code{n_water} to retrieve the 
#' imaginary part of the refractive index of water and by function 
#' \code{a_water} to retrieve the temperature and salinity dependent absorption
#' of pure (saline) water.
#'
#' The sources used in the compilation were:
#' \itemize{
#'   \item Lu (2006);
#'   \item Mason et al. (2016; < 550 nm);
#'   \item Max & Chapados (2009; 1850-2000 nm, > 2500 nm);
#'   \item Pope & Fry (1997);
#'   \item Segelstein (1981; > 700 nm);
#'   \item Wang (2008);
#'   \item Wieliczka et al. (1989).
#' }
#'
#' @format A matrix with 1851 rows and 7 variables:
#' \describe{
#'   \item{wavelength}{wavelength, in nm}
#'   \item{a}{absorption coefficient of water, 1/m}
#'   \item{da/dS}{absorption derivative with respect to salinity, 1/m 1/PSU}
#'   \item{da/dT}{absorption derivative with respect to temperature, 1/m 1/degreeC}
#'   \item{sigma_a}{standard deviation of a}
#'   \item{sigma_da/dS}{standard deviation of da/dS}
#'   \item{sigma_da/dT}{standard deviation of da/dT}
#' }
#'
#' @source WOPP - Water Optical Properties Processor, Version 1.7, 2016.
#'
#' @references
#' Lu, Z. 2006. Optical absorption of pure water in the blue and ultraviolet. 
#' Ph. D. thesis. Texas A&M University.
#' 
#' Mason, J. D.; Cone, M. T.; Fry, E. S. 2016. Ultraviolet (250-550 nm) 
#' absorption spectrum of pure water. Applied Optics 55, 25, 7163-7172. DOI: 
#' 10.1364/AO.55.007163
#'
#' Max, J. J.; Chapados, C. 2009. Isotope effects in liquid water by infrared 
#' spectroscopy. III. H2O and D2O spectra from 6000 to 0 cm-1. Journal of 
#' Chemical Physics 131, 184505. DOI: 10.1063/1.3258646
#'
#' Pope, R. M.; Fry, E. S. 1997. Absorption spectrum (380-700 nm) of pure water: 
#' II. Integrating cavity measurements. Applied Optics 36, 33, 8710-8723. DOI:
#' 10.1364/AO.36.008710
#'
#' Röttgers, R. 2010. Measurements of inherent optical properties of pure water. 
#' Technical note, STSE-WaterRadiance D5, ESA.
#'
#' Röttgers, R.; Doerffer, R.; McKee, D.; Schönfeld, W. 2011. The Water Optical 
#' Properties Processor (WOPP). Pure water spectral absorption, scattering, and 
#' real part of refractive index model. Algorithm Theoretical Basis Document.
#' 
#' Segelstein, D. J. 1975. The complex refractive index of water. Master Thesis,
#' University of Missouri-Kansas City, 175 pp.
#'
#' Wang, L. 2008. Measuring optical absorption coefficient of pure water in UV 
#' using the integrating cavity absorption meter. Ph. D. thesis. Texas A&M 
#' University.
#'
#' Wieliczka, D. M.; Weng, S.; Querry, M. R. 1989. Wedge shaped cell for highly 
#' absorbent liquids: infrared optical constants of water. Applied Optics 28, 9, 
#' 1714-1719. DOI: 10.1364/AO.28.001714
#'
"a_water_wopp"

#' Specific absorption coefficient of cellulose
#'
#' The absorptivity of cellulose was inverted by Jacquemoud et al. (1996) from
#' a large database of leaf reflectance and transmittance using a two stream 
#' model. The retrieved absorptivity was then normalized by the cellulose + 
#' lignin dry mass per leaf area (g / m^2), resulting in the mass specific 
#' absorption coefficient (m^2 / g). Since the results of Jacquemound et al. 
#' (1996) were retrieved with a model based on the Kubelka-Munk theory, the 
#' values are scaled by 2 to be used with radiative transfer theory (Ganapol et 
#' al. 1999).
#'
#' A similar spectral shape was determined by Dawson et al. (1998), but requires 
#' scaling and offset to match the values of Jacquemoud et al. (1996) as used by
#' Fournier et al. (2018) to invert macrophyte properties.
#'
#' @source PROSPEC source code, Version 2.01, 1995.
#'
#' @references
#' Fournier, G. R.; Ardouin, J.-P.; Levesque, M. Modeling Sea Bottom 
#' Hyperspectral Reflectance. Applied Sciences 8,12, 2680. DOI: 
#' 10.3390/app8122680
#'
#' Ganapol, B. D.; Johnson, L. F.; Hlavka, C. A.; Peterson, D. L.; Bond, B. 
#' 1999. LCM2: A coupled leaf/canopy radiative transfer model. Remote Sensing of 
#' Environment 70, 2, 153-166. DOI: 10.1016/S0034-4257(99)00030-9
#'
#' Jacquemoud, S.; Ustin, S. L.; Verdebout, J.; Schmuck, G.; Andreoli, G.; 
#' Hosgood, B.. 1996. Estimating leaf biochemistry using the PROSPECT leaf 
#' optical properties model. Remote Sensing of Environment 56, 3, 194-202. DOI: 
#' 10.1016/0034-4257(95)00238-3 
#'
"a_cellulose"
