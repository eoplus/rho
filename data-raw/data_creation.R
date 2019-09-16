
devtools::load_all()
library(devtools)
#
# Pure (saline) water absorption
#
# From Version 3 of the WOPP:

a_water_wopp <- read.table('data-raw/purewater_abs_coefficients_v3.dat', 
                           skip = 4)
colnames(a_water_wopp) <- c("wavelength", "a", "da/dS", "da/dT", "sigma_a", 
                            "sigma_da/dS", "sigma_da/dt"
devtools::use_data(a_water_wopp, overwrite = TRUE)

#
# Specific absorption coefficient of cellulose + lignin
#
# Data from Jacquemoud et al. 1996, scaled by 2 as required for radiative 
# transfer when compared to Kundela-Munk theory (Ganapol et al. 1999).

a_cellulose <- read.table('data-raw/a_cellulose_jacqemoud96.dat', header = T)
colnames(a_cellulose) <- c("wavelength",  "sa")
a_cellulose[, 2] <- a_cellulose[, 2] * 2
devtools::use_data(a_cellulose, overwrite = TRUE)

#
# Spectral response functions:
#

path <- '/media/alexandre/Alex_data/datasets/spetral_response_functions/distributed_RData'
fls  <- list.files(path, '.RData', full.names = T)

load(fls[grep('l8_iob', fls)])
load(fls[grep('s2a_ib', fls)])
load(fls[grep('s2b_ib', fls)])
load(fls[grep('phra_pld', fls)])
load(fls[grep('phrb_pld', fls)])
load(fls[grep('jss_re', fls)])
load(fls[grep('0c0d_ps', fls)])
load(fls[grep('0e_ps', fls)])
load(fls[grep('0f10_ps', fls)])

rsrf <- list()
rsrf[['l8']] <- list()
rsrf[['l8']][['oli']]   <- srf_oli_l8_iob_mm
rsrf[['s2']] <- list()
rsrf[['s2']][['msia']]  <- srf_msi_s2a_ib_mm
rsrf[['s2']][['msib']]  <- srf_msi_s2b_ib_mm
rsrf[['pld']] <- list()
rsrf[['pld']][['phra']] <- srf_phra_pld_ib_mm
rsrf[['pld']][['phrb']] <- srf_phrb_pld_ib_mm
rsrf[['re']] <- list()
rsrf[['re']][['jss']]   <- srf_jss_re_iob_mm
rsrf[['ps']] <- list()
rsrf[['ps']][['0c0d']]  <- srf_0c0d_ps_iob_mm
rsrf[['ps']][['0e']]    <- srf_0e_ps_iob_mm
rsrf[['ps']][['0f10']]  <- srf_0f10_ps_iob_mm
use_data(rsrf, overwrite = TRUE)

#
# Water NIR similarity spectrum
#

water_nir_sim <- read.csv('ruddick_nir_similarity_spec.csv')
devtools::use_data(water_nir_sim, overwrite = TRUE)
