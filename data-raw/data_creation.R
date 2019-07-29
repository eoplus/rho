
devtools::load_all()

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
