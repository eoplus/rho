
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


