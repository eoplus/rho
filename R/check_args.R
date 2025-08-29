
# Will return TRUE if the condition is NOT met. Condition is that the lengths 
# of all vector arguments must be an integer multiple of the the length of the
# longer vector.

.check_vect_args <- function(arg_ls) {
  arg_lngth <- lengths(arg_ls)
  arg_tab   <- table(arg_lngth)
  n_lngths  <- length(arg_tab)
  x_lnghts  <- as.numeric(names(arg_tab))
  any((x_lnghts[length(x_lnghts)] %% x_lnghts) != 0)
}
