## sanba 0.0.4

* Changed how elapsed times are reported. Also reporting total runtime with multiple runs in VI framework.
* Added functions to compute and display group-specific posterior densities for MCMC estimation methods.
* Fixed some minor bugs and typos.
* Minor edits in the titles of the produced tables.

## R CMD check results

- Running `devtools::check(args = c('--as-cran', '--no-manual'))` and `devtools::check(remote = TRUE, manual = TRUE)` 
  locally, produces no errors, no warnings, and no notes.

- Running `devtools::check_win_devel()`, `devtools::check_win_release()`, and `devtools::check_mac_release()` 
  yields no errors, no warnings, and no notes.

- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.


## Downstream dependencies

There are currently no downstream dependencies for this package.
