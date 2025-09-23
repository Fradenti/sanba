## sanba 0.0.3

* Replaced `arma::is_finite` with `std::isfinite` to ensure compatibility with the latest `RcppArmadillo` release and CRAN policies.
* Implemented small code adjustments to address and resolve compiler warnings.

## R CMD check results

- Running `devtools::check(args = c('--as-cran', '--no-manual'))` and `devtools::check(remote = TRUE, manual = TRUE)` 
  locally, produces the following NOTE:

> checking for future file timestamps ... NOTE
> Unable to verify current time

In contrast, running `devtools::check_win_devel()`, `devtools::check_win_release()`, and `devtools::check_mac_release()` 
  yields no errors, no warnings, and no notes.

- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.


## Downstream dependencies

There are currently no downstream dependencies for this package.
