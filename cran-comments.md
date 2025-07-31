## sanba 0.0.2

We have implemented several revisions to address the feedback received during 
the package's vignette pre-screening review.

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
