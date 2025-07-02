## sanba 0.0.2

We have implemented minor changes in response to comply with the requests of a manuscript pre-screening review.

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces 

> checking for future file timestamps ... NOTE  
  unable to verify current time  
> 0 errors ✔ | 0 warnings ✔ | 1 note ✖  

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

> checking for future file timestamps ... NOTE
  unable to verify current time
> 0 errors ✔ | 0 warnings ✔ | 1 note ✖

- Running `devtools::check_win_devel()` produces

Maintainer: 'Francesco Denti <francescodenti.personal@gmail.com>'

> New submission  
Possibly misspelled words in DESCRIPTION:
  Camerlenghi (20:20)
  Canale (21:25)
  D'Angelo (17:256)
  Denti (17:215, 20:13, 22:25)
  Guindani (20:33, 21:37)
  Variational (3:55)
  Yu (21:33)
  al (17:224, 17:268)
  et (17:221, 17:265)
  
None of the above words are misspelled.  

- Running `devtools::check_mac_release()` produces
  
>* checking installed package size ... NOTE
  installed size is  8.7Mb
  sub-directories of 1Mb or more:
    libs   8.3Mb
  
- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.


## Downstream dependencies

There are currently no downstream dependencies for this package.
