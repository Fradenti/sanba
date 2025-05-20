## sanba 0.0.1

First release of the `sanba` package, which integrates, improves, and downstreams the functionalities of the `R` packages `SANple` and `SANvi`.

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces 

> checking for future file timestamps ... NOTE  
  unable to verify current time  
> 0 errors ✔ | 0 warnings ✔ | 1 note ✖  

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

> checking CRAN incoming feasibility ... [7s/17s] NOTE
> Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
> New submission
> 0 errors ✔ | 0 warnings ✔ | 0 notes ✔

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
