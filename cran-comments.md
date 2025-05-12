## sanba 0.0.1

First release of the `sanba` package, which integrates, improves and downstreams the functionalities of the `R` packages `SANple` and `SANvi`.

## R CMD check results

- Running `devtools::check(args = c('--as-cran','--no-manual'))` locally produces no errors, warnings, or notes.

- Running `devtools::check(remote = TRUE, manual = TRUE)` produces

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Francesco Denti <francescodenti.personal@gmail.com>’
   
   New maintainer:
     Francesco Denti <francescodenti.personal@gmail.com>
   Old maintainer(s):
     Laura D'Angelo <laura.dangelo@live.com>

due to a change in the maintainer of the package.
      
- Running `devtools::check_win_devel()` produces

0 errors | 0 warnings | 1 note

❯ checking CRAN incoming feasibility ... NOTE
Maintainer: 'Francesco Denti <francescodenti.personal@gmail.com>'

New maintainer:
  Francesco Denti <francescodenti.personal@gmail.com>
Old maintainer(s):
  Laura D'Angelo <laura.dangelo@live.com>
  
- Finally, in its current state, this package also passes all the standard checks performed via *GitHub actions*.

## Downstream dependencies

There are currently no downstream dependencies for this package.
