# sanba 0.0.2

* Removed an unused `LinkingTo:` entry from the `DESCRIPTION` file, following a suggested pull request.
* Deleted the unused function `extract_best()`.
* Standardized the title case across documentation files.
* Moved the packages `scales` and `RColorBrewer` from `Depends` to `Imports`.
* Replaced `print()` calls with `message()` in MCMC plotting functions for cleaner output.
* Added various accessor functions for extracting quantities of interest from fitted objects (e.g., runtime, model type, and estimated values), and for interpreting/visualizing the result of posterior inference (e.g., posterior similarity matrix, estimated partition, number of clusters).
* Corrected minor typos and improved clarity in the documentation.
* Added a method for the `estimate_G()` function.
* Created appropriate `summary()` methods; now the partition estimation has its own function, `estimate_cluster()`
* Amended few details (e.g., assignment symbols and logical variables) following the `goodpractice::gp()` suggestions. Also removed `Date` field in the description as suggested.
* Corrected a bug in `estimate_G`; now it works even when there is a single DC detected.
* Made the `plot()` method for the `plot.partition_vi` more robust.

# sanba 0.0.1

* Initial CRAN submission.
* Released to the public the first official version
