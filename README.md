  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/cromanpa94/phruta/branch/main/graph/badge.svg)](https://codecov.io/gh/cromanpa94/phruta?branch=main)
  [![R-CMD-check](https://github.com/cromanpa94/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/cromanpa94/phruta/actions)
  [![](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![](https://travis-ci.com/cromanpa94/phruta.svg?branch=main)](https://travis-ci.com/cromanpa94/phruta?branch=main)
  [![](https://img.shields.io/github/languages/code-size/cromanpa94/phruta.svg)](https://github.com/cromanpa94/phruta)
  [![R build status](https://github.com/cromanpa94/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/cromanpa94/phruta/actions)
  <!-- badges: end -->

## What is `phruta`

The `phruta` R package is designed to simplify the basic phylogenetic pipeline. Specifically, all code is run within the same program and data from intermediate steps are saved in independent folders. Furthrmore, all code is run within the same environment which increases the reproducibility of your analysis. `phruta` retrieves gene sequences, combines newly downloaded and local gene sequences, and performs sequence alignments. 

`phruta` is also able to perform basic phylogenetic inference under `RAxML` on the resulting sequence alignments. The current release allows users to conduct tree dating based on secondary calibrations. `phruta` is essentially a wrapper for alternative R packages and software.


## Installing `phruta`

`phruta` is currently only available through `GitHub`. It can be easily installed using the following code.

```
library(devtools) 
install_github("hadley/dplyr")
```




