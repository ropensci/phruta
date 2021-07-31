  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/cromanpa94/phruta/branch/main/graph/badge.svg)](https://codecov.io/gh/cromanpa94/phruta?branch=main)
  [![R-CMD-check](https://github.com/cromanpa94/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/cromanpa94/phruta/actions)
  [![](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![](https://travis-ci.org/cromanpa94/phruta.svg?branch=main)](https://travis-ci.org/cromanpa94/phruta)
  [![](https://img.shields.io/github/languages/code-size/cromanpa94/phruta.svg)](https://github.com/cromanpa94/phruta)
  [![R build status](https://github.com/cromanpa94/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/cromanpa94/phruta/actions)
  <!-- badges: end -->

## What is `phruta`

The `phruta` R package is designed to simplify the basic phylogenetic pipeline. Furthermore, the fact that all the intermediate steps are saved in independent folder and all the code is run within the same environment increases reproducibility. `phruta` retrieves gene sequences, combines newly downloaded and local gene sequences, and perform sequence alignment. 

`phruta` is also able to perform basic phylogenetic inference under `RAxML` on the resulting sequence alignments. The current release allows users to conduct tree dating based on secondary calibrations. `phruta` is essentially a wrapper for different R packages and additional software.
