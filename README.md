  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/cromanpa94/phruta/branch/main/graph/badge.svg)](https://codecov.io/gh/cromanpa94/phruta?branch=main)
  [![R-CMD-check](https://github.com/cromanpa94/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/cromanpa94/phruta/actions)
  [![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)  
  [![](https://img.shields.io/github/languages/code-size/cromanpa94/phruta.svg)](https://github.com/cromanpa94/phruta)
  [![CodeFactor](https://www.codefactor.io/repository/github/cromanpa94/phruta/badge)](https://www.codefactor.io/repository/github/cromanpa94/phruta)  <!-- badges: end -->

# `phruta` <a href='https://cromanpa94.github.io/phruta'><img src='man/figures/logo.png' align="right" height="300" /></a>

### Assembling phylogenetic trees from taxonomic names.

## What is `phruta`?

The `phruta` R package is designed to simplify the basic phylogenetic pipeline. Specifically, all code is run within the same program and data from intermediate steps are saved in independent folders. Furthrmore, all code is run within the same environment which increases the reproducibility of your analysis. `phruta` retrieves gene sequences, combines newly downloaded and local gene sequences, and performs sequence alignments. 

`phruta` is also able to perform basic phylogenetic inference under `RAxML` on the resulting sequence alignments. The current release allows users to conduct tree dating based on secondary calibrations. `phruta` is essentially a wrapper for alternative R packages and software.


## Why use `phruta`?

`phruta` simplifies the phylogenetic pipeline, increases reproducibility, and helps organizing information used to infer molecular phylogenies.


## Installing `phruta`

`phruta` is currently only available through `GitHub`. It can be easily installed using the following code.

```
library(devtools) 
install_github("cromanpa94/phruta")
```

## Dedication
My package is dedicated to my mom and every single Black woman. I still have lots of things to learn from you. You will always have all my admiration and respect.

The [logo](https://www.flickr.com/photos/gufomusike/3462117620/in/photolist-6NFiPi-xoLbca-FtC6yJ-4nk6wS-x2AZV-b3MUv8-e2B7qj-4uCwwa-e3PJxi-2ePGmUM-b2wBVi-obHf1x-5iP26P-4juoE6-z881E-z88t3-9GmTbQ-dGvrFe-22APdBs-p2t5Zv-8DWQw8-6fAJ2G-7jQhu2-7LEkkL-7vBdyF-jTdXSR-kcntD1-aWGfnx-bk59CK-5JfhKt-6gWfX7-reVehy-bjk7Ki-2xnGjv-dLJbq9-e3VjY3-ugz6U-FGVagm-iqVRuD-YE5pLe-2kPkt84-2kHhswd) features a Palenquera in Cartagena (Colombia). For many folks, Palenqueras are just the Black woman ones who sell fruits in particular Colombian turistic areas. However, palenqueras and Palenque are central to Black identity in Colombia, Latin America, and America. ["Palenque was the first free African town in the Americas"](https://en.wikipedia.org/wiki/San_Basilio_de_Palenque).

_Fruta_ is the Spanish word for _Fruit_. English _ph_ sounds the same as _F_ in Spanish. In `phruta`, _ph_ is relative to phylogenetics. I pronounce `phruta` as _fruta_ in Spanish.

## Additional resources

More details about the functions implemented in `phruta` can be found in the different vignettes associated with the package or in our [website](https://cromanpa94.github.io/phruta/).

## Contributing

Please see our [contributing guide](CONTRIBUTING).

## Contact

Please see the package [DESCRIPTION](DESCRIPTION) for package authors.

