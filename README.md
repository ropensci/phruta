  <!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/ropensci/phruta/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/phruta?branch=main)
  [![R-CMD-check](https://github.com/ropensci/phruta/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/phruta/actions)
  [![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)  
  [![](https://img.shields.io/github/languages/code-size/ropensci/phruta.svg)](https://github.com/ropensci/phruta)
  [![CodeFactor](https://www.codefactor.io/repository/github/ropensci/phruta/badge)](https://www.codefactor.io/repository/github/ropensci/phruta)  
 [![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/458_status.svg)](https://github.com/ropensci/software-review/issues/458)
  <!-- badges: end -->

# The `phruta` `R` package <a href='https://ropensci.github.io/phruta'><img src='man/figures/logo.png' align="right" height="300" /></a>

### Assembling phylogenetic trees from taxonomic names

## What is `phruta`?

The `phruta` R package is designed to simplify the basic phylogenetic pipeline. All the code is run within the same program and data from intermediate steps are saved in independent folders (optional). `phruta` retrieves gene sequences, combines newly downloaded to local gene sequences, performs sequence alignments, and basic phylogenetic inference. 

## Who should consider using `phruta`

The main functions in the `phruta` R package allow for a quick mining and curation of GenBank sequences. This package is designed for students and researchers interested in generating species-level genetic datasets for particular sets of taxa. Specifically, if you have a clade or group of species in mind, `phruta` will help you to assemble a molecular dataset with information available in GenBank.


## Why use `phruta`?

`phruta` simplifies the phylogenetic pipeline, increases reproducibility, and helps organizing information used to infer molecular phylogenies.

## How is `phruta` different from other software?

`phruta` has two core functions. The main applications of these functions is briefly outlined below:

- `sq.retrieve.direct()` and `sq.retrieve.indirect()`: These functions downloads sequences from genbank (nucleotide database) for particular taxa (taxonomic groups or particular species) and a list of genes. 

- `sq.curate()`: After sequences are downloaded from genbank, this function curates sequences within each of the examined genes by detecting sequence outliers and by using taxonomic information. 

In addition to these two main functions, users will be able to align the downloaded sequences, infer phylogenetic trees, and calibrate phylogenies using additional functions in `phruta`.


## Installing `phruta`

`phruta` is currently only available through `GitHub`. It can be easily installed using the following code.

```
library(devtools) 
install_github("ropensci/phruta")
```

Alternatively, you can install `phruta` using:

```
install.packages("phruta", repos = "https://ropensci.r-universe.dev")
```

Please make sure that the `R` packages `msa`, `DECIPHER`, `Biostrings`, and `odseq` are correctly installed. If you are interested in using the development version of `phruta`, please install it using the following code:

```
library(devtools)
install_github("ropensci/phruta", ref = "dev")
```

## Running `phruta` from shiny

I have constructed a shiny app that hosts `phruta` and enables users to run the basic functions in a less-code intensive environment. The app, `salphycon` is currently available in the following [GitHub repo](https://github.com/cromanpa94/salphycon). The shiny app will be live at some point in 2023.


## Installing RAxML <a name="paragraph1"></a>


In `MacOS`, `RAxML` can be easily installed to the `PATH` using one of the two lines below in `conda`:

```{bash eval=FALSE}
conda install -c bioconda/label/cf201901 raxml 
```

```{bash eval=FALSE}
conda install -c bioconda raxml
```

For other `OS` (Windows, Linux), please follow the instructions listed in the official `RAxML` [website](https://cme.h-its.org/exelixis/web/software/raxml/)

Once `RAxML` has been installed to your computer, open `R` and make sure that the following line doesn't throw an error.

```{r eval=FALSE}
system("raxmlHPC")
```

Depending on how `RAxML` was installed, you may want to check if `RAxML` is called from the terminal using `raxmlHPC` or `raxmlHPC`. This string needs to be passed to `tree.raxml` using the argument `raxml_exec`. Please note that this argument corresponds to the `exec` argument in `ips::raxml`. 

Finally, note that `RStudio` sometimes has issues finding *stuff* in the path while using `system()`. If you're using `macOS`, try starting `RStudio` from the command line by running the following line:

```{bash eval=FALSE}
open /Applications/RStudio.app
```

VS code does not suffer of the same issues. In other OS, it might be better to simply avoid using `RStudio` if you're interested in running the phylogenetic functions in `phruta`.

## Installing `PATHd-8` and `treePL` <a name="paragraph2"></a>

There are excellent guides for installing `PATHd-8` and `treePL`. Here, I summarize two potentially relevant options.

First, you can use [Brian O'Meara's](https://github.com/bomeara/phydocker/blob/master/Dockerfile) approach for installing `PATHd-8` in MacOs and linux. I summarize the code in the following [link](https://gist.github.com/cromanpa94/a43bc710a17220f71d796d6590ea7fe4). For Windows users, please use the compiled version of the software provided in the following [link](https://www2.math.su.se/PATHd8/).

Second, you can use homebrew to install `treePL` (Windows, MacOS, and Linux), thanks to Jonathan Chang.

```{bash eval = F}
brew install brewsci/bio/treepl
```

Please check the following [link](https://docs.brew.sh/Homebrew-on-Linux)) if you're interested in running `brew` from Windows and Linux.


## Running `phruta` from `Rstudio` while using `MacOS`?

Only if you're interested in running phylogenetic analyses, please make sure you open `RStudio` using the following code from the terminal:

```{bash eval=FALSE}
open /Applications/RStudio.app
```

## Dedication

My package is dedicated to my mom. I still have lots of things to learn from you. You will always have all my admiration. The [logo](https://www.flickr.com/photos/gufomusike/3462117620/in/photolist-6NFiPi-xoLbca-FtC6yJ-4nk6wS-x2AZV-b3MUv8-e2B7qj-4uCwwa-e3PJxi-2ePGmUM-b2wBVi-obHf1x-5iP26P-4juoE6-z881E-z88t3-9GmTbQ-dGvrFe-22APdBs-p2t5Zv-8DWQw8-6fAJ2G-7jQhu2-7LEkkL-7vBdyF-jTdXSR-kcntD1-aWGfnx-bk59CK-5JfhKt-6gWfX7-reVehy-bjk7Ki-2xnGjv-dLJbq9-e3VjY3-ugz6U-FGVagm-iqVRuD-YE5pLe-2kPkt84-2kHhswd) features a Palenquera in Cartagena (Colombia). For many folks, Palenqueras are just the Black woman ones who sell fruits in particular Colombian turistic areas. However, palenqueras and Palenque are central to Black identity in Colombia, Latin America, and across the America: "Palenque was the first free African town in the Americas"](https://en.wikipedia.org/wiki/San_Basilio_de_Palenque).

## Etymology

_Fruta_ is the Spanish word for _Fruit_. English _ph_ sounds the same as _F_ in Spanish. In `phruta`, _ph_ is relative to phylogenetics. I pronounce `phruta` just as _fruta_ in Spanish.

## Additional resources

More details about the functions implemented in `phruta` can be found in the different vignettes associated with the package or in our [website](https://ropensci.github.io/phruta/).

## Alternatives to `phruta`

Similar functionalities for assembling curated molecular datasets for phylogenetic analyses can be found in [`phylotaR`](https://github.com/ropensci/phylotaR) and [SuperCRUNCH](https://github.com/dportik/SuperCRUNCH). However, note that `phylotaR` is limited to downloading and curating sequences (e.g. doesn't align sequences). Similarly, `SuperCRUNCH` only curates local sequences. `phruta` is closer to the [`SUPERSMART`](https://academic.oup.com/sysbio/article/66/2/152/2418028) and its "new" associated `R` workflow [`SUPERSMARTR`](https://github.com/AntonelliLab/supersmartR). However, most of the applications in the different packages that are part of `SUPERSMARTR` are simplified in `phruta`. 


## Contributing

Please see our [contributing guide](CONTRIBUTING).

## Contact

Please see the package [DESCRIPTION](DESCRIPTION) for package authors.

## Code of conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

