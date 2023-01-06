# Precompiled vignettes

knitr::knit("vignettes/Exporting_data_phruta.Rmd.orig", "vignettes/Exporting_data_phruta.Rmd")
knitr::knit("vignettes/Future_phruta.Rmd.orig", "vignettes/Future_phruta.Rmd")
knitr::knit("vignettes/Phruta_advanced.Rmd.orig", "vignettes/Phruta_advanced.Rmd")
knitr::knit("vignettes/phruta_targetgenes.Rmd.orig", "vignettes/phruta_targetgenes.Rmd")
knitr::knit("vignettes/phruta.Rmd.orig", "vignettes/phruta.Rmd")
knitr::knit("vignettes/Phylogenetics_phruta.Rmd.orig", "vignettes/Phylogenetics_phruta.Rmd")
knitr::knit("vignettes/Future_phruta.Rmd.orig", "vignettes/Future_phruta.Rmd")
knitr::knit("vignettes/usando_phruta.Rmd.orig", "vignettes/usando_phruta.Rmd")
knitr::knit("vignettes/Using_phruta.Rmd.orig", "vignettes/Using_phruta.Rmd")


# Move figures into vignettes/ folder
#figs <- list.files(pattern = ".jp")
#fs::file_move(figs, fs::path("vignettes/", figs))
