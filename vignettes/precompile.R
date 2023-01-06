# Precompiled vignettes

library(knitr)
knit("vignettes/Exporting_data_phruta.Rmd.orig", "vignettes/Exporting_data_phruta.Rmd")
knit("vignettes/Future_phruta.Rmd.orig", "vignettes/Future_phruta.Rmd")
knit("vignettes/Phruta_advanced.Rmd.orig", "vignettes/Phruta_advanced.Rmd")
knit("vignettes/phruta_targetgenes.Rmd.orig", "vignettes/phruta_targetgenes.Rmd")
knit("vignettes/phruta.Rmd.orig", "vignettes/phruta.Rmd")
knit("vignettes/Phylogenetics_phruta.Rmd.orig", "vignettes/Phylogenetics_phruta.Rmd")
knit("vignettes/Future_phruta.Rmd.orig", "vignettes/Future_phruta.Rmd")
knit("vignettes/usando_phruta.Rmd.orig", "vignettes/usando_phruta.Rmd")
knit("vignettes/Using_phruta.Rmd.orig", "vignettes/Using_phruta.Rmd")


# Move figures into vignettes/ folder
#figs <- list.files(pattern = ".jp")
#fs::file_move(figs, fs::path("vignettes/", figs))
