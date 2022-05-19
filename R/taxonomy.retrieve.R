#' Curate sequences from genbank
#'
#' After downloading sequences from genbank, this function curates sequences
#' based on taxonomic information. Note that this function provides two summary
#' datasets. First, the accession numbers. Second, the taxonomic information
#' for each species in the database. The taxonomy strictly follows
#' the gbif taxonomic backbone. The resulting files are saved to
#'\code{"1.CuratedSequences"}. The resulting files also have the most recent
#'curated taxonomy following the gbif taxonomic backbone.
#'
#' @param species_names A vector of species names
#' @param database A name of a database with taxonomic information.
#'                 'gbif' only works for animals and plants. Databases
#'                 follows taxize::classification
#' @param kingdom Optional and only used when database='gbif'. Two possible
#'                options: "animals" or "plants."
#'
#' @import rgbif
#' @import taxize
#'
#' @return data.frame of taxonomic information for the target species
#'         (valid in the database)
#'
#' @examples
#' \dontrun{
#' taxonomy.retrieve(
#'   species_names = c(
#'     "Felis_catus", "PREDICTED:_Vulpes",
#'     "Phoca_largha", "PREDICTED:_Phoca",
#'     "PREDICTED:_Manis", "Felis_silvestris", "Felis_nigripes"
#'   ),
#'   database = "gbif", kingdom = "animals"
#' )
#' }
#' @export

taxonomy.retrieve <-
  function(species_names = NULL,
           database = "gbif",
           kingdom = NULL,
           ranks =
             c("kingdom",
               "phylum",
               "class",
               "order",
               "family",
               "genus",
               "species")) {


    if (database != "gbif") {
      taxo <- classification(species_names, db = database, rows = 1)
      invisible(Taxonomy_species <-
                  as.data.frame(do.call(
                    rbind, lapply(seq_along(taxo), function(x) {
                      if (!is.na(taxo[[x]])) {
                        t(data.frame(taxo[[x]][taxo[[x]][, 2] %in% ranks, 1]))
                      } else {
                        sma <- matrix(nrow = 1, ncol = length(ranks) - 1)
                        cbind(sma, species_names[x])
                      }
                    })
                  )))
      row.names(Taxonomy_species) <- NULL
      colnames(Taxonomy_species) <- ranks
      Taxonomy_species$species <-
        sub(" ", "_", Taxonomy_species$species)
      Taxonomy_species
    } else {
      ## If animals and plants
      gbifkey <-
        lapply(species_names, function(x)
          name_backbone(name = x, kingdom = kingdom))
      keys <- pblapply(seq_along(gbifkey), function(x) {
        if (
          as.character(gbifkey[[x]][which(names(gbifkey[[x]])
                                          == "matchType")]) == "NONE") {
          0
        } else {
          if (length(which(names(gbifkey[[x]]) == "acceptedUsageKey")) == 0) {
            as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "usageKey")])
          } else {
            as.character(gbifkey[[x]][which(names(gbifkey[[x]])
                                            == "acceptedUsageKey")])
          }
        }
      })

      gbif_taxonomy <-
        lapply(unlist(keys), function(x)
          as.data.frame(name_usage(key = x)$data))
      Taxonomy_species <-
        lapply(seq_along(gbif_taxonomy), function(y) {
          sub1 <- gbif_taxonomy[[y]]
          cate <-
            t(data.frame(unlist(lapply(seq_along(ranks), function(x) {
              nu <- which(colnames(sub1) == ranks[x])
              if (length(nu) != 1) {
                NA
              } else {
                sub1[, nu]
              }
            }))))
          colnames(cate) <- ranks
          row.names(cate) <- NULL
          cate
        })
      Taxonomy_species <- do.call(rbind.data.frame, Taxonomy_species)
      Taxonomy_species
    }
  }
