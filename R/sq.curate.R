#' Curate sequences from genbank
#'
#' After downloading sequences from genbank, this function curates sequences based on taxonomic
#' information. Note that this function provides two summary datasets. First, the accession numbers.
#' Second, the taxonomic information for each species in the database. The taxonomy strictly follows
#' the gbif taxonomic backbone. The resulting files are saved to \code{"1.CuratedSequences"}. The
#' resulting files also have the most recent curated taxonomy following the gbif taxonomic backbone.
#'
#' @param filterTaxonomicCriteria A single string of terms (delimited using "|") listing all the
#'                                strings that could be used to identify the species that should
#'                                be in the dataset (character).
#' @param kingdom Two options: animals or plants (character).
#' @param folder The name of the folder where the original sequences are (character).
#'
#' @import stats
#' @import utils
#' @import ape
#' @import rgbif
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' sq.retrieve(
#'   clades = c("Felis", "Vulpes", "Phoca"),
#'   species = "Manis_pentadactyla",
#'   genes = c("ADORA3")
#' )
#' sq.curate(
#'   filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
#'   kingdom = "animals", folder = "0.Sequences"
#' )
#' }
#' @export

sq.curate <- function(filterTaxonomicCriteria = NULL, kingdom = "animals", folder = "0.Sequences") {
  if (is.null(filterTaxonomicCriteria)) stop("Please provide filtering pattern in the filterTaxonomicCriteria argument")
  if (length(filterTaxonomicCriteria) > 1) stop("Use a single string in the filterTaxonomicCriteria argument. \n SOLUTION: Split multiple criteria using |")
  if (is.null(folder)) stop("Folder where curated sequences are saved must be provided")
  match.arg(kingdom, c("animals", "plants"))

  fastaSeqs <- lapply(list.files(folder, full.names = T), read.FASTA)
  names(fastaSeqs) <- list.files(folder, full.names = F)
  seqNames <- sapply(unlist(lapply(fastaSeqs, names)), function(x) paste0(strsplit(x, " ")[[1]][c(2:3)], collapse = "_"))
  seqAccN <- sapply(unlist(lapply(fastaSeqs, names)), function(x) paste0(strsplit(x, " ")[[1]][1], collapse = "_"))
  AccDat <- data.frame("OriginalNames" = unlist(lapply(fastaSeqs, names)), "AccN" = seqAccN, "Species" = seqNames)
  AccDat$file <- rep(list.files(folder, full.names = F), sapply(fastaSeqs, length))

  species_names <- unique(AccDat$Species)

  gbifkey <- lapply(species_names, function(x) name_backbone(name = x, kingdom = kingdom))
  keys <- pblapply(1:length(gbifkey), function(x) {
    if (as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "matchType")]) == "NONE") {
      0
    } else {
      if (length(which(names(gbifkey[[x]]) == "acceptedUsageKey")) == 0) {
        as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "usageKey")])
      } else {
        as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "acceptedUsageKey")])
      }
    }
  })

  gbif_taxonomy <- lapply(unlist(keys), function(x) as.data.frame(name_usage(key = x)$data))
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  Taxonomy_species <- lapply(seq_along(gbif_taxonomy), function(y) {
    sub1 <- gbif_taxonomy[[y]]
    cate <- t(data.frame(unlist(lapply(seq_along(ranks), function(x) {
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

  # Remove PREDICTED species
  if (length(grep("PREDICTED", AccDat$Species)) > 0) {
    dupDel <- AccDat[grep("PREDICTED", AccDat$Species), "OriginalNames"]
    fastaSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% dupDel]
    })
    AccDat <- AccDat[-grep("PREDICTED", AccDat$Species), ]
  }

  # Remove duplicated species
  if (any(duplicated(AccDat[, c(3:4)]))) {
    dupDel <- AccDat[duplicated(AccDat[, c(3:4)]), "OriginalNames"]
    fastaSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% dupDel]
    })
    AccDat <- AccDat[!duplicated(AccDat[, c(3:4)]), ]
  }

  Full_dataset <- cbind.data.frame(Taxonomy_species, species_names)
  Full_dataset$originalSpeciesName <- Full_dataset$species_names
  Full_dataset$species_names <- sub(" ", "_", Full_dataset$species)
  Full_dataset <- na.omit(Full_dataset)

  "%nin%" <- Negate("%in%")

  toDel <- AccDat[which(AccDat$Species %nin% Full_dataset$originalSpeciesName), 1]

  # Remove any non-species species
  if (length(toDel) > 0) {
    fastaSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% toDel]
    })
    AccDat <- AccDat[!AccDat$OriginalNames %in% toDel, ]
  }

  Full_dataset <- Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species, ]
  WrongSpecies <- Full_dataset[!apply(Full_dataset, 1, function(x) any(grepl(filterTaxonomicCriteria, x))), ]
  RightSpecies <- Full_dataset[apply(Full_dataset, 1, function(x) any(grepl(filterTaxonomicCriteria, x))), ]

  if (nrow(WrongSpecies) > 0) {
    seqsToDel <- AccDat[AccDat$Species %in% WrongSpecies, "OriginalNames"]
    AccDat <- AccDat[!AccDat$Species %in% WrongSpecies, ]
    curatedSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% seqsToDel]
    })
    names(curatedSeqs) <- names(fastaSeqs)
  } else {
    curatedSeqs <- fastaSeqs
  }

  ## Rename incorrect synonyms
  toRename <- Full_dataset[Full_dataset$originalSpeciesName != Full_dataset$species_names, ]

  ## Export
  unlink("1.CuratedSequences", recursive = TRUE)
  dir.create("1.CuratedSequences")
  invisible(lapply(seq_along(curatedSeqs), function(y) {
    if (length(names(curatedSeqs[[y]])) > 1) {
      ## Original
      write.FASTA(curatedSeqs[[y]], paste0("1.CuratedSequences/", names(curatedSeqs)[y]))
      ## Renamed
      newNames <- sapply(names(curatedSeqs[[y]]), function(x) paste(strsplit(x, " ")[[1]][c(2:3)], collapse = "_"))
      renamed <- curatedSeqs[[y]]
      if (nrow(toRename) > 0) {
        newNames <- ifelse(newNames %in% toRename$originalSpeciesName, toRename$species_names, newNames)
      }
      names(renamed) <- newNames
      write.FASTA(renamed, paste0("1.CuratedSequences/renamed_", names(curatedSeqs)[y]))
    }
  }))

  AccDat <- AccDat[AccDat$file %in% names(which(table(AccDat$file) > 1)), ]
  Full_dataset <- Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species, ]
  newspp <- sapply(AccDat$Species, function(x) {
    Full_dataset[Full_dataset$originalSpeciesName == x, "species_names"]
  })

  AccDat$OldSpecies <- AccDat$Species
  AccDat$Species <- newspp
  row.names(AccDat) <- NULL

  write.csv(AccDat, "1.CuratedSequences/0.AccessionTable.csv")
  write.csv(Full_dataset, "1.CuratedSequences/1.Taxonomy.csv")
}
