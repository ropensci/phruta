#' Curate sequences from genbank
#'
#' After downloading sequences from genbank, this function
#' curates sequences based on taxonomic
#' information. Note that this function provides two summary datasets.
#' First, the accession numbers.
#' Second, the taxonomic information for each species in the database.
#' The taxonomy strictly follows
#' the gbif taxonomic backbone. The resulting files are saved
#'  to \code{"1.CuratedSequences"}. The
#' resulting files also have the most recent curated taxonomy
#'  following the gbif (or selected database) taxonomic backbone.
#'
#' @param filterTaxonomicCriteria A single string of terms (delimited using "|")
#'                                listing all the strings that could be used to
#'                                identify the species that should
#'                                be in the dataset (character).
#' @param database A name of a database with taxonomic information.
#'                 Although 'gbif' is faster, it only has information for
#'                 animals and plants. Other databases follow
#'                 taxize::classification.
#' @param kingdom Optional and only used when database='gbif'.
#'                Two possible options: "animals" or "plants."
#' @param folder The name of the folder where the original sequences are
#'               located (character).
#' @param removeOutliers Whether odseq:odseq should be used to remove outliers
#' @param minSeqs minimum number of sequences per locus
#'
#' @import stats
#' @import utils
#' @import ape
#' @import rgbif
#' @import taxize
#' @import odseq
#' @import msa
#' @import Biostrings
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' sq.retrieve(
#'   clades = c("Felis", "Vulpes", "Phoca"),
#'   species = "Manis_pentadactyla",
#'   genes = c("ADORA3", "CYTB")
#' )
#' sq.curate(
#'   filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
#'   database = "gbif", kingdom = "animals",
#'   folder = "0.Sequences"
#' )
#' }
#' @export

sq.curate <- function(filterTaxonomicCriteria = NULL,
                      database = "gbif",
                      kingdom = NULL,
                      folder = "0.Sequences",
                      removeOutliers = TRUE,
                      minSeqs=5,
                      threshold=0.05) {
  if (is.null(filterTaxonomicCriteria))
    stop("Please provide filtering pattern in the
         filterTaxonomicCriteria argument")
  if (length(filterTaxonomicCriteria) > 1)
    stop("Use a single string in the filterTaxonomicCriteria argument.
         \n SOLUTION: Split multiple criteria using |")
  if (is.null(folder)) stop("Folder where curated
                            sequences are saved must be provided")

  fastaSeqs <- lapply(list.files(folder, full.names = TRUE), function(x){
    seqs <- read.FASTA(x)
    seqs[!duplicated(names(seqs))]
  })

  names(fastaSeqs) <- list.files(folder, full.names = F)
  seqNames <- unlist(lapply(unlist(lapply(fastaSeqs, names)),
                            function(x){
    paste0(strsplit(x, " ")[[1]][c(2:3)], collapse = "_")
    }))

  seqAccN <- unlist(lapply(unlist(lapply(fastaSeqs, names)), function(x){
    paste0(strsplit(x, " ")[[1]][1], collapse = "_")}))

  AccDat <- data.frame("OriginalNames" = unlist(lapply(fastaSeqs, names)),
                       "AccN" = seqAccN,
                       "Species" = seqNames)
  AccDat$file <- rep(list.files(folder, full.names = F),
                     unlist(lapply(fastaSeqs, length)))

  species_names <- unique(AccDat$Species)

  Taxonomy_species <- if (database == "gbif") {
    taxonomy.retrieve(
      species_names = species_names, database = "gbif",
      kingdom = kingdom
    )
  } else {
    taxonomy.retrieve(species_names = species_names, database = "itis")
  }

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
    AccDat <- AccDat[ !AccDat$OriginalNames %in% dupDel  , ]
  }

  Full_dataset <- cbind.data.frame(Taxonomy_species, species_names)
  Full_dataset$originalSpeciesName <- Full_dataset$species_names
  Full_dataset$species_names <- sub(" ", "_", Full_dataset$species)
  Full_dataset <- na.omit(Full_dataset)

  "%nin%" <- Negate("%in%")

  toDel<-AccDat[which(AccDat$Species %nin% Full_dataset$originalSpeciesName), 1]

  # Remove any "non-species" species
  if (length(toDel) > 0) {
    fastaSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% toDel]
    })
    AccDat <- AccDat[!AccDat$OriginalNames %in% toDel, ]
  }



  ##Detect outliers...
  if(isTRUE(removeOutliers)){
  cat("\n Removing outliers...\n")
  fastaSeqsOutDet <- lapply(fastaSeqs, function(x){
    if(length(x)>minSeqs){
    seqs <- as.list(as.character(x))
    seqs <- unlist(lapply(seqs,paste0,collapse=""))
    seqs <- Biostrings::DNAStringSet(seqs)
    return(seqs)
    }
  })

  fastaSeqsOutDet <- Filter(Negate(is.null), fastaSeqsOutDet)
  fastaSeqsOutDetaln <- lapply(fastaSeqsOutDet, msa::msa)
  resOut <- lapply(fastaSeqsOutDetaln, odseq::odseq, distance_metric = "affine", B = 1000, threshold = threshold)
  names(resOut) <- NULL
  seqsRemove <- names(which(unlist(resOut) ==TRUE))
  AccDel <- gsub("\\..*","",seqsRemove)
  AccDat$AccN <- gsub("\\..*","",AccDat$AccN )
  namesDel <- AccDat[AccDat$AccN  %in%  AccDel,'OriginalNames']

  if (length(toDel) > 0) {
    fastaSeqs <- lapply(fastaSeqs, function(x) {
      x[!names(x) %in% namesDel]
    })
    AccDat <- AccDat[!AccDat$OriginalNames %in% namesDel, ]
  }
  }

  ##
  Full_dataset <-
    Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species,]
  WrongSpecies <-
    Full_dataset[!apply(Full_dataset, 1,
                        function(x) any(grepl(filterTaxonomicCriteria, x))),]
  RightSpecies <-
    Full_dataset[
      apply(Full_dataset, 1,function(x) any(grepl(filterTaxonomicCriteria,x))),]

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
  toRename <-
    Full_dataset[Full_dataset$originalSpeciesName!= Full_dataset$species_names,]

  ## Export
  unlink("1.CuratedSequences", recursive = TRUE)
  dir.create("1.CuratedSequences")
  invisible(lapply(seq_along(curatedSeqs), function(y) {
    if (length(names(curatedSeqs[[y]])) > minSeqs) {
      ## Original
      write.FASTA(curatedSeqs[[y]],
                  paste0("1.CuratedSequences/", names(curatedSeqs)[y]))
      ## Renamed
      newNames <-
        unlist(lapply(names(curatedSeqs[[y]]),
                      function(x) {
                        paste(strsplit(x, " ")[[1]][c(2:3)], collapse = "_")
                        }
                      ))
      renamed <- curatedSeqs[[y]]
      if (nrow(toRename) > 0) {
        newNames <- ifelse(newNames %in% toRename$originalSpeciesName,
                           toRename$species_names, newNames)
      }
      names(renamed) <- newNames
      write.FASTA(renamed, paste0("1.CuratedSequences/renamed_",
                                  names(curatedSeqs)[y]))
    }
  }))

  AccDat <- AccDat[AccDat$file %in% names(which(table(AccDat$file) > minSeqs)), ]
  Full_dataset <-
    Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species, ]
  newspp <- unlist(lapply(AccDat$Species, function(x) {
    Full_dataset[Full_dataset$originalSpeciesName == x, "species_names"]
  }))

  AccDat$OldSpecies <- AccDat$Species
  AccDat$Species <- newspp
  row.names(AccDat) <- NULL

  ##Create a summary of the dataset

  sumTable <-  as.data.frame.matrix(t(table(AccDat$file, AccDat$Species)))
  sumTable$species_names <- row.names(sumTable)
  TableCombined <- merge(Full_dataset,sumTable, by='species_names', all.y =TRUE)

  write.csv(AccDat, "1.CuratedSequences/0.AccessionTable.csv")
  write.csv(Full_dataset, "1.CuratedSequences/1.Taxonomy.csv")
  write.csv(TableCombined, "1.CuratedSequences/2.Taxonomy.Sampling.csv")

}
