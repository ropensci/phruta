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
#' @param mergeGeneFiles A named list, with each element being a character vector
#'                       indicating the names of the files in \code{"0.Sequences"}
#'                       that need to be combined into a single fasta file. For instance,
#'                       you can use this argument to combine CO1 and COI.
#' @param database A name of a database with taxonomic information.
#'                 Although 'gbif' is faster, it only has information for
#'                 animals and plants. Other databases follow
#'                 taxize::classification.
#' @param kingdom Optional and only used when database='gbif'.
#'                Two possible options: "animals" or "plants."
#' @param folder The name of the folder where the original sequences are
#'               located (character).
#' @param sqs.object A list of sequences generated from \code{sq.retrieve.indirect}. Only use if you're
#'                   not interested in download sequences locally.
#' @param removeOutliers Whether  \code{odseq:odseq} should be used to remove outliers
#' @param minSeqs minimum number of sequences per locus
#' @param threshold Relative to \code{odseq::odseq}. Only important if
#'                  \code{removeOutliers = TRUE}
#' @param ranks The taxonomic ranks used to examine the taxonomy of the species
#'              in the \code{0.Sequences} folder.
#'
#' @import stats
#' @import utils
#' @import ape
#' @import rgbif
#' @import taxize
#' @import odseq
#' @import msa
#' @importFrom Biostrings DNAStringSet
#'
#' @name sq.curate
#'
#' @return This function will return an object of class \code{list} with the
#'        following elements. First, the curated sequences with original names.
#'        Second, the curated sequences with species-level names. Third, 
#'        the accession numbers table. Fourth, a summary of taxonomic 
#'        information for all the species sampled in the files.
#'
#' @examples
#' \dontrun{
#' sq.retrieve.direct(
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
                      mergeGeneFiles = NULL,
                      database = "gbif",
                      kingdom = NULL,
                      folder = "0.Sequences",
                      sqs.object = NULL,
                      removeOutliers = TRUE,
                      minSeqs = 5,
                      threshold = 0.05,
                      ranks =
                        c("kingdom",
                          "phylum",
                          "class",
                          "order",
                          "family",
                          "genus",
                          "species")) {
  if (is.null(filterTaxonomicCriteria))
    stop("Please provide filtering pattern in the
         filterTaxonomicCriteria argument")
  if (length(filterTaxonomicCriteria) > 1)
    stop("Use a single string in the filterTaxonomicCriteria argument.
         \n SOLUTION: Split multiple criteria using |")
  if (is.null(folder)) stop("Folder where curated
                            sequences are saved must be provided")


  if(is.null(sqs.object)){

    ##Over-writing?
    if( !isTRUE(pkg.env$.testMode) ) {
      UI <- readline(paste0("This function might overwrite ",
                            "1.CuratedSequences", ". Are you sure you want to continue? (y/n)  "))
      if(UI != 'y') stop('Exiting since you did not press y')
    }


    if (!is.null(mergeGeneFiles)) {

      mergedSeqs <- lapply(mergeGeneFiles, function(x){
        refF <- list.files(folder, pattern = '.fasta')
        targetF <- paste0(unlist(x), '.fasta')
        if (all(targetF %in% refF)) {

          seqs <- lapply(x, function(z) read.FASTA(paste0(folder,"/", z, '.fasta')))

          #unlink(paste0(folder,"/", x, '.fasta')) #remove original files
          file.rename(paste0(folder,"/", x, '.fasta'), paste0(folder,"/", x, '.original.non.merged'))

          comS <- do.call(c,seqs)
          comS[!duplicated(sub(".*? ", "", names(comS)))]
        }else{
          message("\nFiles ", paste(targetF, collapse = " AND "), ", expected to be merged, not found in ", folder,"\n")
        }
      })

      names(mergedSeqs) <- names(mergeGeneFiles)
      invisible(
        lapply(seq_along(mergedSeqs), function(y){
          if(!is.null(mergedSeqs[[y]])){
            write.FASTA(mergedSeqs[[y]],
                        paste0(folder,"/", names(mergedSeqs)[y], ".fasta"))
          }
        })
      )
    }



  fastaSeqs <- lapply(list.files(folder, pattern = '.fasta', full.names = TRUE), function(x){
    seqs <- read.FASTA(x)
    seqs <- seqs[!duplicated(names(seqs))]
    seqs[!duplicated(sub(".*? ", "", names(seqs)))]
  })

  names(fastaSeqs) <- list.files(folder, pattern = '.fasta', full.names = F)

  seqNames <- unlist(lapply(unlist(lapply(fastaSeqs, names)),
                            function(x){
    paste0(strsplit(x, " ")[[1]][2:3], collapse = "_")
    }))

  seqAccN <- unlist(lapply(unlist(lapply(fastaSeqs, names)), function(x){
    paste0(strsplit(x, " ")[[1]][1], collapse = "_")}))

  AccDat <- data.frame("OriginalNames" = unlist(lapply(fastaSeqs, names)),
                       "AccN" = seqAccN,
                       "Species" = seqNames)
  AccDat$file <- rep(list.files(folder,pattern = '.fasta', full.names = F),
                     unlist(lapply(fastaSeqs, length)))

  species_names <- unique(AccDat$Species)

  Taxonomy_species <- if (database == "gbif") {
    taxonomy.retrieve(
      species_names = species_names, database = "gbif",
      kingdom = kingdom, ranks = ranks
    )
  } else {
    taxonomy.retrieve(species_names = species_names, database = "itis", ranks = ranks)
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
    if(length(x) >= minSeqs){
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
  seqsRemove <- names(which(unlist(resOut) == TRUE))
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
    Full_dataset[Full_dataset$originalSpeciesName != Full_dataset$species_names,]


  ## Export
  unlink("1.CuratedSequences", recursive = TRUE)
  dir.create("1.CuratedSequences")
  invisible(lapply(seq_along(curatedSeqs), function(y) {
    if (length(names(curatedSeqs[[y]])) >= minSeqs) {
      ## Original
      ta.cu <- curatedSeqs[[y]]
      ta.cu <- ta.cu[!duplicated(sub(".*? ", "", names(ta.cu)))]

      write.FASTA(ta.cu,
                  paste0("1.CuratedSequences/", names(curatedSeqs)[y]))
      ## Renamed
      newNames <-
        unlist(lapply(names(ta.cu),
                      function(x) {
                        paste(strsplit(x, " ")[[1]][2:3], collapse = "_")
                        }
                      ))

      renamed <- ta.cu
      if (nrow(toRename) > 0) {
        newNames <- ifelse(newNames %in% toRename$originalSpeciesName,
                           toRename$species_names, newNames)
      }
      names(renamed) <- newNames
      renamed <- renamed[!duplicated(names(renamed))]

      write.FASTA(renamed, paste0("1.CuratedSequences/renamed_",
                                  names(curatedSeqs)[y]))
    }
  }))

  perDS <- unlist(lapply(curatedSeqs, function(x){
   spps <- sub(" ", "_", sub(".*? ", "", names(x)))
   spps <- lapply(spps, function(y) strsplit(y, " ")[[1]][1])
   spps <- sub(" ", "", spps)
   spps <- unlist(lapply(spps, function(y){
   Full_dataset[Full_dataset$originalSpeciesName == y, "species_names"]
   }))
   codes <- sub(" .*", "", names(x))
   codes <- gsub("\\..*","",codes)
   codes[!duplicated(spps)]
  } ))

  AccDat$AccN <- gsub("\\..*","",AccDat$AccN)

  AccDat <- AccDat[AccDat$AccN %in% perDS,]

  AccDat <- AccDat[AccDat$file %in% names(which(table(AccDat$file) >= minSeqs)), ]
  Full_dataset <-
    Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species, ]

  newspp <- unlist(lapply(AccDat$Species, function(x) {
    Full_dataset[Full_dataset$originalSpeciesName == x, "species_names"]
  }))

  AccDat$OldSpecies <- AccDat$Species
  AccDat$Species <- newspp
  row.names(AccDat) <- NULL

  Full_dataset <- Full_dataset[!duplicated(Full_dataset$species_names),]


  ##Create a summary of the dataset
  sumTable <-  as.data.frame.matrix(t(table(AccDat$file, AccDat$Species)))
  sumTable$species_names <- row.names(sumTable)
  TableCombined <- merge(Full_dataset,sumTable, by = 'species_names', all.y = TRUE)

  write.csv(AccDat, "1.CuratedSequences/0.AccessionTable.csv")
  write.csv(Full_dataset, "1.CuratedSequences/1.Taxonomy.csv")
  write.csv(TableCombined, "1.CuratedSequences/2.Taxonomy.Sampling.csv")

  }else{

    if (!is.null(mergeGeneFiles)) {
      mergedSeqs <- lapply(mergeGeneFiles, function(x){
        refF <- names(sqs.object)
        targetF <- unlist(x)
        if (all(targetF %in% refF)) {
          seqs <- sqs.object[which(refF %in% targetF)]
          comS <- do.call(c,seqs)
          names(comS) <- gsub("^.*\\.","", names(comS))
          comS[!duplicated(sub(".*? ", "", names(comS)))]
        }else{
          message("\nFiles ", paste(targetF, collapse = " AND "), ", expected to be merged, not found \n")
        }
      })
      names(mergedSeqs) <- names(mergeGeneFiles)

      if(! is.null(mergedSeqs[[1]]) ){
      sqs.object <- sqs.object[-which(names(sqs.object) %in% unlist(mergeGeneFiles))]
      sqs.object <- c(sqs.object, mergedSeqs)
      }else{
        sqs.object
      }

    }



    fastaSeqs <- lapply(sqs.object, function(x){
      seqs <- x[!duplicated(names(x))]
      seqs[!duplicated(sub(".*? ", "", names(seqs)))]
    })


    seqNames <- unlist(lapply(fastaSeqs,
                              function(x){
                                lapply(strsplit(names(x), " "), function(z) paste0(z[2:3], collapse = "_"))
                              }))

    seqAccN <- unlist(lapply(unlist(lapply(fastaSeqs, names)), function(x){
      paste0(strsplit(x, " ")[[1]][1], collapse = "_")}))

    AccDat <- data.frame("OriginalNames" = unlist(lapply(fastaSeqs, names)),
                         "AccN" = seqAccN,
                         "Species" = seqNames)
    AccDat$file <- rep(names(fastaSeqs),
                       unlist(lapply(fastaSeqs, length)))

    species_names <- unique(AccDat$Species)

    Taxonomy_species <- if (database == "gbif") {
      taxonomy.retrieve(
        species_names = species_names, database = "gbif",
        kingdom = kingdom, ranks = ranks
      )
    } else {
      taxonomy.retrieve(species_names = species_names, database = "itis", ranks = ranks)
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
        if(length(x) >= minSeqs){
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
      seqsRemove <- names(which(unlist(resOut) == TRUE))
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
      Full_dataset[Full_dataset$originalSpeciesName != Full_dataset$species_names,]


    ## Export

    seqs.complete.genes <- lapply(seq_along(curatedSeqs), function(y) {
      if (length(names(curatedSeqs[[y]])) >= minSeqs) {
        ## Original
        ta.cu <- curatedSeqs[[y]]
        ta.cu <- ta.cu[!duplicated(sub(".*? ", "", names(ta.cu)))]




        ## Renamed
        newNames <-
          unlist(lapply(names(ta.cu),
                        function(x) {
                          paste(strsplit(x, " ")[[1]][2:3], collapse = "_")
                        }
          ))

        renamed <- ta.cu
        if (nrow(toRename) > 0) {
          newNames <- ifelse(newNames %in% toRename$originalSpeciesName,
                             toRename$species_names, newNames)
        }
        names(renamed) <- newNames
        renamed <- renamed[!duplicated(names(renamed))]

        seqsGene <- list(Original = ta.cu, Renamed = renamed )
        return(seqsGene)
      }
    })

    names(seqs.complete.genes) <- names(curatedSeqs)

    perDS <- unlist(lapply(curatedSeqs, function(x){
      spps <- sub(" ", "_", sub(".*? ", "", names(x)))
      spps <- lapply(spps, function(y) strsplit(y, " ")[[1]][1])
      spps <- sub(" ", "", spps)
      spps <- unlist(lapply(spps, function(y){
        Full_dataset[Full_dataset$originalSpeciesName == y, "species_names"]
      }))
      codes <- sub(" .*", "", names(x))
      codes <- gsub("\\..*","",codes)
      codes[!duplicated(spps)]
    } ))

    AccDat$AccN <- gsub("\\..*","",AccDat$AccN)
    AccDat <- AccDat[AccDat$AccN %in% perDS,]

    AccDat <- AccDat[AccDat$file %in% names(which(table(AccDat$file) >= minSeqs)), ]
    Full_dataset <-
      Full_dataset[Full_dataset$originalSpeciesName %in% AccDat$Species, ]

    newspp <- unlist(lapply(AccDat$Species, function(x) {
      Full_dataset[Full_dataset$originalSpeciesName == x, "species_names"]
    }))

    Full_dataset <- Full_dataset[!duplicated(Full_dataset$species_names),]

    AccDat$OldSpecies <- AccDat$Species
    AccDat$Species <- newspp
    row.names(AccDat) <- NULL

    ##Create a summary of the dataset
    sumTable <-  as.data.frame.matrix(t(table(AccDat$file, AccDat$Species)))
    sumTable$species_names <- row.names(sumTable)
    TableCombined <- merge(Full_dataset,sumTable, by = 'species_names', all.y = TRUE)

    toRet <- list("AccessionTable" = AccDat,
         "Taxonomy" = Full_dataset,
         "Taxonomy.Sampling" = TableCombined,
         "Sequences" = seqs.complete.genes)

    return(toRet)

  }


}
