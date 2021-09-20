#' Run Partitionfinder v.1
#'
#' This function runs partitionfinder v1 within phruta. For now,
#' all analyses are based on genes.
#' Please note that you need at least two gene regions to run partitionfinder.
#'
#' @param folderAlignments Name of the folder where the sequences to align are
#'                         stored (character).
#' @param FilePatterns A string that is common to all the target files in the
#'                     relevant folder (character). Note that
#'                     this argument can be set to \code{"NULL"} if no specific
#'                     pattern wants to be analyzed.
#' @param folderPartitionFinder Name of the new folder where the output files
#'                              are stored (string).
#' @param models Models to run in partitionfinder (string).
#' @param run Run  partitionfinder?
#'
#' @importFrom ape read.FASTA cbind.DNAbin
#' @importFrom ips raxml.partitions write.phy
#' @importFrom utils download.file untar write.csv
#'
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
#'   kingdom = "animals", folder = "0.Sequences"
#' )
#' sq.aln(folder = "1.CuratedSequences")
#' sq.partitionfinderv1(
#'   folderAlignments = "2.Alignments",
#'   FilePatterns = "Masked",
#'   models = "all"
#' )
#' }
#' @export

sq.partitionfinderv1 <- function(folderAlignments = "2.Alignments",
                                 FilePatterns = "Masked",
                                 folderPartitionFinder ="2.1.PartitionFinderv1",
                                 models = "all", run = TRUE) {
  files_fullNames <- list.files(folderAlignments, FilePatterns,
                                full.names = TRUE)
  files <- list.files(folderAlignments, FilePatterns)
  seq <- lapply(lapply(files_fullNames, read.FASTA), as.matrix)
  names(seq) <- files

  concatenated <- do.call(cbind.DNAbin, c(seq,
    fill.with.gaps = TRUE
  ))
  partitions <- do.call(raxml.partitions, seq)

  unlink(folderPartitionFinder, recursive = TRUE)
  dir.create(folderPartitionFinder)
  write.phy(concatenated, paste0(folderPartitionFinder, "/concatenated.phy"))
  write.csv(partitions, paste0(folderPartitionFinder, "/partitions.csv"))

  if (run) {
  ## Set partitions
  if ("PartitionFinder.tar.gz" %in% list.files() == FALSE) {
   download.file("https://github.com/brettc/partitionfinder/archive/v1.1.1.zip",
                 "PartitionFinder.tar.gz")
    untar("PartitionFinder.tar.gz")
  }

  config_file <-
    readLines("partitionfinder-1.1.1/examples/nucleotide/partition_finder.cfg")

  block <- list()
  for (i in seq_along(partitions[,1])) {
    block[[i]] <- paste0("GENE_", i, " ", "=", " ",
                         partitions[i, 3], " ", "-", " ", partitions[i, 4], ";")
  }

  conf_fi_mod <- config_file[-c(16:24)]
  newconf <- append(conf_fi_mod, unlist(block), after = 15)
  newconf[9] <- paste("models =", models, ";")
  newconf[2] <- paste("alignment =", "concatenated.phy;")

  writeLines(newconf, paste0(folderPartitionFinder, "/partition_finder.cfg"))

  tofile <- paste0(folderPartitionFinder, "/partition_finder.cfg")
  system(paste0("python ./partitionfinder-1.1.1/PartitionFinder.py", " ",
                tofile), wait = TRUE)

  unlink("partitionfinder-1.1.1", recursive = TRUE)
  unlink("PartitionFinder.tar.gz", recursive = TRUE)
  }
}
