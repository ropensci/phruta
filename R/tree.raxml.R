#' Tree inference under RAxML
#'
#' Performs tree inference under \code{"RAxML"} for aligned fasta sequences in
#' a given folder (default is \code{"2.Alignments"}).
#'
#' @param folder Name of the folder where the sequences to align are stored (character).
#' @param FilePatterns A string that is common to all the target files in the relevant folder (character). Note that
#'                     this argument can be set to \code{"NULL"} if no specific pattern wants to be analized.
#' @param raxml_exec Where to find \code{"RAxML"} or how to run it from the console? (string).
#' @param Bootstrap Number of bootstrap replicates (numeric).
#' @param outgroup A single string of comma-separated tip labels to be used as outgroup in
#'                 \code{"RAxML"} See \code{"RAxML"} documentation for more details (character).
#' @param ... Arguments passed to \code{"ips::raxml"}.
#'
#' @importFrom ips raxml
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
#' sq.aln(folder = "1.CuratedSequences")
#' tree.raxml(
#'   folder = "2.Alignments", FilePatterns = "Masked",
#'   raxml_exec = "raxmlHPC", Bootstrap = 100,
#'   outgroup = "Manis_pentadactyla"
#' )
#' }
#' @export

tree.raxml <- function(folder = "2.Alignments", FilePatterns = "Masked", raxml_exec = "raxmlHPC", Bootstrap = 100, outgroup, ...) {
  if (is.null(folder)) stop("Please provide folder names")
  if (!is.character(raxml_exec)) stop("Please provide a raxml_exec argument of class character")
  if (!is.numeric(Bootstrap)) stop("Please provide a number for the Bootstrap argument")
  if (Bootstrap == 0) stop("Please indicate more than a single bootstrap replicate")

  files_fullNames <- list.files(folder, FilePatterns, full.names = T)
  files <- list.files(folder, "Masked")
  seq <- lapply(lapply(files_fullNames, read.FASTA), as.matrix)
  names(seq) <- files

  concatenated <- do.call(cbind.DNAbin, c(seq,
    fill.with.gaps = TRUE
  ))
  unlink("3.Phylogeny", recursive = TRUE)
  dir.create("3.Phylogeny")
  mainDir <- getwd()
  setwd(paste0(mainDir, "/", "3.Phylogeny"))
  tr <- raxml(
    DNAbin = concatenated, m = "GTRGAMMA",
    f = "a", N = Bootstrap, p = 1234, x = 1234,
    k = T,
    exec = raxml_exec, threads = 4,
    file = "phruta",
    outgroup = outgroup, ...
  )
  setwd(mainDir)
}
