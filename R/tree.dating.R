#' Tree dating under treePL or
#'
#' Performs tree dating under \code{"treePL"} or \code{"PATHd-8"} based on secondary calibrations.
#' Note that \code{"treePL"} or \code{"PATHd-8"} must be installed in your PATH. How to install
#' \code{"PATHd-8"} in mac \code{"https://gist.github.com/cromanpa94/a43bc710a17220f71d796d6590ea7fe4"} and
#' \code{"treePL"} can be installed using homebrew (brew install brewsci/bio/treepl). Thanks to
#' Brian O'Meara and Jonathan Chang, respectively.
#'
#' @param taxonomyFolder Name of the folder where \code{"1.Taxonomy.csv"}, created duing the
#'                       \code{"sq.curate"} step, is stored (character).
#' @param phylogenyFolder Name of the folder where \code{"RAxML_bipartitions.phruta"}, created duing the
#'                       \code{"tree.raxml"} step, is stored (character).
#' @param ... Arguments passed to \code{"geiger::congruify.phylo"}.
#'
#'
#' @importFrom geiger congruify.phylo
#' @import ape
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
#' tree.dating(
#'   taxonomyFolder = "1.CuratedSequences",
#'   phylogenyFolder = "3.Phylogeny", scale = "treePL"
#' )
#' }
#' @export



tree.dating <- function(taxonomyFolder = "1.CuratedSequences", phylogenyFolder = "3.Phylogeny", ...) {
  if (is.null(taxonomyFolder) | is.null(phylogenyFolder)) stop("Please provide folder names")

  taxonomy <- read.csv(paste0(taxonomyFolder, "/1.Taxonomy.csv"))
  row.names(taxonomy) <- taxonomy$species_names
  TargetTree <- read.tree(paste0(phylogenyFolder, "/RAxML_bipartitions.phruta"))
  resphy <- congruify.phylo(reference = SW.phruta, target = TargetTree, taxonomy = as.matrix(taxonomy), ...)
  names(resphy) <- c("family", "order", "class", "phyla", "kingdom")
  resphy <- Filter(Negate(anyNA), resphy)

  # Export
  if (length(resphy) > 0) {
    unlink("4.Timetree", recursive = TRUE)
    dir.create("4.Timetree")
    invisible(
      lapply(seq_along(resphy), function(x) {
        write.tree(resphy[[x]]$phy, paste0("4.Timetree/", names(resphy)[x], "-levelCalibration.tre"))
        write.csv(resphy[[x]]$calibrations, paste0("4.Timetree/", names(resphy)[x], "-levelCalibration.calibrated.csv"))
      })
    )
  }
}
