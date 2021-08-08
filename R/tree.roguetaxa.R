#' RogueNaRok within phruta
#'
#' Implements the RogueNaRok algorithm for rogue taxon identification within phruta
#'
#' @param folder Name of the folder where the sequences to align are stored (character).
#' @param ... Arguments passed to \code{"Rogue::RogueTaxa"}.
#'
#' @importFrom Rogue RogueTaxa
#' @importFrom ape read.tree
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
#' tree.roguetaxa(folder='3.Phylogeny')
#' }
#' @export

tree.roguetaxa<-function(folder='3.Phylogeny',...){
  trees <- read.tree(paste0(folder,"/RAxML_bootstrap.phruta"))
  BestTree<-read.tree(paste0(folder,"/RAxML_bipartitions.phruta"))
  RT<-RogueTaxa(trees, bestTree=BestTree, ...)
  unlink("3.1.RogueTaxa", recursive = TRUE)
  dir.create('3.1.RogueTaxa')
  write.csv(RT, '3.1.RogueTaxa/RogueTaxa.csv')
}
