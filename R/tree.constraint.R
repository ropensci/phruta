#' Tree inference under RAxML
#'
#' Performs tree inference under \code{"RAxML"} for aligned fasta sequences in
#' a given folder (default is \code{"2.Alignments"}).
#'
#' @param taxonomy_folder Name of the folder where the
#'                        1.Taxonomy file is stored.
#' @param targetColumns Where to find \code{"RAxML"}
#'                      or how to run it from the console? (string).
#' @param Topology A string summarizing the desired topological
#'                 constraint in newick format.
#' @param outgroup Optional and only required when the topology
#'                 argument is not "((ingroup), outgroup);".
#'
#' @importFrom ape read.tree write.tree
#' @importFrom utils read.csv
#'
#' @return None
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
#'   kingdom = "animals", folder = "0.Sequences"
#' )
#'
#' tree.constraint(
#'   taxonomy_folder = "1.CuratedSequences",
#'   targetColumns = c("kingdom", "phylum", "class", "order", "family",
#'   "genus", "species_names"),
#'   Topology = "((ingroup), outgroup);",
#'   outgroup = "Manis_pentadactyla"
#' )
#' tree.constraint(
#'   taxonomy_folder = "1.CuratedSequences",
#'   targetColumns = c("kingdom", "phylum", "class", "order", "family",
#'   "genus", "species_names"),
#'   Topology = "((Felis), (Phoca));"
#' )
#' }
#' @export

tree.constraint <- function(taxonomy_folder = "1.CuratedSequences",
                            targetColumns = c("kingdom", "phylum", "class",
                                              "order", "family", "genus",
                                              "species_names"),
                            Topology = "((ingroup), outgroup);",
                            outgroup = NULL) {
  taxonomy <- read.csv(paste0(taxonomy_folder, "/1.Taxonomy.csv"))

  if (Topology == "((ingroup), outgroup);") {
    Topology1 <- Topology
    ingroup <- taxonomy[!taxonomy$species_names == outgroup, ]
    outgroup <- taxonomy[taxonomy$species_names == outgroup, ]
    clades <- list("ingroup" = getListConstraints(ingroup, targetColumns,
                                                  byClades = FALSE),
                   "outgroup" =
                     getListConstraints(outgroup,
                                        targetColumns, byClades = FALSE))
    for (i in seq_along(clades)) {
      Topology <- gsub(names(clades[i]), clades[[i]], Topology)
    }
  } else {
    cstByClade <- invisible(getListConstraints(taxonomy, targetColumns,
                                               byClades = TRUE))
    Topology1 <- Topology
    TopologyOriginal <- Topology
    Topology <- sub(";", "", Topology, fixed = TRUE)
    Topology <- gsub("[()]", "", Topology)
    Topology <- sub(" ", "", Topology, fixed = TRUE)
    clades <- strsplit(Topology, ",")[[1]]
    for (i in seq_along(clades)) {
      TopologyOriginal <- gsub(clades[i],
                               cstByClade[names(cstByClade) %in% clades[i]],
                               TopologyOriginal)
    }
    Topology <- TopologyOriginal
  }

  unlink("3.2.Phylogeny.constraint", recursive = TRUE)
  dir.create("3.2.Phylogeny.constraint")
  write(Topology1, "3.2.Phylogeny.constraint/OriginalConstraints.cst.tre")
  write(Topology, "3.2.Phylogeny.constraint/phruta.cst.tre")
  tree <- read.tree("3.2.Phylogeny.constraint/phruta.cst.tre")
  write.tree(tree, "3.2.Phylogeny.constraint/phruta.ape.cst.tre")
}
