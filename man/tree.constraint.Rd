% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.constraint.R
\name{tree.constraint}
\alias{tree.constraint}
\title{Tree inference under RAxML}
\usage{
tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus",
    "species_names"),
  Topology = "((ingroup), outgroup);",
  outgroup = NULL
)
}
\arguments{
\item{taxonomy_folder}{Name of the folder where the
1.Taxonomy file is stored.}

\item{targetColumns}{Where to find \code{"RAxML"}
or how to run it from the console? (string).}

\item{Topology}{A string summarizing the desired topological
constraint in newick format.}

\item{outgroup}{Optional and only required when the topology
argument is not "((ingroup), outgroup);".}
}
\value{
None
}
\description{
Performs tree inference under \code{"RAxML"} for aligned fasta sequences in
a given folder (default is \code{"2.Alignments"}).
}
\examples{
\dontrun{
sq.retrieve.direct(
  clades = c("Felis", "Vulpes", "Phoca"),
  species = "Manis_pentadactyla",
  genes = c("ADORA3", "CYTB")
)
sq.curate(
  filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
  kingdom = "animals", folder = "0.Sequences"
)

tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family",
  "genus", "species_names"),
  Topology = "((ingroup), outgroup);",
  outgroup = "Manis_pentadactyla"
)
tree.constraint(
  taxonomy_folder = "1.CuratedSequences",
  targetColumns = c("kingdom", "phylum", "class", "order", "family",
  "genus", "species_names"),
  Topology = "((Felis), (Phoca));"
)
}
}
