% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene.sampling.retrieve.R
\name{gene.sampling.retrieve}
\alias{gene.sampling.retrieve}
\title{Retrieve the distribution of genes for given organism in genbank}
\usage{
gene.sampling.retrieve(
  organism,
  speciesSampling = TRUE,
  npar = 2,
  nSearchesBatch = 499
)
}
\arguments{
\item{organism}{A vector of taxonomic groups (character).}

\item{speciesSampling}{Whether the results should be at the species-level (logical).}

\item{npar}{Number of simultaneous searches (character; optional). The default 2 is
generally fast enough and does not exceed the maximum number of connections
to genbank.}

\item{nSearchesBatch}{Number of searches per batch}
}
\value{
This function returns a \code{data.frame} that comprises the following
columns. First, \code{Gene}, including the name of the relevant
gene sampled in the target taxonomic groups. Second,
\code{Sampled in N species} includes the number of species where the
gene is sampled among the target taxa. Third, \code{PercentOfSampledSpecies}
indicates the percentage of species where the gene is sampled (assuming
GeneBank's taxonomic backbone).
}
\description{
This function retrieves the gene names for given organism in genbank. This function is
useful to explore the genes that are widely sampled for a given taxon
in genbank. However, note that the performance of \code{gene.sampling.retrieve}
is entirely dependant on the quality of submissions to genbank. For instance, if
most of the sequences in genbank (nucleotide) don't have information on the
sequenced region, this function won't provide reliable estimates of the sampling
frequency of each gene in the database. We recommend using this function with caution,
and only when there is no additional information on the genes that are widely sampled
for the target group(s) of interest.
}
\examples{
\dontrun{
test.spp <- gene.sampling.retrieve(organism = "Puma", speciesSampling = TRUE)
test.pop <- gene.sampling.retrieve(organism = "Puma", speciesSampling = FALSE)
}
}
