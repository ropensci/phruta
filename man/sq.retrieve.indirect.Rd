% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sq.retrieve.indirect.R
\name{sq.retrieve.indirect}
\alias{sq.retrieve.indirect}
\title{Retrieve sequences from genbank based on a dataset of accession numbers}
\usage{
sq.retrieve.indirect(acc.table, download.sqs = FALSE)
}
\arguments{
\item{acc.table}{An accession table, ideally generated using \code{acc.table.retrieve}.
The data.frame must have the Species, Acc, and gene column names.}

\item{download.sqs}{Logical indicating whether sequences should be downloaded locally or returned as a list.}
}
\value{
None
}
\description{
Downloads sequences from genbank (nucleotide database) for particular taxa
and genes into a folder called \code{"0.Sequences"}.
}
\examples{
\dontrun{
acc.table.retrieve(
 clades  = c('Felis', 'Vulpes', 'Phoca'),
 species = 'Manis_pentadactyla' ,
 genes   = c("A2AB","ADORA3","ADRB2","APOB",
            "APP","ATP7","BCHE","BDNF",
            "BMI1","BRCA1","BRCA2","CNR1",
            "COI","CREM","CYTB","DMP1",
            "EDG1","ENAM","FBN1","GHR",
            "IRBP","ND1","ND2","PLCB4",
            "PNOC","RAG1a","RAG1b","RAG2",
            "TTN","TYR1","VWF"),
 speciesLevel=TRUE
)

sq.retrieve.indirect(test)

}
}
