#' Retrieve accession numbers and titles for a given combination of species and genes in genbank
#'
#' Retrieve accession numbers and titles for complex searches conducted in genbank.
#' Note that this function depends on \code{acc.retrieve} as it actually uses \code{expand.grid}
#' to find all the relevant organism/gene combinations in genbank.
#'
#' @param clades A vector of clade names (character). Note that this can either be
#'               the name of a clade (e.g. Apis) or the code in NCBI (e.g. txid7459).
#' @param species A vector of species names (logical). Note that this can either be
#'               the name of a clade (e.g. Apis) or the code in NCBI (e.g. txid7459).
#' @param genes A vector of gene names (character; optional).
#' @param speciesLevel Whether the result should be a species-level dataset (logical).
#' @param npar Number of parallel searches (the default is probably the best option).
#' @param nSearchesBatch Number of searches per batch
#'
#' @return data.frame
#'
#' @name acc.table.retrieve
#'
#'
#' @import reutils
#' @import foreach
#' @import doParallel
#' @import doSNOW
#'
#' @examples
#' \dontrun{
#' acc.table.retrieve(
#'  clades  = c('Felis', 'Vulpes', 'Phoca'),
#'  species = 'Manis_pentadactyla' ,
#'  genes   = c("A2AB","ADORA3","ADRB2","APOB",
#'             "APP","ATP7","BCHE","BDNF",
#'             "BMI1","BRCA1","BRCA2","CNR1",
#'             "COI","CREM","CYTB","DMP1",
#'             "EDG1","ENAM","FBN1","GHR",
#'             "IRBP","ND1","ND2","PLCB4",
#'             "PNOC","RAG1a","RAG1b","RAG2",
#'             "TTN","TYR1","VWF"),
#'  speciesLevel = FALSE
#' )
#' }
#' @export


acc.table.retrieve <- function(clades = NULL,
                               species = NULL,
                               genes = NULL,
                               speciesLevel = NULL,
                               npar = 2,
                               nSearchesBatch = 499){

  fullTerms <- expand.grid(c(clades,species), genes, stringsAsFactors = FALSE)
  fullSearch <-
    Map(acc.retrieve,
        organism = fullTerms[,1],
        gene = fullTerms[,2],
        speciesLevel = speciesLevel,
        npar = npar,
        nSearchesBatch = nSearchesBatch)
  fullSearch.df <- do.call(rbind.data.frame, fullSearch)
  row.names(fullSearch.df) <- NULL
  fullSearch.df

}
