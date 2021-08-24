#' Retrieve sequences from genbank
#'
#' Downloads sequences from genbank (nucleotide database) for particular taxa
#' and genes into a folder called \code{"0.Sequences"}.
#'
#' @param clades A vector listing taxonomic groups of interest (character).
#' @param species A vector listing additional species interest (character).
#'                This argument can be used to define additional target species
#'                in the ingroup or species to be sampled in the outgroup
#'                (character).
#' @param genes A vector listing gene names of interest (character).
#' @param maxseqs Maximum number of sequences to retrieve per search
#'                (taxa + gene) (numeric).
#' @param maxlength Maximum lenght of the gene sequence (numeric).
#'
#' @return None
#'
#' @import pbapply
#' @import rentrez
#' @import taxize
#'
#' @examples
#' \dontrun{
#' sq.retrieve(
#'   clades = c("Felis", "Vulpes", "Phoca"),
#'   species = "Manis_pentadactyla",
#'   genes = c("ADORA3", "CYTB")
#' )
#' }
#' @export

sq.retrieve <-
  function(clades = NULL,
           species = NULL,
           genes = NULL,
           maxseqs = 1,
           maxlength = 5000) {
    if (is.null(clades) &
        is.null(species))
      stop("Please provide at least one clade or species")

    if (!is.null(clades)) {
      if (!is.character(clades))
        stop("Please provide at character vector for the clade arguments")
    }
    if (!is.null(species)) {
      if (!is.character(species))
        stop("Please provide at character vector for  the species arguments")
    }

    if (is.null(genes))
      stop("Please provide the name of at least one gene")
    if (!is.numeric(maxseqs) |
        !is.numeric(maxlength))
      stop("maxseqs and maxlength must be numeric")
    if (length(maxseqs) > 1 |
        length(maxlength) > 1)
      stop("Please provide a single number for maxseqs or maxlength")


    taxa <- if (!is.null(clades)) {
      clade.species <-
        downstream(clades,
                   db = "itis",
                   downto = "species",
                   verbose = F)
      clade.species <- do.call(rbind, clade.species)
      c(clade.species$taxonname, species)
    } else {
      species
    }
    unlink("0.Sequences", recursive = TRUE)
    dir.create("0.Sequences")

    singleGene <- function(gene) {
      ret_seqs <- lapply(seq_along(taxa), function(x) {
        tryCatch({
          targetsp <-
            paste0(taxa[x],
                   "[ORGN] AND ",
                   gene,
                   "[Gene] AND 1:",
                   maxlength,
                   "[SLEN]")

          res_rearch <-
            entrez_search(db = "nuccore",
                          term = targetsp,
                          retmax = maxseqs)

          if (length(res_rearch$ids) == 0) {
            return(list(
              taxa = taxa[x],
              gene = gene,
              data = F,
              sequences = NA
            ))
          } else {
            res_seqs <-
              entrez_fetch(db = "nuccore",
                           id = res_rearch$ids,
                           rettype = "fasta")
            return(list(
              taxa = taxa[x],
              gene = gene,
              data = F,
              sequences = res_seqs
            ))
          }
        },
        error = function(e) {
        })
      })

      invisible(lapply(seq_along(taxa), function(x) {
        if (!is.na(ret_seqs[[x]]$sequences)) {
          write(
            ret_seqs[[x]]$sequences,
            paste0("0.Sequences/", gene, ".fasta"),
            sep = "\n",
            append = TRUE
          )
        }
      }))
    }

    invisible(pblapply(genes, singleGene))
  }
