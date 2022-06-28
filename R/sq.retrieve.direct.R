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
#' @param db Follows \code{db} in \code{taxize::downsteam}. Choose from
#'           \code{itis}, \code{gbif}, \code{ncbi}, \code{worms}, or \code{bold}.
#' @param maxseqs Maximum number of sequences to retrieve per search
#'                (taxa + gene) (numeric).
#' @param maxlength Maximum length of the gene sequence (numeric).
#'
#' @return None
#'
#' @import pbapply
#' @import rentrez
#' @import taxize
#' @import parallel
#'
#' @examples
#' \dontrun{
#' sq.retrieve.direct(
#'   clades = c("Felis", "Vulpes", "Phoca"),
#'   species = "Manis_pentadactyla",
#'   genes = c("ADORA3", "CYTB")
#' )
#' }
#' @export

sq.retrieve.direct <-
  function(clades = NULL,
           species = NULL,
           genes = NULL,
           db = 'itis',
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


    ##Over-writing?
    if( !isTRUE(pkg.env$.testMode) ) {
      UI <- readline(paste0("This function might overwrite ",
                            "0.Sequences", ". Are you sure you want to continue? (y/n)  "))
      if(UI != 'y') stop('Exiting since you did not press y')
    }


    if (!is.null(clades)) {
      clade.species <-
        taxize::downstream(clades,
                   db = db,
                   downto = "species",
                   verbose = F,
                   rows = 1)


      clade.species <- do.call(rbind, clade.species)

      if(db == 'gbif'){
        taxa <- c(clade.species$name, species)
      }

      if(db == 'itis'){
        taxa <- c(clade.species$taxonname, species)
      }

      if(db == 'ncbi'){
        taxa <-  c(clade.species$childtaxa_name, species)
      }

      if(db == 'bold'){
        taxa <-  c(clade.species$name, species)
      }

      if(db == 'worms'){
        stop("db = worms is not supported yet")
      }

    } else {
      taxa <- species
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
                   "[TI] AND 1:",
                   maxlength,
                   "[SLEN] NOT Predicted NOT UNVERIFIED NOT sp.")

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

    invisible(
      pblapply(genes, function(x){
      tryCatch({
        singleGene(x)
      }, error = function(e){cat('Skipping...', x, 'try again later...')})
    })
    )
}
