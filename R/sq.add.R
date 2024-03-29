#' Add local sequences to previously downloaded sequences
#'
#' This function adds sequences located in a particular folder
#' (default is \code{"0.AdditionalSequences"}) to
#' fasta files in another folder (default is \code{"0.Sequences"}).
#' Files must be in FASTA format and names
#' of the files should perfectly match the ones in the previously
#' downloaded folder (e.g. \code{"0.Sequences"}).
#' This function creates a third new folder \code{"0.0.OriginalDownloaded"}
#' containing the originally
#' downloaded sequences. The sequences in the \code{"0.Sequences"} folder get
#' replaced with the combined ones.
#' Note that new sequences can be added even for just a fraction of the
#' originally downloaded genes.
#'
#' @param folderDownloaded Name of the folder with downloaded sequences.
#' @param folderNew Name of the folder with local sequences.
#'
#' @importFrom ape write.FASTA
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' sq.add(folderDownloaded = "0.Sequences",
#'        folderNew = "0.AdditionalSequences")
#' }
#' @export

sq.add <- function(folderDownloaded = "0.Sequences",
                   folderNew = "0.AdditionalSequences") {
  if (is.null(folderDownloaded) | is.null(folderNew))
    stop("Please provide folder names")
  if (!is.character(folderDownloaded) | !is.character(folderNew))
    stop("Folder names must be of class character")

  ##Over-writing?
  if( !isTRUE(pkg.env$.testMode) ) {
  UI <- readline(paste0("This function might overwrite ",
                 folderDownloaded, ". Are you sure you want to continue? (y/n)  "))
  if(UI != 'y') stop('Exiting since you did not press y')
  }


  ds.namesFull <- list.files(folderDownloaded, full.names = TRUE)
  ds.names <- list.files(folderDownloaded)
  ls.namesFull <- list.files(folderNew, full.names = TRUE)
  ls.names <- list.files(folderNew)

  ds.sq <- lapply(ds.namesFull, read.FASTA)
  names(ds.sq) <- ds.names
  ls.sq <- lapply(ls.namesFull, read.FASTA)
  names(ls.sq) <- ls.names

  ls.sq <- lapply(ds.names, function(x) {
    tl <- which(ls.names == x)
    if (length(tl) > 0) {
      ls.sq[[tl]]
    }
  })

  combined.sqs <- lapply(seq_along(ds.sq), function(x) {
    if (is.null(ls.sq[[x]])) {
      ds.sq[[x]]
    } else {
      c(ds.sq[[x]], ls.sq[[x]])
    }
  })

  ## Delete sequences in folderDownloaded
  unlink("0.Sequences", recursive = TRUE)
  ## Write combined sequences in folderDownloaded
  dir.create("0.Sequences")
  invisible(
    lapply(seq_along(combined.sqs), function(x) {
      write.FASTA(combined.sqs[[x]], paste0("0.Sequences/", ds.names[x]))
    })
  )

  ## Create and put original downloaded sequences in new folder
  unlink("0.0.OriginalDownloaded", recursive = TRUE)
  dir.create("0.0.OriginalDownloaded")
  invisible(
    lapply(seq_along(ds.sq), function(x) {
      write.FASTA(ds.sq[[x]], paste0("0.0.OriginalDownloaded/", ds.names[x]))
    })
  )
}
