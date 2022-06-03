#' Align sequences
#'
#' Perform multiple sequence alignment on all the fasta files
#' saved in a given folder. Alignment is conducted using
#' \code{"BiocManager::install("DECIPHER")"}. Note that this function first
#' optimizes the direction of the sequences, then aligns
#' using \code{"BDECIPHER"}, and finally masks the resulting
#' alignment (optionally but does it per default). The
#' masking step includes removing common gaps across all the species and removing
#' highly ambiguous positions. The resulting aligned sequences
#' are stored to a new folder
#' \code{"2.Alignments"}.'
#'
#' @param folder Name of the folder where the sequences to align
#'               are stored (character).
#' @param FilePatterns A string that is common to all the target
#'                     files in the relevant folder (character). Note that
#'                     this argument can be set to \code{"NULL"} if no specific
#'                     pattern wants to be analyzed.
#' @param sqs.object A list of sequences generated from sq.curate. Only use if you're
#'                   not interested in download sequences locally.
#' @param mask Removes ambiguous sites (Logical, TRUE or FALSE).
#' @param maxFractionGapsSpecies Maximum fraction of gaps per species (when masked)
#' @param ... Arguments passed to \code{"DECIPHER::AlignSeqs"}.
#'
#' @importFrom methods as
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings readDNAStringSet
#' @import DECIPHER
#' @import msa
#'
#' @return list
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
#' sq.aln(folder = "1.CuratedSequences")
#' }
#' @export

sq.aln <- function(folder = "1.CuratedSequences",
                   FilePatterns = "renamed",
                   sqs.object = NULL,
                   mask = TRUE,
                   maxFractionGapsSpecies = 0.01,
                   ...) {
  if (is.null(folder))
    stop("Folder where curated sequences are saved must be provided")
  if (!is.logical(mask))
    stop("The mask argument must be TRUE or FALSE")

  if(is.null(sqs.object)){

  files <- list.files(folder, FilePatterns)
  files <- sub("renamed_", "", files)
  filesComplete <- list.files(folder, FilePatterns, full.names = TRUE)
  unlink("2.Alignments", recursive = TRUE)
  dir.create("2.Alignments")
  invisible(
    lapply(seq_along(filesComplete), function(x) {
      seqs <- Biostrings::readDNAStringSet(filesComplete[x])
      seqs <- DECIPHER::OrientNucleotides(seqs)
      aligned <- DECIPHER::AlignSeqs(seqs, ...)

      if (mask) {
        alignedNoGaps <- DECIPHER::RemoveGaps(aligned, removeGaps = "common")

        alignedMasked <- DECIPHER::MaskAlignment(alignedNoGaps,
                                       correction = if (length(alignedNoGaps) < 200) {
                                         TRUE
                                         } else {
                                           FALSE
                                           })

        DNAStr <- as(alignedMasked, "DNAStringSet")

        if (max(nchar(as.character(DNAStr))) != 0) {


        ##Species removed while masking the aln
        RemMasking <- !names(aligned) %in% names(DNAStr)

        ##Remove species with not enough data
        NonGaps <- nchar(gsub("-", "", as.character(DNAStr)))
        MaxNumberNonGaps <- as.vector(max(nchar(as.character(DNAStr)))*maxFractionGapsSpecies)

        #if(MaxNumberNonGaps > 100){
        rem <- NonGaps < MaxNumberNonGaps
        SumRem <- cbind.data.frame(Species = names(NonGaps),
                                   NonGaps,
                                   #MaxNumberNonGaps,
                                   removedPerGaps = rem,
                                   removedMasking = RemMasking)
        DNAStr <- DNAStr[!rem]
        DNAStr <- DNAStr[!duplicated(names(DNAStr))]

        #Export everything
        write.csv(SumRem, paste0("2.Alignments/0.Masked.Information_", files[x], ".csv"))
        writeXStringSet(DNAStr,
                        filepath = paste0("2.Alignments/Masked_",
                                          files[x]))
        }else{
          write.csv("Masking failed", paste0("2.Alignments/0.Masked.Information_", files[x], ".csv"))
        }

        }

      writeXStringSet(aligned, filepath = paste0("2.Alignments/Raw_", files[x]))
    })
  )

  } else{
     alns <- lapply(sqs.object$Sequences, function(x) {

        if(!is.null(x)){
        seqs <- Biostrings::DNAStringSet(unlist(lapply(as.character(x$Renamed),paste0, collapse="")))
        seqs <- DECIPHER::OrientNucleotides(seqs)
        aligned <- DECIPHER::AlignSeqs(seqs)

        #if (mask) {
          alignedNoGaps <- DECIPHER::RemoveGaps(aligned, removeGaps = "common")

          alignedMasked <- DECIPHER::MaskAlignment(alignedNoGaps,
                                                   correction = if (length(alignedNoGaps) < 200) {
                                                     TRUE
                                                   } else {
                                                     FALSE
                                                   })

          DNAStr <- as(alignedMasked, "DNAStringSet")

          no.masking <- max(nchar(as.character(DNAStr)))
          if (no.masking != 0) {


            ##Species removed while masking the aln
            RemMasking <- !names(aligned) %in% names(DNAStr)

            ##Remove species with not enough data
            NonGaps <- nchar(gsub("-", "", as.character(DNAStr)))
            MaxNumberNonGaps <- as.vector(max(nchar(as.character(DNAStr)))*maxFractionGapsSpecies)

            #if(MaxNumberNonGaps > 100){
            rem <- NonGaps < MaxNumberNonGaps
            SumRem <- cbind.data.frame(Species = names(NonGaps),
                                       NonGaps,
                                       #MaxNumberNonGaps,
                                       removedPerGaps = rem,
                                       removedMasking = RemMasking)
            DNAStr <- DNAStr[!rem]
            DNAStr <- DNAStr[!duplicated(names(DNAStr))]

          }

        #}

          toRet <- if(no.masking != 0) {
            list("Aln.Original" = as.DNAbin(aligned), "Info.Masked" = SumRem, "Aln.Masked" = as.DNAbin(DNAStr))
          }else{
            list("Aln.Original" = as.DNAbin(aligned))
          }
          return(toRet)
          }
      })
     names(alns) <- names(sqs.object$Sequences)
     return(alns)
  }
}
