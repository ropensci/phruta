#' Align sequences
#'
#' Perform multiple sequence alignment on all the fasta files
#' saved in a given folder. Alignment is conducted using
#' \code{"BiocManager::install("DECIPHER")"}. Note that this function first
#' optimizes the direction of the sequences, then aligns using \code{"BDECIPHER"},
#' and finally masks the resulting alignment (optionally but does it per default). The
#' masking step includes removing common gaps across all the species and removing
#' highly ambiguous positions. The resulting aligned sequences are stored to a new folder
#' \code{"2.Alignments"}.'
#'
#' @param folder Name of the folder where the sequences to align are stored (character).
#' @param FilePatterns A string that is common to all the target files in the relevant folder (character).
#' @param mask Removes ambiguous sites (TRUE or FALSE).
#' @param ... Arguments passed to \code{"DECIPHER::AlignSeqs"}.
#'
#' @import DECIPHER
#' @importFrom methods as
#'
#' @return None
#'
#' @examples
#' \donttest{
#' sq.retrieve(
#'             clades = c('Felis', 'Vulpes', 'Phoca'),
#'             species = 'Manis_pentadactyla' ,
#'             genes = c("ADORA3")
#'            )
#' sq.curate(filterTaxonomicCriteria='Felis|Vulpes|Phoca|Manis',
#'           kingdom='animals', folder='0.Sequences')
#' sq.aln(folder='1.CuratedSequences')
#'}
#' @export

sq.aln<-function(folder='1.CuratedSequences',FilePatterns= 'renamed', mask=T, ...){
  files<-list.files(folder, FilePatterns)
  files<-sub("renamed_",'',files)
  filesComplete<-list.files(folder, FilePatterns, full.names = T)
  unlink("2.Alignments", recursive = TRUE)
  dir.create('2.Alignments')
  invisible(
  lapply(seq_along(filesComplete), function(x){
  seqs <- readDNAStringSet(filesComplete[x])
  seqs <- OrientNucleotides(seqs)
  aligned <- AlignSeqs(seqs, ...)

  if(mask){
    alignedNoGaps<-RemoveGaps(aligned,removeGaps = "common")
    alignedMasked<- MaskAlignment(alignedNoGaps)
    DNAStr = as(alignedMasked, "DNAStringSet")
    writeXStringSet(DNAStr,   file=paste0('2.Alignments/Masked_',files[x]))
  }

  writeXStringSet(aligned,   file=paste0('2.Alignments/',files[x]))
  })
  )
}