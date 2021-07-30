#BiocManager::install("msa")
#BiocManager::install("DECIPHER")

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
