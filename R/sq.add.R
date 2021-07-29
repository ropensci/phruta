sq.add <- function(folderDownloaded='0.Sequences', folderNew='0.AdditionalSequences'){
  ds.namesFull <- list.files(folderDownloaded, full.names = T)
  ds.names <- list.files(folderDownloaded)
  ls.namesFull <- list.files(folderNew, full.names = T)
  ls.names <- list.files(folderNew)

  ds.sq<-lapply(ds.namesFull, read.FASTA)
  names(ds.sq) <- ds.names
  ls.sq<-lapply(ls.namesFull, read.FASTA)
  names(ls.sq)<-ls.names

  ls.sq<- lapply(ds.names, function(x){
    tl<-which(ls.names == x)
    if(length(tl)>0){
      ls.sq[[tl]]
    }
  })

  combined.sqs<-lapply(seq_along(ds.sq), function(x){
    if(is.null(ls.sq[[x]]) ){ ds.sq[[x]]}else{
    c( ds.sq[[x]], ls.sq[[x]])
    }
  } )

  ##Delete sequences in folderDownloaded
  unlink("0.Sequences", recursive = TRUE)
  ##Write combined sequences in folderDownloaded
  dir.create('0.Sequences')
  lapply(seq_along(combined.sqs), function(x){
    write.FASTA(combined.sqs[[x]], paste0('0.Sequences/',ds.names[x]))
  })


  ##Create and put original downloaded sequences in new folder
  dir.create('0.0.OriginalDownloaded')
  lapply(seq_along(ds.sq), function(x){
    write.FASTA(ds.sq[[x]], paste0('0.0.OriginalDownloaded/',ds.names[x]))
  })

}
