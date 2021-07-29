tree.raxml<-function(folder='2.Alignments', FilePatterns= 'Masked', raxml_exec='raxmlHPC', Bootstrap=100, outgroup,...){

    files_fullNames<-list.files(folder, FilePatterns, full.names=T)
    files<-list.files(folder, 'Masked')
    seq<- lapply(lapply(files_fullNames,  read.FASTA), as.matrix)
    names(seq)<-files

    concatenated<-do.call(cbind.DNAbin,c(seq,
                                         fill.with.gaps = TRUE))
    unlink("3.Phylogeny", recursive = TRUE)
    dir.create('3.Phylogeny')
    mainDir<-getwd()
    setwd(paste0(mainDir,'/',"3.Phylogeny" ))
    tr <- raxml(DNAbin=concatenated, m = "GTRGAMMA",
                f = "a", N = Bootstrap, p = 1234, x = 1234,
                k=T,
                exec =raxml_exec ,threads=4,
                file='PhyloPipeR',
                outgroup=outgroup, ...)
    setwd(mainDir)

}



