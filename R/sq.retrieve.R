sq.retrieve <- function(clades=NULL, species=NULL, genes, maxseqs=1, maxlength=5000, export=T) {

  ##Find the species
 taxa <- if(!is.null(clades)){
  clade.species<-downstream(clades, db = 'itis', downto = 'species', verbose = F)
  clade.species<-do.call(rbind,clade.species)
   c(clade.species$taxonname, species)
	}else{
	species
	}

  Full_sequences <- lapply(seq_along(genes), function(y) {
    ret_seqs <- pblapply(seq_along(taxa), function(x) {
      tryCatch({
        targetsp <-
          paste0(taxa[x], '[ORGN] AND ',genes[y], "[Gene] AND 1:" , maxlength, "[SLEN]")

        res_rearch <-
          entrez_search(db = "nuccore",
                        term = targetsp,
                        retmax = maxseqs)

        if (length(res_rearch$ids) == 0) {
          return(list(
            taxa = taxa[x],
            gene = genes[y],
            data = F,
            sequences = NA
          ))

        } else{
          res_seqs <-
            entrez_fetch(db = "nuccore",
                         id = res_rearch$ids,
                         rettype = "fasta")
          return(list(
            taxa = taxa[x],
            gene = genes[y],
            data = F,
            sequences = res_seqs
          ))
        }

      }, error=function(e){})})
  })
  names(Full_sequences) <- genes


  if(export==T){
    unlink("0.Sequences", recursive = TRUE)
    dir.create('0.Sequences')
    lapply(seq_along(genes), function(y) {
      lapply(seq_along(taxa), function(x) {
        if(!is.na( Full_sequences[[y]][[x]]$sequences)){
          write(
            Full_sequences[[y]][[x]]$sequences,
            paste0('0.Sequences/',genes[y],".fasta"),
            sep = "\n",
            append = TRUE
          )
        }
      })
    })

    unsampled_taxa <- lapply(Full_sequences, function(x){
      unsampled<-is.na(unlist(lapply(seq_along(x), function(y) x[[y]][[4]])))
      unlist(lapply(seq_along(x), function(y) x[[y]][[1]]))[unsampled]
    })
    unsampled_taxa<-Reduce(intersect, unsampled_taxa)
    attr(Full_sequences, 'unsampled_taxa') <- unsampled_taxa

    return(Full_sequences)
  }else{
    unsampled_taxa <- lapply(Full_sequences, function(x){
      unsampled<-is.na(unlist(lapply(seq_along(x), function(y) x[[y]][[4]])))
      unlist(lapply(seq_along(x), function(y) x[[y]][[1]]))[unsampled]
    })
    unsampled_taxa<-Reduce(intersect, unsampled_taxa)

    attr(Full_sequences, 'unsampled_taxa') <- unsampled_taxa
    return(Full_sequences)
  }
}
