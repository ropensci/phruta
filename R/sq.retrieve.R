sq.retrieve <- function(clades=NULL, species=NULL, genes, maxseqs=1, maxlength=5000) {

  ##Find the species
 taxa <- if(!is.null(clades)){
  clade.species<-downstream(clades, db = 'itis', downto = 'species', verbose = F)
  clade.species<-do.call(rbind,clade.species)
   c(clade.species$taxonname, species)
	}else{
	species
	}
 unlink("0.Sequences", recursive = TRUE)
 dir.create('0.Sequences')

  singleGene <-  function(gene) {
    ret_seqs <- lapply(seq_along(taxa), function(x) {
      tryCatch({
        targetsp <-
          paste0(taxa[x], '[ORGN] AND ',gene, "[Gene] AND 1:" , maxlength, "[SLEN]")

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

        } else{
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

      }, error=function(e){})})


    invisible(
      lapply(seq_along(taxa), function(x) {
        if(!is.na(ret_seqs[[x]]$sequences)){
          write(
            ret_seqs[[x]]$sequences,
            paste0('0.Sequences/',gene,".fasta"),
            sep = "\n",
            append = TRUE
          )
        }
      })
    )


  }

  Full_sequences = pblapply(genes, singleGene)
  names(Full_sequences) <- genes
#
#   unsampled_taxa <- lapply(Full_sequences, function(x){
#       unsampled<-is.na(unlist(lapply(seq_along(x), function(y) x[[y]][[4]])))
#       unlist(lapply(seq_along(x), function(y) x[[y]][[1]]))[unsampled]
#   })
#   unsampled_taxa<-Reduce(intersect, unsampled_taxa)
#   return(unsampled_taxa)
}
