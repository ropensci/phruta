#' Retrieve accession numbers and titles for a given search in genbank
#'
#' Retrieve accession numbers and titles for searches in genbank. This function
#' is useful for exploring the gene-, population-, and species-level sampling
#' in genbank.
#'
#' @param organism The name of a single taxon (character).
#' @param gene The name of a single gene region (character; optional).
#' @param speciesLevel Whether the result should be a species-level dataset (logical).
#'
#' @return data.frame
#'
#' @import reutils
#' @import foreach
#' @import doParallel
#' @import doSNOW
#'
#' @examples
#' \dontrun{
#' acc.retrieve(organism="Vulpes", gene = 'cytb', species=TRUE)
#' }
#' @export

acc.retrieve <- function(organism, gene=NULL, speciesLevel=FALSE, npar=2){

  if( is.null(gene) & !isTRUE(speciesLevel) ){stop("\nPlease provide the name of a gene region or disable the species-level filtering")}

  get_gene = function(x, search, nObs){
    tryCatch({

      recs_summ <- if(nObs==1){
       reutils::efetch(search,
                                   rettype = "docsum",
                                   retmode = "text")
      }else{
      reutils::efetch(search,
                                   rettype = "docsum",
                                   retmode = "text",
                                   retstart = x, retmax=by)
      }
      doc <- reutils::content(recs_summ)
      xml_data <- XML::xmlToList(doc)
      do.call(rbind.data.frame, lapply(xml_data, function(y){
        cbind.data.frame('Ti'=y[[3]][[1]], "Acc"=y[[2]][[1]])
      }))
    }, error=function(e){})
  }

  base.search <- esearch(term = paste0(organism,"[orgn] ", if(!is.null(gene)){paste0("and " ,gene, "[TITL] ")},
                                       "NOT sp NOT unverified NOT genome NOT aff NOT cf"),
                         db = 'nuccore', usehistory = TRUE)
  xml <- content(base.search, "xml")
  count <- as.numeric(XML::xmlToList(xml)$Count)



  if(count>0){
    message("\nSequences found for gene ", gene, " and organism ", organism)

    myCluster <- makeCluster(npar, type="SOCK")
    registerDoSNOW(myCluster)
    by = 499
    cuts <- seq(1, count, by)
    iterations <- length(cuts)
    invisible(pb <- txtProgressBar(max = iterations, style = 3))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    AccDS <- foreach(x = cuts,
                     .packages = "reutils",
                     .options.snow = opts
                     ,.combine = 'rbind'
    ) %dopar% get_gene(x, search = base.search, nObs=count)

    Species <- sapply(AccDS[,1], function(z) paste(strsplit(z, " ")[[1]][c(1:2)], collapse = " " ))
    AccDS <- cbind.data.frame(Species, AccDS)
    if(!is.null(gene)){AccDS$gene <- gene}
    row.names(AccDS) <- NULL

    if(speciesLevel){
      AccDS[!duplicated(AccDS$Species),]
    }else{
      AccDS
    }

  }else{
    message("\nNo sequences found for gene ", gene, " and organism ", organism)
  }
}


base.search$params$term