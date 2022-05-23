#' Retrieve the distribution of genes for given organism in genbank
#'
#' Retrieve gene names fora given organism in genbank. This function is
#' useful to explore the genes that are widely sampled for a given taxon
#' in genbank.
#'
#' @param organism The name of a single taxon (character).
#' @param npar Number of simultaneous searches (character; optional).
#' @param speciesLevel Whether the results should be at the species-level (logical).
#'
#' @return data.frame
#'
#' @import reutils
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import ape
#'
#' @examples
#' \dontrun{
#' test.spp <- gene.sampling.retrieve(organism = "Puma", speciesSampling=TRUE)
#' test.pop <- gene.sampling.retrieve(organism = "Puma", speciesSampling=FALSE)
#' }
#' @export

gene.sampling.retrieve <- function(organism, speciesSampling=TRUE, npar=2){

  get_gene_list <- function(x, search, nObs, speciesSampling=speciesSampling){
    tryCatch({


    recs_summ <- if(nObs==1){
      reutils::efetch(search,
                      rettype = "ft",
                      retmode = "text")
    }else{
      reutils::efetch(search,
                      rettype = "ft",
                      retmode = "text",
                      retstart = x, retmax=by)
    }

    doc <- reutils::content(recs_summ, 'text')

    newTxtGen <- unlist(strsplit(doc, split = ">Feature gb"))

    newTxt <- lapply(newTxtGen, function(x) unlist(strsplit(x, split = "\n")))
    newTxt <- lapply(newTxt, function(x) unlist(strsplit(x, split = "\t")))

    #genes
    genes <- do.call(rbind, lapply(newTxt, function(x){
      tryCatch({
       xt<- if(length(x)>1){x[[1]]}else{x}
      feat <- sub("|", "", xt , fixed = T)
      feat <- sub("|", "", feat, fixed = T)
      feat <- gsub("\\..*","",feat)
      cbind.data.frame(code = feat, gene = x[which(x == "product")+1])
      }, error=function(e){})
    }))

    if(speciesSampling){

    seqs <- ape::read.GenBank(genes$code, species.names=TRUE, chunk.size=200)
    spp <- attr(seqs, "species")

    spp.genes <- cbind.data.frame(Species=unlist(spp), genes)
    spp.genes[!duplicated(paste0(spp.genes$Species,"_", spp.genes$gene)),]
    }else{
      genes
    }

    }, error=function(e){})
  }

  if(length(organism)>1){
    organism.2 <- paste(paste(organism, "[orgn]"), collapse = " OR ")
  }else{
    organism.2 <- paste0(organism,"[orgn] ")
  }

  base.search <- esearch(term = paste0(organism.2,
                                         " NOT sp NOT unverified NOT genome NOT aff NOT cf NOT predicted NOT TSA NOT EST"),
                           db = 'nuccore', usehistory = TRUE, sort = 'relevance')

  xml <- content(base.search, "xml")
  count <- as.numeric(XML::xmlToList(xml)$Count)

  if(count>0){

    message("\n Genes identified...")

    myCluster <- makeCluster(npar, type="SOCK")
    registerDoSNOW(myCluster)
    by = 499
    cuts <- seq(1, count, by)
    iterations <- length(cuts)
    invisible(pb <- txtProgressBar(max = iterations, style = 3))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    AccDS <- foreach(x = cuts,
                     .packages = c("ape","reutils"),
                     .options.snow = opts
                     ,.combine = 'rbind'
    ) %dopar% get_gene_list(x, search = base.search, nObs=count, speciesSampling=speciesSampling)

    if(speciesSampling){

    tab01 <- ifelse(table(AccDS$gene,AccDS$Species) >1,1,table(AccDS$gene,AccDS$Species))

    resTable <- as.data.frame(sort(rowSums(tab01), decreasing = T))
    resTable$PercentOfSampledSpecies <- 100* resTable[,1] / length(unique(AccDS$Species))
    colnames(resTable)[1] <- "Sampled in N species"
    resTable <- cbind.data.frame(Gene=row.names(resTable),resTable )
    row.names(resTable) <- NULL
    resTable

    }else{
    resTable <- as.data.frame(sort(table(AccDS$gene), decreasing = T))
    resTable$Percent <- 100* resTable[,2] / sum(resTable[,2])
    colnames(resTable)[1] <- "Gene"
    resTable
    }

  }else{
    message("\nNo genes identified...")
  }

}





