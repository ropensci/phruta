#' Retrieve the distribution of genes for given organism in genbank
#'
#' This function retrieves the gene names for given organism in genbank. This function is
#' useful to explore the genes that are widely sampled for a given taxon
#' in genbank. However, note that the performance of \code{gene.sampling.retrieve}
#' is entirely dependant on the quality of submissions to genbank. For instance, if
#' most of the sequences in genbank (nucleotide) don't have information on the
#' sequenced region, this function won't provide reliable estimates of the sampling
#' frequency of each gene in the database. We recommend using this function with caution,
#' and only when there is no additional information on the genes that are widely sampled
#' for the target group(s) of interest.
#'
#'
#' @param organism A vector of taxonomic groups (character).
#' @param speciesSampling Whether the results should be at the species-level (logical).
#' @param npar Number of simultaneous searches (character; optional). The default 2 is
#'             generally fast enough and does not exceed the maximum number of connections
#'             to genbank.
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
#' test.spp <- gene.sampling.retrieve(organism = "Puma", speciesSampling = TRUE)
#' test.pop <- gene.sampling.retrieve(organism = "Puma", speciesSampling = FALSE)
#' }
#' @export

gene.sampling.retrieve <- function(organism, speciesSampling = TRUE, npar = 2){

  x <- NULL
  get_gene_list <- function(x, search, nObs, speciesSampling = speciesSampling){
    tryCatch({

    recs_summ <- if (nObs == 1) {
      reutils::efetch(search,
                      rettype = "ft",
                      retmode = "text")
    }else{
      reutils::efetch(search,
                      rettype = "ft",
                      retmode = "text",
                      retstart = x,
                      retmax = by)
    }

    doc <- reutils::content(recs_summ, 'text')

    newTxtGen <- unlist(strsplit(doc, split = ">Feature gb"))

    newTxt <- lapply(newTxtGen, function(x) unlist(strsplit(x, split = "\n")))
    newTxt <- lapply(newTxt, function(x) unlist(strsplit(x, split = "\t")))

    #genes
    genes <- do.call(rbind, lapply(newTxt, function(x){
      tryCatch({
       xt <- if (length(x) > 1) {x[[1]]} else{x}
      feat <- sub("|", "", xt , fixed = T)
      feat <- sub("|", "", feat, fixed = T)
      feat <- gsub("\\..*","",feat)
      cbind.data.frame(code = feat, gene = x[which(x == "product") + 1])
      }, error = function(e){})
    }))

    if (speciesSampling) {

    seqs <- ape::read.GenBank(genes$code, species.names = TRUE, chunk.size = 200)
    spp <- attr(seqs, "species")

    spp.genes <- cbind.data.frame(Species = unlist(spp), genes)
    spp.genes[!duplicated(paste0(spp.genes$Species,"_", spp.genes$gene)),]
    }else{
      genes
    }

    }, error = function(e){})
  }

  if (length(organism) > 1) {
    organism.2 <- paste(paste(organism, "[orgn]"), collapse = " OR ")
  }else{
    organism.2 <- paste0(organism,"[orgn] ")
  }

  base.search <- esearch(term = paste0(organism.2,
                                         " NOT sp NOT unverified NOT genome NOT aff NOT cf NOT predicted NOT TSA NOT EST"),
                           db = 'nuccore', usehistory = TRUE, sort = 'relevance')

  xml <- content(base.search, "xml")
  count <- as.numeric(XML::xmlToList(xml)$Count)

  if (count > 0) {

    message("\n Genes identified...")

    sys <- Sys.info()["sysname"]
    by = 499
    cuts <- seq(1, count, by)

    if(sys == "Darwin"){
    cl <- makeCluster(npar, type = "SOCK")
    registerDoSNOW(cl)
    iterations <- length(cuts)
    invisible(pb <- txtProgressBar(max = iterations, style = 3))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    AccDS <- foreach(x = cuts,
                     .packages = c("ape","reutils"),
                     .options.snow = opts
                     ,.combine = 'rbind'
    ) %dopar% get_gene_list(x,
                            search = base.search,
                            nObs = count,
                            speciesSampling = speciesSampling)
    stopCluster(cl)
    }

    if(sys == "Linux"){
      AccDS <- pblapply(cuts, function(x){
        get_gene_list(x,
                      search = base.search,
                      nObs = count,
                      speciesSampling = speciesSampling)
      })

      AccDS <- do.call(rbind, AccDS)
    }

    if(sys == "Windows"){

      AccDS <- pblapply(cuts, function(x){
        get_gene_list(x,
                      search = base.search,
                      nObs = count,
                      speciesSampling = speciesSampling)
      })

      AccDS <- do.call(rbind, AccDS)

    }


    if (speciesSampling) {

    tab01 <- ifelse(table(AccDS$gene, AccDS$Species) > 1,1,table(AccDS$gene, AccDS$Species))

    resTable <- as.data.frame(sort(rowSums(tab01), decreasing = T))
    resTable$PercentOfSampledSpecies <- 100 * resTable[,1] / length(unique(AccDS$Species))
    colnames(resTable)[1] <- "Sampled in N species"
    resTable <- cbind.data.frame(Gene = row.names(resTable),resTable )
    row.names(resTable) <- NULL
    resTable

    }else{
    resTable <- as.data.frame(sort(table(AccDS$gene), decreasing = T))
    resTable$Percent <- 100 * resTable[,2] / sum(resTable[,2])
    colnames(resTable)[1] <- "Gene"
    resTable
    }

  }else{
    message("\nNo genes identified...")
  }

}





