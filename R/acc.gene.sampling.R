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
#'
#' @examples
#' \dontrun{
#' test.spp <- gene.sampling.retrieve(organism = "Puma", speciesSampling=TRUE)
#' test.pop <- gene.sampling.retrieve(organism = "Puma", speciesSampling=FALSE)
#' }
#' @export

gene.sampling.retrieve <- function(organism, npar=2, speciesSampling=TRUE){

  get_gene_list <- function(x, search, nObs, speciesSampling=speciesSampling){

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
      feat <- sub("|", "", x[[1]], fixed = T)
      feat <- sub("|", "", feat, fixed = T)
      feat <- gsub("\\..*","",feat)
      cbind.data.frame(code = feat, gene = x[which(x == "product")+1])
      }, error=function(e){})
    }))


    if(speciesSampling){
    chunks <- split(genes$code, ceiling(seq_along(genes$code)/99))

    spp <- do.call(rbind, lapply(chunks, function(z){
      acc.retrieve(organism = z , acc.num = TRUE, speciesLevel = FALSE)
    }))

    spp.genes <- merge(spp, genes, by.x = "Acc", by.y = "code")


    spp.genes[!duplicated(paste0(spp.genes$Species,"_", spp.genes$gene)),]
    }else{
      genes
    }

  }

  base.search <- esearch(term = paste0(organism,"[orgn] ",
                                         "NOT sp NOT unverified NOT genome NOT aff NOT cf"),
                           db = 'nuccore', usehistory = TRUE, sort = 'relevance')

  xml <- content(base.search, "xml")
  count <- as.numeric(XML::xmlToList(xml)$Count)

  if(count>0){

    message("\n Genes identified for ", organism)

    myCluster <- makeCluster(npar, type="SOCK")
    registerDoSNOW(myCluster)
    by = 499
    cuts <- seq(1, count, by)
    iterations <- length(cuts)
    invisible(pb <- txtProgressBar(max = iterations, style = 3))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    AccDS <- foreach(x = cuts,
                     .packages = c("reutils","doSNOW"),
                     .options.snow = opts
                     ,.combine = 'rbind'
    ) %dopar% get_gene_list(x, search = base.search, nObs=count, speciesSampling=speciesSampling)


    resTable <- as.data.frame(sort(table(AccDS$gene), decreasing = T))
    resTable$Percent <- 100* resTable[,2] / sum(resTable[,2])
    colnames(resTable)[1] <- "Gene"
    resTable

  }else{
    message("\nNo sequences found for gene ", gene, " and organism ", organism)
  }

}





