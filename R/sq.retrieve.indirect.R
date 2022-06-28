#' Retrieve sequences from genbank based on a dataset of accession numbers
#'
#' Downloads sequences from genbank (nucleotide database) for particular taxa
#' and genes into a folder called \code{"0.Sequences"}.
#'
#' @param acc.table An accession table, ideally generated using \code{acc.table.retrieve}.
#'                  The data.frame must have the Species, Acc, and gene column names.
#' @param download.sqs Logical indicating whether sequences should be downloaded locally or returned as a list.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' acc.table.retrieve(
#'  clades  = c('Felis', 'Vulpes', 'Phoca'),
#'  species = 'Manis_pentadactyla' ,
#'  genes   = c("A2AB","ADORA3","ADRB2","APOB",
#'             "APP","ATP7","BCHE","BDNF",
#'             "BMI1","BRCA1","BRCA2","CNR1",
#'             "COI","CREM","CYTB","DMP1",
#'             "EDG1","ENAM","FBN1","GHR",
#'             "IRBP","ND1","ND2","PLCB4",
#'             "PNOC","RAG1a","RAG1b","RAG2",
#'             "TTN","TYR1","VWF"),
#'  speciesLevel=TRUE
#' )
#'
#' sq.retrieve.indirect(test)
#'
#' }
#' @export



sq.retrieve.indirect <- function(acc.table, download.sqs = FALSE){


  unlink("0.Sequences", recursive = TRUE)
  dir.create("0.Sequences")

 su <- pblapply(unique(acc.table$gene), function(x){
    acc.table.sub <- acc.table[acc.table$gene == x,]
    seqs <- read.GenBank(acc.table.sub$Acc, species.names = TRUE)

    names(seqs) <- paste0(names(seqs), " ", acc.table.sub$Species)
    if(download.sqs){

      ##Over-writing?
      if( !isTRUE(pkg.env$.testMode) ) {
        UI <- readline(paste0("This function might overwrite ",
                              "0.Sequences", ". Are you sure you want to continue? (y/n)  "))
        if(UI != 'y') stop('Exiting since you did not press y')
      }

    write.dna(
      seqs,
      file=paste0("0.Sequences/", x, ".fasta"),
      format = "fasta"
    )
    }else{
      seqs
    }
  })

 if( !download.sqs ){
   names(su) <- unique(acc.table$gene)
   su
 }

}


