tree.dating <- function(taxonomyFolder="1.CuratedSequences", phylogenyFolder="3.Phylogeny", ...){
  taxonomy <- read.csv(paste0(taxonomyFolder,"/1.Taxonomy.csv"))
  row.names(taxonomy)<-taxonomy$species_names
  TargetTree<-read.tree(paste0(phylogenyFolder,'/RAxML_bipartitions.phruta'))
  resphy=congruify.phylo(reference=TTOL, target=TargetTree,taxonomy=as.matrix(taxonomy), ...)
  names(resphy) <- c('family', 'order', 'class', 'phyla', 'kingdom')
  resphy<-Filter(Negate(anyNA), resphy)

  #Export
  if(length(resphy)>0){
    unlink("4.Timetree", recursive = TRUE)
    dir.create('4.Timetree')
    invisible(
    lapply(seq_along(resphy), function(x){
     write.tree( resphy[[x]]$phy, paste0("4.Timetree/",names(resphy)[x], '-levelCalibration.tre' ))
     write.csv(resphy[[x]]$calibrations, paste0("4.Timetree/",names(resphy)[x], '-levelCalibration.calibrated.csv' ))
    })
    )
  }
}

