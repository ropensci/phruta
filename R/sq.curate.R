sq.curate <- function(filterTaxonomicCriteria=NULL, kingdom='animals', folder='0.Sequences'){

  fastaSeqs<-lapply(list.files(folder, full.names = T), read.FASTA)
  names(fastaSeqs) <- list.files(folder, full.names = F)
  seqNames<-sapply(unlist(lapply(fastaSeqs, names)), function(x) paste0(strsplit(x, " ")[[1]][c(2:3)], collapse = '_'))
  seqAccN<-sapply(unlist(lapply(fastaSeqs, names)), function(x) paste0(strsplit(x, " ")[[1]][1], collapse = '_'))
  AccDat<-data.frame('OriginalNames'= unlist(lapply(fastaSeqs, names)),"AccN"=seqAccN, 'Species'= seqNames)
  AccDat$file <- rep(list.files(folder, full.names = F), sapply(fastaSeqs, length))

    species_names<-unique(AccDat$Species)

    gbifkey <- pblapply(species_names, function(x) name_backbone(name = x, kingdom = kingdom))
    keys<-pblapply(1:length(gbifkey), function(x){

      if(as.character(gbifkey[[x]][which(names(gbifkey[[x]]) == "matchType")]) == "NONE"){0}else{

        if( length(which(names(gbifkey[[x]])=="acceptedUsageKey")) == 0 ){
          as.character(gbifkey[[x]][which(names(gbifkey[[x]])=="usageKey")])

        }else{

          as.character(gbifkey[[x]][which(names(gbifkey[[x]])=="acceptedUsageKey")])

        }

      }

    })
    gbif_taxonomy <- pblapply(unlist(keys), function(x) as.data.frame(name_usage(key = x)$data))
    ranks<-c("kingdom", "phylum", "class","order", "family", "genus", "species")
    Taxonomy_species<-pblapply(seq_along(gbif_taxonomy), function(y){
      sub1<-gbif_taxonomy[[y]]
      cate<-t(data.frame( unlist( lapply(seq_along(ranks), function(x){
        nu<- which(colnames(sub1) == ranks[x])
        if(length(nu) != 1){NA}else{sub1[,nu]}
      }))))

      colnames(cate)<-ranks
      row.names(cate)<-NULL
      cate

    })
    Taxonomy_species<-do.call(rbind.data.frame, Taxonomy_species)


    #Remove PREDICTED species
    if( length(grep('PREDICTED',AccDat$Species )) >0  ){
      dupDel<- AccDat[grep('PREDICTED',AccDat$Species ),'OriginalNames']

      fastaSeqs<-lapply(fastaSeqs, function(x){
        x[!names(x) %in% dupDel]
      })

      AccDat<-AccDat[-grep('PREDICTED',AccDat$Species ),]

    }


    #Remove duplicated species
    if(any(duplicated(AccDat[,c(3:4)]))){
      dupDel<- AccDat[duplicated(AccDat[,c(3:4)]),'OriginalNames']

      fastaSeqs<-lapply(fastaSeqs, function(x){
        x[!names(x) %in% dupDel]
      })

      AccDat<-AccDat[!duplicated(AccDat[,c(3:4)]),]

    }


    Full_dataset<-cbind.data.frame(Taxonomy_species, species_names)
    Full_dataset$originalSpeciesName<-Full_dataset$species_names
    Full_dataset$species_names<-sub(" ", "_", Full_dataset$species)

    Full_dataset<-Full_dataset[Full_dataset$species_names %in% AccDat$Species,]


    WrongSpecies<-Full_dataset[!apply(Full_dataset, 1, function(x) any(grepl(filterTaxonomicCriteria, x))), ]
    RightSpecies<-Full_dataset[apply(Full_dataset, 1, function(x) any(grepl(filterTaxonomicCriteria, x))), ]


    if(nrow(WrongSpecies) > 0){
      cat('At least one incorrect sequence was found... /n')

     seqsToDel<- AccDat[AccDat$Species %in% WrongSpecies,'OriginalNames']
     AccDat<- AccDat[!AccDat$Species %in% WrongSpecies,]
     curatedSeqs<-lapply(fastaSeqs, function(x){
        x[!names(x) %in% seqsToDel]
      })
      names(curatedSeqs)<- names(fastaSeqs)

    }else{
      cat('Everything looks right. Moving the sequences to a new folder (1.CuratedSequences) /n')
      curatedSeqs<-fastaSeqs
    }

    ##Rename incorrect synonyms

    toRename<-Full_dataset[sub(" ","_", Full_dataset$species) != Full_dataset$species_names,]


##Export
    unlink("1.CuratedSequences", recursive = TRUE)
    dir.create('1.CuratedSequences')
    lapply(seq_along(curatedSeqs), function(y) {
      ##Original
        write.FASTA(curatedSeqs[[y]], paste0('1.CuratedSequences/',names(curatedSeqs)[y]))
      ##Renamed
        newNames<-sapply(names(curatedSeqs[[y]]), function(x) paste(strsplit(x,' ')[[1]][c(2:3)], collapse = '_' ))
        renamed<-curatedSeqs[[y]]
        if(nrow(toRename)){
        newNames<-ifelse(newNames %in% toRename$species_names, sub(" ","_" ,toRename$species), newNames  )
        }
        names(renamed)<-newNames
        write.FASTA(renamed, paste0('1.CuratedSequences/renamed_',names(curatedSeqs)[y]) )
    })

    write.csv(AccDat, '1.CuratedSequences/AccessionTable.csv')
    write.csv(Full_dataset, '1.CuratedSequences/Taxonomy.csv')

}

