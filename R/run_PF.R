
run_PF<-function(folder = "ALI"){
  mainDir <- getwd()
  subDir <- if (folder == "ALI") {
    paste0(mainDir, "/Aliscore")
  }  else if (folder == "alig") {
    paste0(mainDir, "/Aligned")
  } else {
    paste0(mainDir, folder)
  }

  setwd(subDir)
  evobiR::SuperMatrix()
  concatenated<-read.FASTA("concatenated.fasta")

##Use PartitionFinder-------

##Set partitions
  setwd(mainDir)
  if( "PartitionFinder.tar.gz" %in% list.files() == FALSE){
  download.file("https://github.com/brettc/partitionfinder/archive/v1.1.1.zip",'PartitionFinder.tar.gz')
  untar("PartitionFinder.tar.gz")} else {}

partitions<-read.csv(paste0(subDir,"/concatenated.partitions.csv"))
config_file<-readLines("partitionfinder-1.1.1/examples/nucleotide/partition_finder.cfg")

block<-list()
for (i in 1:dim(partitions)[1]){
  block[[i]]<-paste0("GENE_",i," ", "="," ", partitions[i,2], " ", "-", " " , partitions[i,3], ";")
}

conf_fi_mod<-config_file[-c(16:24)]
newconf<-append(conf_fi_mod, unlist(block), after=15)
newconf[9]<-paste("models =", "beast;")
newconf[2]<-paste("alignment =", "concatenated.phy;")

subDir_4 <- "PartitionFinder_Analyses"
dir.create(file.path(mainDir, subDir_4))
setwd(file.path(mainDir, subDir_4))

writeLines(newconf, "partition_finder.cfg")
write.dna(concatenated, file = "concatenated.phy", format = "sequential",nbcol=-1,indent = 0,colsep = "",blocksep = 0)


setwd(mainDir)
tofile<-paste0(getwd(), "/PartitionFinder_Analyses/partition_finder.cfg")
system(paste0("python ./partitionfinder-1.1.1/PartitionFinder.py"," ", tofile), wait=F)
cat("You can use R while PartitionFinder is running!
\n Results are written to PartitionFinder_Analyses folder")
}
