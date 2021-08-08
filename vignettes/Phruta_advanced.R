## ----eval=FALSE---------------------------------------------------------------
#  taxonomy.retrieve(species_names=c("Felis_catus", "PREDICTED:_Vulpes",
#                    "Phoca_largha", "PREDICTED:_Phoca" ,
#                    "PREDICTED:_Manis" , "Felis_silvestris" , "Felis_nigripes"),
#                    database='gbif', kingdom=NULL)

## ----eval=FALSE---------------------------------------------------------------
#  taxonomy.retrieve(species_names=c("Felis_catus", "PREDICTED:_Vulpes",
#                    "Phoca_largha", "PREDICTED:_Phoca" ,
#                    "PREDICTED:_Manis" , "Felis_silvestris" , "Felis_nigripes"),
#                    database='gbif', kingdom='animals')

## ----eval=FALSE---------------------------------------------------------------
#  taxonomy.retrieve(species_names=c("Felis_catus", "PREDICTED:_Vulpes",
#                    "Phoca_largha", "PREDICTED:_Phoca" ,
#                    "PREDICTED:_Manis" , "Felis_silvestris" , "Felis_nigripes"),
#                    database='itis')

## ----eval=FALSE---------------------------------------------------------------
#  tree.constraint(
#                  taxonomy_folder = "1.CuratedSequences",
#                  targetColumns = c("kingdom", "phylum", "class", "order",
#                                    "family", "genus", "species_names"),
#                  Topology = "((ingroup), outgroup);",
#                  outgroup = "Manis_pentadactyla"
#   )

## ----eval=FALSE---------------------------------------------------------------
#  tree.constraint( taxonomy_folder = "1.CuratedSequences",
#                   targetColumns = c("kingdom", "phylum", "class",
#                                     "order", "family", "genus", "species_names"),
#                   Topology = "((Felis), (Phoca));"
#   )

## ----eval=FALSE---------------------------------------------------------------
#  sq.partitionfinder1(folderAlignments = "2.Alignments",
#                      FilePatterns = "Masked",
#                      models = "all"
#   )

## ----eval=FALSE---------------------------------------------------------------
#  tree.raxml(folder = "2.Alignments", FilePatterns = "Masked",
#             raxml_exec = "raxmlHPC", Bootstrap = 100,
#             outgroup = "Manis_pentadactyla",
#             partitioned=T
#  )

## ----eval=FALSE---------------------------------------------------------------
#  tree.roguetaxa(folder = "3.Phylogeny")

