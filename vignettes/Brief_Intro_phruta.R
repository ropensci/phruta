## ----eval=FALSE---------------------------------------------------------------
#  sq.retrieve(
#              clades  = c('Felis', 'Vulpes', 'Phoca'),
#              species = 'Manis_pentadactyla' ,
#              genes   = c("A2AB","ADORA3","ADRB2","APOB",
#                        "APP","ATP7","BCHE","BDNF",
#                        "BMI1","BRCA1","BRCA2","CNR1",
#                        "COI","CREM","CYTB","DMP1",
#                        "EDG1","ENAM","FBN1","GHR",
#                        "IRBP","ND1","ND2","PLCB4",
#                        "PNOC","RAG1a","RAG1b","RAG2",
#                        "TTN","TYR1","VWF")
#            )

## ----eval=FALSE---------------------------------------------------------------
#  sq.curate(filterTaxonomicCriteria='Felis|Vulpes|Phoca|Manis',
#            kingdom='animals',
#            folder='0.Sequences')

## ----eval=FALSE---------------------------------------------------------------
#  sq.aln(folder='1.CuratedSequences')

## ----eval=FALSE---------------------------------------------------------------
#  tree.raxml(folder='2.Alignments',
#             FilePatterns= 'Masked',
#             raxml_exec='raxmlHPC',
#             Bootstrap=100,
#             outgroup ="Manis_pentadactyla")

## ----eval=FALSE---------------------------------------------------------------
#  tree.dating(taxonomyFolder="1.CuratedSequences",
#              phylogenyFolder="3.Phylogeny",
#              scale='treePL')

## ----eval=FALSE---------------------------------------------------------------
#  system("raxmlHPC")

