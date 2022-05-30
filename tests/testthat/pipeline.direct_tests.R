unlink(list.dirs("."), recursive = TRUE)


test_that("Retrieve sequences", {
  expect_output(sq.retrieve.direct(
    clades = c("Felis", "Vulpes", "Phoca"),
    species = "Manis_pentadactyla",
    genes = c("ADORA3", "CYTB")
  ))
})

test_that("Curate sequences", {
  expect_invisible(
    sq.curate(
      filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
      kingdom = "animals", folder = "0.Sequences"
    )
  )
})

test_that("Align sequences new", {
  expect_invisible(
    sq.aln(folder = "1.CuratedSequences")
  )
})

## Sequence add

newFas <- read.FASTA("0.Sequences/ADORA3.fasta")
unlink("0.AdditionalSequences", recursive = TRUE)
dir.create("0.AdditionalSequences")
write.FASTA(newFas, "0.AdditionalSequences/ADORA3.fasta")

test_that("Add sequences", {
  expect_snapshot_output(
    sq.add(folderDownloaded = "0.Sequences", folderNew = "0.AdditionalSequences")
  )
})


## Test constraints

taxonomy <- read.csv("1.CuratedSequences/1.Taxonomy.csv")
test_that("Generate list of constraints, not by clade", {
  expect_true(
    class(getListConstraints(dataset = taxonomy,
                             targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                             byClades = FALSE
    )) == "list"
  )
})


test_that("Generate list of constraints, by clade", {
  expect_true(
    class(getListConstraints(taxonomy,
                             targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                             byClades = TRUE
    )) == "list"
  )
})


test_that("Tree constraints non ingroup/outgroup", {
  expect_snapshot_output(
    tree.constraint(
      taxonomy_folder = "1.CuratedSequences",
      targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
      Topology = "((Felis), (Phoca));"
    )
  )
})


test_that("Tree constraints ingroup/outgroup", {
  expect_snapshot_output(
    tree.constraint(
      taxonomy_folder = "1.CuratedSequences",
      targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
      outgroup = "Phoca_largha"
    )
  )
})


test_that("Curate sequences", {
  expect_true(
    class(taxonomy.retrieve(species_names = c("Felis_catus", "PREDICTED:_Vulpes",
                                              "Phoca_largha", "PREDICTED:_Phoca" ,
                                              "PREDICTED:_Manis" , "Felis_silvestris" , "Felis_nigripes"),
                            database = 'gbif', kingdom = 'animals')) == 'data.frame'
  )
})


unlink(list.dirs("."), recursive = TRUE)
