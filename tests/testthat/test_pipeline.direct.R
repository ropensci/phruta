unlink(list.dirs("."), recursive = TRUE)


test_that("Direct pipeline works", {

  expect_output(sq.retrieve.direct(
    clades = c("Felis", "Vulpes", "Phoca"),
    species = "Manis_pentadactyla",
    genes = c("ADORA3", "CYTB")
  ))

  expect_invisible(
    sq.curate(
      filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
      kingdom = "animals", folder = "0.Sequences"
    )
  )

  expect_invisible(
    sq.aln(folder = "1.CuratedSequences")
  )

  newFas <- read.FASTA("0.Sequences/ADORA3.fasta")
  unlink("0.AdditionalSequences", recursive = TRUE)
  dir.create("0.AdditionalSequences")
  write.FASTA(newFas, "0.AdditionalSequences/ADORA3.fasta")

  expect_snapshot_output(
    sq.add(folderDownloaded = "0.Sequences", folderNew = "0.AdditionalSequences")
  )

  taxonomy <- read.csv("1.CuratedSequences/1.Taxonomy.csv")

  expect_true(
    class(getListConstraints(dataset = taxonomy,
                             targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                             byClades = FALSE
    )) == "list"
  )

  expect_true(
    class(getListConstraints(taxonomy,
                             targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                             byClades = TRUE
    )) == "list"
  )

  expect_snapshot_output(
    tree.constraint(
      taxonomy_folder = "1.CuratedSequences",
      targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
      Topology = "((Felis), (Phoca));"
    )
  )

  expect_snapshot_output(
    tree.constraint(
      taxonomy_folder = "1.CuratedSequences",
      targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
      outgroup = "Phoca_largha"
    )
  )

  expect_true(
    class(taxonomy.retrieve(species_names = c("Felis_catus", "PREDICTED:_Vulpes",
                                              "Phoca_largha", "PREDICTED:_Phoca" ,
                                              "PREDICTED:_Manis" , "Felis_silvestris" , "Felis_nigripes"),
                            database = 'gbif', kingdom = 'animals')) == 'data.frame'
  )

})



unlink(list.dirs("."), recursive = TRUE)
