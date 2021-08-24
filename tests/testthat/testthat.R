unlink(list.dirs("."), recursive = TRUE)

# context("Loading sample trees")
test_that("load sample phylogenies", {
  data("SW.phruta")
  expect_equal(length(SW.phruta), 5)
})

test_that("Class of sample phylogenies", {
  data("SW.phruta")
  expect_true(class(SW.phruta) == "multiPhylo")
})

# context("sq.retrieve")

test_that("Error in clades and species", {
  expect_error(sq.retrieve(
    clades = NULL,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("Error in default numeric clades", {
  expect_error(sq.retrieve(
    clades = 1,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})



test_that("Error in default numeric species", {
  expect_error(sq.retrieve(
    clades = NULL,
    species = 1,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("No genes", {
  expect_error(sq.retrieve(
    clades = "Homo",
    species = "Brassica",
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})


test_that("Max sequence as character", {
  expect_error(sq.retrieve(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = "20",
    maxlength = 5000
  ))
})


test_that("Max length as character", {
  expect_error(sq.retrieve(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = 20,
    maxlength = "5000"
  ))
})

test_that("length of maxseqs and maxlength", {
  expect_error(sq.retrieve(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = c(20, 20),
    maxlength = c(5000, 100)
  ))
})


test_that("Silent sq.retrieve", {
  expect_output(sq.retrieve(
    clades = "Psocus",
    species = NULL,
    genes = "HDFJ",
    maxseqs = 1,
    maxlength = 1
  ))
})


# context("tree.dating")

test_that("taxonomyFolder NULL", {
  expect_error(tree.dating(taxonomyFolder = NULL))
})

test_that("phylogenyFolder NULL", {
  expect_error(tree.dating(phylogenyFolder = NULL))
})

# context("sq.curate")

test_that("expect_error default", {
  expect_error(sq.curate())
})

test_that("filterTaxonomicCriteria length >1", {
  expect_error(sq.curate(filterTaxonomicCriteria = c(1, 2)))
})

test_that("folder null", {
  expect_error(sq.curate(filterTaxonomicCriteria = "A", folder = NULL))
})

test_that("Wrong Kingdom", {
  expect_error(sq.curate(filterTaxonomicCriteria = "A", kingdom = "I"))
})

test_that("assuming the function runs", {
  expect_error(sq.curate(filterTaxonomicCriteria = "A", kingdom = "animals"))
})

test_that("assuming the function runs x2", {
  expect_error(sq.curate(filterTaxonomicCriteria = "A", kingdom = "animals"))
})


# context("sq.curate")

test_that("assuming the function runs", {
  expect_invisible(sq.aln(folder = "1.CuratedSequences", FilePatterns = "renamed", mask = T))
})


test_that("NULL folder", {
  expect_error(sq.aln(folder = NULL))
})

test_that("Mask is non-logic", {
  expect_error(sq.aln(mask = "H"))
})

# context("tree.raxml")

test_that("Folder is null tree.raxml", {
  expect_error(tree.raxml(folder = NULL))
})

test_that("raxml_exec is not a character tree.raxml", {
  expect_error(tree.raxml(raxml_exec = 1))
})

test_that("Bootstrap is not a number tree.raxml", {
  expect_error(tree.raxml(Bootstrap = "100"))
})

test_that("Bootstrap is 0 tree.raxml", {
  expect_error(tree.raxml(Bootstrap = 0))
})


# context("sq.add")

test_that("folderDownloaded is null sq.add", {
  expect_error(sq.add(folderDownloaded = NULL))
})

test_that("folderNew is null sq.add", {
  expect_error(sq.add(folderNew = NULL))
})

test_that("sq.add default", {
  expect_silent(sq.add())
})

# Test the pipeline

unlink(list.dirs("."), recursive = TRUE)

test_that("Retrieve sequences", {
  expect_output(sq.retrieve(
    clades = c("Felis", "Vulpes", "Phoca"),
    species = "Manis_pentadactyla",
    genes = c("ADORA3", "CYTB")
  ))
})

# test_that("Curate sequences", {
#  expect_output(
sq.curate(
  filterTaxonomicCriteria = "Felis|Vulpes|Phoca|Manis",
  kingdom = "animals", folder = "0.Sequences"
)
# )
# })

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


##Test constraints

taxonomy <- read.csv("1.CuratedSequences/1.Taxonomy.csv")
test_that("Generate list of constraints, not by clade", {
  expect_true(
    class(getListConstraints(taxonomy,
                             targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                             byClades = F)
    ) == "list")
})


test_that("Generate list of constraints, by clade", {
  expect_true(
    class(getListConstraints(taxonomy,
                       targetColumns = c("kingdom", "phylum", "class", "order", "family", "genus", "species_names"),
                       byClades = T)
  ) == "list")
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

##PartitionFinder

#test_that("Tree constraints ingroup/outgroup", {
#  expect_snapshot_output(
#sq.partitionfinderv1(
#  folderAlignments = "2.Alignments",
#  FilePatterns = "Masked",
#  models = "all"
#)
#)
#})

#Tree rogue

test_that("Error when running tree rogue", {
  expect_error(
tree.roguetaxa(folder = "3.Phylogeny")
)
})


