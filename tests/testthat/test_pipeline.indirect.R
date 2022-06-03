unlink(list.dirs("."), recursive = TRUE)


test_that("Test if the full indirect pipeline works", {
  gs.seqs <- gene.sampling.retrieve(organism = "Phoca", speciesSampling = TRUE)
  targetGenes <- gs.seqs[1,]

  acc.table <- acc.table.retrieve(
    clades  = 'Phoca',
    species = 'Manis_pentadactyla' ,
    genes   = targetGenes$Gene,
    speciesLevel = TRUE
  )

  expect_true(
    class(gs.seqs) == 'data.frame'
  )

  expect_true(
    class(targetGenes) == 'data.frame'
  )

  expect_true(
    class(acc.table) == 'data.frame'
  )

  sq.retrieve.indirect(acc.table, download.sqs = TRUE)

  expect_true(any(grepl('0.Sequences',list.dirs())))

  sq.curate(filterTaxonomicCriteria = 'Felis|Vulpes|Phoca|Manis',
            kingdom = 'animals',
            folder = '0.Sequences',
            removeOutliers = FALSE)

  expect_true(any(grepl('1.CuratedSequences',list.dirs())))

  sq.aln(folder = '1.CuratedSequences', FilePatterns = "renamed_")

  expect_true(any(grepl('2.Alignments',list.dirs())))
}
)


unlink(list.dirs("."), recursive = TRUE)



