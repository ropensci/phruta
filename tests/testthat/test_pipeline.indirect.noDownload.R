
test_that("Test if the full indirect pipeline works without downloading sequences", {

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

  sqs.downloaded <- sq.retrieve.indirect(acc.table = acc.table,
                                         download.sqs = FALSE)


  expect_true(class(sqs.downloaded) == "list"  )


  tb.merged <- list('COI' = c("cytochrome oxidase subunit 1", "cytochrome c oxidase subunit I"))

  sqs.curated <- sq.curate(filterTaxonomicCriteria = 'Felis|Vulpes|Phoca|Manis',
                           mergeGeneFiles = tb.merged,
                           kingdom = 'animals',
                           sqs.object = sqs.downloaded,
                           removeOutliers = FALSE)

  expect_true( class(sqs.curated) == "list" )

  sqs.aln <- sq.aln(sqs.object = sqs.curated)

  expect_true( class(sqs.aln) == "list" )

}
)



