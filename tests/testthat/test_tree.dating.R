unlink(list.dirs("."), recursive = TRUE)

test_that("taxonomyFolder NULL", {
  expect_error(tree.dating(taxonomyFolder = NULL))
})

test_that("phylogenyFolder NULL", {
  expect_error(tree.dating(phylogenyFolder = NULL))
})

unlink(list.dirs("."), recursive = TRUE)
