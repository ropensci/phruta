

unlink(list.dirs("."), recursive = TRUE)

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

unlink(list.dirs("."), recursive = TRUE)
