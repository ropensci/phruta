unlink(list.dirs("."), recursive = TRUE)

test_that("Generate an error when no genes, species or clades are specified", {
  expect_error(sq.retrieve.direct(
    clades = NULL,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("Generate an error when clades are actually numeric and not strings", {
  expect_error(sq.retrieve.direct(
    clades = 1,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})



test_that("Generate an error when species are actually numeric and not strings", {
  expect_error(sq.retrieve.direct(
    clades = NULL,
    species = 1,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("Generate an error when no genes are specified", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})


test_that("Generate an error when the maximum number of sequences is a string instead of numeric", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = "20",
    maxlength = 5000
  ))
})


test_that("Generate an error when the maximum lenght for the sequences retrieved is a string instead of numeric", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = 20,
    maxlength = "5000"
  ))
})

test_that("Generate an error when a vector with more than 1 element is provided for the max* arguments", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = c(20, 20),
    maxlength = c(5000, 100)
  ))
})


test_that("Silent sq.retrieve.direct", {
  expect_output(sq.retrieve.direct(
    clades = "Psocus",
    species = NULL,
    genes = "HDFJ",
    maxseqs = 1,
    maxlength = 1
  ))
})

unlink(list.dirs("."), recursive = TRUE)
