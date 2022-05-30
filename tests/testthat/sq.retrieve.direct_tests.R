unlink(list.dirs("."), recursive = TRUE)

test_that("Error in clades and species", {
  expect_error(sq.retrieve.direct(
    clades = NULL,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("Error in default numeric clades", {
  expect_error(sq.retrieve.direct(
    clades = 1,
    species = NULL,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})



test_that("Error in default numeric species", {
  expect_error(sq.retrieve.direct(
    clades = NULL,
    species = 1,
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})

test_that("No genes", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = NULL,
    maxseqs = 1,
    maxlength = 5000
  ))
})


test_that("Max sequence as character", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = "20",
    maxlength = 5000
  ))
})


test_that("Max length as character", {
  expect_error(sq.retrieve.direct(
    clades = "Homo",
    species = "Brassica",
    genes = "COI",
    maxseqs = 20,
    maxlength = "5000"
  ))
})

test_that("length of maxseqs and maxlength", {
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
