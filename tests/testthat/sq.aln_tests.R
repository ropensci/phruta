unlink(list.dirs("."), recursive = TRUE)


test_that("assuming the function runs", {
  expect_invisible(sq.aln(folder = "1.CuratedSequences", FilePatterns = "renamed", mask = TRUE))
})


test_that("NULL folder", {
  expect_error(sq.aln(folder = NULL))
})

test_that("Mask is non-logic", {
  expect_error(sq.aln(mask = "H"))
})

unlink(list.dirs("."), recursive = TRUE)

