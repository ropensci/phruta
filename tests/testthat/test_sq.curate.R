unlink(list.dirs("."), recursive = TRUE)


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

unlink(list.dirs("."), recursive = TRUE)

