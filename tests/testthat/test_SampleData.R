unlink(list.dirs("."), recursive = TRUE)


test_that("test if the sample phylogenies load", {
  data("SW.phruta")
  expect_equal(length(SW.phruta), 5)
})

test_that("Check that the class of sample phylogenies is multiPhylo", {
  data("SW.phruta")
  expect_true(class(SW.phruta) == "multiPhylo")
})

unlink(list.dirs("."), recursive = TRUE)
