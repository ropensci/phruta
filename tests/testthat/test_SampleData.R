unlink(list.dirs("."), recursive = TRUE)


test_that("load sample phylogenies", {
  data("SW.phruta")
  expect_equal(length(SW.phruta), 5)
})

test_that("Class of sample phylogenies", {
  data("SW.phruta")
  expect_true(class(SW.phruta) == "multiPhylo")
})

unlink(list.dirs("."), recursive = TRUE)
