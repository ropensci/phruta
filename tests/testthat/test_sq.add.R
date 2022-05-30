unlink(list.dirs("."), recursive = TRUE)


test_that("folderDownloaded is null sq.add", {
  expect_error(sq.add(folderDownloaded = NULL))
})

test_that("folderNew is null sq.add", {
  expect_error(sq.add(folderNew = NULL))
})

test_that("sq.add default", {
  expect_silent(sq.add())
})

unlink(list.dirs("."), recursive = TRUE)
