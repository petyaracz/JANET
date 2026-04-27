test_that("generate_natural_classes returns a tibble with required columns", {
  nc <- generate_natural_classes(toy_features)
  expect_s3_class(nc, "tbl_df")
  expect_named(nc, c("members", "definition"))
})

test_that("toy inventory produces exactly 5 natural classes", {
  nc <- generate_natural_classes(toy_features)
  expect_equal(nrow(nc), 5L)
})

test_that("expected natural classes are present", {
  nc <- generate_natural_classes(toy_features)
  expect_true("a"    %in% nc$members)
  expect_true("b, p" %in% nc$members)
  expect_true("a, b" %in% nc$members)
  expect_true("b"    %in% nc$members)
  expect_true("p"    %in% nc$members)
})

test_that("class definitions use +/- notation with square brackets", {
  nc <- generate_natural_classes(toy_features)
  expect_true(all(grepl("^\\[", nc$definition)))
  expect_true(all(grepl("\\]$", nc$definition)))
  expect_true(any(grepl("\\+cons", nc$definition)))
  expect_true(any(grepl("-cons",  nc$definition)))
})

test_that("{b,p} has definition [+cons] only (voice disagrees)", {
  nc <- generate_natural_classes(toy_features)
  row <- nc[nc$members == "b, p", ]
  expect_equal(row$definition, "[+cons]")
})

test_that("{b} has definition [+cons, +voice]", {
  nc <- generate_natural_classes(toy_features)
  row <- nc[nc$members == "b", ]
  expect_equal(row$definition, "[+cons, +voice]")
})

test_that("single-feature inventory works", {
  f <- data.frame(segment = c("a", "b"), cons = c(0L, 1L))
  nc <- generate_natural_classes(f)
  expect_equal(nrow(nc), 2L)
})

test_that("all-NA feature column is ignored", {
  f <- data.frame(
    segment = c("a", "b"),
    cons    = c(0L, 1L),
    dummy   = c(NA_integer_, NA_integer_)
  )
  nc_plain <- generate_natural_classes(
    data.frame(segment = c("a", "b"), cons = c(0L, 1L))
  )
  nc_dummy <- generate_natural_classes(f)
  expect_equal(nrow(nc_dummy), nrow(nc_plain))
})

test_that("missing segment column raises an error", {
  expect_error(generate_natural_classes(data.frame(x = 1L)), "segment")
})
