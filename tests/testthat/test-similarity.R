test_that("compute_segment_similarity returns a tibble with required columns", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  expect_s3_class(sim, "tbl_df")
  expect_named(sim, c("segment1", "segment2", "similarity"))
})

test_that("self-similarities are all 1.0", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  segs <- toy_features$segment
  for (s in segs) {
    row <- sim[sim$segment1 == s & sim$segment2 == s, ]
    expect_equal(row$similarity, 1.0, label = paste0("sim(", s, ",", s, ")"))
  }
})

test_that("similarity matrix is symmetric", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  segs <- toy_features$segment
  for (i in seq_along(segs)) {
    for (j in seq_along(segs)) {
      s_ij <- sim$similarity[sim$segment1 == segs[i] & sim$segment2 == segs[j]]
      s_ji <- sim$similarity[sim$segment1 == segs[j] & sim$segment2 == segs[i]]
      expect_equal(s_ij, s_ji,
        label = paste0("sim(", segs[i], ",", segs[j], ") == sim(", segs[j], ",", segs[i], ")")
      )
    }
  }
})

test_that("sim(a,b) = 0.25 (1 shared class out of 4)", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  val <- sim$similarity[sim$segment1 == "a" & sim$segment2 == "b"]
  expect_equal(val, 0.25)
})

test_that("sim(a,p) = 0.0 (no shared classes)", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  val <- sim$similarity[sim$segment1 == "a" & sim$segment2 == "p"]
  expect_equal(val, 0.0)
})

test_that("sim(b,p) = 0.25 (1 shared class out of 4)", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  val <- sim$similarity[sim$segment1 == "b" & sim$segment2 == "p"]
  expect_equal(val, 0.25)
})

test_that("gap rows are present with similarity 1.0", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  gap_rows <- sim[sim$segment1 == " " | sim$segment2 == " ", ]
  expect_true(nrow(gap_rows) > 0L)
  expect_true(all(gap_rows$similarity == 1.0))
})

test_that("output has n^2 + 2n + 1 rows for n segments", {
  nc  <- generate_natural_classes(toy_features)
  sim <- compute_segment_similarity(toy_features, nc)
  n   <- nrow(toy_features)
  expect_equal(nrow(sim), n^2 + 2L * n + 1L)
})

test_that("NULL natural_classes triggers internal recomputation", {
  sim1 <- compute_segment_similarity(toy_features, NULL)
  nc   <- generate_natural_classes(toy_features)
  sim2 <- compute_segment_similarity(toy_features, nc)
  expect_equal(sim1, sim2)
})
