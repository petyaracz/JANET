test_that("run_janet returns a list with four named tibbles", {
  res <- run_janet(toy_features, c("ab", "ba", "pa"))
  expect_type(res, "list")
  expect_named(res, c("word_distances", "natural_classes", "segment_similarity", "alignment_log"))
  expect_s3_class(res$word_distances,     "tbl_df")
  expect_s3_class(res$natural_classes,    "tbl_df")
  expect_s3_class(res$segment_similarity, "tbl_df")
  expect_s3_class(res$alignment_log,      "tbl_df")
})

test_that("word_distances columns are word1, word2, phon_dist", {
  res <- run_janet(toy_features, c("ab", "ba"))
  expect_named(res$word_distances, c("word1", "word2", "phon_dist"))
})

test_that("natural_classes columns are members, definition", {
  res <- run_janet(toy_features, c("ab", "ba"))
  expect_named(res$natural_classes, c("members", "definition"))
})

test_that("segment_similarity columns are segment1, segment2, similarity", {
  res <- run_janet(toy_features, c("ab", "ba"))
  expect_named(res$segment_similarity, c("segment1", "segment2", "similarity"))
})

test_that("alignment_log columns are correct", {
  res <- run_janet(toy_features, c("ab", "ba"))
  expect_named(res$alignment_log,
    c("word1", "word2", "alignment1", "alignment2", "distances", "summed_distance")
  )
})

test_that("word_distances is symmetric and includes self-distances", {
  words <- c("ab", "ba", "pa")
  res   <- run_janet(toy_features, words)
  d     <- res$word_distances
  # n^2 rows
  expect_equal(nrow(d), length(words)^2L)
  # self-distances == 0
  self <- d[d$word1 == d$word2, ]
  expect_true(all(self$phon_dist == 0.0))
  # symmetry
  for (w1 in words) {
    for (w2 in words) {
      d12 <- d$phon_dist[d$word1 == w1 & d$word2 == w2]
      d21 <- d$phon_dist[d$word1 == w2 & d$word2 == w1]
      expect_equal(d12, d21, label = paste0("sym(", w1, ",", w2, ")"))
    }
  }
})

test_that("duplicate words are silently removed", {
  res1 <- run_janet(toy_features, c("ab", "ba"))
  res2 <- run_janet(toy_features, c("ab", "ba", "ab"))
  expect_equal(nrow(res1$word_distances), nrow(res2$word_distances))
})

test_that("gap_penalty parameter is respected", {
  res1 <- run_janet(toy_features, c("a", "ba"), gap_penalty = 1.0)
  res2 <- run_janet(toy_features, c("a", "ba"), gap_penalty = 2.0)
  d1 <- res1$word_distances$phon_dist[
    res1$word_distances$word1 == "a" & res1$word_distances$word2 == "ba"
  ]
  d2 <- res2$word_distances$phon_dist[
    res2$word_distances$word1 == "a" & res2$word_distances$word2 == "ba"
  ]
  expect_true(d2 > d1)
})

test_that("error on missing segment column", {
  expect_error(run_janet(data.frame(x = 1L), "a"), "segment")
})

test_that("error on non-character words", {
  expect_error(run_janet(toy_features, 1:3), "character")
})

test_that("error on empty word list", {
  expect_error(run_janet(toy_features, character(0L)), "empty")
})

test_that("error on non-scalar gap_penalty", {
  expect_error(run_janet(toy_features, "ab", gap_penalty = c(1, 2)), "gap_penalty")
})
