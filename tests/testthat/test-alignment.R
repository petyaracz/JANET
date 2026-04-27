# For toy_features: dist(a,b) = dist(b,a) = 0.75, dist(a,p) = 1.0,
# dist(b,p) = 0.75, dist(x,x) = 0 for all x.
#
# Alignment of "ab" vs "ba":
#   DP gives total cost 1.5 via a↔b (0.75) + b↔a (0.75).

get_sim <- function() {
  nc <- generate_natural_classes(toy_features)
  compute_segment_similarity(toy_features, nc)
}

test_that("align_all_pairs returns distances and alignment_log tibbles", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  expect_named(res, c("distances", "alignment_log"))
  expect_s3_class(res$distances,     "tbl_df")
  expect_s3_class(res$alignment_log, "tbl_df")
})

test_that("distances tibble has correct columns", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  expect_named(res$distances, c("word1", "word2", "phon_dist"))
})

test_that("alignment_log tibble has correct columns", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  expect_named(res$alignment_log,
    c("word1", "word2", "alignment1", "alignment2", "distances", "summed_distance")
  )
})

test_that("self-distances are 0", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba", "pa"), sim, 1.0)
  self <- res$distances[res$distances$word1 == res$distances$word2, ]
  expect_true(all(self$phon_dist == 0.0))
})

test_that("distance matrix is symmetric", {
  sim  <- get_sim()
  res  <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  d    <- res$distances
  d12  <- d$phon_dist[d$word1 == "ab" & d$word2 == "ba"]
  d21  <- d$phon_dist[d$word1 == "ba" & d$word2 == "ab"]
  expect_equal(d12, d21)
})

test_that("'ab' vs 'ba' distance is 1.5 (a-b + b-a substitutions)", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  d   <- res$distances
  val <- d$phon_dist[d$word1 == "ab" & d$word2 == "ba"]
  expect_equal(val, 1.5)
})

test_that("identical words have distance 0", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ab"), sim, 1.0)
  val <- res$distances$phon_dist[res$distances$word1 == "ab" & res$distances$word2 == "ab"]
  expect_equal(val, 0.0)
})

test_that("gap_penalty changes alignment cost", {
  sim  <- get_sim()
  # "a" vs "ba": optimal with gap_penalty=1 is gap(a) + a↔a + b-del = 1 + 0 + 1 = 2
  # but a↔b + gap(a) = 0.75 + 1 = 1.75. So dist("a","ba") with gap=1 should be 1.75.
  res1 <- janet:::.align_all_pairs(c("a", "ba"), sim, 1.0)
  res2 <- janet:::.align_all_pairs(c("a", "ba"), sim, 2.0)
  d1 <- res1$distances$phon_dist[res1$distances$word1 == "a" & res1$distances$word2 == "ba"]
  d2 <- res2$distances$phon_dist[res2$distances$word1 == "a" & res2$distances$word2 == "ba"]
  expect_true(d2 > d1)
})

test_that("alignment strings use pipe separator", {
  sim <- get_sim()
  res <- janet:::.align_all_pairs(c("ab", "ba"), sim, 1.0)
  log <- res$alignment_log[res$alignment_log$word1 == "ab" & res$alignment_log$word2 == "ba", ]
  expect_true(grepl("\\|", log$alignment1))
  expect_true(grepl("\\|", log$alignment2))
})

test_that("n words yield n^2 distance rows (including self)", {
  sim <- get_sim()
  words <- c("ab", "ba", "pa")
  res   <- janet:::.align_all_pairs(words, sim, 1.0)
  expect_equal(nrow(res$distances), length(words)^2L)
})
