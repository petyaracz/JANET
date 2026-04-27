#' Compute pairwise segment similarity from natural classes
#'
#' Builds a similarity matrix over all segment pairs using the natural-class
#' metric of Frisch, Pierrehumbert & Broe (2004):
#'
#' \deqn{similarity(s_1, s_2) = \frac{|shared|}{|shared| + |non\text{-}shared|}}
#'
#' where *shared* counts natural classes containing both segments and
#' *non-shared* counts classes containing exactly one of the two segments.
#' A gap character (`" "`) is added with similarity `1.0` to every segment
#' (used as the Needleman-Wunsch gap baseline).
#'
#' @param features A data frame with a `segment` column and binary (`0`/`1`/`NA`)
#'   feature columns (same format as for [generate_natural_classes()]).
#' @param natural_classes A tibble as returned by [generate_natural_classes()].
#'   If `NULL` (default), natural classes are computed from `features`.
#'
#' @return A [tibble][tibble::tibble] with columns:
#' \describe{
#'   \item{`segment1`}{First segment.}
#'   \item{`segment2`}{Second segment.}
#'   \item{`similarity`}{Similarity score in \eqn{[0, 1]}.}
#' }
#'
#' @examples
#' features <- data.frame(
#'   segment = c("a", "b", "p"),
#'   cons    = c(0L, 1L, 1L),
#'   voice   = c(1L, 1L, 0L)
#' )
#' nc <- generate_natural_classes(features)
#' compute_segment_similarity(features, nc)
#'
#' @export
compute_segment_similarity <- function(features, natural_classes = NULL) {
  if (!is.data.frame(features)) stop("`features` must be a data frame.")
  if (!"segment" %in% names(features)) stop("`features` must have a 'segment' column.")

  if (is.null(natural_classes)) {
    natural_classes <- generate_natural_classes(features)
  }

  segments <- features[["segment"]]
  n <- length(segments)

  # Parse class membership once: list of character vectors
  class_members <- strsplit(natural_classes[["members"]], ", ")
  n_classes <- length(class_members)

  # Precompute membership matrix: n_classes x n_segments logical matrix
  member_mat <- matrix(FALSE, nrow = n_classes, ncol = n)
  for (k in seq_len(n_classes)) {
    member_mat[k, ] <- segments %in% class_members[[k]]
  }

  # Segment-pair similarity
  n_pairs <- n * n
  s1_out  <- character(n_pairs + 2L * n + 1L)
  s2_out  <- character(n_pairs + 2L * n + 1L)
  sim_out <- double(n_pairs + 2L * n + 1L)

  idx <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      in_i <- member_mat[, i]
      in_j <- member_mat[, j]
      shared     <- sum(in_i & in_j)
      nonshared  <- sum(xor(in_i, in_j))
      denom      <- shared + nonshared
      s1_out[idx]  <- segments[i]
      s2_out[idx]  <- segments[j]
      sim_out[idx] <- if (denom == 0L) 0.0 else shared / denom
      idx <- idx + 1L
    }
  }

  # Gap rows: sim(seg, " ") = sim(" ", seg) = 1.0, sim(" ", " ") = 1.0
  for (s in segments) {
    s1_out[idx]  <- s;   s2_out[idx]  <- " "; sim_out[idx] <- 1.0; idx <- idx + 1L
    s1_out[idx]  <- " "; s2_out[idx]  <- s;   sim_out[idx] <- 1.0; idx <- idx + 1L
  }
  s1_out[idx]  <- " "; s2_out[idx]  <- " "; sim_out[idx] <- 1.0

  tibble::tibble(segment1 = s1_out, segment2 = s2_out, similarity = sim_out)
}
