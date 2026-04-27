#' Run the JANET phonological distance pipeline
#'
#' Computes pairwise phonological distances between words using natural class
#' theory (Frisch, Pierrehumbert & Broe 2004) and Needleman-Wunsch alignment
#' (Dawdy-Hesterberg & Pierrehumbert 2014).
#'
#' The pipeline has two stages:
#' 1. **Segment similarity** — natural classes are derived from the feature
#'    matrix via intersection closure; each segment pair is scored by the ratio
#'    of shared to total (shared + non-shared) natural classes.
#' 2. **Word alignment** — all unique word pairs are aligned with
#'    Needleman-Wunsch dynamic programming using `1 - similarity` as the
#'    substitution cost.
#'
#' @param features A data frame with a `segment` column and binary (`0`/`1`/`NA`)
#'   feature columns. Each row is one phonological segment; each additional
#'   column is a binary distinctive feature (`NA` = not applicable).
#'   Segment labels must be single Unicode characters and must cover every
#'   character that appears in `words`.
#' @param words A character vector of words. Each word is a string whose
#'   individual characters are phonological segments listed in `features`.
#'   Duplicate words are silently removed.
#' @param gap_penalty Numeric scalar. Cost of inserting or deleting a segment
#'   during alignment. Default `1.0`.
#'
#' @return A named list with four [tibble][tibble::tibble]s:
#' \describe{
#'   \item{`word_distances`}{Pairwise phonological distances in long format:
#'     `word1`, `word2`, `phon_dist`. Symmetric (both directions) and includes
#'     self-distances of `0`.}
#'   \item{`natural_classes`}{All natural classes derived from `features`:
#'     `members` (comma-separated sorted segments), `definition` (feature spec
#'     string, e.g. `[+cons, -voice]`).}
#'   \item{`segment_similarity`}{Pairwise segment similarities:
#'     `segment1`, `segment2`, `similarity` ∈ \eqn{[0,1]}.
#'     Includes gap-character (`" "`) rows with similarity `1.0`.}
#'   \item{`alignment_log`}{Detailed alignment for each ordered word pair:
#'     `word1`, `word2`, `alignment1`, `alignment2` (pipe-separated segment
#'     strings), `distances` (pipe-separated per-position costs),
#'     `summed_distance`.}
#' }
#'
#' @references
#' Frisch, S. A., Pierrehumbert, J. B., & Broe, M. B. (2004). Similarity
#' avoidance and the OCP. *Natural Language & Linguistic Theory*, 22(1), 179–228.
#'
#' Dawdy-Hesterberg, L. G., & Pierrehumbert, J. B. (2014). Learnability and
#' generalisation of Arabic broken plural nouns. *Language, Cognition and
#' Neuroscience*, 29(10), 1268–1282.
#'
#' @examples
#' features <- data.frame(
#'   segment = c("a", "b", "p"),
#'   cons    = c(0L, 1L, 1L),
#'   voice   = c(1L, 1L, 0L)
#' )
#' words <- c("ab", "ba", "pa")
#' result <- run_janet(features, words)
#' result$word_distances
#' result$natural_classes
#'
#' @export
run_janet <- function(features, words, gap_penalty = 1.0) {
  if (!is.data.frame(features))
    stop("`features` must be a data frame.")
  if (!"segment" %in% names(features))
    stop("`features` must have a 'segment' column.")
  if (!is.character(words))
    stop("`words` must be a character vector.")
  if (!is.numeric(gap_penalty) || length(gap_penalty) != 1L || is.na(gap_penalty))
    stop("`gap_penalty` must be a single non-NA number.")

  words <- unique(words[nzchar(words)])
  if (length(words) == 0L) stop("`words` is empty after removing blank strings.")

  nat_classes <- generate_natural_classes(features)
  seg_sim     <- compute_segment_similarity(features, nat_classes)
  aln_res     <- .align_all_pairs(words, seg_sim, gap_penalty)

  list(
    word_distances     = aln_res$distances,
    natural_classes    = nat_classes,
    segment_similarity = seg_sim,
    alignment_log      = aln_res$alignment_log
  )
}
