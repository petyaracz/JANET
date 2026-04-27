# -- internal helpers -- #

.basic_cuts <- function(features) {
  feat_names <- setdiff(names(features), "segment")
  cuts <- list()
  for (f in feat_names) {
    col <- features[[f]]
    for (v in c(0L, 1L)) {
      mask <- !is.na(col) & (col == v)
      if (any(mask)) {
        cuts <- c(cuts, list(mask))
      }
    }
  }
  cuts
}

# BFS intersection closure: intersect each discovered set with every basic cut.
# Any n-way intersection is reachable by repeatedly intersecting with one basic
# cut at a time, so we never need to intersect two non-basic-cut sets together.
.intersect_closure <- function(basic_cuts) {
  hash_fn <- function(x) paste(as.integer(x), collapse = "")

  closed <- basic_cuts
  closed_hashes <- vapply(basic_cuts, hash_fn, character(1L))
  queue <- seq_along(basic_cuts)

  while (length(queue) > 0L) {
    ci <- queue[1L]
    queue <- queue[-1L]
    C <- closed[[ci]]
    for (B in basic_cuts) {
      I <- C & B
      if (!any(I)) next
      h <- hash_fn(I)
      if (!h %in% closed_hashes) {
        closed <- c(closed, list(I))
        closed_hashes <- c(closed_hashes, h)
        queue <- c(queue, length(closed))
      }
    }
  }
  closed
}

# Feature spec for a set of segments: include a feature only when every member
# carries the same non-NA value for that feature.
.canonical_spec <- function(mask, features) {
  feat_names <- setdiff(names(features), "segment")
  spec_parts <- character(0L)
  for (f in feat_names) {
    vals <- features[[f]][mask]
    if (!anyNA(vals) && length(unique(vals)) == 1L) {
      sign <- if (vals[1L] == 1) "+" else "-"
      spec_parts <- c(spec_parts, paste0(sign, f))
    }
  }
  if (length(spec_parts) == 0L) return(NULL)
  paste0("[", paste(spec_parts, collapse = ", "), "]")
}

# -- exported function -- #

#' Generate natural classes from a phonological feature matrix
#'
#' Computes all natural classes defined by the intersection closure of binary
#' phonological features, following Frisch, Pierrehumbert & Broe (2004).
#' A natural class is any non-empty set of segments that shares a consistent
#' set of non-missing feature values.
#'
#' @param features A data frame with a `segment` column and binary (`0`/`1`/`NA`)
#'   feature columns. Each row is one segment; each non-`segment` column is a
#'   binary distinctive feature. `NA` means the feature is not applicable to
#'   that segment.
#'
#' @return A [tibble][tibble::tibble] with columns:
#' \describe{
#'   \item{`members`}{Comma-separated, sorted segment names belonging to the class.}
#'   \item{`definition`}{Feature specification string, e.g. `[+cons, -voice]`.}
#' }
#'
#' @examples
#' features <- data.frame(
#'   segment = c("a", "b", "p"),
#'   cons    = c(0L, 1L, 1L),
#'   voice   = c(1L, 1L, 0L)
#' )
#' generate_natural_classes(features)
#'
#' @export
generate_natural_classes <- function(features) {
  if (!is.data.frame(features)) stop("`features` must be a data frame.")
  if (!"segment" %in% names(features)) stop("`features` must have a 'segment' column.")

  segments <- features[["segment"]]
  basic_cuts <- .basic_cuts(features)

  if (length(basic_cuts) == 0L) {
    return(tibble::tibble(members = character(), definition = character()))
  }

  all_masks <- .intersect_closure(basic_cuts)

  results <- lapply(all_masks, function(mask) {
    spec <- .canonical_spec(mask, features)
    if (is.null(spec)) return(NULL)
    list(
      members    = paste(sort(segments[mask]), collapse = ", "),
      definition = spec
    )
  })

  results <- Filter(Negate(is.null), results)

  if (length(results) == 0L) {
    return(tibble::tibble(members = character(), definition = character()))
  }

  tibble::tibble(
    members    = vapply(results, `[[`, character(1L), "members"),
    definition = vapply(results, `[[`, character(1L), "definition")
  )
}
