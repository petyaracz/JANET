# Align two words with Needleman-Wunsch DP.
# distance_dict: named list, key = "seg1|seg2", value = distance (double).
# Returns a list: alignment1, alignment2, distances (pipe-joined strings),
# summed_distance (double), length (int).
.align_pair <- function(word1, word2, distance_dict, gap_penalty) {
  s1 <- strsplit(word1, "")[[1L]]
  s2 <- strsplit(word2, "")[[1L]]
  n  <- length(s1)
  m  <- length(s2)

  M <- matrix(0.0, nrow = n + 1L, ncol = m + 1L)
  M[, 1L] <- seq.int(0L, n) * gap_penalty
  M[1L, ] <- seq.int(0L, m) * gap_penalty

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      key  <- paste0(s1[i], "|", s2[j])
      d    <- distance_dict[[key]]
      if (is.null(d)) d <- gap_penalty
      M[i + 1L, j + 1L] <- min(
        M[i,        j       ] + d,
        M[i,        j + 1L ] + gap_penalty,
        M[i + 1L,   j      ] + gap_penalty
      )
    }
  }

  aln1  <- character(n + m)
  aln2  <- character(n + m)
  dists <- double(n + m)
  pos   <- 0L
  i <- n
  j <- m
  eps <- 1e-9

  while (i > 0L || j > 0L) {
    pos <- pos + 1L
    if (i > 0L && j > 0L) {
      key <- paste0(s1[i], "|", s2[j])
      d   <- distance_dict[[key]]
      if (is.null(d)) d <- gap_penalty
      if (abs(M[i + 1L, j + 1L] - (M[i, j] + d)) < eps) {
        aln1[pos] <- s1[i]; aln2[pos] <- s2[j]; dists[pos] <- d
        i <- i - 1L; j <- j - 1L; next
      }
    }
    if (i > 0L && (j == 0L || abs(M[i + 1L, j + 1L] - (M[i, j + 1L] + gap_penalty)) < eps)) {
      aln1[pos] <- s1[i]; aln2[pos] <- "-"; dists[pos] <- gap_penalty
      i <- i - 1L; next
    }
    aln1[pos] <- "-"; aln2[pos] <- s2[j]; dists[pos] <- gap_penalty
    j <- j - 1L
  }

  len <- pos
  # Reverse in-place (traceback built backwards)
  aln1  <- rev(aln1[seq_len(len)])
  aln2  <- rev(aln2[seq_len(len)])
  dists <- rev(dists[seq_len(len)])

  list(
    alignment1      = paste(aln1,  collapse = "|"),
    alignment2      = paste(aln2,  collapse = "|"),
    distances       = paste(round(dists, 6L), collapse = "|"),
    summed_distance = M[n + 1L, m + 1L],
    length          = len
  )
}

# Align all unique word pairs and return two tibbles: distances and alignment log.
.align_all_pairs <- function(words, similarity_matrix, gap_penalty) {
  words <- unique(words)
  n     <- length(words)

  # Build distance dictionary from similarity matrix (excluding gap rows)
  seg_rows <- similarity_matrix[["segment1"]] != " " & similarity_matrix[["segment2"]] != " "
  keys <- paste0(
    similarity_matrix[["segment1"]][seg_rows], "|",
    similarity_matrix[["segment2"]][seg_rows]
  )
  vals <- 1.0 - similarity_matrix[["similarity"]][seg_rows]
  distance_dict <- setNames(as.list(vals), keys)

  n_pairs  <- n * (n - 1L) / 2L
  log_list  <- vector("list", n_pairs)
  dist_list <- vector("list", n * n)

  log_idx  <- 1L
  dist_idx <- 1L

  for (w in words) {
    dist_list[[dist_idx]] <- list(word1 = w, word2 = w, phon_dist = 0.0)
    dist_idx <- dist_idx + 1L
  }

  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      res <- .align_pair(words[i], words[j], distance_dict, gap_penalty)

      log_list[[log_idx]] <- list(
        word1           = words[i],
        word2           = words[j],
        alignment1      = res$alignment1,
        alignment2      = res$alignment2,
        distances       = res$distances,
        summed_distance = res$summed_distance
      )
      log_idx <- log_idx + 1L

      dist_list[[dist_idx]] <- list(word1 = words[i], word2 = words[j], phon_dist = res$summed_distance)
      dist_idx <- dist_idx + 1L
      dist_list[[dist_idx]] <- list(word1 = words[j], word2 = words[i], phon_dist = res$summed_distance)
      dist_idx <- dist_idx + 1L
    }
  }

  dist_list <- dist_list[!vapply(dist_list, is.null, logical(1L))]
  log_list  <- log_list[!vapply(log_list,  is.null, logical(1L))]

  list(
    distances = tibble::tibble(
      word1     = vapply(dist_list, `[[`, character(1L), "word1"),
      word2     = vapply(dist_list, `[[`, character(1L), "word2"),
      phon_dist = vapply(dist_list, `[[`, double(1L),    "phon_dist")
    ),
    alignment_log = tibble::tibble(
      word1           = vapply(log_list, `[[`, character(1L), "word1"),
      word2           = vapply(log_list, `[[`, character(1L), "word2"),
      alignment1      = vapply(log_list, `[[`, character(1L), "alignment1"),
      alignment2      = vapply(log_list, `[[`, character(1L), "alignment2"),
      distances       = vapply(log_list, `[[`, character(1L), "distances"),
      summed_distance = vapply(log_list, `[[`, double(1L),    "summed_distance")
    )
  )
}
