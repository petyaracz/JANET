# janet

**Joint Alignment and Nonparametric Estimation Toolkit**

[![DOI](https://zenodo.org/badge/1119018640.svg)](https://doi.org/10.5281/zenodo.17980463)

`janet` is an R package for computing phonologically-informed pairwise word
distances. Given a binary phonological feature matrix and a list of words,
it returns segment similarities, natural class definitions, and word-level
phonological distances — the inputs needed for analogical models of
morphological variation.

## How it works

1. **Natural classes** — all non-empty intersections of binary feature cuts are
   enumerated via a BFS closure. Each class is described by the most specific
   feature specification whose members equal that set (Frisch, Pierrehumbert &
   Broe 2004).
2. **Segment similarity** — `sim(s₁, s₂) = shared / (shared + non-shared)`,
   where *shared* counts natural classes containing both segments.
3. **Word alignment** — all unique word pairs are aligned with Needleman-Wunsch
   dynamic programming using `1 − similarity` as the substitution cost
   (Dawdy-Hesterberg & Pierrehumbert 2014).

## Installation

```r
# install.packages("remotes")
remotes::install_github("petyaracz/JANET")
```

## Quick start

```r
library(janet)

# Read your phonological feature matrix (segment col + binary feature cols)
features <- read.delim("features.tsv", na.strings = "")

# Supply a character vector of words
words <- readLines("words.txt")

# Run the full pipeline
result <- run_janet(features, words, gap_penalty = 1.0)
```

`run_janet()` returns a named list of four tibbles:

| Name | Columns | Description |
|------|---------|-------------|
| `word_distances` | `word1`, `word2`, `phon_dist` | Pairwise phonological distances (symmetric, includes self = 0) |
| `natural_classes` | `members`, `definition` | All derived natural classes |
| `segment_similarity` | `segment1`, `segment2`, `similarity` | Pairwise segment similarities |
| `alignment_log` | `word1`, `word2`, `alignment1`, `alignment2`, `distances`, `summed_distance` | Detailed alignment for every pair |

## Example: toy inventory

```r
features <- data.frame(
  segment = c("a", "b", "p"),
  cons    = c(0L, 1L, 1L),
  voice   = c(1L, 1L, 0L)
)

result <- run_janet(features, c("ab", "ba", "pa"))

result$natural_classes
#> # A tibble: 5 × 2
#>   members definition
#>   <chr>   <chr>
#> 1 a       [-cons, +voice]
#> 2 b, p    [+cons]
#> 3 a, b    [+voice]
#> 4 p       [+cons, -voice]
#> 5 b       [+cons, +voice]

result$word_distances
#> # A tibble: 9 × 3
#>   word1 word2 phon_dist
#>   <chr> <chr>     <dbl>
#> 1 ab    ab         0
#> 2 ba    ba         0
#> 3 pa    pa         0
#> 4 ab    ba         1.5
#> ...
```

## Input format

### Feature matrix

A data frame with:
- A `segment` column — one row per phoneme, labels are single Unicode characters.
- One column per binary distinctive feature — values `0`, `1`, or `NA` (= not applicable).

### Word list

A character vector. Each word is a string of characters, every one of which must
appear in the `segment` column of the feature matrix.

> **Digraphs** (e.g. `ng`, `th`) must be mapped to a single character before
> use (e.g. `ng` → `N`).

## Using intermediate outputs

```r
nc  <- generate_natural_classes(features)
sim <- compute_segment_similarity(features, nc)
```

## Hungarian example

```r
vignette("hungarian", package = "janet")
```

Bundled data in `inst/extdata/` includes the Siptár-Törkenczy Hungarian feature
matrix and a 305-word verb-stem list.

## Downstream use

Word distances can feed a kernel ridge regression model for morphological
analogy tasks. See [petyaracz/KRR](https://github.com/petyaracz/KRR).

## Citation

```bibtex
@software{racz_janet_2025,
  author = {Rácz, Péter},
  title  = {JANET: Joint Alignment and Nonparametric Estimation Toolkit},
  year   = {2025},
  url    = {https://github.com/petyaracz/JANET},
  doi    = {10.5281/zenodo.17980463}
}
```

## References

Frisch, S. A., Pierrehumbert, J. B., & Broe, M. B. (2004). Similarity
avoidance and the OCP. *Natural Language & Linguistic Theory*, 22(1), 179–228.

Dawdy-Hesterberg, L. G., & Pierrehumbert, J. B. (2014). Learnability and
generalisation of Arabic broken plural nouns. *Language, Cognition and
Neuroscience*, 29(10), 1268–1282.

Rácz, P., & Lukács, Á. (2024). Variation in the 1sg. indef: More than you
wanted to know. *Acta Linguistica Academica*, 71(1–2), 2–17.

## License

MIT. See [LICENSE](LICENSE).
