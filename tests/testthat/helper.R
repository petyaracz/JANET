# Minimal 3-segment inventory used across test files.
#
# Natural classes (hand-verified):
#   {a}   [-cons, +voice]
#   {b,p} [+cons]
#   {a,b} [+voice]
#   {p}   [+cons, -voice]
#   {b}   [+cons, +voice]
#
# Similarities:
#   sim(a,a) = 1.0   sim(b,b) = 1.0   sim(p,p) = 1.0
#   sim(a,b) = 0.25  sim(a,p) = 0.0   sim(b,p) = 0.25

toy_features <- data.frame(
  segment = c("a", "b", "p"),
  cons    = c(0L, 1L, 1L),
  voice   = c(1L, 1L, 0L)
)
