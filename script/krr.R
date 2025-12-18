setwd('~/Github/phonkrr/')

library(tidyverse)
library(broom)

# Kernel ridge regression function
trainKRR = function(sigma, alpha, train_matrix, test_matrix, target){
  
  train_kernel = exp(-train_matrix^2 / (2 * sigma^2))
  test_kernel = exp(-test_matrix^2 / (2 * sigma^2))
  
  # Fit model: (K + αI)^-1 * y
  n = nrow(train_kernel)
  K_reg = train_kernel + alpha * diag(n)
  coefficients = solve(K_reg, target)
  
  # Predict: K_test * coefficients
  nonword_predictions = test_kernel %*% coefficients
  
  return(as.vector(nonword_predictions))
}

# tune model
tuneModel = function(my_sigma,my_alpha){
  predictions = trainKRR(
    sigma = my_sigma,           # adjust based on your distance scale
    alpha = my_alpha,         # adjust based on overfitting
    train_matrix = train_matrix,
    test_matrix = test_matrix,
    target = real_words$out
  )
  
  # Results
  results = tibble(
    lemma = nonwords,
    pred = predictions
  ) |> 
    left_join(d, by = join_by(lemma))
  
  r = with(results, cor(out,pred, method = 'kendall'))
  return(r)  
}

# -- read -- #

d = read_tsv('source/ikverbs.tsv')  # lemma, p
dist = read_tsv('out/aligned_word_pairs_phonological_distance.tsv.gz')  # word1, word2, distance

# -- pin down outcome var -- #

d$out = plogis(d$log_odds)

# -- setup -- #

# Separate real words from nonwords
real_words = d |>
  filter(type == 'corpus') |> 
  arrange(lemma)

nonwords = d |> 
  filter(type == 'nonword') |> 
  arrange(lemma) |> 
  pull(lemma)

# Get ordered word list for training
word_order = real_words$lemma

# Distance matrices
# Real words × real words
train_dist = dist |> 
  filter(word1 %in% word_order, word2 %in% word_order)

# Nonwords × real words
test_dist = dist |> 
  filter(word1 %in% nonwords, word2 %in% word_order)

# Convert to matrices WITH EXPLICIT ORDERING
train_matrix = train_dist |> 
  pivot_wider(names_from = word2, values_from = phon_dist) |>  # check your distance column name!
  arrange(factor(word1, levels = word_order)) |>
  select(word1, all_of(word_order)) |>
  select(-word1) |> 
  as.matrix()

test_matrix = test_dist |> 
  pivot_wider(names_from = word2, values_from = phon_dist) |>
  arrange(word1) |>
  select(word1, all_of(word_order)) |>
  select(-word1) |> 
  as.matrix()

# -- tune model -- #

tuning = crossing(
  my_sigma = c(1, 2, 3, 4, 8, 10, 20),      # bandwidth for your distance scale
  my_alpha = c(0.01, 0.1, 1, 10, 100)  # regularisation strength
)

tuned = tuning |> 
  mutate(
    r = map2_dbl(my_sigma, my_alpha, ~ tuneModel(.x,.y))
  )

tuned_max = filter(tuned, r == max(r))

predictions = trainKRR(
  sigma = tuned_max$my_sigma,           # adjust based on your distance scale
  alpha = tuned_max$my_alpha,         # adjust based on overfitting
  train_matrix = train_matrix,
  test_matrix = test_matrix,
  target = real_words$out # !!!
)

# Results
results = tibble(
  lemma = nonwords,
  predicted_p = predictions
) |> 
  left_join(d)

results$alpha = tuned_max$my_alpha
results$sigma = tuned_max$my_sigma

results |> 
  write_tsv('out/results.tsv')
