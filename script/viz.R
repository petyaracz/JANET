d = read_tsv('out/results.tsv')

d$predicted_p_scaled = scale(d$predicted_p)

d |> 
  ggplot(aes(out,predicted_p)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  xlab('p(-k)') +
  ylab('predicted p(-k)')

ggsave('ikpreds.png', dpi = 900, width = 3, height = 3)

glm1 = glm(cbind(freq1,freq2) ~ predicted_p_scaled, data = d, family = binomial)

broom::tidy(glm1, conf.int = T)
