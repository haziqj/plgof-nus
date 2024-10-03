library(lavaan.bingof)
library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)
load("R/res_srs_type1.rda")
load("R/res_complex_type1.rda")
load("R/res_complex_ignorewt.RData")

p_type1 <-
  bind_rows(
    mutate(res_srs_type1, sampling = "SRS"),
    filter(res_complex_type1, sampling != "Stratified"),
    res_complex_type1_nowt
  ) |>
  select(name, sim, n, rej_rate = rej_rate5, crit = crit5, sampling) |>
  filter(
    name %in% c("Wald", "WaldVCF", "WaldDiag,MM3", "Pearson,MM3"),
    n == 5000
  ) |>
  mutate(
    name = fct_relabel(name, \(x) gsub(",MM3", "", x)),
    sampling = fct_relevel(sampling, "SRS", "Cluster", "Cluster (ignore wt)", 
                       "Strat-clust", "Strat-clust (ignore wt)"),
    sampling = fct_recode(sampling, "Cluster\n(no wt)" = "Cluster (ignore wt)", 
                           "Strat-clust\n(no wt)" = "Strat-clust (ignore wt)"),
    # name = fct_rev(name)
  ) |>
  ggplot(
    aes(rej_rate, name, col = sampling, shape = sampling, 
             linetype = sampling, group = factor(sampling, 
                                                 levels = rev(c("SRS", 
                                                                "Cluster", 
                                                                "Cluster\n(no wt)",
                                                                "Strat-clust", 
                                                                "Strat-clust\n(no wt)"))))
  ) +
  geom_vline(aes(xintercept = 5 / 100), linetype = "dashed", col = "gray50") +
  geom_pointrange(
    aes(xmin = rej_rate - crit, xmax = rej_rate + crit),
    position = position_dodge(width = 0.5),
    linewidth = 0.8
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  scale_colour_manual(
    values = ggsci::pal_d3()(3)[c(1, 2, 2, 3, 3)],
  ) +
  scale_shape_manual(values = c(19, 19, 4, 19, 4)) +
  scale_linetype_manual(values = c("solid", "solid", "11", "solid", "11")) +
  facet_grid(. ~ sim) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = NULL, x = expression("Rejection rate"), col = "Sampling\ndesign",
       shape = "Sampling\ndesign", linetype = "Sampling\ndesign") +
  guides(col = guide_legend(keyheight = 1.5)
  ); p_type1

save(p_type1, file = "R/p_type1.RData")

