library(tidyverse)
# distributions = [
#   Gamma(16),
#   LogNormal(0, 0.25),
#   Gamma(4),
#   Gamma(2),
#   LogNormal(0, 0.5),
#   Gamma(1),
#   Exponential(),
#   Gamma(0.64),
#   LogNormal(0, 0.75)
# ]

set.seed(0)

# population skewness as response
# put these data into table, three columns pop skewness, sample size, mean sampling skewness
# each row is a different distro
adjusted_skewness <- function(x) {
  n <- length(x)
  sqrt(n * (n - 1)) / (n - 2) * moments::skewness(x)
}
test <- \(distro, r) replicate(r, adjusted_skewness(eval(distro)))

generate_data <- function(n, r, population_skews) {
  raw_distributions <- list(
    test(expr(rgamma(!!n, 16)), r),
    test(expr(rlnorm(!!n, 0, 0.25)), r),
    test(expr(rgamma(!!n, 4)), r),
    test(expr(rgamma(!!n, 2)), r),
    test(expr(rlnorm(!!n, 0, 0.5)), r),
    test(expr(rgamma(!!n, 1)), r),
    test(expr(rf(!!n, 10, 15)), r),
    test(expr(rgamma(!!n, 0.64)), r),
    test(expr(rlnorm(!!n, 0, 0.75)), r)
  )
  mss <- map_dbl(raw_distributions, mean)

  data.frame(
    pop_skewness = population_skews,
    sample_size = n,
    mean_sampling_skewness = mss
  )
}

skew_f <- \(d1, d2)
  ((2 * d1 + d2 - 2) * sqrt(8 * (d2 - 4))) /
    ((d2 - 6) * sqrt(d1 * (d1 + d2 - 2)))
skew_lnorm <- \(sigma) (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
population_skews <- c(
  0.5,
  skew_lnorm(0.25),
  1,
  sqrt(2),
  skew_lnorm(0.5),
  2,
  skew_f(10, 15),
  2.5,
  skew_lnorm(0.75)
)
ns <- c(seq(10, 50, 10), seq(50, 200, 25)) |> unique()
r <- 1e5

skew_data <- map_df(ns, \(n) generate_data(n, r, population_skews))
# skew_data |> write_csv(here::here("skew_data.csv"))
# skew_data <- read_csv(here::here("skew_data"))
skew_data |>
  ggplot(
    aes(
      x = pop_skewness,
      y = mean_sampling_skewness,
      color = fct(as.character(sample_size))
    )
  ) +
  geom_point() +
  geom_line() +
  geom_abline(slope = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Population Skewness vs. Mean Sampling Skewness with Regression Lines",
    x = "Population Skewness",
    y = "Mean Sampling Skewness",
    color = "Sample Size"
  ) +
  theme_minimal()

skew_data |>
  mutate(percent = mean_sampling_skewness / pop_skewness) |>
  ggplot(
    aes(x = sample_size, y = percent, color = fct(as.character(pop_skewness)))
  ) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1) +
  # geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  geom_text(
    data = skew_data |>
      mutate(percent = mean_sampling_skewness / pop_skewness) |>
      group_by(pop_skewness) |>
      slice_tail(n = 1),
    aes(label = round(percent, 2), x = sample_size + 0.5),
    show.legend = FALSE,
    hjust = 0,
    size = 3,
    fontface = "bold"
  )

# TODO: regress pop skewness on mss and n
# TODO: look at each predictor separately
model <- lm(pop_skewness ~ sample_size + mean_sampling_skewness, skew_data)
summary(model)
# fix normal distro
# sampling distribution of means with n=10
# look at percent error in tail probs like in paper

# fix expo distro
# sampling distro with again n=30
# look at percent error in tail probs
