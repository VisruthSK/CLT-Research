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
sampling_distribution <- \(distro, r)
  replicate(r, adjusted_skewness(eval(distro)))

generate_data <- function(n, r, population_skews) {
  raw_distributions <- list(
    sampling_distribution(expr(rgamma(!!n, 16)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.25)), r),
    sampling_distribution(expr(rgamma(!!n, 4)), r),
    sampling_distribution(expr(rgamma(!!n, 2)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.5)), r),
    sampling_distribution(expr(rgamma(!!n, 1)), r),
    sampling_distribution(expr(rf(!!n, 10, 15)), r),
    sampling_distribution(expr(rgamma(!!n, 0.64)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.75)), r)
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

# skew_data <- map_df(ns, \(n) generate_data(n, r, population_skews)) |> write_csv(here::here("skew_data.csv"))
skew_data <- read_csv(here::here("skew_data.csv"))
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

model1 <- lm(
  log(pop_skewness) ~ log(sample_size) + log(mean_sampling_skewness),
  skew_data
)
summary(model1)

# fix normal distro
# sampling distribution of means with n=10
# look at percent error in tail probs like in paper

# fix expo distro
# sampling distro with again n=30
# look at percent error in tail probs

n <- 10
r <- 1e6
bounds <- qnorm(0.975) * c(-1, 1)
percent_error <- \(obs, expe) abs((obs - expe) / expe * 100)

graphing_df <- read_csv(here::here("graphing.csv"))
hist(graphing_df$`Exponential 30 Z-Scores`)
hist(graphing_df$`Exponential 150 Z-Scores`)
hist(graphing_df$`Normal 5 Z-Scores`)
hist(graphing_df$`Normal 10 Z-Scores`)
hist(graphing_df$`Normal 30 Z-Scores`)

tail_percents <- function(sample_distro) {
  tail_values <- c(
    mean(sample_distro < bounds[1]),
    mean(sample_distro > bounds[2])
  )

  list(
    tail_values = tail_values,
    percent_errors = percent_error(tail_values, c(0.025, 0.025))
  )
}

tail_percents(graphing_df$`Exponential 30 Z-Scores`)
tail_percents(graphing_df$`Exponential 150 Z-Scores`)
tail_percents(graphing_df$`Normal 5 Z-Scores`)
tail_percents(graphing_df$`Normal 10 Z-Scores`)
tail_percents(graphing_df$`Normal 30 Z-Scores`)

# Edgeworth Expansions of CDF

edgeworth <- function(z, e, lam, eta) {
  U <- exp(-(z^2 / 2))
  V <- -((lam * (z^2 - 1)) / (6 * sqrt(2 * pi)))
  W <- (3 * eta * (z^4 - 6 * z + 3) - lam^2 * z * (z^4 - 10 * z^2 + 15)) /
    (72 * sqrt(2 * pi))

  suppressWarnings(
    s <- c(
      -sqrt(U) * sqrt(U * V^2 - 4 * W * e) + U * V,
      -sqrt(U) * sqrt(U * V^2 - 4 * W * e) - U * V,
      U * V - sqrt(U) * sqrt(U * V^2 + 4 * W * e),
      U * V + sqrt(U) * sqrt(U * V^2 + 4 * W * e)
    ) /
      (2 * e)
  )
  s <- s[!is.na(s)]
  max(s^2)
}
formula <- \(error, lam) 1.303 * error^-2 * lam^2
standardize <- \(data) (data - mean(data)) / sd(data)

z <- qnorm(0.975)
e <- 0.001
lam <- 2 # expo distro
eta <- 6 # expo distro

err <- (e / (1 - pnorm(z))) |> print()
n_edge <- edgeworth(z, e, lam, eta) |> print()
n_formula <- formula(err, lam) |> print()

samp_dist_formula <- replicate(1e5, mean(rexp(n_form, lam))) |> standardize()
samp_dist_edge <- replicate(1e5, mean(rexp(n_edge))) |> standardize()

hist(samp_dist_formula)
hist(samp_dist_edge)

mean(samp_dist_formula >= z)
mean(samp_dist_edge >= z)
