library(tidyverse)
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
  sampling_distributions <- list(
    sampling_distribution(expr(rgamma(!!n, 16)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.25)), r),
    sampling_distribution(expr(rgamma(!!n, 4)), r),
    sampling_distribution(expr(rgamma(!!n, 2)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.5)), r),
    sampling_distribution(expr(rgamma(!!n, 1)), r),
    sampling_distribution(expr(rgamma(!!n, 0.64)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.75)), r)
  )
  distribution_names <- c(
    "Gamma (α = 16, θ = 1)",
    "Log-normal (μ = 0, σ = 0.25)",
    "Gamma (α = 4, θ = 1)",
    "Gamma (α = 2, θ = 1)",
    "Log-normal (μ = 0, σ = 0.5)",
    "Gamma (α = 1, θ = 1)",
    "Gamma (α = 0.64, θ = 1)",
    "Log-normal (μ = 0, σ = 0.75)"
  )

  data.frame(
    distribution = distribution_names,
    pop_skewness = population_skews,
    sample_size = n,
    mean_sampling_skewness = map_dbl(sampling_distributions, mean),
    median_sampling_skewness = map_dbl(sampling_distributions, median),
    lower_quantile = map_dbl(
      sampling_distributions,
      \(x) quantile(abs(x), 0.025)
    ),
    upper_quantile = map_dbl(
      sampling_distributions,
      \(x) quantile(abs(x), 0.975)
    )
  )
}

skew_lnorm <- \(sigma) (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
population_skews <- c(
  0.5,
  skew_lnorm(0.25),
  1,
  sqrt(2),
  skew_lnorm(0.5),
  2,
  2.5,
  skew_lnorm(0.75)
)
ns <- c(seq(10, 50, 10), seq(50, 200, 25)) |> unique()
r <- 1e5

skew_data <- map_df(ns, \(n) generate_data(n, r, population_skews)) |>
  write_csv(here::here("skew_data.csv"))
skew_data <- read_csv(here::here("skew_data.csv"))
# skew_data |> mutate(Median = (lower_quantile + upper_quantile) / 2) |> View()

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
  theme_bw()

skew_data |>
  mutate(
    percent = mean_sampling_skewness / pop_skewness,
    distribution = fct(
      distribution,
      levels = unique(distribution[order(pop_skewness)])
    )
  ) |>
  ggplot(
    aes(x = sample_size, y = percent, color = distribution)
  ) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1) +
  theme_bw() +
  ggrepel::geom_label_repel(
    data = skew_data |>
      mutate(percent = mean_sampling_skewness / pop_skewness) |>
      group_by(pop_skewness) |>
      slice_tail(n = 1),
    aes(label = round(percent, 2), x = sample_size + 0.5),
    show.legend = FALSE,
    size = 3.5,
    fontface = "bold",
    segment.color = "grey",
    min.segment.length = 0,
    force = 2,
    nudge_x = 10,
    direction = "y",
    max.overlaps = 17
  ) +
  labs(
    title = "Convergence Rate of Sample Skewness to Population Skewness",
    x = "Sample Size",
    y = "Percent of Population Skewness",
    color = "Distribution"
  )

ggsave(
  here::here("Figures", "skewness_convergence.png"),
  width = 12,
  height = 6.75,
  dpi = 1000
)

exp_data <- skew_data |>
  filter(distribution == "Gamma (α = 1, θ = 1)") |>
  mutate(
    log_correction = log(pop_skewness / mean_sampling_skewness - 1),
    log_sample_size = log(sample_size)
  )
exp_data |>
  ggplot(aes(x = log_sample_size, y = log_correction)) +
  geom_point() +
  theme_bw() +
  geom_smooth(se = FALSE) +
  labs(
    x = "Log Sample Size",
    y = "Log Correction Factor", # (log(pop_skewness/mean_sampling_skewness - 1)
    title = "Linear Relationship Between Log Sample Size and Log Correction Factor"
  )

model <- lm(log_correction ~ log_sample_size, data = exp_data)
summary(model)

coeffs <- coef(model)
exp_corrected_skewness <- function(sample_skew, n) {
  unname(sample_skew * (exp(coeffs[1] + coeffs[2] * log(n)) + 1))
}

skew_data |>
  filter(pop_skewness < 3) |>
  mutate(
    corrected_skew = exp_corrected_skewness(
      mean_sampling_skewness,
      sample_size
    ),
    percent = corrected_skew / pop_skewness,
    distribution = fct(
      distribution,
      levels = unique(distribution[order(pop_skewness)])
    )
  ) |>
  ggplot(aes(x = sample_size, y = percent, color = distribution)) +
  geom_point() +
  geom_hline(yintercept = 1)

# REMOVE FROM HERE ON
perskewgamma <- lapply(
  ns,
  \(n) sampling_distribution(expr(rgamma(!!n, 4)), r)
) |>
  map_dbl(\(x) ecdf(x)(1))

plot(ns, perskewgamma)

perskewlnorm <- lapply(
  ns,
  \(n) sampling_distribution(expr(rlnorm(!!n, 0, 0.5)), r)
) |>
  map_dbl(\(x) ecdf(x)(1.75))

plot(ns, perskewlnorm)

modele <- lm(c(perskewgamma, perskewlnorm) ~ rep(ns, 2))

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

# TODO: look at abs skewness and see if it converges to pop skewness
# TODO: bootstrap sample skewness and look at CI; see how often wrong
# TODO: Try KS

theta <- adjusted_skewness(x)
y <- replicate(1e5, adjusted_skewness(sample(x, 10, TRUE)))
thetas <- quantile(y, c(0.025, 0.975))
c(2 * theta - thetas[2], 2 * theta - thetas[1])

library(boot)

boot_skew <- \(data, indices) adjusted_skewness(data[indices])

r <- 1e4
n <- 10
# x <- rexp(n)

# withr::with_options(
#   options(boot.ncpus = 10),
#   upperbounds <- replicate(
#     r / 10,
#     boot.ci(
#       boot(x, boot_skew, r, parallel = "multicore"),
#       conf = 0.6,
#       type = "bca"
#     )$bca[5]
#   ) |>
#     as.data.frame()
# )

conf_levels <- seq(0.4, 0.95, 0.05)
results_df <- data.frame(matrix(ncol = length(conf_levels), nrow = r / 10))
colnames(results_df) <- paste0("conf_", conf_levels)

for (i in seq_along(conf_levels)) {
  conf_level <- conf_levels[i]

  withr::with_options(
    options(boot.ncpus = 10),
    results_df[, i] <- replicate(
      r / 10,
      {
        x <- rexp(n)
        boot.ci(
          boot(x, boot_skew, r, parallel = "multicore"),
          conf = conf_level,
          type = "bca"
        )$bca[5]
      }
    )
  )
}

results_df |> write_csv(here::here("upperbounds.csv"))

results_df |> summarise(across(everything(), \(x) mean(x)))
