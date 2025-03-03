library(tidyverse)
library(flextable)
set.seed(0)

# Figures
df <- read_csv("means.csv") |>
  group_by(Distribution) |>
  mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
  filter(
    `Lower Tail` >= 0.02 & `Lower Tail` <= 0.03,
    `Upper Tail` >= 0.02 & `Upper Tail` <= 0.03
  ) |>
  filter(`Sample Size` == min(`Sample Size`)) |> # get the smallest sample size for each distribution
  select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
  arrange(Skewness, Distribution)

ft <- flextable(df) |>
  colformat_double(j = c("Skewness", "Sampling Skewness"), digits = 3) |>
  font(fontname = "Trebuchet MS") |>
  autofit()

save_as_image(ft, path = here::here("Figures", "table.png"), res = 2000)

graphing <- read_csv("graphing.csv")
distro_plot <- function(col_name, x, y, title) {
  sample_size <- str_extract(col_name, "\\d+")
  upper <- sum(graphing[[col_name]] >= 1.96) / length(graphing[[col_name]])
  lower <- sum(graphing[[col_name]] <= -1.96) / length(graphing[[col_name]])

  sub <- ggplot() +
    geom_line(aes(x, y)) +
    labs(title = title) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

  main <- ggplot(graphing, aes(x = .data[[col_name]])) +
    geom_histogram(bins = 100, aes(y = after_stat(density))) +
    stat_function(
      fun = dnorm,
      args = list(
        mean = mean(graphing[[col_name]]),
        sd = sd(graphing[[col_name]])
      ),
      color = "blue",
      linewidth = 1.25
    ) +
    geom_vline(xintercept = 1.96, linetype = "dashed", color = "blue") +
    annotate(
      "text",
      x = 2.5,
      y = 0.1,
      label = "1.96 SDs",
      color = "blue",
      size = 5
    ) +
    annotate(
      "text",
      x = -2.5,
      y = 0.1,
      label = "-1.96 SDs",
      color = "blue",
      size = 5
    ) +
    annotate(
      "text",
      x = 3,
      y = 0.05,
      label = round(upper, 4),
      color = "black",
      size = 5
    ) +
    annotate(
      "text",
      x = -3,
      y = 0.05,
      label = round(lower, 4),
      color = "black",
      size = 5
    ) +
    geom_vline(xintercept = -1.96, linetype = "dashed", color = "blue") +
    labs(
      title = paste(
        "Sampling Distribution of Standardized Means for n =",
        sample_size
      ),
      x = "Z Score",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 12)
    )

  combined <- main + patchwork::inset_element(sub, 0.6, 0.6, 1, 1)

  ggsave(
    here::here(
      "Figures",
      paste0(
        str_replace(str_sub(col_name, 1, -10), " ", "_"),
        ".png"
      )
    ),
    plot = combined,
    width = 12,
    height = 6.75,
    dpi = 200
  )
}

x <- seq(0, 5, length.out = 100)
y <- dexp(x)
c("Exponential 30 Z-Scores", "Exponential 150 Z-Scores") |>
  walk(
    distro_plot,
    x = x,
    y = y,
    title = "Exponential Population with Rate = 1"
  )

x <- seq(-5, 5, length.out = 100)
y <- dnorm(x)
c("Normal 10 Z-Scores") |>
  walk(
    distro_plot,
    x = x,
    y = y,
    title = "Standard Normal Population"
  )

read_csv("means.csv") |>
  filter(`Sample Size` <= 500) |>
  mutate(
    Distribution = paste(
      str_replace(Distribution, "\\{.*\\}", " "),
      round(Skewness, 2)
    ),
    Distribution = fct(
      Distribution,
      levels = unique(Distribution[order(-Skewness)])
    )
  ) |>
  arrange(Skewness) |>
  ggplot(aes(x = `Sample Size`, color = Distribution)) +
  geom_rect(
    aes(xmin = 0, xmax = Inf, ymin = 0.02, ymax = 0.03),
    fill = "grey",
    linewidth = 0
  ) +
  geom_hline(yintercept = 0.025, linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = `Upper Tail`), linewidth = 1) +
  geom_line(aes(y = `Lower Tail`), linewidth = 1) +
  geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.75) +
  annotate(
    "text",
    x = 40,
    y = 0.04,
    label = "n = 30",
    hjust = 0,
    size = 5,
    color = "black"
  ) +
  labs(
    title = "Upper and Lower Tail Weights of Sampling Distributions",
    x = "Sample Size",
    y = "Tail Weight"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Figures", "Tail_Weights.png"),
  width = 12,
  height = 6.75,
  dpi = 1000
)

# Analysis
model <- lm(Skewness ~ sqrt(`Sample Size`), df)
summary(model)

df |>
  ggplot(aes(x = sqrt(`Sample Size`), y = Skewness, label = Distribution)) +
  geom_abline(
    slope = coef(model)[2],
    intercept = coef(model)[1],
    color = "blue",
    linetype = "dashed",
    linewidth = 1.5
  ) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  ggrepel::geom_label_repel(seed = 0, nudge_x = 1, size = 5) +
  labs(
    title = "Linear Relationship between Skewness and Empirical Minimum Square Root Sample Size",
    x = "Square Root of Sample Size",
    y = "Skewness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 23),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Figures", "Skew_Sample_Size.png"),
  width = 15,
  height = 8.5,
  dpi = 1000
)

n_per_bound <- function(percent) {
  bounds <- 0.025 * percent * c(-1, 1) + 0.025
  temp <- read_csv("means.csv", show_col_types = FALSE) |>
    group_by(Distribution) |>
    mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
    filter(
      `Lower Tail` >= bounds[1] & `Lower Tail` <= bounds[2],
      `Upper Tail` >= bounds[1] & `Upper Tail` <= bounds[2]
    ) |>
    filter(`Sample Size` == min(`Sample Size`)) |> # get the smallest sample size for each distribution
    select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
    arrange(Skewness, Distribution)

  stopifnot(nrow(temp) == 9)

  coef(lm(Skewness ~ sqrt(`Sample Size`), temp))["sqrt(`Sample Size`)"]
}

x <- 30:6 / 100
y <- unname(sapply(x, n_per_bound))

model <- lm(y ~ x)
summary(model)

data.frame(x, y) |>
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Relationship between Normality Error and Beta_sqrt(sample size)",
    x = "Error Rate",
    y = "Coefficient for sqrt(`Sample Size`) in lm"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Figures", "20_percent.png"),
  width = 15,
  height = 8.5,
  dpi = 1000
)

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
    title = "Population Skewness vs Mean Sampling Skewness with Regression Lines",
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
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  ) +
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
    log_correction = log((pop_skewness / mean_sampling_skewness) - 1),
    log_sample_size = log(sample_size)
  )

exp_data |>
  ggplot(aes(x = log_sample_size, y = log_correction)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  labs(
    x = "Log Sample Size",
    y = "Log Correction Factor", # (log(pop_skewness/mean_sampling_skewness - 1)
    title = "Linear Relationship Between Log Sample Size and Log Correction Factor for Exponential Distribution"
  )

ggsave(
  here::here("Figures", "exponential_correction.png"),
  width = 12,
  height = 6.75,
  dpi = 1000
)

model <- lm(log_correction ~ log_sample_size, data = exp_data)
summary(model)

coeffs <- coef(model)
exp_corrected_skewness <- \(skew, n)
  unname(skew * (exp(coeffs[1] + coeffs[2] * log(n)) + 1))

skew_data |>
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
  geom_line() +
  geom_hline(yintercept = 1) +
  theme_bw() +
  labs(
    x = "Sample Size",
    y = "Percentage of Sample Size",
    title = "Corrected Sample Skewnesses",
    color = "Distribution"
  )

ggsave(
  here::here("Figures", "corrected_sample_skewness.png"),
  width = 12,
  height = 6.75,
  dpi = 1000
)
