library(tidyverse)
library(flextable)
set.seed(0)

poster_plot_dpi <- 400
poster_table_res <- 600

# Figures
df <- read_csv("means.csv.gz") |>
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

save_as_image(
  ft,
  path = here::here("Figures", "table.png"),
  res = poster_table_res
)
expected_n_distros <- read_csv("means.csv.gz", show_col_types = FALSE) |>
  summarize(n = n_distinct(Distribution)) |>
  pull(n)

graphing <- read_csv("graphing.csv.gz") |>
  bind_cols(read_csv("gamma_graphing.csv.gz"))

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
    dpi = poster_plot_dpi
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

x <- seq(0, 40, length.out = 100)
y <- dgamma(x, shape = 16, rate = 1)
distro_plot("Gamma 10 Z-Scores", x, y, "Gamma(16, 1) Population")

read_csv("means.csv.gz") |>
  filter(`Sample Size` <= 500) |>
  mutate(
    Total = `Upper Tail` + `Lower Tail`,
    Distribution = paste(
      str_replace(Distribution, "\\{.*\\}", " "),
      round(Skewness, 2)
    ),
    Distribution = fct(
      as.character(Distribution),
      levels = as.character(unique(Distribution[order(-Skewness)]))
    )
  ) |>
  arrange(Skewness) |>
  ggplot(aes(x = `Sample Size`, color = Distribution)) +
  geom_rect(
    aes(xmin = 0, xmax = Inf, ymin = 0.02, ymax = 0.03),
    fill = "grey",
    linewidth = 0,
    show.legend = FALSE
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
  dpi = poster_plot_dpi
)

# read_csv("means.csv.gz") |>
#   filter(`Sample Size` <= 500) |>
#   mutate(
#     Total = `Upper Tail` + `Lower Tail`,
#     Distribution = paste(
#       str_replace(Distribution, "\\{.*\\}", " "),
#       round(Skewness, 2)
#     ),
#     Distribution = fct(
#       as.character(Distribution),
#       levels = as.character(unique(Distribution[order(-Skewness)]))
#     )
#   ) |>
#   arrange(Skewness) |>
#   ggplot(aes(x = `Sample Size`, color = Distribution)) +
#   geom_rect(
#     aes(xmin = 0, xmax = Inf, ymin = 0.02, ymax = 0.03),
#     fill = "grey",
#     linewidth = 0,
#     show.legend = FALSE
#   ) +
#   geom_hline(yintercept = 0.025, linetype = "dashed", linewidth = 1) +
#   geom_line(aes(y = `Total`), linewidth = 1) +
#   geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.75) +
#   annotate(
#     "text",
#     x = 40,
#     y = 0.04,
#     label = "n = 30",
#     hjust = 0,
#     size = 5,
#     color = "black"
#   ) +
#   labs(
#     title = "Total Tail Weights of Sampling Distributions",
#     x = "Sample Size",
#     y = "Tail Weight"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = 20),
#     axis.title = element_text(size = 17),
#     axis.text = element_text(size = 12)
#   )

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
  ggrepel::geom_label_repel(seed = 0, nudge_x = 1, size = 3) +
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
  dpi = poster_plot_dpi
)

n_per_bound <- function(percent) {
  bounds <- 0.025 * percent * c(-1, 1) + 0.025
  temp <- read_csv("means.csv.gz", show_col_types = FALSE) |>
    group_by(Distribution) |>
    mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
    filter(
      `Lower Tail` >= bounds[1] & `Lower Tail` <= bounds[2],
      `Upper Tail` >= bounds[1] & `Upper Tail` <= bounds[2]
    ) |>
    filter(`Sample Size` == min(`Sample Size`)) |> # get the smallest sample size for each distribution
    select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
    arrange(Skewness, Distribution)

  stopifnot(nrow(temp) == expected_n_distros)

  coef(lm(Skewness ~ sqrt(`Sample Size`), temp))["sqrt(`Sample Size`)"]
}

x <- 50:16 / 100
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
  dpi = poster_plot_dpi
)

adjusted_skewness <- function(x) {
  n <- length(x)
  d <- x - mean(x)
  sqrt(n * (n - 1)) / (n - 2) * sum(d^3) / sum(d^2)^1.5 * length(x)^0.5
}
sampling_distribution <- \(distro, r) {
  replicate(r, adjusted_skewness(eval(distro)))
}

generate_data <- function(n, r, population_skews) {
  sampling_distributions <- list(
    sampling_distribution(expr(rgamma(!!n, 16)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.25)), r),
    sampling_distribution(expr(rgamma(!!n, 4)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.35)), r),
    sampling_distribution(expr(rgamma(!!n, 2.56)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.4)), r),
    sampling_distribution(expr(rgamma(!!n, 2)), r),
    sampling_distribution(expr(rgamma(!!n, 1.78)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.45)), r),
    sampling_distribution(expr(rgamma(!!n, 1.31)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.5)), r),
    sampling_distribution(expr(rgamma(!!n, 1)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.55)), r),
    sampling_distribution(expr(rgamma(!!n, 0.79)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.6)), r),
    sampling_distribution(expr(rgamma(!!n, 0.64)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.65)), r),
    sampling_distribution(expr(rgamma(!!n, 0.53)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.7)), r),
    sampling_distribution(expr(rgamma(!!n, 0.44)), r),
    sampling_distribution(expr(rlnorm(!!n, 0, 0.75)), r)
  )
  distribution_names <- c(
    "Gamma (α = 16, θ = 1)",
    "Log-normal (μ = 0, σ = 0.25)",
    "Gamma (α = 4, θ = 1)",
    "Log-normal (μ = 0, σ = 0.35)",
    "Gamma (α = 2.56, θ = 1)",
    "Log-normal (μ = 0, σ = 0.4)",
    "Gamma (α = 2, θ = 1)",
    "Gamma (α = 1.78, θ = 1)",
    "Log-normal (μ = 0, σ = 0.45)",
    "Gamma (α = 1.31, θ = 1)",
    "Log-normal (μ = 0, σ = 0.5)",
    "Gamma (α = 1, θ = 1)",
    "Log-normal (μ = 0, σ = 0.55)",
    "Gamma (α = 0.79, θ = 1)",
    "Log-normal (μ = 0, σ = 0.6)",
    "Gamma (α = 0.64, θ = 1)",
    "Log-normal (μ = 0, σ = 0.65)",
    "Gamma (α = 0.53, θ = 1)",
    "Log-normal (μ = 0, σ = 0.7)",
    "Gamma (α = 0.44, θ = 1)",
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
skew_gamma <- \(alpha) 2 / sqrt(alpha)
population_skews <- c(
  skew_gamma(16),
  skew_lnorm(0.25),
  1,
  skew_lnorm(0.35),
  skew_gamma(2.56),
  skew_lnorm(0.4),
  sqrt(2),
  skew_gamma(1.78),
  skew_lnorm(0.45),
  skew_gamma(1.31),
  skew_lnorm(0.5),
  skew_gamma(1),
  skew_lnorm(0.55),
  skew_gamma(0.79),
  skew_lnorm(0.6),
  skew_gamma(0.64),
  skew_lnorm(0.65),
  skew_gamma(0.53),
  skew_lnorm(0.7),
  skew_gamma(0.44),
  skew_lnorm(0.75)
)
ns <- c(
  seq(10, 100, 10),
  seq(100, 200, 25),
  seq(200, 500, 50),
  seq(500, 1000, 200)
) |>
  unique()
r <- 1e5

# TODO: extend this to lower skewnesses, missing small skew Gamma and log normal
# skew_data <- map_df(ns, \(n) generate_data(n, r, population_skews)) |>
#   write_csv(here::here("skew_data.csv"))
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
    title = "Population Skewness vs Mean Sampling Skewness with Regression Lines",
    x = "Population Skewness",
    y = "Mean Sampling Skewness",
    color = "Sample Size"
  ) +
  theme_bw()

skew_data |>
  mutate(
    percent = mean_sampling_skewness / pop_skewness,
    distribution = round(pop_skewness, 2), # Only keep the skewness digits
    distribution = fct(
      as.character(distribution),
      levels = as.character(unique(distribution[order(pop_skewness)]))
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
      mutate(
        percent = mean_sampling_skewness / pop_skewness,
        distribution = round(pop_skewness, 2), # Only keep the skewness digits
        distribution = fct(
          as.character(distribution),
          levels = as.character(unique(distribution[order(pop_skewness)]))
        )
      ) |>
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
    max.overlaps = Inf
  ) +
  labs(
    title = "Convergence Rate of Sample Skewness to Population Skewness",
    x = "Sample Size",
    y = "Percent of Population Skewness",
    color = "Skewness"
  )

ggsave(
  here::here("Figures", "skewness_convergence.png"),
  width = 12,
  height = 6.75,
  dpi = poster_plot_dpi
)

plot_skew_conv <- function(df, legend_position = "right") {
  legend_breaks <- levels(droplevels(df$distribution))

  hline_data <- df |>
    distinct(distribution, pop_skewness)

  df |>
    ggplot(
      aes(
        x = sample_size,
        y = mean_sampling_skewness,
        color = distribution
      )
    ) +
    geom_point() +
    geom_line() +
    # geom_hline(
    #   data = hline_data,
    #   aes(
    #     yintercept = pop_skewness,
    #     color = distribution
    #   ),
    #   linetype = "dashed",
    #   show.legend = FALSE
    # ) +
    scale_color_manual(
      values = full_dist_colors,
      limits = full_dist_levels,
      breaks = legend_breaks
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 12),
      legend.position = legend_position
    ) +
    labs(
      title = "Convergence Rate of Sample Skewness to Population Skewness",
      x = "Sample Size",
      y = "Mean Sampling Skewness",
      color = "Population Skewness"
    )
}

target_dists <- c("0.5", "2", "3.26")

full_dist_levels <- plot_data |>
  pull(distribution) |>
  levels()

full_dist_colors <- setNames(
  scales::hue_pal()(length(full_dist_levels)),
  full_dist_levels
)

plot_data <- skew_data |>
  mutate(
    distribution = round(pop_skewness, 2),
    distribution = fct(
      as.character(distribution),
      levels = as.character(unique(distribution[order(pop_skewness)]))
    )
  ) |>
  filter(distribution %in% target_dists)

target_dists |>
  walk(\(dist) {
    dist_file <- gsub("\\.", "_", dist)

    plot_data |>
      filter(
        distribution %in% target_dists[seq_len(match(dist, target_dists))]
      ) |>
      plot_skew_conv() |>
      ggsave(
        filename = glue::glue("Figures/skew_convergence_{dist_file}.png"),
        plot = _,
        width = 12,
        height = 6.75,
        dpi = poster_plot_dpi
      )
  })

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
  dpi = poster_plot_dpi
)

model <- lm(log_correction ~ log_sample_size, data = exp_data)
summary(model)

coeffs <- coef(model)
exp_corrected_skewness <- \(skew, n) {
  unname(skew * (exp(coeffs[1] + coeffs[2] * log(n)) + 1))
}

skew_data |>
  mutate(
    percent = exp_corrected_skewness(
      mean_sampling_skewness,
      sample_size
    ) /
      pop_skewness,
    distribution = round(pop_skewness, 2),
    distribution = fct(
      as.character(distribution),
      levels = as.character(unique(distribution[order(pop_skewness)]))
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
      mutate(
        percent = exp_corrected_skewness(
          mean_sampling_skewness,
          sample_size
        ) /
          pop_skewness,
        distribution = round(pop_skewness, 2),
        distribution = fct(
          as.character(distribution),
          levels = as.character(unique(distribution[order(pop_skewness)]))
        )
      ) |>
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
    max.overlaps = Inf
  ) +
  labs(
    x = "Sample Size",
    y = "Percentage of Population Skewness",
    title = "Corrected Sample Skewnesses",
    color = "Skewness"
  )

ggsave(
  here::here("Figures", "corrected_sample_skewness.png"),
  width = 12,
  height = 6.75,
  dpi = poster_plot_dpi
)

correction_factor <- \(n) exp_corrected_skewness(1, n)

corrected_sample_skewness_table <- tibble(
  sample_skewness = df$Skewness / correction_factor(df$`Sample Size`),
  min_n = df$`Sample Size`
) |>
  arrange(sample_skewness)

corrected_sample_skewness_model <- lm(
  sqrt(min_n) ~ sample_skewness,
  data = corrected_sample_skewness_table
)
corrected_sample_skewness_model |> summary()

corrected_sample_skewness_model_pred <- predict(
  corrected_sample_skewness_model
)^2
corrected_sample_skewness_model_rmse <- sqrt(
  mean(
    (corrected_sample_skewness_table$min_n -
      corrected_sample_skewness_model_pred)^2
  )
) |>
  print()

corrected_sample_skewness_formula <- \(x) (2 + 5.2 * x)^2
corrected_sample_skewness_formula_rmse <- sqrt(
  mean(
    (corrected_sample_skewness_table$min_n -
      corrected_sample_skewness_formula(
        corrected_sample_skewness_table$sample_skewness
      ))^2
  )
) |>
  print()

corrected_sample_skewness_model_line <- tibble(
  sample_skewness = seq(
    min(corrected_sample_skewness_table$sample_skewness),
    max(corrected_sample_skewness_table$sample_skewness),
    length.out = 100
  )
) |>
  mutate(
    sample_size = predict(
      corrected_sample_skewness_model,
      newdata = pick(everything())
    )^2
  )

corrected_sample_skewness_model_line |>
  ggplot(aes(sample_skewness, sample_size)) +
  geom_line(
    linewidth = 1.5,
    color = "#08519C"
  ) +
  labs(
    title = "Sample Size vs Sample Skewness (1 sample means)",
    x = "Sample Skewness",
    y = "Sample Size"
  ) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Figures", "Sample_Skewness_Sample_Size.png"),
  width = 15,
  height = 8.5,
  dpi = poster_plot_dpi
)

# t stats two sided
tibble(
  sample_skewness = seq(
    min(t_corrected_sample_skewness_table$sample_skewness),
    max(t_corrected_sample_skewness_table$sample_skewness),
    length.out = 100
  )
) |>
  mutate(
    sample_size = predict(
      t_corrected_sample_skewness_model,
      newdata = pick(sample_skewness)
    )^2
  ) |>
  ggplot(aes(sample_skewness, sample_size)) +
  geom_line(
    linewidth = 1.5,
    color = "#08519C"
  ) +
  geom_point(
    data = tibble(
      sample_skewness = 0.595,
      required_n = 81
    ),
    aes(sample_skewness, required_n),
    inherit.aes = FALSE,
    size = 4,
    color = "#D62728"
  ) +
  geom_text(
    data = tibble(
      sample_skewness = 0.595,
      required_n = 81
    ),
    aes(
      sample_skewness,
      required_n,
      label = "Party-arrival data"
    ),
    inherit.aes = FALSE,
    nudge_x = 0.035,
    nudge_y = 5,
    hjust = 0,
    color = "#D62728",
    size = 4
  ) +
  labs(
    title = "Minimum Sample Size vs Sample Skewness (Two-Tailed <i>t</i>-Statistics)",
    x = "Sample Skewness",
    y = "Sample Size"
  ) +
  coord_cartesian(expand = FALSE) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = ggtext::element_markdown(size = 24),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = here::here(
    "Figures",
    paste0(
      "Corrected_Sample_Skewness_Sample_Size_t_statistics_",
      on,
      ".png"
    )
  ),
  width = 15,
  height = 8.5,
  dpi = plot_dpi
)

# ------------------------ two sample ------------------------

diff_means <- read_csv("difference_means.csv.gz", show_col_types = FALSE) |>
  mutate(
    Ratio = case_when(
      `Sample Size 2` == `Sample Size 1` ~ "1:1",
      `Sample Size 2` == ceiling(`Sample Size 1` / 2) ~ "1:2",
      TRUE ~ "1:3"
    ),
    Distribution = str_replace(Distribution, "\\{.*\\}", " ")
  )

c("1:1", "1:2", "1:3") |>
  walk(\(ratio_val) {
    diff_means |>
      filter(Ratio == ratio_val) |>
      group_by(Distribution, Ratio) |>
      filter(
        `Lower Tail` >= 0.02 & `Lower Tail` <= 0.03,
        `Upper Tail` >= 0.02 & `Upper Tail` <= 0.03
      ) |>
      filter(`Sample Size 1` == min(`Sample Size 1`)) |>
      ungroup() |>
      select(
        Distribution,
        Ratio,
        Skewness,
        `Sampling Skewness`,
        `Sample Size 1`,
        `Sample Size 2`
      ) |>
      mutate(`Combined n` = `Sample Size 1` + `Sample Size 2`) |>
      arrange(Ratio, Skewness, Distribution) |>
      flextable() |>
      colformat_double(j = c("Skewness", "Sampling Skewness"), digits = 3) |>
      font(fontname = "Trebuchet MS") |>
      autofit() |>
      save_as_image(
        path = here::here(
          "Figures",
          paste0(
            "table_diff_means_",
            str_replace_all(ratio_val, ":", "_"),
            ".png"
          )
        ),
        res = poster_table_res
      )
  })

c("1:1", "1:2", "1:3") |>
  walk(\(ratio_val) {
    p <- diff_means |>
      filter(`Sample Size 1` <= 500, Ratio == ratio_val) |>
      mutate(
        Distribution = paste(Distribution, round(Skewness, 2)),
        Distribution = fct(
          as.character(Distribution),
          levels = as.character(unique(Distribution[order(-Skewness)]))
        )
      ) |>
      arrange(Skewness) |>
      ggplot(aes(x = `Sample Size 1`, color = Distribution)) +
      geom_rect(
        aes(xmin = 0, xmax = Inf, ymin = 0.02, ymax = 0.03),
        fill = "grey",
        linewidth = 0,
        show.legend = FALSE
      ) +
      geom_hline(yintercept = 0.025, linetype = "dashed", linewidth = 1) +
      geom_line(aes(y = `Upper Tail`), linewidth = 1) +
      geom_line(aes(y = `Lower Tail`), linewidth = 1) +
      labs(
        title = paste(
          "Upper and Lower Tail Weights for Difference in Means (",
          ratio_val,
          ")",
          sep = ""
        ),
        x = "Sample Size 1",
        y = "Tail Weight"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 12)
      )

    ggsave(
      here::here(
        "Figures",
        paste0(
          "Tail_Weights_Diff_Means_",
          str_replace_all(ratio_val, ":", "_"),
          ".png"
        )
      ),
      plot = p,
      width = 12,
      height = 6.75,
      dpi = poster_plot_dpi
    )
  })

diff_df <- diff_means |>
  group_by(Distribution, Ratio) |>
  filter(
    `Lower Tail` >= 0.02 & `Lower Tail` <= 0.03,
    `Upper Tail` >= 0.02 & `Upper Tail` <= 0.03
  ) |>
  filter(`Sample Size 1` == min(`Sample Size 1`)) |>
  ungroup() |>
  arrange(Ratio, Skewness, Distribution)

diff_df |>
  ggplot(aes(
    x = sqrt(`Sample Size 1`),
    y = Skewness,
    color = Ratio,
    label = Distribution
  )) +
  geom_point(size = 2.5) +
  geom_smooth(
    data = diff_df,
    aes(x = sqrt(`Sample Size 1`), y = Skewness, color = Ratio),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggrepel::geom_label_repel(seed = 0, size = 3.5, show.legend = FALSE) +
  labs(
    title = "Skewness vs Minimum sqrt(Sample Size 1) by Ratio",
    x = "Square Root of Sample Size 1",
    y = "Skewness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Figures", "Skew_Sample_Size_Diff_Means.png"),
  width = 15,
  height = 8.5,
  dpi = poster_plot_dpi
)
