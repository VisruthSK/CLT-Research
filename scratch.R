library(tidyverse)
library(flextable)
set.seed(0)

t_stats <- read_csv("t_statistics.csv.gz", show_col_types = FALSE) |>
  mutate(Total = `Upper Tail` + `Lower Tail`)

expected_n_distros <- t_stats |>
  summarize(n = n_distinct(Distribution)) |>
  pull(n)

minimum_n_table <- function(bounds = NULL, on = c("total", "tails")) {
  on <- match.arg(on)
  if (is.null(bounds)) {
    bounds <- if (on == "total") {
      relative_error_bounds(0.05, 0.20)
    } else {
      relative_error_bounds(0.025, 0.20)
    }
  }

  t_stats |>
    filter(
      if (on == "total") {
        between(Total, bounds[1], bounds[2])
      } else {
        between(`Lower Tail`, bounds[1], bounds[2]) &
          between(`Upper Tail`, bounds[1], bounds[2])
      }
    ) |>
    group_by(Distribution) |>
    slice_min(`Sample Size`, n = 1, with_ties = TRUE) |>
    ungroup() |>
    mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
    select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
    arrange(Skewness, Distribution)
}

ordered_distribution_labels <- function(df) {
  df |>
    mutate(
      Distribution = paste(
        str_replace(Distribution, "\\{.*\\}", " "),
        round(Skewness, 2)
      ),
      Distribution = fct(
        as.character(Distribution),
        levels = as.character(unique(Distribution[order(-Skewness)]))
      )
    ) |>
    arrange(Skewness)
}

relative_error_bounds <- function(center, percent) {
  center * (1 + percent * c(-1, 1))
}

predict_n_from_skewness <- function(model_coeffs, skewness) {
  as.integer(unname(ceiling(
    ((skewness - model_coeffs[1]) / model_coeffs[2])^2
  )))
}

on <- "tails"
bounds <- if (on == "total") {
  relative_error_bounds(0.05, 0.20)
} else {
  relative_error_bounds(0.025, 0.20)
}
table_res <- 600
plot_dpi <- 450

table_df <- minimum_n_table(bounds = bounds, on = on)
stopifnot(nrow(table_df) == expected_n_distros)

ft <- flextable(table_df) |>
  colformat_double(
    j = c("Skewness", "Sampling Skewness"),
    digits = 3
  ) |>
  font(fontname = "Trebuchet MS") |>
  autofit()

save_as_image(
  ft,
  path = here::here("Figures", paste0("table_t_statistics_", on, ".png")),
  res = table_res
)

tail_plot <- t_stats |>
  ordered_distribution_labels() |>
  ggplot(aes(x = `Sample Size`, color = Distribution)) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = Inf,
      ymin = relative_error_bounds(0.025, 0.20)[1],
      ymax = relative_error_bounds(0.025, 0.20)[2]
    ),
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
    title = "Upper and Lower Tail Weights of t-Statistic Sampling Distributions",
    x = "Sample Size",
    y = "Tail Weight"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

total_tail_plot <- t_stats |>
  ordered_distribution_labels() |>
  ggplot(aes(x = `Sample Size`, color = Distribution)) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = Inf,
      ymin = relative_error_bounds(0.05, 0.10)[1],
      ymax = relative_error_bounds(0.05, 0.10)[2]
    ),
    fill = "grey",
    linewidth = 0,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = Total), linewidth = 1) +
  geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.75) +
  annotate(
    "text",
    x = 40,
    y = 0.06,
    label = "n = 30",
    hjust = 0,
    size = 5,
    color = "black"
  ) +
  labs(
    title = "Total Tail Weights of t-Statistic Sampling Distributions",
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
  here::here("Figures", "Tail_Weights_t_statistics.png"),
  plot = tail_plot,
  width = 12,
  height = 6.75,
  dpi = plot_dpi
)

ggsave(
  here::here("Figures", "Total_Tail_Weights_t_statistics.png"),
  plot = total_tail_plot,
  width = 12,
  height = 6.75,
  dpi = plot_dpi
)

avg_sample_skewness_df <- t_stats |>
  filter(`Sample Size` <= 500) |>
  mutate(
    `Average Sample Skewness (%)` = 100 * `Average Sample Skewness` / Skewness
  ) |>
  ordered_distribution_labels()

avg_sample_skewness_plot <- avg_sample_skewness_df |>
  ggplot(aes(
    x = `Sample Size`,
    y = `Average Sample Skewness (%)`,
    color = Distribution
  )) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 30, linetype = "dashed", linewidth = 0.75) +
  annotate(
    "text",
    x = 40,
    y = max(avg_sample_skewness_df$`Average Sample Skewness (%)`) * 0.95,
    label = "n = 30",
    hjust = 0,
    size = 5,
    color = "black"
  ) +
  labs(
    title = "Average Sample Skewness as Percentage of Total Skewness",
    x = "Sample Size",
    y = "Average Sample Skewness (%)"
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
    paste0("Avg_Sample_Skewness_t_statistics_", on, ".png")
  ),
  plot = avg_sample_skewness_plot,
  width = 12,
  height = 6.75,
  dpi = plot_dpi
)

pop_model <- lm(Skewness ~ sqrt(`Sample Size`), table_df)
summary(pop_model)

pop_skew_plot <- table_df |>
  ggplot(aes(
    x = sqrt(`Sample Size`),
    y = Skewness,
    label = Distribution
  )) +
  geom_abline(
    slope = coef(pop_model)[2],
    intercept = coef(pop_model)[1],
    color = "blue",
    linetype = "dashed",
    linewidth = 1.5
  ) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  labs(
    title = paste(
      "Linear Relationship between Population Skewness and",
      "Empirical Minimum Square Root Sample Size for t Statistics"
    ),
    x = "Square Root of Sample Size",
    y = "Population Skewness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 23),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  ) +
  ggrepel::geom_label_repel(seed = 0, nudge_x = 1, size = 5)

ggsave(
  here::here(
    "Figures",
    paste0("Skew_Sample_Size_t_statistics_pop_", on, ".png")
  ),
  plot = pop_skew_plot,
  width = 15,
  height = 8.5,
  dpi = plot_dpi
)

# 10% relative error tolerance for t statistics
ten_percent_error <- 0.10
ten_percent_tail_bounds <- relative_error_bounds(0.025, ten_percent_error)
ten_percent_table_df <- minimum_n_table(
  bounds = ten_percent_tail_bounds,
  on = "tails"
)

ten_percent_pop_model <- lm(
  Skewness ~ sqrt(`Sample Size`),
  ten_percent_table_df
)
summary(ten_percent_pop_model)

expo_empirical_10_percent <- t_stats |>
  filter(
    str_detect(Distribution, "^Exponential"),
    between(
      `Lower Tail`,
      ten_percent_tail_bounds[1],
      ten_percent_tail_bounds[2]
    ),
    between(
      `Upper Tail`,
      ten_percent_tail_bounds[1],
      ten_percent_tail_bounds[2]
    )
  ) |>
  slice_min(`Sample Size`, n = 1, with_ties = TRUE)

expo_required_n_10_percent <- predict_n_from_skewness(
  coef(ten_percent_pop_model),
  skewness = 2
)

cat(
  paste(
    "Estimated n required for an exponential population with 10% relative",
    "tail-error tolerance for t statistics:",
    expo_required_n_10_percent
  ),
  "\n"
)

exp_data <- t_stats |>
  filter(
    str_detect(Distribution, "^Exponential") |
      Distribution == "Gamma{Float64}(α=1.0, θ=1.0)"
  ) |>
  mutate(
    log_correction = log((Skewness / `Average Sample Skewness`) - 1),
    log_sample_size = log(`Sample Size`)
  )

model <- lm(log_correction ~ log_sample_size, data = exp_data)
summary(model)

coeffs <- coef(model)
exp_corrected_skewness <- \(skew, n) {
  unname(skew * (exp(coeffs[1] + coeffs[2] * log(n)) + 1))
}

corrected_table_df <- table_df |>
  left_join(
    t_stats |>
      mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
      select(
        Distribution,
        `Sample Size`,
        `Average Sample Skewness`
      ),
    by = c("Distribution", "Sample Size")
  ) |>
  mutate(
    `Corrected Skewness` = exp_corrected_skewness(
      `Average Sample Skewness`,
      `Sample Size`
    )
  )

t_correction_factor <- \(n) exp_corrected_skewness(1, n)

t_corrected_sample_skewness_table <- tibble(
  sample_skewness = table_df$Skewness /
    t_correction_factor(table_df$`Sample Size`),
  min_n = table_df$`Sample Size`
) |>
  arrange(sample_skewness)

t_corrected_sample_skewness_model <- lm(
  sqrt(min_n) ~ sample_skewness,
  data = t_corrected_sample_skewness_table
)
t_corrected_sample_skewness_model |> summary()

t_corrected_sample_skewness_model_pred <- predict(
  t_corrected_sample_skewness_model
)^2
t_corrected_sample_skewness_model_rmse <- sqrt(
  mean(
    (t_corrected_sample_skewness_table$min_n -
      t_corrected_sample_skewness_model_pred)^2
  )
) |>
  print()

t_corrected_sample_skewness_formula_coeffs <- round(
  coef(t_corrected_sample_skewness_model),
  1
) |>
  print()

t_corrected_sample_skewness_formula <- \(x) {
  (t_corrected_sample_skewness_formula_coeffs[1] +
    t_corrected_sample_skewness_formula_coeffs[2] * x)^2
}
t_corrected_sample_skewness_formula_rmse <- sqrt(
  mean(
    (t_corrected_sample_skewness_table$min_n -
      t_corrected_sample_skewness_formula(
        t_corrected_sample_skewness_table$sample_skewness
      ))^2
  )
) |>
  print()

sample_skewness_breaks <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
pop_model_coeffs <- coef(pop_model)

sample_skewness_lookup <- map_int(
  sample_skewness_breaks,
  \(sample_skewness) {
    n_grid <- 5:10000
    n_grid[
      which(
        (pop_model_coeffs[1] +
          pop_model_coeffs[2] * sqrt(n_grid)) /
          exp_corrected_skewness(1, n_grid) >=
          sample_skewness
      )[1]
    ]
  }
)

tibble(
  Metric = c("Sample skewness", "Min sample size"),
  `.25` = c(".25", as.character(sample_skewness_lookup[1])),
  `.5` = c(".5", as.character(sample_skewness_lookup[2])),
  `.75` = c(".75", as.character(sample_skewness_lookup[3])),
  `1.0` = c("1.0", as.character(sample_skewness_lookup[4])),
  `1.5` = c("1.5", as.character(sample_skewness_lookup[5])),
  `2.0` = c("2.0", as.character(sample_skewness_lookup[6]))
) |>
  print() |>
  flextable() |>
  font(fontname = "Trebuchet MS") |>
  autofit() |>
  save_as_image(
    path = here::here(
      "Figures",
      paste0("table_t_statistics_sample_skewness_", on, ".png")
    ),
    res = table_res
  )

# Ensuring table/correction matches means
means_sample_skewness_breaks <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
slide_table <- c(10, 20, 35, 50, 100, 165)

means_table_df <- read_csv("means.csv.gz", show_col_types = FALSE) |>
  filter(
    between(
      `Lower Tail`,
      relative_error_bounds(0.025, 0.20)[1],
      relative_error_bounds(0.025, 0.20)[2]
    ),
    between(
      `Upper Tail`,
      relative_error_bounds(0.025, 0.20)[1],
      relative_error_bounds(0.025, 0.20)[2]
    )
  ) |>
  group_by(Distribution) |>
  slice_min(`Sample Size`, n = 1, with_ties = TRUE) |>
  ungroup() |>
  select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
  arrange(Skewness, Distribution)

means_pop_model <- lm(Skewness ~ sqrt(`Sample Size`), means_table_df)

skew_data <- read_csv("skew_data.csv", show_col_types = FALSE)

means_exp_data <- skew_data |>
  filter(distribution == "Gamma (α = 1, θ = 1)") |>
  mutate(
    log_correction = log((pop_skewness / mean_sampling_skewness) - 1),
    log_sample_size = log(sample_size)
  )

means_correction_model <- lm(
  log_correction ~ log_sample_size,
  data = means_exp_data
)
means_coeffs <- coef(means_correction_model)

means_exp_corrected_skewness <- \(skew, n) {
  unname(skew * (exp(means_coeffs[1] + means_coeffs[2] * log(n)) + 1))
}

means_pop_model_coeffs <- coef(means_pop_model)

means_sample_skewness_lookup <- map_int(
  means_sample_skewness_breaks,
  \(sample_skewness) {
    n_grid <- 5:10000
    n_grid[
      which(
        (means_pop_model_coeffs[1] +
          means_pop_model_coeffs[2] * sqrt(n_grid)) /
          means_exp_corrected_skewness(1, n_grid) >=
          sample_skewness
      )[1]
    ]
  }
)

means_slide_formula_lookup <- map_int(
  means_sample_skewness_breaks,
  \(sample_skewness) {
    n_grid <- 5:10000
    n_grid[
      which(
        n_grid^1.365 /
          (6 * (n_grid^0.865 + exp(1.696))) >=
          sample_skewness
      )[1]
    ]
  }
)

tibble(
  sample_skewness = means_sample_skewness_breaks,
  clt_r_lookup = means_sample_skewness_lookup,
  slide_formula_lookup = means_slide_formula_lookup,
  slide_table = slide_table
)

# Two-sample t statistics
two_sample_t <- read_csv(
  "two_sample_t_statistics.csv.gz",
  show_col_types = FALSE
) |>
  mutate(
    Ratio = case_when(
      `Sample Size 2` == `Sample Size 1` ~ "1:1",
      `Sample Size 2` == ceiling(`Sample Size 1` / 2) ~ "1:2",
      TRUE ~ "1:3"
    ),
    Distribution = str_replace(Distribution, "\\{.*\\}", " "),
    Total = `Upper Tail` + `Lower Tail`
  )

minimum_two_sample_t_table <- function(
  bounds = relative_error_bounds(0.025, 0.20),
  on = c("total", "tails")
) {
  on <- match.arg(on)

  two_sample_t |>
    filter(
      if (on == "total") {
        between(Total, bounds[1], bounds[2])
      } else {
        between(`Lower Tail`, bounds[1], bounds[2]) &
          between(`Upper Tail`, bounds[1], bounds[2])
      }
    ) |>
    group_by(Distribution, Ratio) |>
    slice_min(`Sample Size 1`, n = 1, with_ties = TRUE) |>
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
    arrange(Ratio, Skewness, Distribution)
}

two_sample_t_min_df <- minimum_two_sample_t_table(on = "tails")

c("1:1", "1:2", "1:3") |>
  walk(\(ratio_val) {
    two_sample_t_min_df |>
      filter(Ratio == ratio_val) |>
      flextable() |>
      colformat_double(j = c("Skewness", "Sampling Skewness"), digits = 3) |>
      font(fontname = "Trebuchet MS") |>
      autofit() |>
      save_as_image(
        path = here::here(
          "Figures",
          paste0(
            "table_two_sample_t_statistics_",
            str_replace_all(ratio_val, ":", "_"),
            ".png"
          )
        ),
        res = table_res
      )
  })

c("1:1", "1:2", "1:3") |>
  walk(\(ratio_val) {
    p <- two_sample_t |>
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
        aes(
          xmin = 0,
          xmax = Inf,
          ymin = relative_error_bounds(0.025, 0.20)[1],
          ymax = relative_error_bounds(0.025, 0.20)[2]
        ),
        fill = "grey",
        linewidth = 0,
        show.legend = FALSE
      ) +
      geom_hline(yintercept = 0.025, linetype = "dashed", linewidth = 1) +
      geom_line(aes(y = `Upper Tail`), linewidth = 1) +
      geom_line(aes(y = `Lower Tail`), linewidth = 1) +
      labs(
        title = paste(
          "Upper and Lower Tail Weights for Two-Sample t Statistics (",
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
          "Tail_Weights_Two_Sample_t_",
          str_replace_all(ratio_val, ":", "_"),
          ".png"
        )
      ),
      plot = p,
      width = 12,
      height = 6.75,
      dpi = plot_dpi
    )
  })

two_sample_t_min_df |>
  ggplot(aes(
    x = sqrt(`Sample Size 1`),
    y = Skewness,
    color = Ratio,
    label = Distribution
  )) +
  geom_point(size = 2.5) +
  geom_smooth(
    data = two_sample_t_min_df,
    aes(x = sqrt(`Sample Size 1`), y = Skewness, color = Ratio),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggrepel::geom_label_repel(seed = 0, size = 3.5, show.legend = FALSE) +
  labs(
    title = "Skewness vs Minimum sqrt(Sample Size 1) for Two-Sample t Statistics",
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
  here::here("Figures", "Skew_Sample_Size_Two_Sample_t.png"),
  width = 15,
  height = 8.5,
  dpi = plot_dpi
)
