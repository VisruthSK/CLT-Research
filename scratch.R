library(tidyverse)
library(flextable)
set.seed(0)
t_stats <- read_csv("t_statistics.csv.gz", show_col_types = FALSE) |>
  mutate(Total = `Upper Tail` + `Lower Tail`)

expected_n_distros <- t_stats |>
  summarize(n = n_distinct(Distribution)) |>
  pull(n)

minimum_n_table <- function(bounds = c(0.04, 0.06), on = c("total", "tails")) {
  on <- match.arg(on)
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

on <- "tails"
bounds <- if (on == "total") c(0.04, 0.06) else c(0.02, 0.03)
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
    aes(xmin = 0, xmax = Inf, ymin = 0.045, ymax = 0.055),
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

model <- lm(Skewness ~ sqrt(`Sample Size`), table_df)
summary(model)

skew_plot <- table_df |>
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
  labs(
    title = "Linear Relationship between Skewness and Empirical Minimum Square Root Sample Size for t Statistics",
    x = "Square Root of Sample Size",
    y = "Skewness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 23),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  ) +
  ggrepel::geom_label_repel(seed = 0, nudge_x = 1, size = 5)

ggsave(
  here::here("Figures", paste0("Skew_Sample_Size_t_statistics_", on, ".png")),
  plot = skew_plot,
  width = 15,
  height = 8.5,
  dpi = plot_dpi
)
