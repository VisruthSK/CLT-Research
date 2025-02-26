library(tidyverse)
library(flextable)

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

save_as_image(ft, path = here::here("Poster", "table.png"), res = 2000)

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
      "Poster",
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
    plot.title = element_text(size = 23),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12)
  )

ggsave(
  here::here("Poster", "Tail_Weights.png"),
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
  here::here("Poster", "Skew_Sample_Size.png"),
  width = 15,
  height = 8.5,
  dpi = 1000
)
