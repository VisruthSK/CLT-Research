library(tidyverse)
library(flextable)
set.seed(0)

poster_plot_dpi <- 400
poster_table_res <- 600

# Ensure Figures directory exists
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# 1. Read datasets
means_df <- read_csv("means.csv.gz", show_col_types = FALSE) |>
  select(Distribution, Skewness) |>
  distinct()

boot_df <- read_csv("bootstrap_comparison.csv.gz", show_col_types = FALSE) |>
  left_join(means_df, by = "Distribution") |>
  mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " "))

# ==============================================================================
# ANALYSIS 1: Bootstrap Tails compared to Asymptotic Nominal Value (0.025)
# ==============================================================================
# We find when bootstrap tails are within 20% of 0.025 (i.e., [0.02, 0.03])
nominal_value <- 0.025
error_margin <- 0.20
lower_bound <- nominal_value * (1 - error_margin)
upper_bound <- nominal_value * (1 + error_margin)

cat(sprintf("Finding convergence to nominal value [%.3f, %.3f]...\n", lower_bound, upper_bound))

df_nominal <- boot_df |>
  group_by(Distribution) |>
  filter(
    `Bootstrap Lower Tail` >= lower_bound & `Bootstrap Lower Tail` <= upper_bound,
    `Bootstrap Upper Tail` >= lower_bound & `Bootstrap Upper Tail` <= upper_bound
  ) |>
  filter(`Sample Size` == min(`Sample Size`)) |>
  select(Distribution, Skewness, `Sampling Skewness`, `Sample Size`) |>
  arrange(Skewness, Distribution)

# Fit linear model of Skewness on sqrt(n)
model_nominal <- lm(Skewness ~ sqrt(`Sample Size`), data = df_nominal)
print(summary(model_nominal))

# Plot skewness vs sqrt(n) for nominal convergence
df_nominal |>
  ggplot(aes(x = sqrt(`Sample Size`), y = Skewness, label = Distribution)) +
  geom_abline(
    slope = coef(model_nominal)[2],
    intercept = coef(model_nominal)[1],
    color = "red",
    linetype = "dashed",
    linewidth = 1.5
  ) +
  geom_point(size = 3, color = "darkred") +
  ggrepel::geom_label_repel(seed = 0, nudge_x = 1, size = 3) +
  labs(
    title = "Bootstrap Tail Convergence to Nominal Value (0.025 ± 20%)",
    subtitle = paste("R-squared:", round(summary(model_nominal)$r.squared, 3)),
    x = "Square Root of Minimum Sample Size",
    y = "Population Skewness"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

ggsave(
  filename = "Figures/bootstrap_nominal_convergence.png",
  width = 12,
  height = 6.75,
  dpi = poster_plot_dpi
)




# ==============================================================================
# COMPARISON TABLE: Sampling vs. Bootstrap Convergence
# ==============================================================================
# Compare minimum sample size for Sampling distro vs Bootstrap distro to reach nominal bounds [0.02, 0.03]
sampling_nominal <- read_csv("means.csv.gz", show_col_types = FALSE) |>
  group_by(Distribution) |>
  mutate(Distribution = str_replace(Distribution, "\\{.*\\}", " ")) |>
  filter(
    `Lower Tail` >= lower_bound & `Lower Tail` <= upper_bound,
    `Upper Tail` >= lower_bound & `Upper Tail` <= upper_bound
  ) |>
  filter(`Sample Size` == min(`Sample Size`)) |>
  select(Distribution, Skewness, Sampling_N = `Sample Size`)

comparison_df <- df_nominal |>
  select(Distribution, Skewness, Bootstrap_N = `Sample Size`) |>
  left_join(sampling_nominal, by = c("Distribution", "Skewness")) |>
  select(Distribution, Skewness, `Sampling Min N` = Sampling_N, `Bootstrap Min N` = Bootstrap_N) |>
  arrange(Skewness, Distribution)

# Print comparison
print(comparison_df)

# Save Comparison Table as image
ft_comparison <- flextable(comparison_df) |>
  colformat_double(j = "Skewness", digits = 3) |>
  font(fontname = "Trebuchet MS") |>
  autofit()

save_as_image(
  ft_comparison,
  path = here::here("Figures", "bootstrap_vs_sampling_table.png"),
  res = poster_table_res
)

cat("Analysis complete. Visualizations saved to Figures/ directory.\n")
