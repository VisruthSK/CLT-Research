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

# ==============================================================================
# COMBINED MODEL PLOT: CLT vs. Bootstrap
# ==============================================================================
clt_points <- sampling_nominal |>
  select(Distribution, Skewness, N = Sampling_N) |>
  mutate(Type = "CLT (Sampling)")

boot_points <- df_nominal |>
  select(Distribution, Skewness, N = `Sample Size`) |>
  mutate(Type = "Bootstrap")

combined_df <- bind_rows(clt_points, boot_points)

model_clt <- lm(Skewness ~ sqrt(N), data = clt_points)
model_boot <- lm(Skewness ~ sqrt(N), data = boot_points)

# Write summary of models
cat("\nCLT Model summary:\n")
print(summary(model_clt))
cat("\nBootstrap Model summary:\n")
print(summary(model_boot))

combined_df |>
  ggplot(aes(x = sqrt(N), y = Skewness, color = Type, shape = Type)) +
  geom_point(size = 3) +
  geom_abline(
    slope = coef(model_clt)[2], intercept = coef(model_clt)[1],
    color = "#1f77b4", linetype = "dashed", linewidth = 1.2
  ) +
  geom_abline(
    slope = coef(model_boot)[2], intercept = coef(model_boot)[1],
    color = "#d62728", linetype = "dashed", linewidth = 1.2
  ) +
  labs(
    title = "CLT Sampling vs. Bootstrap Tail Convergence Model",
    subtitle = sprintf(
      "CLT: Skewness = %.4f * sqrt(N) + %.4f (R² = %.3f)\nBootstrap: Skewness = %.4f * sqrt(N) + %.4f (R² = %.3f)",
      coef(model_clt)[2], coef(model_clt)[1], summary(model_clt)$r.squared,
      coef(model_boot)[2], coef(model_boot)[1], summary(model_boot)$r.squared
    ),
    x = "Square Root of Minimum Sample Size",
    y = "Population Skewness"
  ) +
  theme_bw() +
  scale_color_manual(values = c("CLT (Sampling)" = "#1f77b4", "Bootstrap" = "#d62728")) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

ggsave(
  filename = "Figures/clt_vs_bootstrap_model.png",
  width = 12,
  height = 6.75,
  dpi = poster_plot_dpi
)

cat("Analysis complete. Visualizations saved to Figures/ directory.\n")
