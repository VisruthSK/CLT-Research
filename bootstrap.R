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

# ==============================================================================
# PART 2: Expected Sample Skewness Bias Correction and Modeling
# ==============================================================================
# Load skew_data.csv to model the sample skewness downward bias
skew_data <- read_csv("skew_data.csv", show_col_types = FALSE)

exp_data <- skew_data |>
  filter(distribution == "Gamma (α = 1, θ = 1)") |>
  mutate(
    log_correction = log((pop_skewness / mean_sampling_skewness) - 1),
    log_sample_size = log(sample_size)
  )
model_bias <- lm(log_correction ~ log_sample_size, data = exp_data)
coeffs <- coef(model_bias)

exp_corrected_skewness <- \(skew, n) {
  unname(skew * (exp(coeffs[1] + coeffs[2] * log(n)) + 1))
}
correction_factor <- \(n) exp_corrected_skewness(1, n)

# Apply correction factor to compute expected sample skewness at minimum N
corrected_clt_table <- sampling_nominal |>
  select(Distribution, Skewness, min_n = Sampling_N) |>
  mutate(sample_skewness = Skewness / correction_factor(min_n))

corrected_bootstrap_table <- df_nominal |>
  select(Distribution, Skewness, min_n = `Sample Size`) |>
  mutate(sample_skewness = Skewness / correction_factor(min_n))

# Fit the models: sqrt(min_n) ~ sample_skewness
corrected_sample_skewness_model <- lm(
  sqrt(min_n) ~ sample_skewness,
  data = corrected_clt_table
)
corrected_bootstrap_model <- lm(
  sqrt(min_n) ~ sample_skewness,
  data = corrected_bootstrap_table
)

cat("\nCorrected Sample Skewness CLT Model:\n")
print(summary(corrected_sample_skewness_model))
cat("\nCorrected Sample Skewness Bootstrap Model:\n")
print(summary(corrected_bootstrap_model))

# Create comparison plot
clt_sample_points <- corrected_clt_table |>
  select(sample_skewness, min_n) |>
  mutate(Type = "CLT (Sampling)")

boot_sample_points <- corrected_bootstrap_table |>
  select(sample_skewness, min_n) |>
  mutate(Type = "Bootstrap")

combined_sample_df <- bind_rows(clt_sample_points, boot_sample_points)

sample_skewness_seq <- seq(
  min(combined_sample_df$sample_skewness),
  max(combined_sample_df$sample_skewness),
  length.out = 100
)

line_clt <- tibble(
  sample_skewness = sample_skewness_seq,
  Type = "CLT (Sampling)"
) |>
  mutate(sample_size = predict(corrected_sample_skewness_model, newdata = pick(sample_skewness))^2)

line_boot <- tibble(
  sample_skewness = sample_skewness_seq,
  Type = "Bootstrap"
) |>
  mutate(sample_size = predict(corrected_bootstrap_model, newdata = pick(sample_skewness))^2)

combined_lines <- bind_rows(line_clt, line_boot)

combined_sample_df |>
  ggplot(aes(x = sample_skewness, y = min_n, color = Type, shape = Type)) +
  geom_point(size = 3) +
  geom_line(data = combined_lines, aes(y = sample_size), linewidth = 1.5) +
  labs(
    title = "Sample Size vs. Expected Sample Skewness (CLT vs. Bootstrap)",
    subtitle = sprintf(
      "CLT: N = (%.4f * S + %.4f)² (RMSE = %.3f)\nBootstrap: N = (%.4f * S + %.4f)² (RMSE = %.3f)",
      coef(corrected_sample_skewness_model)[2], coef(corrected_sample_skewness_model)[1],
      sqrt(mean((clt_sample_points$min_n - predict(corrected_sample_skewness_model)^2)^2)),
      coef(corrected_bootstrap_model)[2], coef(corrected_bootstrap_model)[1],
      sqrt(mean((boot_sample_points$min_n - predict(corrected_bootstrap_model)^2)^2))
    ),
    x = "Expected Sample Skewness (S)",
    y = "Minimum Sample Size (N)"
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
  filename = "Figures/clt_vs_bootstrap_sample_skewness.png",
  width = 12,
  height = 6.75,
  dpi = poster_plot_dpi
)

cat("Analysis complete. Visualizations saved to Figures/ directory.\n")
