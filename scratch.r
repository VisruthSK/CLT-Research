library(tidyverse)
library(moments)
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
r <- 100000

# looking at moments::skewness we can see it calculates "g1", being m3/m2^3/2;
# we can get G1 as defined in Joanes, Gill by the following adjustment.
# They claim it performs better as an estimator for non-normal distributions
adjusted_skewness <- function(x) {
  n <- length(x)
  sqrt(n * (n - 1)) / (n - 2) * skewness(x)
}

temp <- rgamma(30, 16)
adjusted_skewness(temp)

# adjskew <- replicate(r, adjusted_skew(rgamma(30, 7.5)))
# hist(adjskew)
# mean(adjskew)

# gamma30 <- replicate(r, skewness(rgamma(30, 16)))
# hist(gamma30)
# mean(gamma30)

test <- function(distro, r) {
  temp <- replicate(r, adjusted_skewness(eval(distro)))
  # hist(temp)
  print(mean(temp))
  temp
}

# Gamma 16 distro, skew = 0.5
ns <- seq(10, 250, 20)
datans <- lapply(ns, function(n) test(expr(rgamma(!!n, 16)), r))
means <- datans |>
  lapply(mean) |>
  unlist()

plot(ns, means)
michaelis_menten <- function(n, Vmax, Km) {
  Vmax * n / (Km + n)
}
model <- nls(
  means ~ michaelis_menten(ns, Vmax, Km),
  start = list(Vmax = 0.5, Km = median(ns))
)
summary(model)
Vmax_est <- coef(model)["Vmax"]
Km_est <- coef(model)["Km"]
plot(ns, means, xlab = "n", ylab = "Mean", main = "Michaelis-Menten Fit")
curve(michaelis_menten(x, Vmax_est, Km_est), add = TRUE, col = "red")

# Expo distro, skew = 2
ns <- seq(10, 250, 20)
datans <- lapply(ns, function(n) test(expr(rexp(!!n)), r))
means <- datans |>
  lapply(mean) |>
  unlist()

plot(ns, means)
model <- nls(
  means ~ michaelis_menten(ns, Vmax, Km),
  start = list(Vmax = 2, Km = median(ns))
)
summary(model)
Vmax_est <- coef(model)["Vmax"]
Km_est <- coef(model)["Km"]
plot(ns, means, xlab = "n", ylab = "Mean", main = "Michaelis-Menten Fit")
curve(michaelis_menten(x, Vmax_est, Km_est), add = TRUE, col = "red")

# Fix n, population skewness vs mean sampling distro skewnwess
n <- 10
sampling_distro16 <- test(expr(rgamma(!!n, 16)), r)
population_skew16 <- 0.5
sampling_distro4 <- test(expr(rgamma(!!n, 4)), r)
population_skew4 <- 1
sampling_distro2 <- test(expr(rgamma(!!n, 2)), r)
population_skew2 <- sqrt(2)
sampling_distro1 <- test(expr(rgamma(!!n, 1)), r)
population_skew1 <- 2
sampling_distro64 <- test(expr(rgamma(!!n, 0.64)), r)
population_skew64 <- 2.5










# population skewness as response
# put these data into table, three columns pop skewness, sample size, mean sampling skewness
# each row is a different distro
adjusted_skewness <- function(x) {
  n <- length(x)
  sqrt(n * (n - 1)) / (n - 2) * moments::skewness(x)
}
test <- \(distro, r) replicate(r, adjusted_skewness(eval(distro)))

generate_data <- function(n, r, population_skews) {
  raw_distributions <- list(
    test(expr(rgamma(!!n, 16)), r),
    test(expr(rlnorm(!!n, 0, 0.25)), r),
    test(expr(rgamma(!!n, 4)), r),
    test(expr(rgamma(!!n, 2)), r),
    test(expr(rlnorm(!!n, 0, 0.5)), r),
    test(expr(rgamma(!!n, 1)), r),
    test(expr(rgamma(!!n, 0.64)), r),
    test(expr(rlnorm(!!n, 0, 0.75)), r)
  )
  mss <- map_dbl(raw_distributions, mean)

  data.frame(
    pop_skewness = population_skews,
    sample_size = n,
    mean_sampling_skewness = mss
  )
}

population_skews <- c(
  0.5,
  (exp(0.25^2) + 2) * sqrt(exp(0.25^2) - 1),
  1,
  sqrt(2),
  (exp(0.5^2) + 2) * sqrt(exp(0.5^2) - 1),
  2,
  2.5,
  (exp(0.75^2) + 2) * sqrt(exp(0.75^2) - 1)
)
ns <- c(seq(10, 50, 10), seq(50, 100, 25)) |> unique()
r <- 1e6

skew_data <- map_df(ns, \(n) generate_data(n, r, population_skews))
skew_data |> write_csv(here::here("skew_data"))
# skew_data <- read_csv(here::here("skew_data"))
ggplot(skew_data, aes(x = pop_skewness, y = mean_sampling_skewness, color = fct(as.character(sample_size)))) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Population Skewness vs. Mean Sampling Skewness with Regression Lines",
    x = "Population Skewness",
    y = "Mean Sampling Skewness",
    color = "Sample Size"
  ) +
  theme_minimal()











data_func <- function(distro, n, r, ...) {
  replicate(r, adjusted_skewness(distro(n, ...)))

  # distro(n, ...)
  # sampling_distro <- test(expr(!!distro(!!n, 1)), r)
  # sampling_distro
}

data_func <- function(distro, n, r, ...) {
  replicate(r, adjusted_skewness(distro(n, ...)))
}

data_func(rgamma, n, 10, shape = 16)

popskews <- c(
  population_skew16,
  population_skew4,
  population_skew2,
  population_skew1,
  population_skew64
)
samp_skews <- c(
  mean(sampling_distro16),
  mean(sampling_distro4),
  mean(sampling_distro2),
  mean(sampling_distro1),
  mean(sampling_distro64)
)

plot(popskews, samp_skews)

# find pathways to help students better connect with the material, talk about both in context of teaching stats.

gamma1016 <- test(expr(rgamma(10, 16)), r) # skew 0.5
gamma3016 <- test(expr(rgamma(30, 16)), r) # skew 0.5
exp301 <- test(expr(rexp(30)), r) # skew 2

boostrapped_skewness <- function(x) {
  adjusted_skewness(sample(x, length(x), TRUE))
}

tempdata <- rgamma(30, 16) # actual skewness is 0.5

adjusted_skewness(tempdata)

boot_skew <- replicate(r, boostrapped_skewness(tempdata))
hist(boot_skew)
mean(boot_skew)
mean(boot_skew <= 0.5)

sample_size <- function(x) {
  -5.19 * x^4 - 1.44 * x^3 + 5.65 * x^2 + 14.75 * x + 70.07
}
skews <- c(0.5, 0.778, 1, 1.414, 1.75, 2, 2.5, 3.263)
sample_size_data <- sapply(skews, function(x) sample_size(x))
names(sample_size_data) <- skews

wrap_test <- function(distro, r) mean(test(distro, r))

ns <- seq(5, 105, 10)
skew_over_n <- map_dbl(ns, function(x) wrap_test(expr(rgamma(!!x, 16)), r))
names(skew_over_n) <- ns
exp_skew_over_n <- map_dbl(ns, function(x) wrap_test(expr(rgamma(!!x, 1)), r))
names(exp_skew_over_n) <- ns
lnorm_skew_over_n <- map_dbl(
  ns,
  function(x) wrap_test(expr(rlnorm(!!x, 0, .75)), r)
)
names(lnorm_skew_over_n) <- ns

plot(ns, skew_over_n)
plot(ns, exp_skew_over_n)
plot(ns, lnorm_skew_over_n)

plot(ns, skew_over_n / 0.5, ylim = c(0, 1))
plot(ns, exp_skew_over_n / 2, ylim = c(0, 1))
plot(ns, lnorm_skew_over_n / 3.263, ylim = c(0, 1))

skew_over_n / 0.5
exp_skew_over_n / 2
lnorm_skew_over_n / 3.263
