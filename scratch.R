library(tidyverse)
skewnesses <- seq(2, 10)
ns <- c(seq(100, 1500, 50))

expand.grid(skewnesses, ns) |>
  rename(skew = Var1, n = Var2) |>
  mutate(
    required_n = 36 * skew^2,
    ok = n >= required_n,
    c_skew = skew * (exp(-0.865 * log(n) + 1.696) + 1),
    c_required_n = ceiling(36 * c_skew^2),
    c_ok = n >= c_required_n,
    c_c_skew = skew * (exp(-0.865 * log(c_required_n) + 1.696) + 1),
    c_c_required_n = ceiling(36 * c_c_skew^2),
    c_c_ok = n >= c_c_required_n,
  )
# |>
# group_by(skew) |>
# filter(c_ok) |>
# slice(1)

acceptable_ns <- c(9, 19, 33, 51, 99, 164)

tibble(skew = skewnesses, n = acceptable_ns) |>
  mutate(
    c_skew = skew * (exp(-0.865 * log(n) + 1.696) + 1),
    c_required_n = ceiling(36 * c_skew^2),
    c_ok = n >= c_required_n
  )

tibble(skew = skewnesses, n = acceptable_ns - 1) |>
  mutate(
    c_skew = skew * (exp(-0.865 * log(n) + 1.696) + 1),
    c_required_n = ceiling(36 * c_skew^2),
    c_ok = n >= c_required_n
  )

# q-p/sqrt(npq)
bin_skew <- \(n, p, q = 1 - p) (q - p) / sqrt(n * p * q)
bern_skew <- \(p, q = 1 - p) (q - p) / sqrt(p * q)

p <- 0.05
bin_skew(bern_skew(p)^2 * 36, p)

set.seed(0)
n <- 165
t_stats <- replicate(1e5, t.test(rexp(n))$statistic)
normed_t_stats <- (t_stats - mean(t_stats)) / sd(t_stats)
true <- qt(0.975, n - 1)
observed <- ecdf(normed_t_stats) |> quantile(c(0.025, 0.975))
1 - (observed[1] + true) / true * 100
