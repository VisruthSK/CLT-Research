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
r <- 1000000

adjusted_skewness <- function(x) {
  n <- length(x)
  sqrt(n * (n - 1)) / (n - 2) * moments::skewness(x)
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

gamma1016 <- test(expr(rgamma(10, 16)), r)

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

sample_size <- function(x) -5.19 * x^4 - 1.44 * x^3 + 5.65 * x^2 + 14.75 * x + 70.07
skews <- c(0.5, 0.778, 1, 1.414, 1.75, 2, 2.5, 3.263)
sample_size_data <- sapply(skews, function(x) sample_size(x))
names(sample_size_data) <- skews


wrap_test <- function(distro, r) mean(test(distro, r))


ns <- seq(5, 105, 10)
skew_over_n <- map_dbl(ns, function(x) wrap_test(expr(rgamma(!!x, 16)), r))
names(skew_over_n) <- ns
exp_skew_over_n <- map_dbl(ns, function(x) wrap_test(expr(rgamma(!!x, 1)), r))
names(exp_skew_over_n) <- ns
lnorm_skew_over_n <- map_dbl(ns, function(x) wrap_test(expr(rlnorm(!!x, 0, .75)), r))
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
