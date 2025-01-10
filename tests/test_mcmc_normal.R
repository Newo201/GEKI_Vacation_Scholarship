pacman::p_load(pacman, testthat, purrr)

source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/mcmc_normal.R')

x <- rep(1, 50)
prior_params = list(alpha.mean = 0, alpha.sd = 5, sigma2.mean = 0, sigma2.sd = 2)

lupost_normal_test <- partial(lupost_normal, theta.init = c(1, 2), x = x, prior_params = prior_params)

test_that('log_posterior_density is a scalar', {
  expect_equal(is.atomic(lupost_normal_test(rep(1, 50))), TRUE)
})

test_that('Extreme values have probability 0', {
  expect_equal(exp(lupost_normal_test(rep(-10, 50))), 0)
  expect_equal(exp(lupost_normal_test(rep(10, 50))), 0)
})

test_that('Density peaks at the true value of alpha', {
  expect_gt(lupost_normal_test(rep(1, 50)), lupost_normal_test(rep(2, 50)))
  expect_gt(lupost_normal_test(rep(1, 50)), lupost_normal_test(rep(0, 50)))
})

# y_seq <- seq(-5, 5, length.out = 11)
# res <- c()
# 
# for (y in y_seq) {
#   res <- c(res, lupost_normal(rep(y, 50), c(1, 2), x, prior_params))
# }
# 
# plot(y_seq, res)