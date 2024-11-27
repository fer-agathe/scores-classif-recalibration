library(tidyverse)
# Figure motivation----
predicted_risk <- seq(0, 1, length.out = 1000)
f_inverse_logistic <- function(x, x0 = .5, k=.1) {
  x0 - k * log(1/x-1)
}

observed <- f_inverse_logistic(predicted_risk, x0 =.45, k=.1)

pdf("../figs/autocalib.pdf", width = 4, height=4)
par(mar = c(4.1, 4.1, 1.1, 1.1))
plot(
  predicted_risk, observed, t= "l",
  xlim = c(0, 1), ylim = c(0,1),
  xlab = "Predicted risk", ylab = "Observed proportion",
  col = "#56B4E9", lwd=2
)
abline(a = 0, b = 1, lty = 2, col = "black", lwd = 2)
text(x = .7, y = .9, substitute(paste(bold("Perfect"))), col = "black")
text(x = .2, y = .5, substitute(paste(bold("Underestimation\n of lower risks"))), col = "#56B4E9")
text(x = .8, y = .3, substitute(paste(bold("Overestimation\n of higher risks"))), col = "#56B4E9")
dev.off()

png("../figs/thumbnail.png", width = 320, height=256, res = 75, bg = "transparent")
par(mar = c(4.1, 4.1, 1.1, 1.1))
plot(
  predicted_risk, observed, t= "l",
  xlim = c(0, 1), ylim = c(0,1),
  xlab = "Predicted risk", ylab = "Observed proportion",
  col = "#56B4E9", lwd=2
)
abline(a = 0, b = 1, lty = 2, col = "black", lwd = 2)
text(x = .7, y = .9, substitute(paste(bold("Perfect"))), col = "black")
text(x = .2, y = .5, substitute(paste(bold("Underestimation\n of lower risks"))), col = "#56B4E9")
text(x = .8, y = .3, substitute(paste(bold("Overestimation\n of higher risks"))), col = "#56B4E9")
dev.off()


# Figure simulated probabilities
source("../scripts/functions/simul-data.R")

nb_obs <- 10000

coefficients <- list(
  # First category (baseline, 2 covariates)
  c(0.5, 1),  # scenario 1, 0 noise variable
  c(0.5, 1),  # scenario 2, 10 noise variables
  c(0.5, 1),  # scenario 3, 50 noise variables
  c(0.5, 1),  # scenario 4, 100 noise variables
  # Second category (same as baseline, with lower number of 1s)
  c(0.5, 1),  # scenario 5, 0 noise variable
  c(0.5, 1),  # scenario 6, 10 noise variables
  c(0.5, 1),  # scenario 7, 50 noise variables
  c(0.5, 1),  # scenario 8, 100 noise variables
  # Third category (same as baseline but with 5 num. and 5 categ. covariates)
  c(0.1, 0.2, 0.3, 0.4, 0.5, 0.01, 0.02, 0.03, 0.04, 0.05),
  c(0.1, 0.2, 0.3, 0.4, 0.5, 0.01, 0.02, 0.03, 0.04, 0.05),
  c(0.1, 0.2, 0.3, 0.4, 0.5, 0.01, 0.02, 0.03, 0.04, 0.05),
  c(0.1, 0.2, 0.3, 0.4, 0.5, 0.01, 0.02, 0.03, 0.04, 0.05),
  # Fourth category (nonlinear predictor, 3 covariates)
  c(0.5, 1, .3),  # scenario 5, 0 noise variable
  c(0.5, 1, .3),  # scenario 6, 10 noise variables
  c(0.5, 1, .3),  # scenario 7, 50 noise variables
  c(0.5, 1, .3)  # scenario 8, 100 noise variables
)

# Mean parameter for the normal distribution to draw from to draw num covariates
mean_num <- list(
  # First category (baseline, 2 covariates)
  rep(0, 2),  # scenario 1, 0 noise variable
  rep(0, 2),  # scenario 2, 10 noise variables
  rep(0, 2),  # scenario 3, 50 noise variables
  rep(0, 2),  # scenario 4, 100 noise variables
  # Second category (same as baseline, with lower number of 1s)
  rep(0, 2),  # scenario 5, 0 noise variable
  rep(0, 2),  # scenario 6, 10 noise variables
  rep(0, 2),  # scenario 7, 50 noise variables
  rep(0, 2),  # scenario 8, 100 noise variables
  # Third category (same as baseline but with 5 num. and 5 categ. covariates)
  rep(0, 5),
  rep(0, 5),
  rep(0, 5),
  rep(0, 5),
  # Fourth category (nonlinear predictor, 3 covariates)
  rep(0, 3),
  rep(0, 3),
  rep(0, 3),
  rep(0, 3)
)
# Sd parameter for the normal distribution to draw from to draw num covariates
sd_num <- list(
  # First category (baseline, 2 covariates)
  rep(1, 2),  # scenario 1, 0 noise variable
  rep(1, 2),  # scenario 2, 10 noise variables
  rep(1, 2),  # scenario 3, 50 noise variables
  rep(1, 2),  # scenario 4, 100 noise variables
  # Second category (same as baseline, with lower number of 1s)
  rep(1, 2),  # scenario 5, 0 noise variable
  rep(1, 2),  # scenario 6, 10 noise variables
  rep(1, 2),  # scenario 7, 50 noise variables
  rep(1, 2),  # scenario 8, 100 noise variables
  # Third category (same as baseline but with 5 num. and 5 categ. covariates)
  rep(1, 5),
  rep(1, 5),
  rep(1, 5),
  rep(1, 5),
  # Fourth category (nonlinear predictor, 3 covariates)
  rep(1, 3),
  rep(1, 3),
  rep(1, 3),
  rep(1, 3)
)

params_df <- tibble(
  scenario = 1:16,
  coefficients = coefficients,
  n_num = c(rep(2, 8), rep(5, 4), rep(3, 4)),
  add_categ = c(rep(FALSE, 8), rep(TRUE, 4), rep(FALSE, 4)),
  n_noise = rep(c(0, 10, 50, 100), 4),
  mean_num = mean_num,
  sd_num = sd_num,
  size_train = rep(nb_obs, 16),
  size_valid = rep(nb_obs, 16),
  size_calib = rep(nb_obs, 16),
  size_test = rep(nb_obs, 16),
  transform_probs = c(rep(FALSE, 4), rep(TRUE, 4), rep(FALSE, 4), rep(FALSE, 4)),
  linear_predictor = c(rep(TRUE, 12), rep(FALSE, 4)),
  seed = 202105
)
rm(coefficients, mean_num, sd_num)


simu_data <- simulate_data_wrapper(
  scenario = 1,
  params_df = params_df,
  repn = 1
)

simu_data$data

true_prob <- simu_data$data$probs_train

pdf("../figs/true-prob-dgp1.pdf", width = 5, height=4)
par(mar = c(2.1, 2.1, 1.1, 1.1))
hist(
  true_prob,
  breaks = seq(0, 1, by = .05),
  col = "gray",
  border = "white",
  xlab = "p", ylab = "",
  xlim = c(0, 1),
  main = "True Probabilities"
)
dev.off()
