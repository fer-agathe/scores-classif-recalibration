library(tidyverse)
library(ggh4x)
library(ggrepel)
library(rpart)
library(locfit)
library(philentropy)

# Colours for train/validation/test
colour_samples <- c(
  "Train" = "#0072B2",
  "Validation" = "#009E73",
  "Calibration" = "#CC79A7",
  "Test" = "#D55E00"
)

colour_recalib <- c(
  "None" = "#88CCEE",
  "Platt" = "#44AA99",
  "Isotonic" = "#882255"
)

source("../scripts/functions/simul-data.R")
library(ks)
source("../scripts/functions/subsample_target_distribution.R")

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

source("../scripts/functions/metrics.R")

library(xgboost)

recalibrate <- function(obs_calib,
                        obs_test,
                        pred_calib,
                        pred_test,
                        method = c("platt", "isotonic")) {
  data_calib <- tibble(d = obs_calib, scores = pred_calib)
  data_test <- tibble(d = obs_test, scores = pred_test)

  if (method == "platt") {
    lr <- glm(d ~ scores, family = binomial(link = 'logit'), data = data_calib)
    score_c_calib <- predict(lr, newdata = data_calib, type = "response")
    score_c_test <- predict(lr, newdata = data_test, type = "response")
  } else if (method == "isotonic") {
    iso <- isoreg(x = data_calib$scores, y = data_calib$d)
    fit_iso <- as.stepfun(iso)
    score_c_calib <- fit_iso(data_calib$scores)
    score_c_test <- fit_iso(data_test$scores)

  } else {
    stop("Unrecognized method: platt or isotonic only")
  }
  # Format results in tibbles:
  # For calibration set
  tb_score_c_calib <- tibble(
    d = obs_calib,
    p_u = pred_calib,
    p_c = score_c_calib
  )
  # For test set
  tb_score_c_test <- tibble(
    d = obs_test,
    p_u = pred_test,
    p_c = score_c_test
  )

  list(
    tb_score_c_calib = tb_score_c_calib,
    tb_score_c_test = tb_score_c_test
  )

}

#' Computes the performance and calibration metrics for an xgb model,
#' depending on the number of iterations kept.
#'
#' @param nb_iter number of boosting iterations to keep
#' @param params hyperparameters of the current model
#' @param fitted_xgb xgb estimated model
#' @param tb_train_xgb train data (in xgb.DMatrix format)
#' @param tb_valid_xgb validation data (in xgb.DMatrix format)
#' @param tb_calib_xgb calibration data (in xgb.DMatrix format)
#' @param tb_test_xgb test data (in xgb.DMatrix format)
#' @param simu_data simulated dataset
#' @param true_prob list with true probabilities on train, calibration,
#'  validation and test sets
get_metrics_nb_iter <- function(nb_iter,
                                params,
                                fitted_xgb,
                                tb_train_xgb,
                                tb_valid_xgb,
                                tb_calib_xgb,
                                tb_test_xgb,
                                simu_data,
                                true_prob) {

  ind <- params$ind
  max_depth <- params$max_depth
  tb_train <- simu_data$data$train |> rename(d = y)
  tb_valid <- simu_data$data$valid |> rename(d = y)
  tb_calib <- simu_data$data$calib |> rename(d = y)
  tb_test <- simu_data$data$test |> rename(d = y)

  # Predicted scores
  scores_train <- predict(fitted_xgb, tb_train_xgb, iterationrange = c(1, nb_iter))
  scores_valid <- predict(fitted_xgb, tb_valid_xgb, iterationrange = c(1, nb_iter))
  scores_calib <- predict(fitted_xgb, tb_calib_xgb, iterationrange = c(1, nb_iter))
  scores_test <- predict(fitted_xgb, tb_test_xgb, iterationrange = c(1, nb_iter))

  # Recalibration
  # Platt scaling
  res_recalibration_platt <- recalibrate(
    obs_calib = tb_calib$d,
    obs_test = tb_test$d,
    pred_calib = scores_calib,
    pred_test = scores_test,
    method = "platt"
  )
  scores_c_platt_calib <- res_recalibration_platt$tb_score_c_calib$p_c
  scores_c_platt_test <- res_recalibration_platt$tb_score_c_test$p_c

  # Isotonic regression
  res_recalibration_iso <- recalibrate(
    obs_calib = tb_calib$d,
    obs_test = tb_test$d,
    pred_calib = scores_calib,
    pred_test = scores_test,
    method = "isotonic"
  )
  scores_c_iso_calib <- res_recalibration_iso$tb_score_c_calib$p_c
  scores_c_iso_test <- res_recalibration_iso$tb_score_c_test$p_c

  ## Histogram of scores----
  breaks <- seq(0, 1, by = .05)
  scores_train_hist <- hist(scores_train, breaks = breaks, plot = FALSE)
  scores_calib_hist <- hist(scores_calib, breaks = breaks, plot = FALSE)
  scores_valid_hist <- hist(scores_valid, breaks = breaks, plot = FALSE)
  scores_test_hist <- hist(scores_test, breaks = breaks, plot = FALSE)
  scores_c_platt_calib_hist <- hist(scores_c_platt_calib, breaks = breaks, plot = FALSE)
  scores_c_platt_test_hist <- hist(scores_c_platt_test, breaks = breaks, plot = FALSE)
  scores_c_iso_calib_hist <- hist(scores_c_iso_calib, breaks = breaks, plot = FALSE)
  scores_c_iso_test_hist <- hist(scores_c_iso_test, breaks = breaks, plot = FALSE)

  scores_hist <- list(
    train = scores_train_hist,
    valid = scores_valid_hist,
    calib = scores_calib_hist,
    test = scores_test_hist,
    calib_c_platt = scores_c_platt_calib_hist,
    test_c_platt = scores_c_platt_test_hist,
    calib_c_iso = scores_c_iso_calib_hist,
    test_c_iso = scores_c_iso_test_hist,
    scenario = simu_data$scenario,
    ind = ind,
    repn = simu_data$repn,
    max_depth = params$max_depth,
    nb_iter = nb_iter
  )

  ## Estimation of P(q1 < score < q2)----
  prop_btw_q_h <- function(s, sample_name, recalib_name) {
    map(
      c(.1, .2, .3, .4),
      ~prop_btw_quantiles(s = s, q1 = .x)
    ) |>
      list_rbind() |>
      mutate(sample = sample_name, recalib = recalib_name)
  }

  proq_scores_train <- prop_btw_q_h(
    scores_train, sample_name = "train", recalib_name = "none"
  )
  proq_scores_valid <- prop_btw_q_h(
    scores_valid, sample_name = "valid", recalib_name = "none"
  )
  proq_scores_calib <- prop_btw_q_h(
    scores_calib, sample_name = "calib", recalib_name = "none"
  )
  proq_scores_test <- prop_btw_q_h(
    scores_test, sample_name = "test", recalib_name = "none"
  )
  proq_scores_c_platt_calib <- prop_btw_q_h(
    scores_c_platt_calib, sample_name = "calib", recalib_name = "platt"
  )
  proq_scores_c_platt_test <- prop_btw_q_h(
    scores_c_platt_test, sample_name = "test", recalib_name = "platt"
  )
  proq_scores_c_iso_calib <- prop_btw_q_h(
    scores_c_iso_calib, sample_name = "calib", recalib_name = "isotonic"
  )
  proq_scores_c_iso_test <- prop_btw_q_h(
    scores_c_iso_test, sample_name = "test", recalib_name = "isotonic"
  )


  ## Dispersion Metrics----
  disp_train <- dispersion_metrics(
    true_probas = true_prob$train, scores = scores_train
  ) |>
    mutate(sample = "train", recalib = "none")
  disp_valid <- dispersion_metrics(
    true_probas = true_prob$valid, scores = scores_valid
  ) |>
    mutate(sample = "valid", recalib = "none")

  disp_calib <- dispersion_metrics(
    true_probas = true_prob$calib, scores = scores_calib
  ) |>
    mutate(sample = "calib", recalib = "none")

  disp_test <- dispersion_metrics(
    true_probas = true_prob$test, scores = scores_test
  ) |>
    mutate(sample = "test", recalib = "none")


  disp_c_platt_calib <- dispersion_metrics(
    true_probas = true_prob$calib, scores = scores_c_platt_calib
  ) |>
    mutate(sample = "calib", recalib = "platt")

  disp_c_platt_test <- dispersion_metrics(
    true_probas = true_prob$test, scores = scores_c_platt_test
  ) |>
    mutate(sample = "test", recalib = "platt")

  disp_c_iso_calib <- dispersion_metrics(
    true_probas = true_prob$calib, scores = scores_c_iso_calib
  ) |>
    mutate(sample = "calib", recalib = "isotonic")

  disp_c_iso_test <- dispersion_metrics(
    true_probas = true_prob$test, scores = scores_c_iso_test
  ) |>
    mutate(sample = "test", recalib = "isotonic")

  # Performance and Calibration Metrics
  # We add very small noise to predicted scores
  # otherwise the local regression may crash
  scores_train_noise <- scores_train +
    runif(n = length(scores_train), min = 0, max = 0.01)
  scores_train_noise[scores_train_noise > 1] <- 1
  metrics_train <- compute_metrics(
    obs = tb_train$d, scores = scores_train_noise, true_probas = true_prob$train
  ) |> mutate(sample = "train", recalib = "none")

  scores_valid_noise <- scores_valid +
    runif(n = length(scores_valid), min = 0, max = 0.01)
  scores_valid_noise[scores_valid_noise > 1] <- 1
  metrics_valid <- compute_metrics(
    obs = tb_valid$d, scores = scores_valid_noise, true_probas = true_prob$valid
  ) |> mutate(sample = "valid", recalib = "none")

  scores_calib_noise <- scores_calib +
    runif(n = length(scores_calib), min = 0, max = 0.01)
  scores_calib_noise[scores_calib_noise > 1] <- 1
  metrics_calib <- compute_metrics(
    obs = tb_calib$d, scores = scores_calib_noise, true_probas = true_prob$calib
  ) |> mutate(sample = "calib", recalib = "none")

  scores_test_noise <- scores_test +
    runif(n = length(scores_test), min = 0, max = 0.01)
  scores_test_noise[scores_test_noise > 1] <- 1
  metrics_test <- compute_metrics(
    obs = tb_test$d, scores = scores_test_noise, true_probas = true_prob$test
  ) |> mutate(sample = "test", recalib = "none")

  # With recalibrated scores (platt)
  scores_c_platt_calib_noise <- scores_c_platt_calib +
    runif(n = length(scores_c_platt_calib), min = 0, max = 0.01)
  scores_c_platt_calib_noise[scores_c_platt_calib_noise > 1] <- 1
  metrics_c_platt_calib <- compute_metrics(
    obs = tb_calib$d, scores = scores_c_platt_calib_noise,
    true_probas = true_prob$calib
  ) |> mutate(sample = "calib", recalib = "platt")

  scores_c_platt_test_noise <- scores_c_platt_test +
    runif(n = length(scores_c_platt_test), min = 0, max = 0.01)
  scores_c_platt_test_noise[scores_c_platt_test_noise > 1] <- 1
  metrics_c_platt_test <- compute_metrics(
    obs = tb_test$d, scores = scores_c_platt_test_noise,
    true_probas = true_prob$test
  ) |> mutate(sample = "test", recalib = "platt")

  # With recalibrated scores (isotonic)
  scores_c_iso_calib_noise <- scores_c_iso_calib +
    runif(n = length(scores_c_iso_calib), min = 0, max = 0.01)
  scores_c_iso_calib_noise[scores_c_iso_calib_noise > 1] <- 1
  metrics_c_iso_calib <- compute_metrics(
    obs = tb_calib$d, scores = scores_c_iso_calib_noise,
    true_probas = true_prob$calib
  ) |> mutate(sample = "calib", recalib = "isotonic")

  scores_c_iso_test_noise <- scores_c_iso_test +
    runif(n = length(scores_c_iso_test), min = 0, max = 0.01)
  scores_c_iso_test_noise[scores_c_iso_test_noise > 1] <- 1
  metrics_c_iso_test <- compute_metrics(
    obs = tb_test$d, scores = scores_c_iso_test_noise,
    true_probas = true_prob$test
  ) |> mutate(sample = "test", recalib = "isotonic")

  tb_metrics <- metrics_train |>
    bind_rows(metrics_valid) |>
    bind_rows(metrics_calib) |>
    bind_rows(metrics_test) |>
    bind_rows(metrics_c_platt_calib) |>
    bind_rows(metrics_c_platt_test) |>
    bind_rows(metrics_c_iso_calib) |>
    bind_rows(metrics_c_iso_test) |>
    left_join(
      disp_train |>
        bind_rows(disp_valid) |>
        bind_rows(disp_calib) |>
        bind_rows(disp_test) |>
        bind_rows(disp_c_platt_calib) |>
        bind_rows(disp_c_platt_test) |>
        bind_rows(disp_c_iso_calib) |>
        bind_rows(disp_c_iso_test),
      by = c("sample", "recalib")
    ) |>
    mutate(
      scenario = simu_data$scenario,
      ind = ind,
      repn = simu_data$repn,
      max_depth = params$max_depth,
      nb_iter = nb_iter
    )

  tb_prop_scores <- proq_scores_train |>
    bind_rows(proq_scores_valid) |>
    bind_rows(proq_scores_calib) |>
    bind_rows(proq_scores_test) |>
    bind_rows(proq_scores_c_platt_calib) |>
    bind_rows(proq_scores_c_platt_test) |>
    bind_rows(proq_scores_c_iso_calib) |>
    bind_rows(proq_scores_c_iso_test) |>
    mutate(
      scenario = simu_data$scenario,
      ind = ind,
      repn = simu_data$repn,
      max_depth = params$max_depth,
      nb_iter = nb_iter
    )

  list(
    scenario = simu_data$scenario,     # data scenario
    ind = ind,                         # index for grid
    repn = simu_data$repn,             # data replication ID
    nb_iter = nb_iter,                 # number of boosting iterations
    tb_metrics = tb_metrics,           # table with performance/calib/divergence
    #  metrics
    tb_prop_scores = tb_prop_scores,   # table with P(q1 < score < q2)
    scores_hist = scores_hist          # histogram of scores
  )
}

#' Train an xgboost model and compute performance, calibration, and dispersion
#' metrics
#'
#' @param params tibble with hyperparameters for the simulation
#' @param ind index of the grid (numerical ID)
#' @param simu_data simulated data obtained with `simulate_data_wrapper()`
simul_xgb <- function(params,
                      ind,
                      simu_data) {
  tb_train <- simu_data$data$train |> rename(d = y)
  tb_valid <- simu_data$data$valid |> rename(d = y)
  tb_calib <- simu_data$data$calib |> rename(d = y)
  tb_test <- simu_data$data$test |> rename(d = y)
  true_prob <-
    list(
      train = simu_data$data$probs_train,
      valid = simu_data$data$probs_valid,
      calib = simu_data$data$probs_calib,
      test = simu_data$data$probs_test
    )

  ## Format data for xgboost----
  tb_train_xgb <- xgb.DMatrix(
    data = model.matrix(d ~ -1 + ., tb_train), label = tb_train$d
  )
  tb_valid_xgb <- xgb.DMatrix(
    data = model.matrix(d ~ -1 + ., tb_valid), label = tb_valid$d
  )
  tb_calib_xgb <- xgb.DMatrix(
    data = model.matrix(d ~ -1 + ., tb_calib), label = tb_calib$d
  )
  tb_test_xgb <- xgb.DMatrix(
    data = model.matrix(d ~ -1 + ., tb_test), label = tb_test$d
  )
  # Parameters for the algorithm
  param <- list(
    max_depth = params$max_depth, #Note: root node is indexed 0
    eta = params$eta,
    nthread = 1,
    objective = "binary:logistic",
    eval_metric = "auc"
  )
  watchlist <- list(train = tb_train_xgb, eval = tb_valid_xgb)

  ## Estimation----
  xgb_fit <- xgb.train(
    param, tb_train_xgb,
    nrounds = params$nb_iter_total,
    watchlist,
    verbose = 0
  )

  # Then, for each boosting iteration number up to params$nb_iter_total
  # compute the predicted scores and evaluate the metrics
  resul <- map(
    seq(2, params$nb_iter_total),
    ~get_metrics_nb_iter(
      nb_iter = .x,
      params = params,
      fitted_xgb = xgb_fit,
      tb_train_xgb = tb_train_xgb,
      tb_valid_xgb = tb_valid_xgb,
      tb_calib_xgb = tb_calib_xgb,
      tb_test_xgb = tb_test_xgb,
      simu_data = simu_data,
      true_prob = true_prob
    ),
  )
  resul
}

simulate_xgb_scenario <- function(scenario, params_df, repn) {
  # Generate Data
  simu_data <- simulate_data_wrapper(
    scenario = scenario,
    params_df = params_df,
    repn = repn
  )

  # Looping over the grid hyperparameters for the scenario
  res_simul <- vector(mode = "list", length = nrow(grid))
  cli::cli_progress_bar("Iteration grid", total = nrow(grid), type = "tasks")
  for (j in 1:nrow(grid)) {
    curent_params <- grid |> dplyr::slice(!!j)
    res_simul[[j]] <- simul_xgb(
      params = curent_params,
      ind = curent_params$ind,
      simu_data = simu_data
    )
    cli::cli_progress_update()
  }


  # The metrics computed for all set of hyperparameters (identified with `ind`)
  # and for each number of boosting iterations (`nb_iter`), for the current
  # scenario (`scenario`) and current replication number (`repn`)
  metrics_simul <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "tb_metrics") |> list_rbind()
  ) |>
    list_rbind()

  # P(q_1<s(x)<q_2)
  prop_scores_simul <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "tb_prop_scores") |> list_rbind()
  ) |>
    list_rbind()

  # Histogram of estimated scores
  scores_hist <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "scores_hist")
  )

  list(
    metrics_simul = metrics_simul,
    scores_hist = scores_hist,
    prop_scores_simul = prop_scores_simul
  )
}

grid <- expand_grid(
  # max_depth = c(2, 4, 6),
  max_depth = c(2,4,6),
  nb_iter_total = 400,
  eta = 0.3
) |>
  mutate(ind = row_number())

repns_vector <- 1:100

simulate_xgb_scenario <- function(scenario, params_df, repn) {
  # Generate Data
  simu_data <- simulate_data_wrapper(
    scenario = scenario,
    params_df = params_df,
    repn = repn
  )

  # Looping over the grid hyperparameters for the scenario
  res_simul <- vector(mode = "list", length = nrow(grid))
  cli::cli_progress_bar("Iteration grid", total = nrow(grid), type = "tasks")
  for (j in 1:nrow(grid)) {
    curent_params <- grid |> dplyr::slice(!!j)
    res_simul[[j]] <- simul_xgb(
      params = curent_params,
      ind = curent_params$ind,
      simu_data = simu_data
    )
    cli::cli_progress_update()
  }


  # The metrics computed for all set of hyperparameters (identified with `ind`)
  # and for each number of boosting iterations (`nb_iter`), for the current
  # scenario (`scenario`) and current replication number (`repn`)
  metrics_simul <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "tb_metrics") |> list_rbind()
  ) |>
    list_rbind()

  # Sanity check
  # metrics_simul |> count(scenario, repn, ind, sample, nb_iter) |>
  #   filter(n > 1)

  # P(q_1<s(x)<q_2)
  prop_scores_simul <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "tb_prop_scores") |> list_rbind()
  ) |>
    list_rbind()

  # Sanity check
  # prop_scores_simul |> count(scenario, repn, ind, sample, nb_iter)

  # Histogram of estimated scores
  scores_hist <- map(
    res_simul,
    function(simul_grid_j) map(simul_grid_j, "scores_hist")
  )

  list(
    metrics_simul = metrics_simul,
    scores_hist = scores_hist,
    prop_scores_simul = prop_scores_simul
  )
}

