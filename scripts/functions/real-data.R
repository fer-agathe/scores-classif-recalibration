# Dependencies:
# tidyverse, gam, gamsel, recipes, MASS, latex2exp
# functions/metrics.R

#' Split dataset into train and test set
#'
#' @param data dataset
#' @param prop_train proportion in the train test (default to .8)
#' @param seed desired seed (default to `NULL`)
#'
#' @returns a list with two elements: the train set, the test set
split_train_test <- function(data,
                             prop_train = .8,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  size_train <- round(prop_train * nrow(data))
  ind_sample <- sample(1:nrow(data), replace = FALSE, size = size_train)

  list(
    train = data |> dplyr::slice(ind_sample),
    test = data |> dplyr::slice(-ind_sample)
  )
}

#' Split prior scores into train and test set
#'
#' @param prior scores
#' @param prop_train proportion in the train test (default to .8)
#' @param seed desired seed (default to `NULL`)
#'
#' @returns a list with two elements: the train set, the test set of prior scores
split_train_test_prior <- function(prior,
                                   prop_train = .8,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  size_train <- round(prop_train * length(prior))
  ind_sample <- sample(1:length(prior), replace = FALSE, size = size_train)

  list(
    train = prior[ind_sample],
    test = prior[-ind_sample]
  )
}

#' One-hot encoding, and renaming variables to avoid naming that do not respect
#' r old naming conventions
#'
#' @param data_train train set
#' @param data_test test set
#' @param target_name name of the target (response) variable
#' @param intercept should a column for an intercept be added? Default to
#'  `FALSE`
#'
#' @returns list with five elements:
#'  - `train`: train set
#'  - `test`: test set
#'  - `initial_covariate_names`: vector of names of all explanatory variables
#'  - `categ_names`: vector of new names of categorical variables (if any)
#'  - `covariate_names`: vector of new names of all explanatory variables (including
#'     categorical ones).
encode_dataset <- function(data_train,
                           data_test,
                           target_name,
                           intercept = FALSE) {

  col_names <- colnames(data_train)
  col_names_covariates <- col_names[-which(col_names == target_name)]
  new_names_covariates <- str_c("X_", 1:length(col_names_covariates))
  data_train <- data_train |>
    rename_with(.cols = all_of(col_names_covariates), .fn = ~new_names_covariates)
  data_test <- data_test |>
    rename_with(.cols = all_of(col_names_covariates), .fn = ~new_names_covariates)

  data_rec <- recipes::recipe(
    formula(str_c(target_name, " ~ .")),
    data = data_train
  )

  ref_cell <- data_rec |> recipes::step_dummy(
    recipes::all_nominal(), -recipes::all_outcomes(),
    one_hot = TRUE
  ) |>
    recipes::prep(training = data_train)

  X_train_dmy <- recipes::bake(ref_cell, new_data = data_train)
  X_test_dmy  <- recipes::bake(ref_cell, new_data = data_test)

  # Identify categorical variables
  # Bake the recipe to apply the transformation
  df_transformed <- recipes::bake(ref_cell, new_data = NULL)
  # Get the names of the transformed data
  new_names <- names(X_train_dmy)
  original_vars <- names(data_train)
  categ_names <- setdiff(new_names, original_vars)
  covariate_names <- colnames(X_train_dmy)
  covariate_names <- covariate_names[!covariate_names == target_name]
  categ_names <- categ_names[!categ_names == target_name]
  list(
    train = X_train_dmy,
    test = X_test_dmy,
    initial_covariate_names = col_names_covariates,
    categ_names = categ_names,
    covariate_names = covariate_names
  )
}

## Estimation Functions----
#' Train a GAMSEL model
#'
#' @param data_train train set
#' @param data_test test set
#' @param target_name name of the target (response) variable
#' @param degrees degree for the splines
#' @param return_model if TRUE, the estimated model is returned
#'
#' @returns list with estimated scores on train set (`scores_train`) and on
#'  test set (`scores_test`)
train_gamsel <- function(data_train,
                         data_test,
                         target_name,
                         degrees = 6,
                         return_model = FALSE) {
  # Encode dataset so that categorical variables become dummy variables
  data_dmy <- encode_dataset(
    data_train = data_train,
    data_test = data_test,
    target_name = target_name,
    intercept = FALSE
  )
  # Estimation
  X_dmy_train <- data_dmy$train |> dplyr::select(-!!target_name)
  X_dmy_train <- X_dmy_train |> mutate(across(everything(), as.numeric))
  X_dmy_test <- data_dmy$test |> dplyr::select(-!!target_name)
  X_dmy_test <- X_dmy_test |> mutate(across(everything(), as.numeric))
  y_train <- data_dmy$train |> dplyr::pull(!!target_name)
  y_test <- data_dmy$test |> dplyr::pull(!!target_name)

  deg <- rep(NA, ncol(X_dmy_train))
  col_names_X <- colnames(X_dmy_train)
  nb_val <- map_dbl(
    col_names_X, ~X_dmy_train |> pull(.x) |> unique() |> length()
  )
  for (i_var_name in 1:ncol(X_dmy_train)) {
    var_name <- col_names_X[i_var_name]
    if (var_name %in% data_dmy$categ_names) {
      deg[i_var_name] <- 1
    } else {
      deg[i_var_name] <- min(nb_val[i_var_name]-1, degrees)
    }
  }
  gamsel_cv <- gamsel::cv.gamsel(
    x = as.data.frame(X_dmy_train), y = y_train, family = "binomial",
    degrees = deg
  )
  gamsel_out <- gamsel::gamsel(
    x = as.data.frame(X_dmy_train), y = y_train, family = "binomial",
    degrees = deg,
    lambda = gamsel_cv$lambda.min
  )
  # Scores on train and test set
  scores_train <- predict(
    gamsel_out, newdata = as.data.frame(X_dmy_train), type = "response")[, 1]
  scores_test <- predict(
    gamsel_out, newdata = as.data.frame(X_dmy_test), type = "response")[, 1]
  scores_train[which(is.na(scores_train))] <-
    1/(1 + exp(-predict(gamsel_out,
                        newdata = as.data.frame(X_dmy_train[which(is.na(scores_train)),]))
               [, 1]))
  scores_test[which(is.na(scores_test))] <-
    1/(1 + exp(-predict(gamsel_out,
                        newdata = as.data.frame(X_dmy_test[which(is.na(scores_test)),]))
               [, 1]))

  if (return_model == TRUE) {
    res <- list(
      scores_train = scores_train,
      scores_test = scores_test,
      fit = fit)
  } else {
    list(scores_train = scores_train, scores_test = scores_test, fit = NULL)
  }
}

## Distribution of scores----

#' Maximum-likelihood fitting of Beta distribution on scores
#'
#' @param scores vector of estimated scores
#' @param shape1 non-negative first parameter of the Beta distribution
#' @param shape1 non-negative second parameter of the Beta distribution
#'
#' @returns An object of class `fitdistr`, a list with four components
#'  (see: MASS::fitdistr())
#'  - `estimate`: the parameter estimates
#'  - `sd`: the estimated standard errors
#'  - `vcov`: the estimated variance-covariance matrix
#'  - `loglik`: the log-likelihood
fit_beta_scores <- function(scores, shape1 = 1, shape2 = 1) {
  # Fit a beta distribution
  mle_fit <- MASS::fitdistr(
    scores, "beta", start = list(shape1 = 1, shape2 = 1)
  )
  mle_fit
}

#' Estimation of a GLM-logistic, a GAM and a GAMSEL model on a classification
#' task. Then, on estimated scores from the test set, fits a Beta distribution.
#'
#' @param dataset dataset with response variable and predictors
#' @param target_name name of the target (response) variable
#' @param seed desired seed (default to `NULL`)
#'
#' @returns A list with the following elements:
#'  - `scores_glm`: scores on train and test set (in a list) from the GLM
#'  - `scores_gam`: scores on train and test set (in a list) from the GAM
#'  - `scores_gamsel`: scores on train and test set (in a list) from the GAMSEL
#'  - `mle_glm`: An object of class "fitdistr" for the GLM model
#'    (see fit_beta_scores())
#'  - `mle_gamsel`: An object of class "fitdistr" for the GAM model
#'    (see fit_beta_scores())
#'  - `mle_gamsel`: An object of class "fitdistr" for the GAMSEL model
#'    (see fit_beta_scores())
get_beta_fit <- function(dataset,
                         target_name,
                         seed = NULL) {
  # Split data into train/test
  data <- split_train_test(data = dataset, prop_train = .7, seed = seed)

  # Train a GAMSEL model
  scores_gamsel <- train_gamsel(
    data_train = data$train, data_test = data$test, target_name = target_name,
    degrees = 6
  )
  # Add a little noise to the estimated scores to avoid being in [0,1] and be
  # in (0,1) instead.
  x_gamsel <- (scores_gamsel$scores_test * (1 - 1e-6)) + 1e-6 / 2
  # Fit a Beta distribution on these scores
  mle_gamsel <- fit_beta_scores(scores = x_gamsel[!is.nan(x_gamsel)])

  list(
    scores_gamsel = scores_gamsel,
    mle_gamsel = mle_gamsel
  )
}

## Plots----

#' Plots the histogram of scores estimated with GAMSEL
#' add densities of Beta distribution for whith the parameters have been
#' estimated using scores from the GLM, the GAM, or the GAMSEL model
#'
#' @param fit_resul results obtained from get_beta_fit()
#' @param title title of the graph (e.g.: dataset name)
plot_hist_scores_beta <- function(fit_resul, title = NULL) {
  val_u <- seq(0, 1, length = 651)
  layout(mat = matrix(1:2), heights = c(3,1))
  type <- "mle_gamsel"
  dbeta_val <- dbeta(
    val_u,
    fit_resul[[type]]$estimate[1],
    fit_resul[[type]]$estimate[2]
  )
  y_lim <- c(
    0, range(dbeta_val[!is.infinite(dbeta_val)], na.rm = TRUE)
    |> max(na.rm = TRUE)
  )

  # Histogram of scores obtained with the GAMSEL, on test set
  par(mar = c(4.1, 4.1, 1, 2.1))
  hist(
    fit_resul$scores_gamsel$scores_test,
    breaks = seq(0, 1, by = .05), probability = TRUE,
    main = title, xlab = latex2exp::TeX("$\\hat{s}(x)$"),
    ylim = y_lim
  )

  legend_name <- "GAMSEL"
  colour <- "#E69F00"
  type <- "mle_gamsel"
  lines(
    val_u,
    dbeta_val,
    col = colour,lwd = 1.5
  )
  par(mar = c(0, 4.1, 0, 2.1))
  plot.new()
  legend(
    xpd = TRUE, ncol = 1,
    "center",
    title = "Model",
    lwd = 1.5,
    col = colour,
    legend = legend_name
  )
}
