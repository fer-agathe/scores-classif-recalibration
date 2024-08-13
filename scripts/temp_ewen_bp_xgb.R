
scenario <- 1
repn <- 1
# max_depth <- 2

plot_bp_interest <- function(metrics_interest,
                             scores_hist_interest,
                             label,
                             recalib_method) {
  subtitle <- str_c(
    "Depth = ", metrics_interest$max_depth, ", ",
    "MSE = ", round(metrics_interest$mse, 2), ", ",
    "AUC = ", round(metrics_interest$AUC, 2), ", \n",
    "Brier = ", round(metrics_interest$brier, 2), ",",
    "ICI = ", round(metrics_interest$ici, 2), ", ",
    "KL = ", round(metrics_interest$KL_20_true_probas, 2)
  )

  if (recalib_method == "none") {
    plot(
      main = str_c(label, " (iter = ", metrics_interest$nb_iter,")"),
      scores_hist_interest$test,
      xlab = latex2exp::TeX("$\\hat{s}(x)$"),
      ylab = ""
    )
  } else if (recalib_method == "platt") {
    plot(
      main = str_c(label, " (iter = ", metrics_interest$nb_iter,")"),
      scores_hist_interest$test_c_platt,
      xlab = latex2exp::TeX("$\\hat{s}(x)$"),
      ylab = "",
      col = colour_recalib[["Platt"]]
    )
  } else if (recalib_method == "iso") {
    plot(
      main = str_c(label, " (iter = ", metrics_interest$nb_iter,")"),
      scores_hist_interest$test_c_iso,
      xlab = latex2exp::TeX("$\\hat{s}(x)$"),
      ylab = "",
      col = colour_recalib[["Isotonic"]]
    )
  }
  mtext(side = 3, line = -0.25, adj = .5, subtitle, cex = .5)
}





plot_bp_xgb <- function(scenario,
                        repn,
                        paper_version = FALSE) {
  # Focus on current scenario
  scores_hist_scenario <- scores_hist_all[[scenario]]
  # Focus on a particular replication
  scores_hist_repn <- scores_hist_scenario[[repn]]
  # # Focus on a value for max_depth
  max_depth_val <- map_dbl(scores_hist_repn, ~.x[[1]]$max_depth)
  # i_max_depth <- which(max_depth_val == max_depth)
  # scores_hist <- scores_hist_repn[[i_max_depth]]

  # True Probabilities
  simu_data <- simulate_data_wrapper(
    scenario = scenario,
    params_df = params_df,
    repn = repn # only one replication here
  )
  true_prob <- simu_data$data$probs_train

  for (recalib_method in c("none", "platt", "iso")) {

    i_method <- match(recalib_method, c("none", "platt", "iso"))
    recalib_method_lab <- c("None", "Platt", "Isotonic")[i_method]

    # The metrics for the corresponding simulations, on the validation set
    metrics_xgb_current_valid <-
      metrics_xgb_all |>
      filter(
        scenario == !!scenario,
        repn == !!repn,
        sample == "Validation",
        recalib == "None"
      )
    # and on the test set
    metrics_xgb_current_test <-
      metrics_xgb_all |>
      filter(
        scenario == !!scenario,
        repn == !!repn,
        sample == "Test",
        recalib == recalib_method_lab
      )

    if (paper_version == FALSE) {
      hist(
        true_prob,
        breaks = seq(0, 1, by = .05),
        xlab = "p", ylab = "",
        main = "True Probabilities",
        xlim = c(0, 1)
      )
      mtext(
        side = 2, recalib_method_lab, line = 2.5, cex = 1,
        font.lab = 2
      )
      # Iterations of interest----
      ## Start of iterations
      scores_hist_start <- scores_hist_repn[[1]][[1]]
      metrics_start <- metrics_xgb_current_test |>
        filter(
          nb_iter == scores_hist_start$nb_iter,
          max_depth == scores_hist_start$max_depth
        )

      plot_bp_interest(
        metrics_interest = metrics_start,
        scores_hist_interest = scores_hist_start,
        label = "Start",
        recalib_method = recalib_method
      )

      ## End of iterations
      scores_hist_end <- scores_hist_repn[[1]][[length(scores_hist_repn[[1]])]]
      metrics_end <- metrics_xgb_current_test |>
        filter(
          nb_iter == scores_hist_end$nb_iter,
          max_depth == scores_hist_start$max_depth
        )
      plot_bp_interest(
        metrics_interest = metrics_end,
        scores_hist_interest = scores_hist_end,
        label = "End",
        recalib_method = recalib_method
      )

      ## Iteration with min MSE on validation set
      metrics_valid_mse_star <- metrics_xgb_current_valid |> arrange(mse) |>
        dplyr::slice(1)
      nb_iter_mse <- metrics_valid_mse_star$nb_iter
      max_depth_mse_star <- metrics_valid_mse_star$max_depth
      i_max_depth_mse_star <- which(max_depth_val == max_depth_mse_star)
      # Metrics at the same iteration on the test set
      metrics_min_mse <-
        metrics_xgb_current_test |>
        filter(
          nb_iter == !!nb_iter_mse,
          max_depth == max_depth_mse_star
        )
      # Note: indexing at 0 here...
      ind_mse <- which(map_dbl(scores_hist_repn[[i_max_depth_mse_star]], "nb_iter") == nb_iter_mse)
      scores_hist_min_mse <- scores_hist_repn[[i_max_depth_mse_star]][[ind_mse]]
      plot_bp_interest(
        metrics_interest = metrics_min_mse,
        scores_hist_interest = scores_hist_min_mse,
        label = "MSE*",
        recalib_method = recalib_method
      )
    }
    ## Iteration with max AUC on validation set
    metrics_valid_auc_star <-
      metrics_xgb_current_valid |> arrange(desc(AUC)) |> dplyr::slice(1)
    nb_iter_auc <- metrics_valid_auc_star$nb_iter
    max_depth_auc_star <- metrics_valid_auc_star$max_depth
    i_max_depth_auc_star <- which(max_depth_val == max_depth_auc_star)

    metrics_max_auc <-
      metrics_xgb_current_test |>
      filter(nb_iter == !!nb_iter_auc, max_depth == max_depth_auc_star)
    # Note: indexing at 0 here...
    ind_auc <- which(map_dbl(scores_hist_repn[[i_max_depth_auc_star]], "nb_iter") == nb_iter_auc)
    scores_hist_max_auc <- scores_hist_repn[[i_max_depth_auc_star]][[ind_auc]]
    plot_bp_interest(
      metrics_interest = metrics_max_auc,
      scores_hist_interest = scores_hist_max_auc,
      label = "AUC*",
      recalib_method = recalib_method
    )
    if (paper_version == TRUE) {
      mtext(
        side = 2, recalib_method_lab, line = 2.5, cex = 1,
        font.lab = 2
      )
    }

    ## Min Brier on validation set
    metrics_valid_brier_star <-
      metrics_xgb_current_valid |> arrange(brier) |> dplyr::slice(1)
    nb_iter_brier <- metrics_valid_brier_star$nb_iter
    max_depth_brier_star <- metrics_valid_brier_star$max_depth
    i_max_depth_brier_star <- which(max_depth_val == max_depth_brier_star)

    metrics_min_brier <-
      metrics_xgb_current_test |>
      filter(nb_iter == !!nb_iter_brier, max_depth == max_depth_brier_star)
    ind_brier <- which(map_dbl(scores_hist_repn[[i_max_depth_brier_star]], "nb_iter") == nb_iter_brier)
    scores_hist_min_brier <- scores_hist_repn[[i_max_depth_brier_star]][[ind_brier]]
    plot_bp_interest(
      metrics_interest = metrics_min_brier,
      scores_hist_interest = scores_hist_min_brier,
      label = "Brier*",
      recalib_method = recalib_method
    )

    ## Min ICI on validation set
    metrics_valid_ici_star <-
      metrics_xgb_current_valid |> arrange(ici) |> dplyr::slice(1)
    nb_iter_ici <-   metrics_valid_ici_star$nb_iter
    max_depth_ici_star <- metrics_valid_ici_star$max_depth
    i_max_depth_ici_star <- which(max_depth_val == max_depth_ici_star)

    metrics_min_ici <-
      metrics_xgb_current_test |>
      filter(nb_iter == !!nb_iter_ici, max_depth == max_depth_ici_star)
    ind_ici <- which(map_dbl(scores_hist_repn[[i_max_depth_ici_star]], "nb_iter") == nb_iter_ici)
    scores_hist_min_ici <- scores_hist_repn[[i_max_depth_ici_star]][[ind_ici]]
    plot_bp_interest(
      metrics_interest = metrics_min_ici,
      scores_hist_interest = scores_hist_min_ici,
      label = "ICI*",
      recalib_method = recalib_method
    )

    ## Min KL on validation set
    metrics_valid_kl_star <-
      metrics_xgb_current_valid |> arrange(KL_20_true_probas) |> dplyr::slice(1)
    nb_iter_kl <-   metrics_valid_kl_star$nb_iter
    max_depth_kl_star <- metrics_valid_kl_star$max_depth
    i_max_depth_kl_star <- which(max_depth_val == max_depth_kl_star)

    metrics_min_kl <-
      metrics_xgb_current_test |>
      filter(nb_iter == !!nb_iter_kl, max_depth == max_depth_kl_star)
    ind_kl <- which(map_dbl(scores_hist_repn[[i_max_depth_kl_star]], "nb_iter") == nb_iter_kl)
    scores_hist_min_kl <- scores_hist_repn[[i_max_depth_kl_star]][[ind_kl]]
    plot_bp_interest(
      metrics_interest = metrics_min_kl,
      scores_hist_interest = scores_hist_min_kl,
      label = "KL*",
      recalib_method = recalib_method
    )
  }
}


for (scenario in 1:16) {
  pdf(
    file = str_c("../figs/bp_synthetic_xbg_", scenario, ".pdf"),
    height = 4.5, width = 10
  )
  par(mfrow = c(3,4), mar = c(4.1, 4, 3.5, 1.5))
  plot_bp_xgb(scenario = scenario, repn = 1, paper_version = TRUE)
  dev.off()
}


# Suite avec graphiques ----

df_plot <-
  metrics_xgb_all |>
  mutate(
    dgp = case_when(
      scenario %in% 1:4 ~ 1,
      scenario %in% 5:8 ~ 2,
      scenario %in% 9:12 ~ 3,
      scenario %in% 13:16 ~ 4
    ),
    dgp = factor(dgp, levels = 1:4, labels = str_c("DGP ", 1:4)),
    no_noise = c(0, 10, 50, 100)[(scenario-1)%%4 + 1],
    no_noise = factor(
      no_noise, levels = c(no_noise),
      labels = str_c(no_noise, " noise variables")
    )
  ) |>
  select(
    dgp, no_noise, scenario, recalib, ind, sample, nb_iter, max_depth,
    brier, ici, KL_20_true_probas
  ) |>
  group_by(dgp, no_noise, scenario, recalib, ind, sample, nb_iter, max_depth) |>
  summarise(
    brier = mean(brier),
    ici = mean(ici),
    KL_20_true_probas = mean(KL_20_true_probas),
    .groups = "drop"
  ) |>
  mutate(
    max_depth = factor(
      max_depth,
      levels = c(2, 4, 6)
    )
  )

formatter1000 <- function(x) x*1000


p_brier <- ggplot(
  data = df_plot |> arrange(nb_iter) |> filter(max_depth == 2),
  mapping = aes(x = brier, y = KL_20_true_probas)
) +
  geom_path(
    mapping = aes(colour = sample, linetype = recalib),
    arrow = arrow(type = "closed", ends = "last",
                  length = unit(0.08, "inches"))
  ) +
  # facet_wrap(~scenario) +
  ggh4x::facet_grid2(dgp~no_noise, scales = "free_y", independent = "y") +
  labs(
    x = latex2exp::TeX("Calibration (Brier), $\\times 10^{3}$, log scale"),
    y = "KL Divergence"
  ) +
  scale_x_log10(labels = formatter1000) + scale_y_log10() +
  scale_colour_manual("Sample", values = colour_samples) +
  scale_linetype_discrete("Recalibration") +
  theme_paper() +
  theme(legend.key.width = unit(1.5, "cm"))


ggsave(p_brier, file = "../figs/xgb-kl-calib-brier-leaves-all.pdf",
       width = 10, height = 8)

p_brier


p_ici <- ggplot(
  data = df_plot |> arrange(nb_iter) |> filter(max_depth == 2),
  mapping = aes(x = ici, y = KL_20_true_probas)
) +
  geom_path(
    mapping = aes(colour = sample, linetype = recalib),
    arrow = arrow(type = "closed", ends = "last",
                  length = unit(0.08, "inches"))
  ) +
  # facet_wrap(~scenario) +
  ggh4x::facet_grid2(dgp~no_noise, scales = "free_y", independent = "y") +
  labs(
    x = latex2exp::TeX("Calibration (ICI), $\\times 10^{3}$, log scale"),
    y = "KL Divergence"
  ) +
  scale_x_log10(labels = formatter1000) + scale_y_log10() +
  scale_colour_manual("Sample", values = colour_samples) +
  scale_linetype_discrete("Recalibration") +
  theme_paper() +
  theme(legend.key.width = unit(1.5, "cm"))


ggsave(p_ici, file = "../figs/xgb-kl-calib-ici-leaves-all.pdf",
       width = 10, height = 8)

p_ici



# Tableaux----

# Identify the smallest tree on the validation set, when the scores are not
# recalibrated
smallest_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(nb_iter) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "smallest") |>
  ungroup()

# Identify the largest tree
largest_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(desc(nb_iter)) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "largest") |>
  ungroup()

# Identify tree with highest AUC on test set
highest_auc_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(desc(AUC)) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "largest_auc") |>
  ungroup()

# Identify tree with lowest MSE
lowest_mse_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(mse) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "lowest_mse") |>
  ungroup()

# Identify tree with lowest brier
lowest_brier_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(brier) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "lowest_brier") |>
  ungroup()

# Identify tree with lowest ICI
lowest_ici_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(ici) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "lowest_ici") |>
  ungroup()

# Identify tree with lowest KL
lowest_kl_xgb <-
  metrics_xgb_all |>
  filter(sample == "Validation", recalib == "None") |>
  group_by(scenario, repn) |>
  arrange(KL_20_true_probas) |>
  slice_head(n = 1) |>
  select(scenario, repn, ind, nb_iter, recalib) |>
  mutate(result_type = "lowest_kl") |>
  ungroup()

# Merge these
models_of_interest_xgb <-
  smallest_xgb |>
  bind_rows(largest_xgb) |>
  bind_rows(highest_auc_xgb) |>
  bind_rows(lowest_mse_xgb) |>
  bind_rows(lowest_brier_xgb) |>
  bind_rows(lowest_ici_xgb) |>
  bind_rows(lowest_kl_xgb)

models_of_interest_metrics <- NULL
for (recalibration_method in c("None", "Platt", "Isotonic")) {
  # Add metrics now
  models_of_interest_metrics <-
    models_of_interest_metrics |>
    bind_rows(
      models_of_interest_xgb |> select(-recalib) |>
        left_join(
          metrics_xgb_all |>
            filter(
              recalib == recalibration_method,
              sample %in% c("Validation", "Test")
            ),
          by = c("scenario", "repn", "ind", "nb_iter"),
          relationship = "many-to-many" # (calib, test)
        )
    )
}


models_of_interest_metrics <-
  models_of_interest_metrics |>
  mutate(
    result_type = factor(
      result_type,
      levels = c(
        "smallest", "largest", "lowest_mse", "largest_auc",
        "lowest_brier", "lowest_ici", "lowest_kl"),
      labels = c(
        "Smallest", "Largest", "MSE*", "AUC*",
        "Brier*", "ICI*", "KL*"
      )
    )
  )
