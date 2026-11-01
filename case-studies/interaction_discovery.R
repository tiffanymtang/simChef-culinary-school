# simChef case study: interaction discovery with iRF
# Run from the project root with source("case-studies/interaction_discovery.R")

rm(list = ls())
library(here)
library(simChef)
library(randomForest)
library(iRF)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyselect)
library(yardstick)
library(rlang)
library(future)
set.seed(331)

## Initial DGP setup (only need to run once, rather than for every simulation replicate)

# Load in enhancer data
load(here::here("case-studies", "data", "enhancer.Rdata")) #> X, Y, varnames.all, (train.id, test.id)
tf_feature_idxs <- 46:80 # keep TF features only
X <- X[, tf_feature_idxs]

# Get the variable groups for the TF features to use in iRF
var_groups <- gsub("_", ".", varnames.all[tf_feature_idxs, 2])

# Helper function to split observations within each class so train/test retain the same class balance.
stratified_split <- function(y, train_prop) {
  train <- lapply(
    split(seq_along(y), y), 
    function(idx) {
      sample(idx, floor(length(idx) * train_prop))
    }
  ) |> 
    unlist(use.names = FALSE)
  list(train = sort(train), test = setdiff(seq_along(y), train))
}

# Estimate RF-noise probabilities (to generate realistic noise in the response variable)
estimate_rf_probs <- function(X, y, n_splits = 20, ...) {
  probabilities <- matrix(NA_real_, nrow(X), n_splits)
  for (b in seq_len(n_splits)) {
    split <- stratified_split(y, train_prop = 0.5)
    train_ids <- split$train
    test_ids <- split$test
    fit <- randomForest::randomForest(
      x = X[train_ids, , drop = FALSE],
      y = factor(y[train_ids], levels = c(0, 1)),
      ...
    )
    probabilities[test_ids, b] <- predict(
      fit, X[test_ids, , drop = FALSE], type = "prob"
    )[, "1"]
  }
  return(rowMeans(probabilities, na.rm = TRUE))
}
rf_probs <- estimate_rf_probs(X, Y)

## Define DGP(s) for simulation study ------------------------------------------

#' Simulate X and y data
#' 
#' @param X A data frame of covariates
#' @param rf_probs A vector of RF-noise probabilities for each observation
#' @param n The number of training observations to sample
#' 
#' @return A list containing
#' - X_train: training covariates
#' - y_train: training response
#' - X_test: test covariates
#' - y_test: test response
#' - var_groups: variable groups
#' - true_interaction: the true interaction terms
generate_xy_data <- function(X, rf_probs, n) {
  
  # simulate f(x) from LSS function
  fx <- (X[, "kr1"] > 1.25 | X[, "kr2"] > 1.25) &
    (X[, "twi1"] > 1.25 | X[, "twi2"] > 1.25) & 
    (X[, "D1"] > 1.25) &
    (X[, "hb1"] > 1.25 | X[, "hb2"] > 1.25) & 
    (X[, "wt_ZLD"] > 75)
  
  # incorporate RF-noise into the response variable
  y <- as.integer(fx | (runif(nrow(X)) < rf_probs))
  
  # split the data into train and test sets
  split <- stratified_split(y, train_prop = n / nrow(X))
  train_ids <- split$train
  test_ids <- split$test

  return(list(
    X_train = X[train_ids, , drop = FALSE],
    y_train = factor(y[train_ids], levels = c(0, 1)),
    X_test = X[test_ids, , drop = FALSE],
    y_test = factor(y[test_ids], levels = c(0, 1)),
    true_interaction = c("D1", "hb", "kr", "twi", "wt.ZLD")
  ))
}

# Create the DGP object for the simulation study
dgp <- simChef::create_dgp(
  .dgp_fun = generate_xy_data,
  .name = "Enhancer LSS with RF Noise",
  # additional arguments to pass to the DGP function
  X = X, rf_probs = rf_probs, n = 400
)

## Define Method(s) for simulation study ---------------------------------------

#' Fit iRF
#' 
#' @param X_train Training covariates
#' @param y_train Training response
#' @param X_test Test covariates
#' @param y_test Test response
#' @param true_interaction The true interaction terms
#' @param niter The number of iterations for iRF
#' @param ... Additional arguments to pass to the iRF function
#' 
#' @return A list containing the fitted iRF model and the test response
#' 
#' @details
#' The inputs should include named parameters that match the outputs of the 
#' DGP function. 
fit_irf <- function(X_train, y_train, X_test, y_test, 
                    true_interaction, niter, ...) {
  fit <- iRF::iRF(
    x = X_train, y = y_train, 
    xtest = X_test, ytest = y_test,
    iter.return = niter, int.return = niter,
    signed = FALSE, verbose = FALSE, n.core = 1,
    ...
  )
  return(list(
    rf_fit = fit$rf.list,
    interactions = fit$interaction,
    # also output objects necessary for evaluation
    y_test = y_test,
    true_interaction = true_interaction
  ))
}

# Create the method object for the simulation study
irf <- simChef::create_method(
  .method_fun = fit_irf,
  .name = "iRF",
  # additional arguments to pass to the method function
  varnames.grp = var_groups,
  niter = 5, 
  n.bootstrap = 20, 
  rit.param = list(depth = 5, ntree = 100, nchild = 2, class.id = 1, min.nd = 1)
)

## Define Evaluators for simulation study --------------------------------------

#' Evaluator wrapper function
#' 
#' @param fit_results A data frame of fit results (automatically inputted by `simChef`)
#' @param vary_params A character vector of parameters to vary across (automatically inputted by `simChef`)
#' @param eval_id The name of the evaluation metric to compute
#' @param eval_fun A function that computes the evaluation metric
eval_summary <- function(fit_results, vary_params = NULL, eval_id, eval_fun) {
  group_ids <- c(".dgp_name", ".method_name", vary_params, "interaction_order", "active")
  eval_out <- fit_results |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      {{eval_id}} := eval_fun(
        rf_fit = rf_fit, 
        interactions = interactions, 
        y_test = y_test, 
        true_interaction = true_interaction
      )
    ) |>
    tidyr::unnest(cols = {{eval_id}}) |>
    dplyr::group_by(dplyr::across(dplyr::any_of(group_ids))) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(eval_id), 
        list(mean = mean, raw = list), 
        .names = "{.fn}_{.col}"
      ),
      .groups = "drop"
    )
  if ("n" %in% colnames(eval_out)) {
    eval_out <- eval_out |> 
      dplyr::mutate(n = factor(n, levels = sort(unique(n))))
  }
  return(eval_out)
}

#' Evaluate prediction performance using AUPRC
eval_prediction <- function(rf_fit, y_test, ...) {
  yardstick::pr_auc_vec(
    y_test, rf_fit$test$votes[, "1"], event_level = "second"
  )
}

#' Evaluate interaction recovery by interaction order
eval_interaction_recovery <- function(interactions, true_interaction, ...) {
  if (is.null(interactions) || !nrow(interactions)) {
    return(list(data.frame(
      interaction_order = 2:length(true_interaction),
      recovery_rate = 0
    )))
  }
  lapply(
    2:length(true_interaction),
    function(order) {
      # get all possible true interactions of the given order
      targets <- combn(
        true_interaction, order,
        FUN = function(x) paste(sort(x), collapse = "_")
      )
      # check if any of the found interactions match the true interactions
      recovered <- sapply(
        targets,
        function(target) {
          target_features <- strsplit(target, "_", fixed = TRUE)[[1]]
          any(
            sapply(
              interactions$int, 
              function(candidate) {
                all(target_features %in% strsplit(candidate, "_", fixed = TRUE)[[1]])
              }
            )
          )
        }
      )
      data.frame(
        interaction_order = order,
        recovery_rate = mean(recovered)
      )
    }
  ) |> 
    dplyr::bind_rows() |> 
    list()
}

#' Evaluate stability of active (true) and inactive (false) interactions
eval_stability <- function(interactions, true_interaction, ...) {
  if (is.null(interactions) || !nrow(interactions)) {
    return(
      list(data.frame(
        interaction_order = integer(), 
        active = logical(),
        stability = numeric()
      ))
    )
  }
  features <- strsplit(interactions$int, "_", fixed = TRUE)
  list(data.frame(
    interaction_order = lengths(features),
    active = sapply(features, function(x) all(x %in% true_interaction)),
    stability = interactions$stability
  ))
}

# Create the evaluator objects for the simulation study
prediction <- simChef::create_evaluator(
  .eval_fun = eval_summary,
  .name = "AUPRC",
  eval_id = "auprc",
  eval_fun = eval_prediction
)
interaction_recovery <- simChef::create_evaluator(
  .eval_fun = eval_summary,
  .name = "Interaction Recovery",
  eval_id = "recovery_rate",
  eval_fun = eval_interaction_recovery
)
stability <- simChef::create_evaluator(
  .eval_fun = eval_summary,
  .name = "Interaction Stability",
  eval_id = "stability",
  eval_fun = eval_stability
)

## Define Visualizers for simulation study -------------------------------------

#' Plot prediction performance (AUPRC) across varying parameters
plot_prediction <- function(fit_results = NULL, eval_results, vary_params = NULL,
                            eval_name = "AUPRC") {
  plt_df <- eval_results[[eval_name]] |> 
    tidyr::unnest(cols = raw_auprc) |> 
    dplyr::select(
      tidyselect::all_of(
        c(".dgp_name", ".method_name", vary_params, "raw_auprc")
      )
    )
  plt_df |> 
    ggplot2::ggplot() +
    ggplot2::aes(
      x = !!rlang::sym(vary_params[[1]]), 
      y = raw_auprc,
      fill = !!rlang::sym(vary_params[[1]])
    ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      x = vary_params[[1]], y = eval_name
    ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(fill = "none")
}

#' Plot interaction recovery across varying parameters
plot_interaction_recovery <- function(fit_results = NULL, eval_results, vary_params = NULL,
                                     eval_name = "Interaction Recovery") {
  plt_df <- eval_results[[eval_name]] |> 
    dplyr::select(
      tidyselect::all_of(
        c(".dgp_name", ".method_name", vary_params, "interaction_order", "mean_recovery_rate")
      )
    )
  plt_df |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = as.factor(interaction_order), 
      y = mean_recovery_rate, 
      fill = !!rlang::sym(vary_params[[1]])
    ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      x = "Interaction Order", y = eval_name, fill = vary_params[[1]]
    ) +
    ggplot2::theme_minimal()
}

#' Plot interaction stability across varying parameters
plot_stability <- function(fit_results = NULL, eval_results, vary_params = NULL,
                           eval_name = "Interaction Stability") {
  plt_df <- eval_results[[eval_name]] |> 
    tidyr::unnest(cols = raw_stability) |>
    dplyr::select(
      tidyselect::all_of(
        c(".dgp_name", ".method_name", vary_params, "interaction_order", "active", "raw_stability")
      )
    ) |> 
    dplyr::mutate(
      active = dplyr::case_when(
        active ~ "True Interaction",
        !active ~ "False Interaction"
      )
    )
  plt_df |> 
    ggplot2::ggplot() +
    ggplot2::aes(
      x = as.factor(interaction_order), 
      y = raw_stability, 
      fill = active
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(!!rlang::sym(vary_params[[1]])),
      labeller = purrr::partial(ggplot2::label_both, sep = " = ")
    ) +
    ggplot2::geom_boxplot(
      width = 0.7,
      position = ggplot2::position_dodge2(
        width = 0.8, preserve = "single"
      )
    ) +
    # ggplot2::scale_x_discrete(
    #   limits = as.character(2:5), drop = FALSE
    # ) +
    ggplot2::labs(
      x = "Interaction Order", y = eval_name, fill = ""
    ) +
    ggplot2::theme_minimal()
}

# Create the visualizer objects for the simulation study
prediction_plot <- create_visualizer(
  .viz_fun = plot_prediction,
  .name = "AUPRC Plot"
)
interaction_recovery_plot <- create_visualizer(
  .viz_fun = plot_interaction_recovery,
  .name = "Interaction Recovery Plot"
)
stability_plot <- create_visualizer(
  .viz_fun = plot_stability,
  .name = "Interaction Stability Plot"
)

## Create and run the simulations ----------------------------------------------
# create simulation recipe
experiment <- create_experiment(
  name = "Example iRF Simulation",
  save_dir = here::here("case-studies", "results", "Example iRF Simulation")
) |>
  add_dgp(dgp) |> 
  add_method(irf) |> 
  add_evaluator(prediction) |>
  add_evaluator(interaction_recovery) |>
  add_evaluator(stability) |>
  add_visualizer(prediction_plot) |>
  add_visualizer(interaction_recovery_plot) |>
  add_visualizer(stability_plot) |>
  add_vary_across(.dgp = "Enhancer LSS with RF Noise", n = c(400, 1200, 2000))

# Run the simulation study sequentially (this will take a while)
results <- run_experiment(
  experiment, n_reps = 20, save = TRUE, use_cached = TRUE
)

# # Run the simulation study in parallel
# future::plan(future::multisession, workers = 5) # adjust number of workers as needed
# results <- run_experiment(
#   experiment, n_reps = 20, save = TRUE,
#   future.globals = c("stratified_split")
# )

## Document simulation study ---------------------------------------------------

# Initialize documentation > Then need to fill out .md descriptions (can use coding agents to help with this)
init_docs(experiment)
# Render the documentation to HTML
render_docs(experiment)
