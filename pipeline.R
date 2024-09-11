library(readr)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(lubridate)
library(furrr)
library(progressr)
library(AdaptiveConformal)
library(opera)

source("R/learners.R")
source("R/loss.R")

apply_learners <- function(prevdata, newdata, learners, taus, prefix = "arrivals", progress = NULL) {
  id <- newdata$id[[1]]
  my_id <- newdata$my_id[[1]]
  date <- newdata$date_numeric[[1]]
  print(paste0("Starting: ", id[1]))
  
  results <- tibble(
    learner = names(learners)
  ) %>%
    mutate(preds = map(learner, function(learner) {
      filename <- glue::glue("raw_output_original/{learner}/{prefix}_{id}_{date}.rds")
      if(file.exists(filename)) {
        if(!is.null(progress)) progress()
        preds <- read_rds(filename)
      }
      else {
        tryCatch({
          preds <- learners[[learner]](prevdata, newdata)
          write_rds(preds, filename)
        }, error = function(cond) {
          message(cond)
          preds <- matrix(NA, ncol = length(taus), nrow = nrow(newdata))
        })
      }
      newdata %>% select(id, my_id, date, hospitalized, arrivals) %>% mutate(q0.1 = preds[,1], q0.5 = preds[,2], q0.9 = preds[,3])
    })) %>%
    unnest(c(preds)) %>%
    mutate(my_id = my_id)
  
  return(results)
}

setup_opera <- function(id, data, tau, learners, outcome, model, progress = NULL) {
  filename <- glue::glue("raw_output_original/{model}/{outcome}_{id}_{tau}.rds")
  if(file.exists(filename)) {
    if(!is.null(progress)) {
      progress()
    }
    return(read_rds(filename))
  }
  else {
    Y <- data[[outcome]]
    X <- as.matrix(data[, learners])
    
    print(glue::glue("Starting {model}"))
    
    mix <- opera::mixture(Y = Y, experts = X, model = model, loss.type = list(name = "pinball", tau = tau))
    
    write_rds(mix, filename)
    
    if(!is.null(progress)) {
      progress()
    }
    return(mix)
  }
}

get_mix_weights <- function(mix) {
  as_tibble(mix$weights)
}

online_aggregation <- function(test_set_preds, outcome_column) {
  setup <- test_set_preds %>%
    select(id, my_id, date, learner, !!outcome_column, q0.1, q0.5, q0.9) %>%
    pivot_longer(q0.1:q0.9, names_to = "quantile") %>%
    mutate(quantile = as.numeric(str_sub(quantile, 2, 4))) %>%
    mutate(value = inv_transform(value)) %>%
    group_by(id, quantile) %>%
    nest() %>%
    mutate(data = map(data, function(data) pivot_wider(data, names_from = "learner"))) %>%
    ungroup()
  
  learners <- unique(test_set_preds$learner)
  globals <- c("setup_opera")
  packages <- c("dplyr", "tidyr", "lubridate", "opera", "progressr", "readr")
  
  p <- NULL
  opera_results <- setup %>%
    mutate(
      ewa = furrr::future_pmap(list(id, data, quantile), setup_opera, outcome = outcome_column, learners = learners, model = "EWA", progress = p, .options = furrr_options(seed = TRUE, globals = globals, packages = packages)),
      boa = furrr::future_pmap(list(id, data, quantile), setup_opera, outcome = outcome_column, learners = learners, model = "BOA", progress = p, .options = furrr_options(seed = TRUE, globals = globals, packages = packages)),
      ewa_weights = map(ewa, get_mix_weights),
      boa_weights = map(boa, get_mix_weights),
      ewa_predictions = map(ewa, function(mix) tibble(ewa = mix$prediction[,1])),
      boa_predictions = map(boa, function(mix) tibble(boa = mix$prediction[,1])),
    )
  
  options <- furrr_options(seed = TRUE, globals = globals, packages = packages)
    opera_results <- opera_results %>%
      mutate(
        sl = furrr::future_pmap(list(id, data, quantile), setup_opera, outcome = outcome_column, learners = learners, model = "SuperLearner", progress = p, .options = options),
        sl_weights  = map(sl, get_mix_weights),
        sl_predictions  = map(sl, function(mix) tibble(sl = mix$prediction[,1]))
      )
  opera_results
}

plan(multisession, workers = 8)

print("Reading data...")
X_grouped <- read_rds("output/X_grouped.rds")
print("Finished reading data.")

options(scipen=999)

taus <- c(0.1, 0.5, 0.9)

#
# Step 1: base prediction algorithms
#

form_arrivals <- arrivals ~ arrivals_lag1 + arrivals_lag7 + arrivals_diff + 
  christmas + nye + ny + july14 + 
  year + month + wday + holiday + holiday_lag1 + 
  vacances_zone_c + flu_lagged + diarrhea_lagged + chickenpox_lagged + max_temp_lag1 + min_temp_lag1 + prcp_lag1

form_hospitalized <- hospitalized ~ hospitalized_lag1 + hospitalized_lag7 + hospitalized_diff + 
  christmas + nye + ny + july14 + 
  year + month + wday + holiday + holiday_lag1 + 
  vacances_zone_c + flu_lagged + diarrhea_lagged + chickenpox_lagged + max_temp_lag1 + min_temp_lag1 + prcp_lag1

form_arrivals_arimax     <- arrivals     ~ holiday + holiday_lag1 + vacances_zone_c + christmas + nye + ny + july14 + flu_lagged + diarrhea_lagged + chickenpox_lagged + max_temp_lag1 + min_temp_lag1 + prcp_lag1
form_hospitalized_arimax <- hospitalized ~ holiday + holiday_lag1 + vacances_zone_c + christmas + nye + ny + july14 + flu_lagged + diarrhea_lagged + chickenpox_lagged + max_temp_lag1 + min_temp_lag1 + prcp_lag1

form_arrivals_qgam     <- arrivals ~ arrivals_lag1 + arrivals_lag7 + arrivals_diff +
  christmas + nye + ny + july14 +
  year + month + wday + holiday + holiday_lag1 + vacances_zone_c + s(flu_lagged) + s(diarrhea_lagged) + s(chickenpox_lagged) + s(max_temp_lag1) + s(min_temp_lag1) + s(prcp_lag1)
form_hospitalized_qgam <- hospitalized ~ hospitalized_lag1 + hospitalized_lag7 + hospitalized_diff + 
  christmas + nye + ny + july14 +
  year + month + wday + holiday + holiday_lag1 + vacances_zone_c + s(flu_lagged) + s(diarrhea_lagged) + s(chickenpox_lagged) + s(max_temp_lag1) + s(min_temp_lag1) + s(prcp_lag1)

arrivals_learners = list(
  bjml = learner_bjml(form_arrivals, taus),
  qreg = learner_quantreg(form_arrivals, taus),
  grf = learner_grf(form_arrivals, taus),
  drf = learner_drf(form_arrivals, taus),
  gbm = learner_gbm(form_arrivals, taus),
  arima = learner_arima(form_arrivals, taus),
  arimax = learner_arimax(form_arrivals_arimax, taus),
  qgam = learner_qgam(form_arrivals_qgam, taus),
  one_ahead = learner_one_ahead(form_arrivals, taus, "date_numeric"),
  week_ahead = learner_week_ahead(form_arrivals, taus, "date_numeric")
)

hospitalized_learners = list(
  bjml = learner_bjml(form_hospitalized, taus),
  qreg = learner_quantreg(form_hospitalized, taus),
  grf = learner_grf(form_hospitalized, taus),
  drf = learner_drf(form_hospitalized, taus),
  gbm = learner_gbm(form_hospitalized, taus),
  arima = learner_arima(form_hospitalized, taus),
  arimax = learner_arimax(form_hospitalized_arimax, taus),
  qgam = learner_qgam(form_hospitalized_qgam, taus),
  one_ahead = learner_one_ahead(form_hospitalized, taus, "date_numeric"),
  week_ahead = learner_week_ahead(form_hospitalized, taus, "date_numeric")
)

print("Starting computation...")
globals <- c("form_arrivals", "form_hospitalized", "form_arrivals_qgam", "form_hospitalized_qgam", "form_arrivals_arimax", "form_hospitalized_arimax", "taus", "apply_learners", "form_hospitalized2", "form_arrivals2")
packages <- c("dplyr", "readr", "tibble", "purrr", "tidyr", "forecast", "grf", "quantreg", "bsts", "qgam", "gbm", "drf", "prophet")
options <- furrr_options(seed = TRUE, globals = globals, packages = packages, scheduling = 2)

arrivals_base_predictions <- X_grouped %>%
  ungroup() %>%
  mutate(pred_arrivals = furrr::future_map2(prevdata, data, apply_learners, prefix = "arrivals", learners = arrivals_learners, taus = taus, .options = options)) %>%
  select(pred_arrivals) %>% 
  unnest(pred_arrivals)

hospitalizations_base_predictions <- X_grouped %>%
  ungroup() %>%
  mutate(pred_arrivals = furrr::future_map2(prevdata, data, apply_learners, prefix = "hospitalized", learners = hospitalized_learners, taus = taus, .options = options)) %>%
  select(pred_arrivals) %>% 
  unnest(pred_arrivals)

#
# Step 2: ensemble methods
#
loss_lower <- loss_quantile(taus[1])
loss_upper <- loss_quantile(taus[2])
loss_int <- loss_interval(1 - (taus[2] - taus[1]))
loss_median <- loss_lt()
mae <- loss_mae()
mape <- loss_mape()

learners <- unique(arrivals_base_predictions$learner)

transform <- function(x) x
inv_transform <- function(x) ifelse(x < 0, 0, x)

arrival_ensembles <- arrivals_base_predictions %>%
  filter(!is.na(arrivals)) %>%
  online_aggregation("arrivals")

arrival_ensemble_preds <- arrival_ensembles %>% 
  select(-ewa, -boa, -sl) %>%
  unnest(c(data, ewa_predictions, boa_predictions, sl_predictions)) %>%
  pivot_longer(cols = all_of(c(learners, "ewa", "boa", "sl"))) %>%
  select(-ewa_weights, -boa_weights, -sl_weights)

hospitalization_ensembles <- hospitalizations_base_predictions %>%
  filter(!is.na(hospitalized)) %>%
  online_aggregation("hospitalized")

hospitalization_ensemble_preds <- hospitalization_ensembles %>%
  select(-ewa, -boa, -sl) %>%
  unnest(c(data, ewa_predictions, boa_predictions, sl_predictions)) %>%
  pivot_longer(cols = all_of(c(learners, "ewa", "boa", "sl"))) %>%
  select(-ewa_weights, -boa_weights, -sl_weights)

#
# Step 3: Conformalization
#

apply_conformal <- function(data, outcome, method, alpha = 0.8, ...) {
  x <- as.matrix(data[, c("0.1", "0.9")])
  y <- data[[outcome]]
  
  aci(y, x, alpha = alpha, method = method, parameters = list(gamma_grid = gammas, interval_constructor = "linear"))
}

conformalize <- function(data, outcome_column, gammas) {
  data %>%
    filter(quantile %in% c(0.1, 0.9), !(name %in% c("week_ahead", "one_ahead"))) %>%
    pivot_wider(names_from = "quantile", values_from = "value") %>%
    filter(!is.na(`0.1`), !is.na(`0.9`)) %>%
    group_by(id, my_id, name) %>%
    nest() %>%
    mutate(results = map(data, apply_conformal, outcome = outcome_column, method = "AgACI", gammas = gammas),
           preds = map2(data, results, function(data, results) {
             data %>%
               mutate(
                 lower_conformal = results$interval[,1], 
                 upper_conformal = results$interval[,2]
               )
           })
    ) 
}

gammas <- seq(0, 10, 0.5)
conformalized_arrivals <- arrival_ensemble_preds %>%
  conformalize("arrivals", gammas)
conformalized_hospitalized <- hospitalization_ensemble_preds %>%
  conformalize("hospitalized", gammas)

write_rds(conformalized_arrivals, "output/conformalized_arrivals.rds")
write_rds(conformalized_hospitalized, "output/conformalized_hospitalized.rds")

write_rds(arrival_ensemble_preds, "output/arrivals_ensemble_preds.rds")
write_rds(hospitalization_ensemble_preds, "output/hospitalizations_ensemble_preds.rds")


#
# Variable importance
#
var_imp <- function(x, y, pval = FALSE, B = 100) {
  btstrp <- rep(NA, B)
  ## type?
  vals <- unique(x)
  nvals <- length(vals)
  if (nvals == 1) {
    type <- "constant"
    out <- NA
  } else if (nvals <= 5) {
    type <- "categorical"
  } else {
    x <- as.numeric(x)
    type <- "numeric"
  }
  if (type == "numeric") {
    out <- abs(cor(x, y, method = "spearman", use = "complete.obs"))
    if (pval) {
      for (bb in 1:B) {
        btstrp[bb] <- abs(cor(sample(x), y, method = "spearman", use = "complete.obs"))
      }
      out <- mean(out <= btstrp)
    }
  } else if (type == "categorical") {
    ybar <- mean(y, na.rm = TRUE)
    denominator <- sum((y - ybar)^2, na.rm = TRUE)
    numerator <- tibble(x = x, y = y) %>%
      filter(!is.na(x), !is.na(y)) %>%
      group_by(x) %>%
      summarize(diff_of_means = mean(y) - ybar,
                count = n()) %>%
      summarize(numerator = sum(count * diff_of_means^2)) %>%
      pull(numerator)
    out <- sqrt(numerator / denominator)
    if (pval) {
      for (bb in 1:B) {
        numerator <- tibble(x = sample(x), y = y) %>%
          group_by(x) %>%
          summarize(diff_of_means = mean(y) - ybar,
                    count = n()) %>%
          summarize(numerator = sum(count * diff_of_means^2)) %>%
          pull(numerator)
        btstrp[bb] <- sqrt(numerator / denominator)
      }
      out <- mean(out <= btstrp)
    }
  }
  return(out)
}

get_type <- function(x) {
  ## type?
  vals <- unique(x)
  nvals <- length(vals)
  if (nvals == 1) {
    type <- NA
  } else if (nvals <= 5) {
    type <- nvals + 1
  } else {
    type <- 6
  }
  return(type)
}

pval <- function(x, type, PVALS, na.rm = TRUE) {
  pval <- mean(x <= PVALS[[type]], na.rm = na.rm)
  return(pval)
}

pval.vectorize <- Vectorize(pval, vectorize.args = c("x", "type"))

arrival_preds_with_x <- arrival_ensemble_preds %>%
  filter(quantile == 0.5) %>%
  left_join(X_grouped %>% ungroup() %>% select(data) %>% unnest(data), by = c("date", "id", "my_id"))

hospitalization_preds_with_x <- arrival_ensemble_preds %>%
  filter(quantile == 0.5) %>%
  left_join(X_grouped %>% ungroup() %>% select(data) %>% unnest(data), by = c("date", "id", "my_id"))

arrival_preds_with_x

arrivals_vars <- labels(terms(form_arrivals))
hospitalization_vars <- labels(terms(form_hospitalized))

arrivals_varimp <- expand_grid(
  variable = arrivals_vars,
  name = c("boa", "ewa", "sl")
) %>%
  mutate(importance = map2_dbl(variable, name, \(variable, name) {
    print(name)
    var_imp(arrival_preds_with_x[arrival_preds_with_x$name == name,][[variable]], 
            arrival_preds_with_x[arrival_preds_with_x$name == name,]$value)
  }),
  pval = map2_dbl(variable, name, \(variable, name) {
    print(name)
    var_imp(arrival_preds_with_x[arrival_preds_with_x$name == name,][[variable]], 
            arrival_preds_with_x[arrival_preds_with_x$name == name,]$value, pval = TRUE, B = 1e3)
  })
  ) %>%
  pivot_wider(names_from = "name", values_from = c("importance", "pval")) %>%
  arrange(-importance_boa)

hospitalization_varimp <- expand_grid(
  variable = hospitalization_vars,
  name = c("boa", "ewa", "sl")
) %>%
  mutate(importance = map2_dbl(variable, name, \(variable, name) {
    print(name)
    var_imp(hospitalization_preds_with_x[hospitalization_preds_with_x$name == name, ][[variable]], 
            hospitalization_preds_with_x[hospitalization_preds_with_x$name == name, ]$value)
  }),
  pval = map2_dbl(variable, name, \(variable, name) {
    print(name)
    var_imp(hospitalization_preds_with_x[hospitalization_preds_with_x$name == name, ][[variable]], 
            hospitalization_preds_with_x[hospitalization_preds_with_x$name == name, ]$value,
            pval = TRUE, B = 1e3)
  }),
  ) %>%
  pivot_wider(names_from = "name", values_from = c("importance", "pval")) %>%
  arrange(-importance_boa)

arrivals_varimp %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  select(variable, importance_boa, pval_boa, importance_ewa, pval_ewa, importance_sl, pval_sl) %>%
  knitr::kable(format = "latex")

hospitalization_varimp %>%
  select(variable, importance_boa, pval_boa, importance_ewa, pval_ewa, importance_sl, pval_sl) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  knitr::kable(format = "latex")
