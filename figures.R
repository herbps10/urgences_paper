library(tidyverse)
library(patchwork)

source("R/learners.R")
source("R/loss.R")

evaluate_conformal <- function(conformalized, outcome_column, groups = c("id", "name")) {
  conformalized %>%
    pivot_longer(c(`0.1`, `0.9`, lower_conformal, upper_conformal), names_to = "quantile") %>%
    mutate(name = ifelse(quantile %in% c("0.1", "0.9"), name, paste0(name, " conformal"))) %>%
    mutate(quantile = case_when(
      quantile %in% c("0.1", "0.9") ~ quantile,
      quantile == "lower_conformal" ~ "0.1",
      quantile == "upper_conformal" ~ "0.9"
    )) %>%
    pivot_wider(names_from = "quantile", values_from = "value") %>%
    ungroup() %>%
    mutate(
      loss_interval = loss_int(`0.1`, `0.9`, .[[outcome_column]]),
      coverage_below = .[[outcome_column]] < `0.1`,
      coverage_above = .[[outcome_column]] > `0.9`,
      coverage = `0.1` <= .[[outcome_column]] & `0.9` >= .[[outcome_column]],
      interval_width = `0.9` - `0.1`,
      wday = wday(date),
    ) %>%
    group_by(.[,groups]) %>%
    summarize_at(vars(c(loss_interval, coverage_below, coverage_above, coverage, interval_width)), mean)
}

point_performance <- function(results, outcome_column) {
  results %>%
    ungroup() %>%
    filter(quantile == 0.5) %>%
    mutate(mae = mae(inv_transform(value), .data[[outcome_column]]),
           mape = mape(inv_transform(value), .data[[outcome_column]]),
           loss_median = loss_median(inv_transform(value), .data[[outcome_column]])) %>%
    group_by(id, name) %>%
    filter(!is.na(mae), !is.na(mape), !is.na(loss_median)) %>%
    summarize_at(vars(mae, mape, loss_median), list(mean = mean, median = median)) %>%
    select(-loss_median_median)
}

point_performance_by_month <- function(results, outcome_column) {
  results %>%
    ungroup() %>%
    filter(quantile == 0.5) %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
    mutate(mae = mae(inv_transform(value), .data[[outcome_column]]),
           mape = mape(inv_transform(value), .data[[outcome_column]]),
           loss_median = loss_median(inv_transform(value), .data[[outcome_column]])) %>%
    group_by(name, month) %>%
    filter(!is.na(mae), !is.na(mape), !is.na(loss_median)) %>%
    summarize_at(vars(mae, mape, loss_median), list(mean = mean, median = median)) %>%
    select(-loss_median_median)
}

interval_performance <- function(results, outcome_column) {
  results %>%
    filter(quantile %in% c(0.1, 0.9), name != "one_ahead") %>%
    ungroup() %>%
    pivot_wider(names_from = "quantile", values_from = "value") %>%
    mutate(
      loss_interval = loss_int(inv_transform(`0.1`), inv_transform(`0.9`), .[[outcome_column]]),
      coverage_below = .[[outcome_column]] < inv_transform(`0.1`),
      coverage_above = .[[outcome_column]] > inv_transform(`0.9`),
      coverage = inv_transform(`0.1`) <= .[[outcome_column]] & inv_transform(`0.9`) >= .[[outcome_column]],
      interval_width = inv_transform(`0.9`) - inv_transform(`0.1`)
    ) %>%
    group_by(id, name) %>%
    summarize_at(vars(loss_interval, coverage_below, coverage, coverage_above, interval_width), mean)
}

mae <- loss_mae()
mape <- loss_mape()
loss_lower <- loss_quantile(taus[1])
loss_upper <- loss_quantile(taus[2])
loss_int <- loss_interval(1 - (taus[2] - taus[1]))
loss_median <- loss_lt()
taus <- c(0.1, 0.5, 0.8)

transform <- function(x) x
inv_transform <- function(x) ifelse(x < 0, 0, x)

arrival_ensemble_preds <- read_rds("output/arrivals_ensemble_preds.rds")
hospitalization_ensemble_preds <- read_rds("output/hospitalizations_ensemble_preds.rds")

arrival_conformal <- read_rds("output/conformalized_arrivals.rds") %>%
  select(my_id, name, preds) %>%
  unnest(preds)
hospitalization_conformal <- read_rds("output/conformalized_hospitalized.rds") %>%
  select(my_id, name, preds) %>%
  unnest(preds)

######## Table 1 #######

analysis_dataset <- read_rds("analysis_data/analysis_data.rds")

n_missing = function(x, ...) sum(is.na(x))

tableone <- analysis_dataset %>%
  group_by(id) %>%
  summarize_at(vars(arrivals, hospitalized), .funs = list(missing = n_missing, mean = mean, min = min, max = max), na.rm = TRUE) %>%
  pivot_longer(cols = c(starts_with("arrivals"), starts_with("hospitalized"))) %>%
  tidyr::separate(name, into = c("type", "stat")) %>%
  group_by(type, stat) %>%
  summarize_at(vars(value), .funs = list(mean = mean, min = min, max = max)) %>%
  mutate(type = ifelse(type == "hospitalized", "hospitalizations", type)) %>%
  mutate(stat = case_when(
    stat == "mean" ~ glue::glue("Average {type} in one day"),
    stat == "max" ~  glue::glue("Most {type} in one day"),
    stat == "min" ~  glue::glue("Fewest {type} in one day"),
    stat == "missing" ~  glue::glue("Days of missing {type}"),
  )) %>%
  mutate_if(is.numeric, scales::label_number(0.1)) %>%
  mutate(range = glue::glue("({min}-{max})")) %>%
  ungroup()
  #tidyr::separate(name, into = c("type", "x1", "x2"), sep = "_") %>%
  #pivot_wider(names_from = c("type", "x2"), values_from = value) %>%
  #select(x1, arrivals_min, arrivals_mean, arrivals_max, hospitalized_min, hospitalized_mean, hospitalized_max)

tableone %>%
  select(stat, mean, range) %>%
  knitr::kable(format = "latex")

######## Figure 1 #######

plot_ed <- function(ids) {
  p1 <- analysis_dataset %>%
    filter(id %in% ids) %>%
    select(id, date, arrivals) %>%
    ggplot(aes(x = date, y = arrivals)) +
    geom_point(size = 0.5) +
    labs(y = "Arrivals", x = "Date") +
    ggtitle("(A) Arrivals") +
    ggpubr::theme_pubclean()
  
  p2 <- analysis_dataset %>%
    filter(id %in% ids) %>%
    select(id, date, hospitalized) %>%
    ggplot(aes(x = date, y = hospitalized)) +
    geom_point(size = 0.5) +
    labs(y = "Hospitalizations", x = "Date") +
    ggtitle("(B) Hospitalizations") +
    ggpubr::theme_pubclean()
  
  p1 / p2
}
# Redacted plot_ed()
ggsave("plots/example_timeseries.pdf", width = 9, height = 6)


###### Table 3 ######
method_names <- tribble(
  ~name, ~title,
  "boa", "Bernstein Online Aggregation",
  "ewa", "Exponentially Weighted Average",
  "sl", "Super Learner",
  "bjml", "BJML",
  "qreg", "Quantile Regression",
  "qgam", "Quantile Generalized Additive Model",
  "bart", "Bayesian Additive Regression Trees",
  "gbm", "Gradient Boosted Machine",
  "grf", "Generalized Random Forest",
  #"prophet", "Prophet",
  "qreg_basic", "Quantile Regression\n(month, day)",
  "qreg_basic2", "Quantile Regression\n(all calendar + lagged outcome)",
  "qreg_basic3", "Quantile Regression\n(all calendar)",
  "drf", "Distributional Random Forest",
  "arima", "ARIMA",
  "arimax", "ARIMAX",
  "week_ahead", "Y_{t-7}",
  "one_ahead", "Y_{t-1}"
)

#
# Point performance tables
#
tab_arrivals <- point_performance(arrival_ensemble_preds, "arrivals") %>%
  group_by(name) %>%
  summarize(mae_median = median(mae_mean, na.rm = TRUE), mae_mean = mean(mae_mean, na.rm = TRUE), mape_mean = mean(mape_mean, na.rm = TRUE)) %>%
  arrange(-mae_mean) %>%
  left_join(method_names) %>%
  select(title, mae_mean, mape_mean) %>%
  mutate(mae_mean = scales::number(mae_mean, accuracy = 0.1),
         mape_mean = scales::percent(mape_mean, accuracy = 0.1))

tab_hospitalizations <- point_performance(hospitalization_ensemble_preds, "hospitalized") %>%
  group_by(name) %>%
  summarize(mae_mean = mean(mae_mean, na.rm = TRUE), mape_mean = mean(mape_mean, na.rm = TRUE)) %>%
  arrange(-mae_mean) %>%
  left_join(method_names) %>%
  select(title, mae_mean, mape_mean) %>%
  mutate(mae_mean = scales::number(mae_mean, accuracy = 0.1),
         mape_mean = scales::percent(mape_mean, accuracy = 0.1))

tab_arrivals %>% left_join(tab_hospitalizations, by = c("title")) %>%
  rename(
    Algorithm = title, 
    `MAE arrivals` = mae_mean.x,
    `MAPE arrivals` = mape_mean.x,
    `MAE hospitalizations` = mae_mean.y,
    `MAPE hospitalizations` = mape_mean.y
    ) %>%
  knitr::kable(format = "latex")

###### Table 4 ######
tab_arrivals <- point_performance(arrival_ensemble_preds, "arrivals") %>%
  filter(!(name %in% c("bjml", "one_ahead", "week_ahead", "ewa", "sl", "boa"))) %>%
  group_by(id) %>%
  mutate(rank_mae = rank(mae_mean),
         rank_mape = rank(mape_mean)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(number1_mae = mean(rank_mae == 1),
            number1_mape = mean(rank_mape == 1)) 

tab_hospitalized <- point_performance(hospitalization_ensemble_preds, "hospitalized") %>%
  filter(!(name %in% c("bjml", "one_ahead", "week_ahead", "ewa", "sl", "boa"))) %>%
  group_by(id) %>%
  mutate(rank_mae = rank(mae_mean),
         rank_mape = rank(mape_mean)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(number1_mae = mean(rank_mae == 1),
            number1_mape = mean(rank_mape == 1)) 

tab_arrivals %>%
  bind_rows(tab_hospitalized) %>%
  left_join(method_names) %>%
  mutate_if(is.numeric, function(x) scales::label_percent()(x)) %>%
  select(title, number1_mae, number1_mape) %>%
  knitr::kable(format = "latex")

##### Figure 4 #####

average_demand <- X_grouped %>% 
  ungroup() %>%
  filter(date == max(date)) %>% 
  mutate(all_data = map2(data, prevdata, bind_rows)) %>% 
  select(all_data) %>% 
  unnest(cols = c(all_data)) %>%
  group_by(id) %>%
  summarize(mean_arrivals = mean(arrivals, na.rm = TRUE),
            mean_hospitalized = mean(hospitalized, na.rm = TRUE))

p1 <- point_performance(arrival_ensemble_preds, "arrivals") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  ggplot(aes(x = mean_arrivals, y = mae_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "Mean daily arrivals", y = "Mean Absolute Error", title = "(A) Arrivals")

p2 <- point_performance(hospitalization_ensemble_preds, "hospitalized") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  ggplot(aes(x = mean_hospitalized, y = mae_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  labs(x = "Mean daily hospitalizations", y = "Mean Absolute Error", title = "(B) Hospitalizations") +
  theme(legend.position = "bottom")
p1 / p2

ggsave("plots/performance_comparison_by_ed.pdf", width = 8, height = 6)

##### Figure 8 #####
p1 <- point_performance(arrival_ensemble_preds, "arrivals") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  ggplot(aes(x = mean_arrivals, y = mape_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "Mean daily arrivals", y = "Mean Absolute\nPercentage Error", title = "(A) Arrivals")

p2 <- point_performance(hospitalization_ensemble_preds, "hospitalized") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  ggplot(aes(x = mean_hospitalized, y = mape_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  labs(x = "Mean daily hospitalizations", y = "Mean Absolute\nPercentage Error", title = "(B) Hospitalizations") +
  theme(legend.position = "bottom")
p1 / p2 & plot_annotation(title = "MAPE by Emergency Department")
ggsave("plots/performance_comparison_by_ed_mape.pdf", width = 8, height = 6)

###### Figure 9 #####
p1 <- point_performance(arrival_ensemble_preds, "arrivals") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  filter(mean_arrivals > 100) %>%
  ggplot(aes(x = mean_arrivals, y = mape_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "Mean daily arrivals", y = "Mean Absolute\nPercentage Error", title = "(A) Arrivals")

p2 <- point_performance(hospitalization_ensemble_preds, "hospitalized") %>%
  filter(name %in% c("boa", "qreg", "gbm")) %>%
  left_join(method_names) %>%
  left_join(average_demand) %>%
  rename(Algorithm = title) %>%
  filter(mean_hospitalized > 10) %>%
  ggplot(aes(x = mean_hospitalized, y = mape_mean, color = Algorithm)) +
  geom_line(aes(group = id), color = "black", alpha = 0.5) +
  geom_point() +
  ggpubr::theme_pubclean() +
  labs(x = "Mean daily hospitalizations", y = "Mean Absolute\nPercentage Error", title = "(B) Hospitalizations") +
  theme(legend.position = "bottom")
p1 / p2 & plot_annotation(title = "MAPE by Emergency Department (Zoom)")
ggsave("plots/performance_comparison_by_ed_mape_zoom.pdf", width = 8, height = 6)

bold_best_percent <- function(p, target = 80) {
  x <- as.numeric(str_replace(p, "%", ""))
  lowest <- abs((x - target)) == min(abs(x - target))
  p[lowest] <- str_c("\\textbf{", p[lowest], "}")
  p
}

bold_lowest <- function(p) {
  x <- as.numeric(p)
  lowest <- x == min(x, na.rm = TRUE)
  p[lowest] <- str_c("\\textbf{", p[lowest], "}")
  p
}

##### Table 5 #####
conformal_table <- function(data) {
  data %>%
    mutate(conformal = ifelse(str_detect(name, "conformal"), "conformal", "original"),
         name = str_replace(name, " conformal", "")) %>%
    group_by(name, conformal) %>%
    summarize(coverage = mean(coverage), interval_width = mean(interval_width), loss_interval = mean(loss_interval)) %>%
    pivot_wider(names_from = conformal, values_from = c(coverage, interval_width, loss_interval)) %>%
    left_join(method_names) %>%
    ungroup() %>%
    select(title, starts_with("coverage"), starts_with("interval_width"), starts_with("loss_interval")) %>%
    arrange(-interval_width_conformal) %>%
    mutate_at(vars(starts_with("coverage")), scales::percent_format(accuracy = 0.1)) %>%
    mutate_at(vars(starts_with("interval_width")), scales::number_format(accuracy = 0.1)) %>%
    mutate_at(vars(starts_with("loss_interval")), scales::number_format(accuracy = 0.1)) %>%
    mutate_at(vars(starts_with("coverage")), bold_best_percent) %>%
    mutate_at(vars(c(starts_with("interval_width"), starts_with("loss_interval"))), bold_lowest) %>%
    mutate_if(is.character, ~str_replace(., "\\%", "\\\\%"))
}

tab_arrivals <- evaluate_conformal(arrival_conformal, "arrivals") %>%
  conformal_table()

tab_hospitalizations <- evaluate_conformal(hospitalization_conformal, "hospitalized") %>%
  conformal_table()
  
tab_arrivals %>%
  select(title, coverage_original, coverage_conformal, interval_width_original, interval_width_conformal) %>%
  knitr::kable(format = "latex", escape = FALSE)

tab_hospitalizations %>%
  select(title, coverage_original, coverage_conformal, interval_width_original, interval_width_conformal) %>%
  knitr::kable(format = "latex", escape = FALSE)

#
# Interval performance figures
#
comparison_plot <- function(data) {
  data %>% 
    ungroup() %>%
    mutate(conformal = ifelse(str_detect(name, "conformal"), "Conformalized", "Original"),
           name = str_replace(name, " conformal", "")) %>%
    filter(name %in% c("ewa", "boa", "sl")) %>%
    #filter(name %in% c("qreg_basic", "qreg_basic2", "qreg_basic3")) %>%
    left_join(method_names) %>%
    ggplot(aes(x = coverage, y = interval_width)) +
    geom_line(aes(group = id), alpha = 0.1) +
    geom_point(aes(color = conformal), alpha = 0.7) +
    geom_vline(xintercept = 0.8, lty = 2) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(x = "Empirical coverage", y = "Interval width") +
    ggpubr::theme_pubclean() +
    theme(legend.title = element_blank()) +
    facet_wrap(~title) +
    theme(legend.position = "right")
}

p1 <- evaluate_conformal(arrival_conformal, "arrivals") %>%
  comparison_plot() +
  ggtitle(label = "(A) Arrivals")

p2 <- evaluate_conformal(hospitalization_conformal, "hospitalized") %>%
  comparison_plot() +
  ggtitle(label = "(B) Hospitalizations")

p1 / p2

ggsave("plots/conformalized_comparisons.pdf", width = 8, height = 6)


#
# Examples
#
target_id <- 18

arrivals_conformal_preds <- arrival_conformal %>%
  filter(my_id == target_id, name == "boa")

hospitalization_conformal_preds <- hospitalization_conformal %>%
  filter(my_id == target_id, name == "boa")

arrivals_preds <- arrival_ensemble_preds %>%
  filter(my_id == target_id, name == "boa", quantile == 0.5)

hospitalized_preds <- hospitalization_ensemble_preds %>%
  filter(my_id == target_id, name == "boa", quantile == 0.5)


p1 <- arrivals_conformal_preds %>%
  ggplot(aes(x = date, y = arrivals)) +
  geom_point(size = 0.5, color = "#619CFF") +
  geom_line(aes(y = lower_conformal), alpha = 0.5) +
  geom_line(data = arrivals_preds, aes(y = value)) +
  geom_line(aes(y = upper_conformal), alpha = 0.5) +
  labs(x = "Date", y = "Arrivals") +
  ggpubr::theme_pubclean() +
  ggtitle("(A) Arrivals")

p2 <- hospitalization_conformal_preds %>%
  ggplot(aes(x = date, y = hospitalized)) +
  geom_point(size = 0.5, color = "#619CFF") +
  geom_line(aes(y = lower_conformal), alpha = 0.5) +
  geom_line(data = hospitalized_preds, aes(y = value)) +
  geom_line(aes(y = upper_conformal), alpha = 0.5) +
  labs(x = "Date", y = "Hospitalizations") +
  ggpubr::theme_pubclean() +
  ggtitle("(B) Hospitalizations")

#(p1 + p1_zoom) / (p2 + p2_zoom)
p1 / p2

ggsave("plots/example_forecasts.pdf", width = 8, height = 6)

#
# Example
#
target_id <- 7
context <- arrival_conformal %>%
  filter(my_id == target_id, name == "boa") %>%
  filter(month(date) %in% 4:5)

conformal_preds <- arrival_conformal %>%
  filter(my_id == target_id, name == "boa") %>%
  filter(month(date) == 6, day(date) == 1)

preds <- arrival_ensemble_preds %>%
  filter(name == "boa", my_id == target_id, quantile == 0.5) %>%
  filter(month(date) == 6, day(date) == 1)

context %>%
  ggplot(aes(x = date, y = arrivals)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(data = conformal_preds, aes(ymin = lower_conformal, ymax = upper_conformal), width = 0, size = 1) +
  geom_point(data = preds, aes(y = value), size = 3) +
  labs(x = "Date", y = "Arrivals")
 
ggsave("plots/example_goal.pdf", width = 8, height = 4)

conformal_preds <- arrival_conformal %>%
  filter(my_id == target_id, name == "boa")

preds <- arrival_ensemble_preds %>%
  filter(name == "boa", my_id == target_id, quantile == 0.5) 

ggplot(conformal_preds, aes(x = date, y = arrivals)) +
  geom_point(alpha = 0.5) +
  geom_line(data = preds, aes(y = value)) +
  geom_line(aes(y = lower_conformal), alpha = 0.5) +
  geom_line(aes(y = upper_conformal), alpha = 0.5)


#
# Pipeline example
#

target_id <- 7

base_preds <- arrival_ensemble_preds  %>%
  left_join(method_names) %>%
  mutate(title = ifelse(title == "Quantile Generalized Additive Model", "Quantile Generalized\nAdditive Model", title)) %>%
  filter(my_id == target_id, name %in% c("qreg", "grf", "qgam")) %>%
  filter(month(date) %in% 5:9) %>%
  pivot_wider(names_from = quantile, values_from = value)

base_pred_metrics <- base_preds %>%
  group_by(title) %>%
  summarize(coverage = mean(`0.1` < arrivals & `0.9` > arrivals),
            mae = mean(abs(`0.5` - arrivals))) %>%
  mutate(coverage = glue::glue("EmpCov: {scales::label_percent(0.1)(coverage)}"),
         mae = glue::glue("MAE: {scales::label_number(0.1)(mae)}"))

base_preds %>%
  ggplot(aes(x = date, y = arrivals)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = `0.5`, color = name)) +
  geom_line(aes(y = `0.1`, color = name), alpha = 0.5) +
  geom_line(aes(y = `0.9`, color = name), alpha = 0.5) +
  geom_text(data = base_pred_metrics, aes(y = 223, x = ymd("2018-07-25"), label = coverage), hjust = 0, size = 3) +
  geom_text(data = base_pred_metrics, aes(y = 207, x = ymd("2018-07-25"), label = mae), hjust = 0, size = 3) +
  facet_wrap(~title) +
  labs(x = "", y = "Arrivals") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none")

ggsave("plots/pipeline_example/step1.pdf", width = 7, height = 2.5)

ensemble_preds <- arrival_ensemble_preds %>%
  filter(my_id == target_id, name %in% c("boa")) %>%
  filter(month(date) %in% 5:9)  %>%
  pivot_wider(names_from = quantile, values_from = value) %>%
  mutate(title2 = "Bernstein Online Aggregation")

ensemble_pred_metrics <- ensemble_preds %>%
  group_by(title2) %>%
  summarize(coverage = mean(`0.1` < arrivals & `0.9` > arrivals),
            mae = mean(abs(`0.5` - arrivals))) %>%
  mutate(coverage = glue::glue("EmpCov: {scales::label_percent(0.1)(coverage)}"),
         mae = glue::glue("MAE: {scales::label_number(0.1)(mae)}"))
  
base_preds  %>%
  ggplot(aes(x = date, y = arrivals)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = `0.5`, color = name), alpha = 0.5) +
  geom_line(aes(y = `0.1`, color = name), alpha = 0.2) +
  geom_line(aes(y = `0.9`, color = name), alpha = 0.2) +
  geom_line(data = ensemble_preds, aes(y = `0.5`), alpha = 1) +
  geom_line(data = ensemble_preds, aes(y = `0.1`), alpha = 0.5) +
  geom_line(data = ensemble_preds, aes(y = `0.9`), alpha = 0.5) +
  geom_text(data = ensemble_pred_metrics, aes(y = 223, x = ymd("2018-09-10"), label = coverage), hjust = 0, size = 3) +
  geom_text(data = ensemble_pred_metrics, aes(y = 207, x = ymd("2018-09-10"), label = mae), hjust = 0, size = 3) +
  facet_wrap(~title2) +
  labs(x = "Date", y = "Arrivals") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none") +
  labs(x = "")

ggsave("plots/pipeline_example/step2.pdf", width = 7, height = 2.5)

conformal_preds <- arrival_conformal %>%
  filter(my_id == target_id, name == "boa") %>%
  filter(month(date) %in% 5:9)

conformal_pred_metrics <- conformal_preds %>%
  mutate(title = "Conformalized Bernstein Online Aggregation") %>%
  group_by(title) %>%
  summarize(coverage = mean(lower_conformal < arrivals & upper_conformal > arrivals)) %>%
  mutate(coverage = glue::glue("EmpCov: {scales::label_percent(0.1)(coverage)}"))

arrival_ensemble_preds %>%
  left_join(method_names) %>%
  filter(my_id == target_id, name %in% c("boa")) %>%
  mutate(title = "Conformalized Bernstein Online Aggregation") %>%
  filter(month(date) %in% 5:9) %>%
  pivot_wider(names_from = quantile, values_from = value) %>%
  ggplot(aes(x = date, y = arrivals)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = `0.5`)) +
  geom_line(aes(y = `0.1`), alpha = 0.4) +
  geom_line(aes(y = `0.9`), alpha = 0.4) +
  geom_line(data = conformal_preds, aes(y = lower_conformal), color = "red", alpha = 0.5) +
  geom_line(data = conformal_preds, aes(y = upper_conformal), color = "red", alpha = 0.5) +
  geom_text(data = conformal_pred_metrics, aes(y = 213, x = ymd("2018-09-10"), label = coverage), hjust = 0, size = 3) +
  geom_text(data = ensemble_pred_metrics, aes(y = 197, x = ymd("2018-09-10"), label = mae), hjust = 0, size = 3) +
  facet_wrap(~title) +
  labs(x = "", y = "Arrivals") +
  theme(legend.position = "none") +
  ggpubr::theme_pubclean() +
  labs(x = "")

ggsave("plots/pipeline_example/step3.pdf", width = 7, height = 2.5)


hospitalizations_by_day_codes <- read_rds("analysis_data/hospitalized_by_day_codes.rds")

code_descriptions <- tribble(
  ~code, ~description,
  "J21", "J21: acute bronchiolitis",
  "J45", "J45: asthma",
  "S72", "S72: fracture of femur"
)

day_codes <- hospitalizations_by_day_codes %>%
  group_by(code, DAT_ENT) %>% 
  summarize(hospitalized = sum(hospitalized))  %>%
  group_by(DAT_ENT) %>%
  mutate(total = sum(hospitalized)) %>%
  mutate(prop = hospitalized / total)

codes <- c("J21", "J45", "S72")

day_codes %>%
  filter(code %in% codes) %>%
  ungroup() %>%
  left_join(code_descriptions) %>%
  ggplot(aes(x = DAT_ENT, y = prop)) +
  geom_point(size = 0.1) +
  scale_y_continuous(labels = scales::label_percent()) +
  facet_wrap(~description, ncol = 3) +
  ggpubr::theme_pubclean() +
  labs(x = "Date", y = "Hospitalizations\n(% of daily total)")

ggsave("plots/hospitalizations_examples.pdf", width = 9, height = 3)


p1 <- point_performance_by_month(arrival_ensemble_preds %>% filter(name %in% c("ewa", "boa", "sl")), "arrivals") %>% 
  mutate(name = case_when(name == "ewa" ~ "Exponentially Weighted Average", name == "boa" ~ "Bernstein Online Aggregation", name == "sl" ~ "Super Learner")) %>%
  ggplot(aes(x = month, y = mae_mean)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~name) + 
  labs(x = "Month", y = "MAE", title = "Arrivals") + 
  scale_x_continuous(breaks = 1:12)

p2 <- point_performance_by_month(hospitalization_ensemble_preds %>% filter(name %in% c("ewa", "boa", "sl")), "hospitalized") %>% 
  mutate(name = case_when(name == "ewa" ~ "Exponentially Weighted Average", name == "boa" ~ "Bernstein Online Aggregation", name == "sl" ~ "Super Learner")) %>%
  ggplot(aes(x = month, y = mae_mean)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~name) + 
  labs(x = "Month", y = "MAE", title = "Hospitalizations") + 
  scale_x_continuous(breaks = 1:12)

p1 / p2
ggsave(filename = "plots/training_quality_over_time.pdf")
