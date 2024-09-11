#source("R/drf_boosted.R")

learner_one_ahead <- function(formula, taus, date, ...) {
  outcome <- all.vars(formula)[1]
  function(data, newdata, return_fit = FALSE) {
    matrix(newdata[[paste0(outcome, "_lag1")]], byrow = FALSE, ncol = length(taus), nrow = nrow(newdata))
  }
}

learner_week_ahead <- function(formula, taus, date, ...) {
  outcome <- all.vars(formula)[1]
  function(data, newdata, return_fit = FALSE) {
    matrix(newdata[[paste0(outcome, "_lag7")]], byrow = FALSE, ncol = length(taus), nrow = nrow(newdata))
  }
}

learner_bjml <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  function(data, newdata, return_fit = FALSE) {
    preds <- matrix(ncol = length(taus), nrow = nrow(newdata))

    for(index in 1:nrow(newdata)) {
      # Get week number
      wnumber <- lubridate::week(newdata$date2[index])
      
      # Filter for same week number
      if(wnumber == 53) {
        wdata <- data[lubridate::week(data$date2) >= 52,]
      }
      else {
        wdata <- data[lubridate::week(data$date2) == wnumber, ]
      }
      
      preds[index,] <- quantile(wdata[[outcome]], taus, na.rm = TRUE)
    }
    
    return(preds)
  }
}

learner_quantreg <- function(formula, taus, scale = FALSE, ...) {
  covars <- all.vars(formula)[-1]
  
  function(data, newdata, return_fit = FALSE) {
    start <- Sys.time()
    
    if(length(scale) > 1 || !(scale == FALSE)) {
      for(covar in scale) {
        newdata[[covar]] <- (newdata[[covar]] - mean(data[[covar]], na.rm = TRUE)) / sd(data[[covar]], na.rm = TRUE)
        data[[covar]] <- (data[[covar]] - mean(data[[covar]], na.rm = TRUE)) / sd(data[[covar]], na.rm = TRUE)
      }
    }
    
    fits <- list()
    for(i in 1:length(taus)) {
      fits[[i]] <- quantreg::rq(formula, tau = taus[i], data = data, ...)
    }
    time <- Sys.time() - start
    
    preds <- matrix(ncol = length(taus), nrow = nrow(newdata))
    for(i in 1:length(taus)) {
      preds[, i] <- predict(fits[[i]], newdata = newdata)
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fits))
    }
    else {
      return(preds)
    }
  }
}

learner_prophet <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  
  function(data, newdata, return_fit = FALSE) {
    data$ds <- as.Date(data$date2)
    data$y <- data[[outcome]]
    
    start <- Sys.time()
    fit <- prophet::prophet(data)
    time <- Sys.time() - start
    
    future <- prophet::make_future_dataframe(fit, periods = nrow(newdata), include_history = FALSE)
    pred <- predict(fit, future)
    
    preds <- matrix(c(
      pred$yhat_lower,
      pred$yhat,
      pred$yhat_upper
    ), ncol = 3, nrow = nrow(newdata), byrow = FALSE)
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}

learner_qgam <- function(formula, taus, ...) {
  function(data, newdata, return_fit = FALSE) {
    start <- Sys.time()
    fit <- qgam::mqgam(formula, data = data, qu = taus, ...)
    time <- Sys.time() - start
    
    preds <- matrix(ncol = length(taus), nrow = nrow(newdata))
    for(i in 1:length(taus)) {
      preds[, i] <- predict(fit$fit[[i]], newdata = newdata)
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}

learner_grf <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  function(data, newdata, return_fit = FALSE) {
    nas <- is.na(data[[outcome]])
    start <- Sys.time()
    fit <- grf::quantile_forest(
      X = data[!nas, covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric),
      Y = data[[outcome]][!nas],
      quantiles = taus,
      ...
    )
    time <- Sys.time() - start
    
    preds <- predict(fit, newdata = newdata[, covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric))$predictions
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}


learner_bart <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  function(data, newdata, return_fit = FALSE) {
    nas <- is.na(data[[outcome]])
    start <- Sys.time()
    fit <- BART::wbart(
      x.train = as.matrix(data[!nas, covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric)),
      y.train = data[[outcome]][!nas],
      x.test = as.matrix(newdata[, covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric)),
      ...
    )
    time <- Sys.time() - start
    
    browser()
    preds <- t(apply(fit$yhat.test, 2, quantile, taus))
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}


learner_drf <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  function(data, newdata, return_fit = FALSE) {
    data <- tidyr::drop_na(data[,c(outcome, covars)])
    start <- Sys.time()
    fit <- drf::drf(
      X = data[,covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric),
      Y = data[[outcome]],
      ...
    )
    time <- Sys.time() - start
    
    preds <- matrix(unlist(
      predict(fit, newdata = newdata[, covars] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric),
              functional = "quantile",
              quantiles = c(0.1, 0.5, 0.9)
              )
    ), ncol = 3, nrow = nrow(newdata), byrow = FALSE)
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}

learner_drf_boosted <- function(formula, taus, iterations = 5, ...) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  function(data, newdata, return_fit = FALSE) {
    data <- tidyr::drop_na(data[,c(outcome, covars)])
    start <- Sys.time()
    fit <- drf_boosted(
      X = data[,covars] %>% mutate_if(is.character, as.factor),
      Y = data[[outcome]],
      iterations = iterations,
      ...
    )
    time <- Sys.time() - start
    
    preds <- matrix(unlist(
      predict_drf_boosted(fit, newdata = newdata[, covars] %>% mutate_if(is.character, as.factor))
    ), ncol = 3, nrow = nrow(newdata), byrow = FALSE)
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit, time = time))
    }
    else {
      return(preds)
    }
  }
}

learner_gbm <- function(formula, taus, ...) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  
  function(data, newdata, return_fit = FALSE) {
    nas <- is.na(data[[outcome]])
    data <- data[!nas,] %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric)
    start <- Sys.time()
    
    fits <- list()
    for(i in 1:length(taus)) {
      fits[[i]] <- gbm::gbm(formula, data = data, n.trees = 2e3, distribution = list(name = "quantile", alpha = taus[i]), ...)
    }
    time <- Sys.time() - start
    
    newdata <- newdata %>% mutate_if(is.character, as.factor) %>% mutate_all(as.numeric)
    
    preds <- matrix(ncol = length(taus), nrow = nrow(newdata))
    for(i in 1:length(taus)) {
      preds[, i] <- predict(fits[[i]], newdata)
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fits), time = time)
    }
    else {
      return(preds)
    }
  }
}

learner_arimax <- function(formula, taus) {
  outcome <- all.vars(formula)[1]
  covars <- all.vars(formula)[-1]
  
  function(data, newdata, return_fit = FALSE) {
    y <- ts(c(data[[outcome]], newdata[[outcome]]), frequency = 365)
    x <- as.matrix(bind_rows(data[,covars], newdata[,covars]))
    fit <- forecast::auto.arima(window(y, time(y)[1], time(y)[nrow(data)]), seasonal = FALSE, xreg = x[1:nrow(data),], allowdrift = FALSE)
    
    preds <- matrix(ncol = 3, nrow = nrow(newdata))
    for(index in 1:nrow(newdata)) {
      pred <- forecast::forecast(fit, xreg = x[(nrow(data) + index),,drop=FALSE], level = 80)
      preds[index, ] <- c(pred$lower, pred$mean, pred$upper)
      if(nrow(newdata) > 1) {
      	fit <- forecast::Arima(y = window(y, time(y)[1], time(y)[nrow(data) + index]), xreg =  x[1:(nrow(data) + index),,drop=FALSE], model = fit)
      }
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit))
    }
    else {
      return(preds)
    }
  }
}

learner_arima <- function(formula, taus) {
  outcome <- all.vars(formula)[1]
  
  function(data, newdata, return_fit = FALSE) {
    y <- ts(c(data[[outcome]], newdata[[outcome]]), frequency = 365)
    fit <- forecast::auto.arima(window(y, time(y)[1], time(y)[nrow(data)]), seasonal = FALSE, allowdrift = FALSE)
    
    preds <- matrix(ncol = 3, nrow = nrow(newdata))
    for(index in 1:nrow(newdata)) {
      pred <- forecast::forecast(fit, level = 80, h = 1)
      preds[index, ] <- c(pred$lower, pred$mean, pred$upper)
      if(nrow(newdata) > 1) {
      	fit <- forecast::Arima(window(y, time(y)[1], time(y)[nrow(data) + index]), model = fit)
      }
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit))
    }
    else {
      return(preds)
    }
  }
}

learner_tbats <- function(formula, taus) {
  outcome <- all.vars(formula)[1]
  
  function(data, newdata, return_fit = FALSE) {
    y <- ts(c(data[[outcome]], newdata[[outcome]]), frequency = 365)
    fit <- forecast::tbats(window(y, time(y)[1], time(y)[nrow(data)]))
    
    preds <- matrix(ncol = 3, nrow = nrow(newdata))
    for(index in 1:nrow(newdata)) {
      pred <- forecast::forecast(fit, level = 80, h = 1)
      preds[index, ] <- c(pred$lower, pred$mean, pred$upper)
      if(nrow(newdata) > 1) {
      	fit <- forecast::tbats(window(y, time(y)[1], time(y)[nrow(data) + index]), model = fit)
      }
    }
    
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit))
    }
    else {
      return(preds)
    }
  }
}

learner_bsts <- function(formula, taus, date, retry = 5) {
  outcome <- all.vars(formula)[1]
  
  function(data, newdata, return_fit = FALSE) {
    data <- data[!is.na(data[[outcome]]),]
    
    y <- zoo::zoo(data[[outcome]], data[[date]])
    setup <- AddSemilocalLinearTrend(list(), y)
    setup <- AddAr(setup, y, 3)
    setup <- AddMonthlyAnnualCycle(setup, y, data[[date]][1])
    
    attempt <- 0
    fit <- NULL
    preds <- matrix(NA, ncol = 3, nrow = nrow(newdata))
    while(attempt < retry && is.null(fit)) {
      try({
        print(attempt)
        fit <- bsts(data[[outcome]], state.specification=setup, niter=1000)
	pred <- predict(fit, newdata = newdata, quantiles = taus)
        preds <- matrix(c(
          pred$interval 
        ), ncol = 3, nrow = nrow(newdata), byrow = FALSE)
      })
      attempt <- attempt + 1
    }
        
    if(return_fit == TRUE) {
      return(list(preds = preds, fit = fit))
    }
    else {
      return(preds)
    }
  }
}
