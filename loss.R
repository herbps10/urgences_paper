loss_quantile <- function(tau) {
  function(pred, observed) {
    #pmax(tau * (pred - observed), (tau - 1) * (pred - observed))
    ifelse(pred <= observed, tau * abs(pred - observed), (1 - tau) * abs(pred - observed))
  }
}

loss_interval <- function(alpha) {
  function(lower, upper, observed) {
    na <- is.na(observed)
    lower <- lower[!na]
    upper <- upper[!na]
    observed <- observed[!na]
    (upper - lower) + 2 / alpha * (lower - observed) * (observed < lower) + 2 / alpha * (observed - upper) * (observed > upper)
  }
}

loss_mape <- function() {
  function(pred, observed) {
    na <- is.na(observed)
    Metrics::ape(ifelse(observed[!na] == 0, 1, observed[!na]), pred[!na])
  }
}

loss_smape <- function() {
  function(pred, observed) {
    na <- is.na(observed)
    abs(observed[!na] - pred[!na]) / ((observed[!na] + pred[!na]))
  }
}

loss_mae <- function() {
  function(pred, observed) {
    na <- is.na(observed)
    abs(pred[!na] - observed[!na])
  }
}

loss_lt <- function() {
  function(pred, observed) {
    na <- is.na(observed)
    pred[!na] < observed[!na]
  }
}
