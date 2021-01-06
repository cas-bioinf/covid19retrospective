inv_logit <- function(x) {  1/(1+exp(-x)) }

calibration_intercept <- function(lp, y, family = binomial, ...) {
  #lp <- log(p / (1 - p))
  resp = y == "dead"
  fit <- glm(resp ~ 1, family = family, offset = lp, ...)
  variance <- diag(vcov(fit))
  c(estimate = coef(fit), se = sqrt(variance), var = variance)
}

calibration_slope <- function(lp, y, family = binomial, ...) {
  #lp <- log(p / (1 - p))
  resp = y == "dead"
  fit <- glm(resp ~ lp + 1, family = family, ...)
  variance <- diag(vcov(fit))[2]
  c(estimate = coef(fit)[2], se = sqrt(variance), var = variance)
}

auc <- function(lp, y, method = "delong", direction = "<", ...) {
  fit <- pROC::auc(response = y, predictor = lp, direction = direction, method = method,  ...)
  variance <- var(fit)
  ci <- confint(fit)
  c(estimate = fit[[1]][[1]], se = sqrt(variance), var = variance, ci.lb = ci$ci.lb, ci.ub = ci$ci.ub)
}

net_benefit <- function(lp, y, threshold = .50, test_harm = 0, case_label = "dead", ...) {
  p = inv_logit(lp)
  tp <- sum((p >= threshold) & (y == case_label))
  fp <- sum((p >= threshold) & (y != case_label))
  n <- length(p)

  if (threshold == 1) { # by definition
    return(0)
  } else{return(tp/n - fp/n * threshold / (1-threshold) - test_harm)}
}

net_benefit_bootstrap <- function(lp, y, threshold = .50, test_harm = 0, case_label = "dead", I = 100, ...) {
  p = inv_logit(lp)
  nb <- rep(NA, I)
  for (i in seq_len(I)) {
    s <- sample.int(length(p), length(p), replace = TRUE)
    nb[i] <- net_benefit(p[s], y[s], threshold, test_harm, case_label, ...)
  }
  se <- sd(nb)
  est <- net_benefit(p, y, threshold, test_harm, case_label, ...)

  c(estimate = est, se = se, var = se^2)
}

oe_ratio <- function(lp, y, case_label = "dead", log = TRUE, ...) {
  p = inv_logit(lp)
  O <- sum(y == case_label)
  E <- sum(p)
  N <- length(p)
  g <- if(log) "log(OE)" else NULL

  fit <- metamisc::oecalc(O = O, E = E, N = N, g = g, ...)

  c(estimate = fit$theta,
    se = fit$theta.se,
    var = fit$theta.se^2,
    ci.lb = fit$theta.cilb,
    ci.ub = fit$theta.ciub,
    O = O,
    E = E,
    N = N,
    log = log)
}

validate.mids <- function(data, model, measure, alpha = .05, ...) {
  # From mice
  call <- match.call()
  if (!is.mids(data)) {
    stop("The data must have class mids")
  }

  # Make sure that these two are functions:
  model   <- match.fun(model)
  measure <- match.fun(measure)

  # Here we have two steps instead of one.
  # 1. Predict the outcome
  predicted_values <- lapply(seq_len(data$m),
                             function(i) model(data = complete(data, i)))

  # 2. Estimate performance
  analyses <- sapply(seq_len(data$m),
                     function(i) measure(
                       lp = predicted_values[[i]],
                       y = complete(data, i)$mortality, ...))
  analyses <- as.data.frame(t(analyses))
  analyses$ci.lb <- analyses$estimate + qnorm(alpha/2)   * analyses$se
  analyses$ci.ub <- analyses$estimate + qnorm(1-alpha/2) * analyses$se

  # From mice
  object <- list(call = call,
                 call1 = data$call,
                 nmis = data$nmis,
                 analyses = analyses)
  oldClass(object) <- c("mira")
  object
}

pool.validate.mids <- function(object, alpha = .05, ...) {
  fit <- pool.scalar(Q = object$analyses$estimate,
                     U = object$analyses$var,
                     n = nrow(object$analyses))

  list("invididual" = data.frame(estimate = fit$qhat,
                                 var = fit$u,
                                 se = sqrt(fit$u),
                                 ci.lb = object$analyses$ci.lb,
                                 ci.ub = object$analyses$ci.ub),
       "pooled" = data.frame(estimate = fit$qbar,
                             variance = fit$t,
                             within = fit$ubar,
                             between = fit$b,
                             se = sqrt(fit$t),
                             ci.lb = fit$qbar + qnorm(alpha/2) * sqrt(fit$t),
                             ci.ub = fit$qbar + qnorm(1 - alpha/2) * sqrt(fit$t),
                             df = fit$df,
                             r = fit$r,
                             fmi = fit$fmi,
                             m = fit$m))
}

pool.auc.mids <- function(object, alpha = .05, ...) {
  m <- length(object$analyses$se)

  logit <- function(x) log(x / (1-x))
  inv_logit <- function(x) {1/(1+exp(-x))}

  # Estimates from individual studies
  ind <- list(auc = data.frame(est = object$analyses$estimate,
                               se = object$analyses$se,
                               ci.lb = object$analyses$ci.lb,
                               ci.ub = object$analyses$ci.ub))
  ind$logit.auc = data.frame(est = logit(ind$auc$est),
                             se = ind$auc$se / (ind$auc$est * (1 - ind$auc$est)),
                             ci.lb = logit(ind$auc$ci.lb),
                             ci.ub = logit(ind$auc$ci.ub))

  # Pooled across imputed data sets
  logit.auc <- list()
  logit.auc$est <- mean(ind$logit.auc$est)
  logit.auc$within <- mean(ind$logit.auc$se^2)
  logit.auc$between <- (1 + (1/m)) * var(ind$logit.auc$est)
  logit.auc$var <- logit.auc$within + logit.auc$between
  logit.auc$se <- sqrt(logit.auc$var)
  logit.auc$ci.lb <- logit.auc$est + qnorm(alpha/2)     * logit.auc$se
  logit.auc$ci.ub <- logit.auc$est + qnorm(1 - alpha/2) * logit.auc$se
  logit.auc$m <- m

  auc <- list()
  auc$est <- inv_logit(logit.auc$est)
  auc$ci.lb <- inv_logit(logit.auc$ci.lb)
  auc$ci.ub <- inv_logit(logit.auc$ci.ub)
  auc$m <- m

  return(list(individual = ind,
              pooled = list(logit.auc = as.data.frame(logit.auc),
                            auc = as.data.frame(auc))))
}

pool.oe.mids <- function(object, alpha = .05, ...) {
  fit <- pool.validate.mids(object, alpha = alpha, ...)

  if (!object$analyses$log[[1]]) return(list(oe = fit))

  oe <- list(individual = exp(fit$invididual))
  oe$individual$var <- NULL
  oe$individual$se <- NULL

  oe$pooled <- data.frame(estimate = exp(fit$pooled$estimate),
                          ci.lb = exp(fit$pooled$ci.lb),
                          ci.ub = exp(fit$pooled$ci.ub))

  list(log_oe = fit,
       oe = oe)
}
