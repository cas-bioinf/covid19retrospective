brmhmm <- function(brmshmmdata) {
  d <- validate_brmshmmdata(brmshmmdata)

  prepdata <- make_data_hmm(d)

  bfit <- brms::brm(
    formula = d$formula,
    family = rate_hmm_family,
    data = prepdata$brmsdata,
    prior = d$prior,
    stanvars = rate_hmm_stanvars(prepdata$standata)
    )

  structure(list(
    brmsfit = bfit,
    data = brmshmmdata,
    data_processed = prepdata
  ), class = "brmshmmfit")
}


validate_brmshmmfit <- function(fit) {
  if(!inherits(fit, "brmshmmfit")) {
    stop("Not of class brmshmmfit")
  }
  fit
}


