brmhmm <- function(brmshmmdata, cache_file = NULL, ...) {
  d <- validate_brmshmmdata(brmshmmdata)

  prepdata <- make_data_hmm(d)

  bfit <- brm_with_cache(
    cache_file = cache_file,
    formula = brms::brmsformula(d$formula, family = rate_hmm_family),
    data = prepdata$brmsdata,
    prior = d$prior,
    stanvars = rate_hmm_stanvars(prepdata$standata),
    ...
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

summary.brmshmmfit <- function(fit) {
  validate_brmshmmfit(fit)
  structure(list(
    brmssummary = summary(fit$brmsfit)#,
    #hmmsummary = summary(fit$brmsfit$fit, pars =  c("sensitivity", "other_observations_probs"))
    ),
    class = "summary.brmshmmfit"
  )
}

print.summary.brmshmmfit <- function(s) {
  print(s$brmssummary)
  #print(s$hmmsummary$summary)
}
