stan_jm_with_cache <- function(formulaLong, dataLong, formulaEvent, dataEvent, cache_file = NULL, ...) {
  fit <- NULL
  if(!is.null(cache_file) && file.exists(cache_file)) {
    cache_contents <- readRDS(cache_file)
    if(!is.list(cache_contents) || is.null(cache_contents$formulaLong) ||
       is.null(cache_contents$dataLong) ||
       is.null(cache_contents$formulaEvent) ||
       is.null(cache_contents$dataEvent) ||
       is.null(cache_contents$fit)) {
      message("Invalid cache content")
    } else {
      if(identical(formulaLong, cache_contents$formulaLong) &&
         identical(dataLong, cache_contents$dataLong) &&
         identical(formulaEvent, cache_contents$formulaEvent) &&
         identical(dataEvent, cache_contents$dataEvent)
         ) {
        fit <- cache_contents$fit
      } else {
        message("Cache file out of date, refitting")
      }
    }
  }

  if(is.null(fit)) {
    fit <- rstanarm::stan_jm(formulaLong = formulaLong, dataLong = dataLong,
                              formulaEvent = formulaEvent, dataEvent = dataEvent,
                              ...)
    if(!is.null(cache_file)) {
      saveRDS(loo::nlist(formulaLong, dataLong, formulaEvent, dataEvent, fit), file = cache_file)
    }
  }

  fit
}
