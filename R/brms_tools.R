brm_with_cache <- function(cache_file = NULL, ...) {
  fit <- NULL
  stancode <- brms::make_stancode(...)
  standata <- brms::make_standata(...)
  if(!is.null(cache_file) && file.exists(cache_file)) {
    cache_contents <- readRDS(cache_file)
    if(!is.list(cache_contents) || is.null(cache_contents$stancode) ||
       is.null(cache_contents$standata) || is.null(cache_contents$fit)) {
      message("Invalid cache content")
    } else if(!identical(normalize_stancode(stancode), normalize_stancode(cache_contents$stancode))) {
      message("Model code out of date, refitting")
    } else if(!identical(standata, cache_contents$standata)) {
      message("Model data out of date, refitting")
    } else {
      fit <- cache_contents$fit
    }
  }

  if(is.null(fit)) {
    fit <- brms::brm(...)
    if(!is.null(cache_file) && brms:::contains_samples(fit)) {
      saveRDS(loo::nlist(stancode, standata, fit), file = cache_file)
    }
  }

  fit
}

normalize_stancode <- function(x) {
  # Remove single-line comments
  x <- gsub("//[^\n\r]*[\n\r]", " ", x)
  x <- gsub("//[^\n\r]*$", " ", x)
  # Remove multi-line comments
  # x <- gsub("/\\*([^*]*[^/])*\\*/", " ", x)
  x <- gsub("/\\*([^*]*(\\*[^/])?)*\\*/", " ", x)
  # Standardize whitespace (including newlines)
  x <- gsub("[[:space:]]+"," ", x)

  trimws(x)
}
