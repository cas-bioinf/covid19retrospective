sum_percent <- function(x) {
  paste0(sum(x, na.rm = TRUE), " (", scales::percent(mean(x, na.rm = TRUE), accuracy = 1), ")")
}

mean_iqr <- function(x) {
  paste0(round(mean(x, na.rm = TRUE)), " (", round(quantile(x, 0.25, na.rm = TRUE)), " - ", round(quantile(x, 0.75, na.rm = TRUE)), ")")
}
