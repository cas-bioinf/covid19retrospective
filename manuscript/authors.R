get_authors_string <- function() {
  authors_csv <- read.csv(here::here("private_data", "authors.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")

  authors_csv$Affiliation = factor(authors_csv$Affiliation, levels = unique(authors_csv$Affiliation))

  max_affiliation <- 0

  authors_strings <- character( nrow(authors_csv))
  for(i in 1:nrow(authors_csv)) {
    if(as.integer(authors_csv$Affiliation[i]) > max_affiliation) {
      affiliation <- paste0("^[", authors_csv$Affiliation[i], "]")
      max_affiliation <- as.integer(authors_csv$Affiliation[i])
    } else {
      affiliation <- paste0("^", as.integer(authors_csv$Affiliation[i]), "^")
    }
    authors_strings[i] <- paste0(gsub(" ", "&nbsp;", authors_csv$Name[i]), affiliation)
  }

  paste0(authors_strings, collapse = ", ")
}
