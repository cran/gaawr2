## ----setup, include=FALSE-----------------------------------------------------
set.seed(0)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.height = 8,
  fig.path = "web/",
  fig.width = 8,
  collapse = TRUE,
  comment = "#>",
  dev = "CairoPNG")

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
pkgs <- c("httr","httpuv","jsonlite","plumber","seqminer")
for (p in pkgs) if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
    if (!requireNamespace(p)) warning(paste0("This vignette needs package `", p, "'; please install"))
}
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

## ----url----------------------------------------------------------------------
check_url_availability <- function(url) {
# Send a GET request to the URL
  response <- tryCatch({
    httr::GET(url)
  }, error = function(e) {
# If there was an error, return FALSE
    return(NULL)
  })
# If response is NULL or status code is not 200, return FALSE
  if (is.null(response) || httr::status_code(response) != 200) {
    message(paste("Warning: The URL", url, "is not accessible. Please check the link."))
    return(FALSE)
  }
# If status code is 200, the URL is available
  message(paste("The URL", url, "is accessible."))
  return(TRUE)
}

url_to_check <- "http://www.zotero.org/styles/nature-genetics"
is_available <- check_url_availability(url_to_check)

if (is_available) {
  message("Using CSL as usual.")
} else {
  message("Using fallback to local CSL file instead.")
}

