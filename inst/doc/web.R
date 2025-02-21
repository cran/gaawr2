## ----setup, include=FALSE-----------------------------------------------------
set.seed(0)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.path = "plumber/",
  collapse = TRUE,
  comment = "#>",
  dev = "png")

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
pkgs <- c("httr","httpuv","jsonlite","plumber","seqminer")
for (p in pkgs) if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
    if (!requireNamespace(p)) warning(paste0("This vignette needs package `", p, "'; please install"))
}
invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))

