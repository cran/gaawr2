#!/usr/bin/bash

# convert -density 300 logo.svg logo.png
# https://www.freeconvert.com/image-converter

Rscript -e '
      knitr::knit("README.Rmd")
      library(pkgdown)
    # init_site()
    # roxygen2::roxygenise()
      devtools::document()
      build_site()
  '
