## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----gaawr2_incl, echo=FALSE, warning=FALSE-----------------------------------
htmltools::tagList(
  htmltools::tags$iframe(
    id = "gapframe",
    src = "https://jinghuazhao.github.io/R/vignettes/gaawr2.html",
    style = "
      width: 100%;
      border: none;
      display: block;
      min-height: 200px;
    "
  ),
  htmltools::tags$script(htmltools::HTML("
window.addEventListener('message', function(event) {
  if (!event.data || event.data.type !== 'setHeight') return;
  var iframe = document.getElementById('gapframe');
  if (!iframe) return;
  iframe.style.height = Math.max(event.data.height, 200) + 'px';
});
  "))
)

