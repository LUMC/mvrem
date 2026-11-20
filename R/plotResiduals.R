#' Residual plots of a fitted model
#'
#' One panel per equation with the residuals of each observation, observed value
#' minus fitted mean (see [residuals.mvremModel()]), against the time variable
#' (`type = "time"`, the points of a subject connected) or against the fitted
#' mean (`type = "fitted"`), with a smooth curve of the mean residual and a line
#' at zero. Residuals are on the scale of the response's mean: proportions for
#' `trials()` responses, rates for `exposure()` responses. Observations with a
#' missing response have no residual and are not drawn.
#'
#' By default the plot is an HTML widget (see [htmlPlot()]): the subject under
#' the mouse is identified and highlighted in all panels, as in
#' [plot.mvremModel()].
#'
#' @inheritParams plot.mvremModel
#' @param type `"time"` (default): residuals against the time variable;
#'   `"fitted"`: against the fitted mean.
#' @return The widget, or `x` invisibly for base graphics.
#' @examples
#' \dontrun{
#' plotResiduals(m)
#' plotResiduals(m, type = "fitted")
#' plotResiduals(m, "x1", interactive = FALSE)   # base graphics
#' }
#' @export
plotResiduals <- function(x, ...) UseMethod("plotResiduals")

#' @rdname plotResiduals
#' @export
plotResiduals.mvremModel <- function(x, response = NULL, type = c("time", "fitted"), time = NULL,
                                     group = NULL, n = NULL, ids = NULL, ncol = NULL,
                                     interactive = NULL, colorBy = NULL, ...) {
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  d <- plotResidualsData(x, response, type, time, group, defaultSubjectCount(n, interactive), ids,
                         colorBy)
  ncol <- autoNcol(length(d$panels), ncol)
  kinds <- c(obs = "res", obsMean = "resMean")
  labels <- c(obs = "residuals (subjects)", obsMean = "mean residual")
  if (isTRUE(interactive)) return(htmlLines(d, ncol, kinds = kinds, labels = labels))
  drawLinesPanels(d, ncol, linesLegend(d$haveFit, kinds, labels))
  invisible(x)
}

# The panels of plotResiduals(), in the format of plotLinesData(): per equation the
# subjects drawn (x = time or fitted mean, y = residual), the limits over all rows
# and the smooth mean residual.
plotResidualsData <- function(x, response, type, time, group, n, ids, colorBy = NULL) {
  needFit(x)
  prep <- bothParts(x$data)
  colors <- colorSpec(prep, colorBy)
  eqs <- bothParts(x$forms)
  response <- residualResponses(x, response)
  time <- if (type == "time") plotTimeVariable(eqs, prep, time, timeVariable(x)) else NA_character_
  group <- plotGroupVariable(eqs, prep, group)
  shown <- shownSubjects(prep, group, n, ids)

  res <- stats::residuals(x, response)
  panels <- list()
  for (r in response) {
    p <- eqs[[r]]
    df <- res[res$response == r, , drop = FALSE]
    if (type == "time" && !time %in% names(prep[[r]])) {   # `res` has the columns of ALL responses
      panels[[r]] <- list(response = r, missing = TRUE,
                          title = paste0(r, ": no variable `", time, "`"))
      next
    }
    xx <- if (type == "time") df[[time]] else df$mean
    yy <- df$residual
    gg <- subjectOfRows(df, group)
    draw <- if (is.na(group)) unique(gg) else intersect(shown, unique(gg))
    ok <- is.finite(xx) & is.finite(yy)
    w <- observationWeights(x$inla, r)
    size <- if (!is.null(w) && length(w) == nrow(df)) weightSizes(w)   # rows as in fittedValues()
    # (the colouring variables of this response only: `res` has the columns of all responses)
    subjects <- subjectSeries(draw, gg, ok, xx, list(y = yy, size = size),
                              cv = colorColumns(df[intersect(names(df), names(prep[[r]]))], colors))
    subjects <- Filter(function(s) length(s$x) > 0L, subjects)
    panels[[r]] <- list(
      response = r, missing = FALSE, xlab = if (type == "time") time else "fitted mean",
      ylab = paste(c(paste0("residual of ", p$lhs), rateDenominators(p)), collapse = " / "),
      # points of a subject are joined only along time (against the fitted mean the order
      # of the points means nothing)
      grouped = !is.na(group), connect = type == "time" && !is.na(group), hline = 0,
      title = equationTitle(p), note = if (!is.null(size)) weightNote(p),
      xlim = range(xx[ok]), ylim = range(c(0, yy[ok])),
      subjects = subjects, useFit = FALSE, smoothObs = smoothCurve(xx, yy), smoothFit = NULL)
  }
  panels <- unname(panels)
  if (type == "time") panels <- shareXlim(panels)   # against the fitted mean the scales differ
  list(panels = panels, haveFit = FALSE, time = time, group = group, colors = colors)
}
