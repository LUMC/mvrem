#' Interactive HTML plots
#'
#' The plotting functions ([plot.mvremModel()], [plotResiduals()], [plot.mvremPredictions()],
#' [plotRandomEffects()], [plotVarianceComponents()], [plotFixedEffects()], [plotAssociations()]) return
#' an interactive HTML widget with `interactive = TRUE`. The widget is
#' self-contained (own JavaScript, no external libraries or network access), so it works
#' in the RStudio viewer, in R Markdown / Quarto HTML output (also with embedded
#' resources) and as a stand-alone file (`htmlwidgets::saveWidget()`).
#'
#' * **Hover**: a tooltip names the subject (or term) under the mouse with its values
#'   (time, observed and fitted value; posterior mean and interval; ...), and the
#'   subject is highlighted while the others fade. The highlight is **linked**: the
#'   same subject is highlighted in all panels of the widget (e.g. the marker, the
#'   outcome, the intercept and the slope panels). For a fitted model a subject has
#'   two lines per panel, both highlighted in red: its observations (solid) and
#'   its fitted values (dotted).
#' * **Click** pins a subject (several can be pinned); a click on a pinned
#'   subject, double click on the background or Escape clears the pins.
#' * **Search box**: type (part of) an id to highlight the matching subjects.
#' * **Legend**: click an entry to switch observed / fitted elements on and off.
#'
#' The widgets are the default when the suggested package htmlwidgets is
#' installed (`interactive = NULL`); without it, or with `interactive = FALSE`
#' or `options(mvrem.interactive = FALSE)`, the plots are drawn with base
#' graphics. `interactive = TRUE` without htmlwidgets is an error.
#'
#' @param panels List of panel specifications (internal format, see the code of
#'   the plotting functions).
#' @param ncol Number of panel columns.
#' @param legend Legend entries (`kind`, `label`, `color`, `dash`) or `NULL`.
#' @param linked If `TRUE`, ids are shared by the panels (linked highlighting).
#' @param groupBy The colour / shape grouping of the points (`groupSpec()`), or `NULL`.
#' @param search Label of the search box (`NULL`: none). The box highlights the
#'   elements whose id contains the typed text, or whose subject-level `groupBy`
#'   value does (a variable constant within every subject), and lists them in a
#'   dropdown (groups first); a pick (click, or arrow keys + Enter) highlights
#'   exactly that subject or group.
#' @return An `htmlwidget` (class `mvremPlot`).
#' @name htmlPlot
#' @keywords internal
NULL

# `interactive = NULL` (the default of the plotting functions): HTML widgets if htmlwidgets is
# installed (and `options(mvrem.interactive = FALSE)` is not set), base graphics otherwise.
# An explicit TRUE needs htmlwidgets (error), FALSE always draws with base graphics.
useInteractive <- function(interactive) {
  if (is.null(interactive))
    return(isTRUE(getOption("mvrem.interactive", TRUE)) && hasHtmlwidgets())
  if (isTRUE(interactive)) requireHtmlwidgets()
  isTRUE(interactive)
}

hasHtmlwidgets <- function() requireNamespace("htmlwidgets", quietly = TRUE)

requireHtmlwidgets <- function() {
  if (!hasHtmlwidgets())
    stop("`interactive = TRUE` needs the htmlwidgets package: install.packages(\"htmlwidgets\").",
         call. = FALSE)
}

#' @rdname htmlPlot
htmlPlot <- function(panels, ncol = 2L, legend = NULL, linked = TRUE, search = "subject",
                     groupBy = NULL) {
  requireHtmlwidgets()
  payload <- list(panels = unname(panels), ncol = as.integer(ncol), legend = unname(legend),
                  linked = linked, search = search)
  # colour / shape by a categorical variable (groupSpec()): the dropdowns list `vars` /
  # `shapeVars`, lines and points carry the values `cv`
  if (length(groupBy$vars)) payload$groupBy <- groupBy
  htmlwidgets::createWidget(
    "mvremPlot", payload, package = "mvrem",
    sizingPolicy = htmlwidgets::sizingPolicy(padding = 0, browser.fill = FALSE,
                                             knitr.figure = FALSE, knitr.defaultWidth = "100%"))
}

# Numbers of the widget spec: 5 significant digits keep the JSON (and the page) small.
r5 <- function(v) signif(as.numeric(v), 5)

# One panel of the widget spec from a panel of the plot data: title, axis labels and limits,
# `connect` (join a subject's points), its lines and curves, and the optional reference line /
# note. Shared by the line plots, the survival panel and the predictions.
widgetPanel <- function(pn, lines, curves = NULL, connect = TRUE, hline = NULL) {
  c(list(title = pn$title, xlab = pn$xlab, ylab = pn$ylab, xlim = r5(pn$xlim), ylim = r5(pn$ylim),
         connect = connect, lines = lines),
    if (!is.null(curves)) list(curves = curves),
    if (!is.null(hline)) list(hline = hline), if (!is.null(pn$note)) list(note = pn$note))
}

# The widget of subject panels (linked across panels, with the subject search box and the
# colour / shape dropdowns): the line plots and the predictions.
subjectWidget <- function(panels, ncol, legend, group, groups)
  htmlPlot(panels, ncol, widgetLegend(legend), linked = TRUE, search = subjectSearchLabel(group),
           groupBy = groups)

# panels of plot.mvremModel() (data from plotLinesData(); a survival panel goes to
# htmlSurvival()) and plotResiduals() (plotResidualsData(); `kinds` renames the observed
# elements, `labels` the legend); the base-graphics counterpart is drawLinesPanels(). Subject lines carry no colour: the widget
# draws them in its default tones per kind (observed / fitted), which the highlight overrides.
htmlLines <- function(d, ncol, kinds = c(obs = "obs", obsMean = "obsMean"),
                      labels = c(obs = "observations (subjects)", obsMean = "observed mean")) {
  panels <- lapply(d$panels, function(pn) {
    if (isTRUE(pn$missing)) return(list(title = pn$title, empty = TRUE))
    if (isTRUE(pn$survival)) return(htmlSurvival(pn))
    lines <- list()
    for (s in pn$subjects) {
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = kinds[["obs"]], x = r5(s$x), y = r5(s$y)),
                                       if (!is.null(s$size)) list(size = r5(s$size)),   # weights()
                                       if (length(s$cv)) list(cv = s$cv))               # colorBy / shapeBy
      if (!is.null(s$fit))
        lines[[length(lines) + 1L]] <- list(id = s$id, kind = "fit", x = r5(s$x), y = r5(s$fit))
    }
    # predictions of missing observations: open circles, not joined (`hollow`)
    for (s in pn$predicted)
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = "fitMissing", x = r5(s$x), y = r5(s$fit),
                                            color = plotColors$accent, hollow = TRUE),
                                       if (length(s$cv)) list(cv = s$cv))
    curve <- function(cv, kind, label, color, dash)
      if (!is.null(cv)) list(list(kind = kind, label = label, color = color, dash = dash,
                                  width = 3, x = r5(cv$x), y = r5(cv$y)))
    curves <- c(curve(pn$smoothObs, kinds[["obsMean"]], labels[["obsMean"]], plotColors$widget$mean, ""),
                curve(pn$smoothFit, "fitMean", "fitted mean", plotColors$accent, "6 3"))
    widgetPanel(pn, lines, if (length(curves)) curves else list(), panelConnect(pn), pn$hline)
  })
  subjectWidget(panels, ncol, c(linesLegend(d$haveFit, kinds, labels, hasPredicted(d)), survivalEntries(d)),
                d$group, d$groups)
}

# Label of the search box of the subject plots: the grouping variable ("id"), "subject" without one;
# the placeholder reads "highlight id" (+ " or clinic" for searchable subject-level variables).
subjectSearchLabel <- function(group) if (is.na(group)) "subject" else group

# panels with one interval per row: posterior mean and 95% interval (caterpillar / forest)
htmlIntervals <- function(title, id, label, mean, lo, hi, xlab = "posterior mean and 95% interval",
                          ylab = "", ylabels = FALSE) {
  tip <- sprintf("%s\n%s [%s, %s]", label, signif(mean, 4), signif(lo, 4), signif(hi, 4))
  p <- list(title = title, xlab = xlab, ylab = ylab, yReverse = TRUE, vline = 0,
            intervals = list(id = id, pos = seq_along(id), mean = r5(mean), lo = r5(lo), hi = r5(hi),
                             tip = tip))
  if (ylabels) p$ylabels <- list(pos = seq_along(id), text = label)
  p
}

# scatter panel (e.g. random intercept vs slope), one point per subject
htmlPoints <- function(title, id, x, y, xlab, ylab, cv = NULL) {
  list(title = title, xlab = xlab, ylab = ylab, vline = 0, hline = 0,
       points = c(list(id = id, x = r5(x), y = r5(y),
                       tip = sprintf("%s\n%s = %s, %s = %s", id, xlab, signif(x, 4), ylab, signif(y, 4))),
                  if (length(cv)) list(cv = cv)))     # colorBy / shapeBy: values per point
}
