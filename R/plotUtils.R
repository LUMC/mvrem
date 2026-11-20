# Pieces shared by the plots (R/plot.R, R/plotResiduals.R, R/plotResults.R, R/htmlPlot.R and
# plot.mvremPredictions() in R/predict.R): the colours, the panel grid, the choice of the
# subjects drawn, the per-subject series, panel titles and the legends of the line plots.

# ---- colours ---------------------------------------------------------------------------------

# Base graphics and the widget use different tones for the same element ON PURPOSE: base
# graphics draws the subject lines semi-transparent (where many overlap they darken), the
# widget draws them opaque and light (its highlight fades the other subjects instead). The
# blue accent (fitted / predicted values, estimates) and the ellipse colour are shared. The
# widget's JavaScript (inst/htmlwidgets/mvremPlot.js) has the same widget tones as defaults
# for lines without a `color`, and the accent for intervals and points.
plotAccent <- "#1b6ca8"
plotColors <- list(
  accent = plotAccent,
  ellipse = "#e08a00",                                   # pairs plot: 95% ellipses
  densities = c(plotAccent, "#d95f02", "#1b9e77", "#7570b3", "#e6ab02", "#666666"),
  base = list(
    subject = grDevices::adjustcolor("grey40", 0.45),   # observed values / residuals per subject
    subjectFaint = grDevices::adjustcolor("grey40", 0.35),   # observations behind predictions
    fit = grDevices::adjustcolor(plotAccent, 0.55),     # fitted values per subject
    mean = "grey15",                                    # smooth observed mean / mean residual
    interval = "grey60",                                # caterpillar intervals
    forestInterval = "grey50",                          # forest intervals (drawn thicker)
    reference = "grey50"),                              # zero lines
  widget = list(
    subject = "#c9c9c9",
    mean = "#222222",
    accentLight = "#9cc3e3"),   # fitted values per subject, prediction bounds
  # colouring by a grouping variable (colorBy): one colour per level (at most maxColorLevels)
  groups = c("#1b6ca8", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
             "#a6761d", "#17becf", "#bcbd22", "#8c564b", "#555555"))

# ---- panel grid ------------------------------------------------------------------------------

# Default number of panel columns: one row for up to two panels, then two or three columns.
autoNcol <- function(np, ncol = NULL) {
  if (!is.null(ncol)) return(ncol)
  if (np <= 2) np else if (np <= 6) 2L else 3L
}

# Base graphics: a grid of `np` panels (`ncol` columns, default autoNcol()); returns the old
# `par()` settings for on.exit().
panelLayout <- function(np, ncol = NULL) {
  ncol <- autoNcol(np, ncol)
  graphics::par(mfrow = c(ceiling(np / ncol), ncol), mar = c(4, 4, 2.5, 1), mgp = c(2.4, 0.8, 0))
}

# Base graphics: an empty panel with only a title (e.g. a survival response in plot()).
emptyPanel <- function(title) {
  graphics::plot.new()
  graphics::title(title, cex.main = 0.9)
}

# ---- subjects --------------------------------------------------------------------------------

# Default number of subjects drawn: the widget stays readable with many more lines, as the
# highlight picks one out.
defaultSubjectCount <- function(n, interactive) {
  if (!is.null(n)) return(n)
  if (interactive) 300L else 30L
}

# The subjects drawn: `ids` if given, else at most `n` levels of `group` (over all
# responses' data), evenly spaced over the sorted levels so that a subset still covers the
# whole range of ids (and the same subjects are drawn in every panel).
shownSubjects <- function(prep, group, n, ids = NULL) {
  if (!is.null(ids)) return(as.character(ids))
  if (is.na(group)) return(character())
  levs <- sort(unique(unlist(lapply(prep, function(df) as.character(df[[group]])))))
  levs[unique(round(seq(1, max(length(levs), 1), length.out = min(n, length(levs)))))]
}

# The subject of each row of `df` (as character); all rows one pseudo-subject `fill` without
# a grouping variable.
subjectOfRows <- function(df, group, fill = "all") {
  if (is.na(group) || is.null(df[[group]])) rep(fill, nrow(df)) else as.character(df[[group]])
}

# One series per subject in `ids`: the rows of the subject (`gg`) that are `keep`, ordered by
# `x`; `list(id, x, <values>)` with each element of the named list `values` subset alike
# (a NULL element stays NULL).
subjectSeries <- function(ids, gg, keep, x, values, cv = NULL) {
  lapply(ids, function(g) {
    k <- which(gg == g & keep)
    k <- k[order(x[k])]
    c(list(id = g, x = x[k]), lapply(values, function(v) if (!is.null(v)) v[k]),
      # values of the colouring variables at the points (colorBy), named by variable
      if (length(cv)) list(cv = lapply(cv, function(v) v[k])))
  })
}

# ---- labels ----------------------------------------------------------------------------------

# Panel title of an equation: `lhs  (outcome, beta)`.
equationTitle <- function(p) {
  paste0(p$lhs, "  (", if (p$part == "out") "outcome" else "longitudinal", ", ", p$family$family, ")")
}

# The denominators of a response shown as a proportion / rate (`trials()`, then
# `exposure()`), deparsed; empty for other responses.
rateDenominators <- function(p) {
  vapply(Filter(Negate(is.null), list(p$aux$trials, p$aux$exposure)), deparse1, "")
}

# ---- legends ---------------------------------------------------------------------------------

# A legend entry in both renderers: `kind` and `color` / `dash` for the widget (the kind is
# what a click on the entry switches off), `col` / `lwd` / `lty` for base graphics.
legendEntry <- function(kind, label, color, dash, col, lwd, lty, pch = NA) {
  list(kind = kind, label = label, color = color, dash = dash, col = col, lwd = lwd, lty = lty,
       pch = pch, shape = if (is.na(pch)) "line" else "circle")
}

# Whether some panel of plotLinesData() has predictions of missing observations.
hasPredicted <- function(d) any(vapply(d$panels, function(pn) length(pn$predicted) > 0L, NA))

# Legend of the line plots (plot(), plotResiduals()): the subjects and their smooth mean,
# with a fit also the fitted values and their mean. `kinds` / `labels` rename the observed
# elements (residual plots).
linesLegend <- function(haveFit, kinds = c(obs = "obs", obsMean = "obsMean"),
                        labels = c(obs = "observations (subjects)", obsMean = "observed mean"),
                        predicted = FALSE) {
  b <- plotColors$base
  w <- plotColors$widget
  c(list(legendEntry(kinds[["obs"]], labels[["obs"]], w$subject, "", b$subject, 1, 1),
         legendEntry(kinds[["obsMean"]], labels[["obsMean"]], w$mean, "", b$mean, 2.5, 1)),
    if (haveFit)
      list(legendEntry("fit", "fitted (subjects)", w$accentLight, "2 3", b$fit, 1.2, 3),
           legendEntry("fitMean", "fitted mean", plotColors$accent, "6 3", plotColors$accent, 2.5, 2)),
    if (haveFit && predicted)
      list(legendEntry("fitMissing", "predicted (missing observation)", plotColors$accent, "",
                       plotColors$accent, 1, 0, pch = 1)))
}

# The legend in the widget's format (htmlPlot(legend = )).
widgetLegend <- function(entries) lapply(entries, function(e) e[c("kind", "label", "color", "dash", "shape")])

# The legend in base graphics (top left of the first panel).
baseLegend <- function(entries) {
  graphics::legend("topleft", bty = "n", cex = 0.75, seg.len = 1.5,
                   legend = vapply(entries, `[[`, "", "label"),
                   col = vapply(entries, `[[`, "", "col"),
                   lwd = vapply(entries, `[[`, 0, "lwd"),
                   lty = vapply(entries, `[[`, 0, "lty"),
                   pch = vapply(entries, function(e) if (is.null(e$pch)) NA_real_ else e$pch, 0))
}

# ---- line panels (base graphics) -------------------------------------------------------------

# Whether the points of a subject are joined by lines: panels may say so (`connect`, residuals
# against the fitted mean are not), else whenever there are subjects (`grouped`).
panelConnect <- function(pn) if (!is.null(pn$connect)) pn$connect else pn$grouped

# Base-graphics renderer of the line panels of plotLinesData() / plotResidualsData(), the
# counterpart of htmlLines(): subject lines (observed, fitted), the smooth means, a reference
# line `hline` and the legend `entries` (linesLegend()) in the first panel.
drawLinesPanels <- function(d, ncol, entries) {
  op <- panelLayout(length(d$panels), ncol)
  on.exit(graphics::par(op), add = TRUE)
  b <- plotColors$base
  for (i in seq_along(d$panels)) {
    pn <- d$panels[[i]]
    if (isTRUE(pn$missing)) {
      emptyPanel(pn$title)
      next
    }
    graphics::plot(pn$xlim, pn$ylim, type = "n", xlab = pn$xlab, ylab = pn$ylab,
                   main = pn$title, cex.main = 0.95)
    if (!is.null(pn$hline)) graphics::abline(h = pn$hline, lty = 3, col = b$reference)
    if (!is.null(pn$note)) graphics::mtext(pn$note, side = 3, adj = 1, line = 0.1, cex = 0.6)
    connect <- panelConnect(pn)
    spec <- d$colors
    for (s in pn$subjects) {
      # colorBy: points by their level; the line too when the level is the same at all points
      sc <- seriesColors(s, spec, b$subject); pc <- sc$points; lc <- sc$line
      if (connect) graphics::lines(s$x, s$y, col = lc, lwd = 0.8)
      graphics::points(s$x, s$y, col = pc, pch = 16, cex = 0.5 * if (is.null(s$size)) 1 else s$size)
      if (pn$useFit && pn$grouped) graphics::lines(s$x, s$fit, col = b$fit, lwd = 1.2, lty = 3)
      if (pn$useFit && !pn$grouped) graphics::points(s$x, s$fit, col = b$fit, pch = 16, cex = 0.4)
    }
    if (!is.null(pn$smoothObs)) graphics::lines(pn$smoothObs, col = b$mean, lwd = 2.5)
    if (!is.null(pn$smoothFit)) graphics::lines(pn$smoothFit, col = plotColors$accent, lwd = 2.5, lty = 2)
    for (s in pn$predicted) graphics::points(s$x, s$fit, col = plotColors$accent, pch = 1, cex = 0.8)
    if (i == 1L) baseLegend(entries)
    if (i == 1L) colorLegend(d$colors)
  }
}

# Smooth curve (lowess) of `y` against `x` over the finite pairs; NULL when there are too few
# points or distinct x values for lowess() to be meaningful.
smoothCurve <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L || length(unique(x[ok])) < 3L) return(NULL)
  stats::lowess(x[ok], y[ok])
}

# One x range for all panels of a plot against time: the panels show the same variable, so
# they share its range (and therefore its tick marks, which both renderers derive from the
# range). Panels without data (`missing`) or without a finite range are left as they are.
shareXlim <- function(panels) {
  has <- vapply(panels, function(p) !isTRUE(p$missing) && length(p$xlim) == 2L &&
                  all(is.finite(p$xlim)), NA)
  if (sum(has) < 2L) return(panels)
  r <- range(unlist(lapply(panels[has], `[[`, "xlim")))
  panels[has] <- lapply(panels[has], function(p) { p$xlim <- r; p })
  panels
}

# ---- observation weights ---------------------------------------------------------------------

# The `weights()` of the observation rows of response `r` (the stack's `scale` column, one value
# per prepared data row), or NULL when the response has none.
observationWeights <- function(inl, r) {
  inp <- inl$inputs[[r]]
  if (!"scale" %in% names(inp)) return(NULL)
  inl$stacks[[r]]$main[[inp[["scale"]]]]
}

# Relative point radius for weights `w`: the point AREA is proportional to the weight, relative
# to the median weight (a point of median weight keeps the usual size); clamped to 0.3 - 3 times
# the usual radius so that extreme weights stay visible and do not cover the panel.
weightSizes <- function(w) {
  if (is.null(w)) return(NULL)
  m <- stats::median(w[is.finite(w) & w > 0])
  if (!is.finite(m)) return(NULL)
  pmin(pmax(sqrt(w / m), 0.3), 3)
}

# Panel note of a weighted response.
weightNote <- function(p) paste0("point area proportional to ", deparse1(p$aux$weights))

# ---- colouring by a grouping variable (colorBy) ------------------------------------------------

maxColorLevels <- 12L

# The variables a line plot can be coloured by: categorical columns (factor, character, logical)
# of the model's data with 2 to maxColorLevels distinct values over all responses (so not the
# subject id), in the order of appearance; `.key` excluded.
colorCandidates <- function(prep) {
  cand <- unique(unlist(lapply(prep, function(df)
    names(df)[vapply(df, function(v) is.factor(v) || is.character(v) || is.logical(v), NA)])))
  cand <- setdiff(cand, keyColumn)
  ok <- vapply(cand, function(v) {
    n <- length(colorLevelsOf(prep, v))
    n >= 2L && n <= maxColorLevels
  }, NA)
  cand[ok]
}

# The levels of variable `v` over all data frames: factor levels in their order, else sorted.
colorLevelsOf <- function(prep, v) {
  vals <- lapply(prep, `[[`, v)
  vals <- Filter(Negate(is.null), vals)
  levs <- unique(unlist(lapply(vals, function(x) if (is.factor(x)) levels(droplevels(x)) else character())))
  other <- sort(unique(unlist(lapply(vals, function(x) if (!is.factor(x)) as.character(x)))))
  setdiff(unique(c(levs, other)), NA_character_)
}

# The colouring of a line plot: the candidate variables, their levels, the palette and the
# selected variable (`colorBy`, NULL = grey), checked.
colorSpec <- function(prep, colorBy = NULL) {
  vars <- colorCandidates(prep)
  if (!is.null(colorBy)) {
    if (!isString(colorBy) || !colorBy %in% vars)
      stop("`colorBy` must be one of the categorical variables of the data with 2 to ",
           maxColorLevels, " levels: ",
           if (length(vars)) paste(vars, collapse = ", ") else "(none in this model)", ".", call. = FALSE)
  }
  list(vars = vars, levels = stats::setNames(lapply(vars, colorLevelsOf, prep = prep), vars),
       palette = plotColors$groups, selected = colorBy)
}

# The values of the colouring variables in the rows of `df` (as character; NULL when `df` lacks
# the variable), for subjectSeries(cv = ).
colorColumns <- function(df, spec) {
  Filter(Negate(is.null), stats::setNames(lapply(spec$vars, function(v)
    if (!is.null(df[[v]])) as.character(df[[v]])), spec$vars))
}

# Base graphics: the colours of a subject series under colorBy: `points` by their level, `line`
# = the level's colour when it is the same at all points, else `default`.
seriesColors <- function(s, spec, default, alpha = 0.8) {
  if (is.null(spec$selected)) return(list(points = default, line = default))
  pc <- grDevices::adjustcolor(colorOf(s$cv[[spec$selected]], spec, default), alpha)
  list(points = pc, line = if (length(unique(pc)) == 1L) pc[1L] else default)
}

# Base graphics: the legend of the colorBy levels (top right; nothing when no variable is selected).
colorLegend <- function(spec) {
  if (is.null(spec$selected)) return(invisible())
  levs <- spec$levels[[spec$selected]]
  graphics::legend("topright", bty = "n", cex = 0.7, pch = 16, title = spec$selected,
                   legend = levs, col = colorOf(levs, spec, plotColors$base$subject))
}

# Base graphics: the colour of each value `vals` of the selected variable (`default` when no
# variable is selected, the value is missing or the series has no such variable).
colorOf <- function(vals, spec, default) {
  v <- spec$selected
  if (is.null(v) || is.null(vals)) return(rep(default, max(1L, length(vals))))
  k <- match(vals, spec$levels[[v]])
  out <- spec$palette[(k - 1L) %% length(spec$palette) + 1L]
  out[is.na(k)] <- default
  out
}
