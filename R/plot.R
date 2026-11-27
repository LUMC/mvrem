#' Plot the data of a model, and its fit
#'
#' One panel per equation (in the order the formulae were passed to [jsm()]): the
#' observations of each subject as a thin line with points against a time
#' variable, and a smooth curve of the observed mean. It works on a model that
#' was only constructed (`jsm(..., fit = FALSE)`), to look at the input data
#' before fitting. For a fitted model, the fitted values (INLA's posterior mean
#' of the mean of each observation, per subject) and their smooth average are added;
#' observations are drawn with solid lines, fitted values with dotted (subjects) or
#' dashed (mean) lines.
#'
#' Responses with `trials()` are shown as proportions (successes / trials) and
#' responses with `exposure()` as rates (counts / exposure), so that they are
#' comparable with the fitted values, which are on this scale. Transformed
#' responses (e.g. `log(w)`) are shown on the transformed scale.
#'
#' A survival response (one row per subject) is drawn against its own time: the
#' Kaplan-Meier curve with the censored times as open circles and, for a fitted
#' model, the model's survival curve averaged over the subjects, with a 95% band
#' when the fit keeps INLA's configurations (see [survivalCurves()]). With
#' `colorBy` a variable of the survival data, there is one curve per level.
#'
#' By default (`interactive = NULL`, when the suggested package htmlwidgets is
#' installed) the plot is an HTML widget (see [htmlPlot()]): the
#' subject under the mouse is identified (id, time, observed and fitted value)
#' and highlighted in all panels at once, a click pins a subject, a search box
#' finds subjects while typing (by id, or by the value of a subject-level
#' categorical variable such as a clinic) and lists them in a dropdown (a pick
#' highlights that subject or all subjects of that clinic), and the legend
#' switches elements on and off. The widget works offline. `interactive =
#' FALSE` draws with base graphics.
#'
#' @param x A model (class `mvremModel`), e.g. from [jsm()].
#' @param response Names of the responses to plot (default: all).
#' @param time Name of the variable on the horizontal axis. Default: the numeric
#'   covariate used by most equations (typically the time).
#' @param group Name of the grouping variable (subjects). Default: the grouping
#'   variable of most random-effect terms. `NA` draws no subject lines.
#' @param n Maximum number of subjects drawn (evenly spaced over the levels of
#'   the grouping variable), so that dense data stay readable. Default: 30, and
#'   300 for the interactive widget.
#' @param ids Subjects to draw (values of the grouping variable); overrides `n`.
#' @param fitted If `TRUE` (default) and the model is fitted, add the fitted values.
#' @param ncol Number of panel columns (default: automatic).
#' @param interactive `NULL` (default): an HTML widget with mouse-over
#'   identification of the subjects if the htmlwidgets package is installed,
#'   else base graphics; `TRUE` / `FALSE` force one of them (see [htmlPlot()]).
#' @param ... Ignored.
#' @param colorBy,shapeBy `NULL` (default: grey circles) or the name of a categorical variable
#'   of the data (factor, character or logical; 2 to 12 levels for `colorBy`, 2 to 6 for
#'   `shapeBy`; e.g. a clinic or a device, also a column of the data that no formula uses) by
#'   which the observed points are coloured / shaped; a subject's line takes the colour too when
#'   the variable is constant within the subject. The widget offers all such variables in two
#'   dropdowns ("color:", "shape:"); these arguments only set the initial choice. In the
#'   widget, clicking a level in the legend of the chosen variable hides that subset of the data
#'   (its points, and the subjects whose value is constant), clicking again shows it.
#' @return The widget, or `x` invisibly for base graphics.
#' @examples
#' \dontrun{
#' m <- jsm(list(y ~ 1 + t + scaled(x1)), list(x1 ~ 1 + t + (1 + t | id)),
#'          data = d, family = list(y = betaFamily()), fit = FALSE)
#' plot(m)               # the input data (HTML widget: hover identifies a subject)
#' m <- refit(m)
#' plot(m)               # with the fitted values
#' plot(m, interactive = FALSE)   # base graphics
#' }
#' @export
plot.mvremModel <- function(x, response = NULL, time = NULL, group = NULL, n = NULL,
                            ids = NULL, fitted = TRUE, ncol = NULL, interactive = NULL,
                            colorBy = NULL, shapeBy = NULL, ...) {
  interactive <- useInteractive(interactive)
  d <- plotLinesData(x, response, time, group, defaultSubjectCount(n, interactive), ids, fitted,
                     colorBy, shapeBy)
  ncol <- autoNcol(length(d$panels), ncol)
  if (isTRUE(interactive)) return(htmlLines(d, ncol))
  drawLinesPanels(d, ncol, linesLegend(d$haveFit, predicted = hasPredicted(d)))
  invisible(x)
}

# The data of the panels of plot.mvremModel(): per equation the subjects drawn
# (observed values, fitted values), the axis limits (all rows, not only the
# subjects drawn) and the smooth means. Shared by the base and HTML renderers.
plotLinesData <- function(x, response, time, group, n, ids, fitted, colorBy = NULL, shapeBy = NULL) {
  inl <- x$inla
  prep <- plotData(x)
  groups <- groupSpec(prep, colorBy, shapeBy)
  eqs <- bothParts(x$forms)
  response <- selectResponses(x, response, observed = TRUE)   # in the order of the model

  time <- plotTimeVariable(eqs, prep, time, timeVariable(x))
  group <- plotGroupVariable(eqs, prep, group)
  shown <- shownSubjects(prep, group, n, ids)

  haveFit <- isTRUE(fitted) && !is.null(x$fit) && !is.null(x$fit$summary.fitted.values)
  fv <- if (haveFit) x$fit$summary.fitted.values
  block <- if (haveFit) as.character(stackRows(inl)$.block)

  panels <- list()
  for (r in response) {
    p <- eqs[[r]]
    df <- prep[[r]]
    if (isTRUE(inl$survival[[r]])) {
      # Kaplan-Meier + the model's survival curve (E31), against the survival time
      panels[[r]] <- survivalPanel(x, r, groups, group, haveFit)
      next
    }
    if (!time %in% names(df)) {
      panels[[r]] <- list(response = r, missing = TRUE,
                          title = paste0(r, ": no variable `", time, "`"))
      next
    }
    # the observed values on the scale of the fitted values: proportions / rates
    st <- inl$stacks[[r]]$main
    y <- st[[inl$responses[[r]]]]
    aux <- inl$inputs[[r]]
    if ("Ntrials" %in% names(aux)) y <- y / st[[aux[["Ntrials"]]]]
    if ("E" %in% names(aux)) y <- y / st[[aux[["E"]]]]
    # the label names every denominator (both only possible through inlaFamily()), as plotResiduals()
    ylab <- paste(c(p$lhs, rateDenominators(p)), collapse = " / ")
    tt <- df[[time]]
    gg <- subjectOfRows(df, group)
    f <- if (haveFit) fv[observationRows(block, r), , drop = FALSE]
    useFit <- haveFit && nrow(f) == nrow(df)   # fitted rows must match the data rows one to one

    draw <- if (is.na(group)) unique(gg) else intersect(shown, unique(gg))
    # weights(): point area proportional to the weight of the observation
    size <- weightSizes(observationWeights(inl, r))
    subjects <- subjectSeries(draw, gg, !is.na(y), tt, list(y = y, fit = if (useFit) f$mean, size = size),
                              cv = groupColumns(df, groups))
    # INLA's predictions of the missing observations (response NA in a kept row): drawn as
    # open circles, not joined to the subject's lines
    predicted <- if (useFit && anyNA(y))
      Filter(function(s) length(s$x) > 0L, subjectSeries(draw, gg, is.na(y), tt, list(fit = f$mean),
                                                         cv = groupColumns(df, groups)))
    panels[[r]] <- list(
      response = r, missing = FALSE, xlab = time, ylab = ylab, grouped = !is.na(group),
      title = equationTitle(p),
      xlim = range(tt, finite = TRUE), ylim = range(c(y, if (useFit) f$mean), finite = TRUE),
      subjects = subjects, predicted = predicted, useFit = useFit,
      note = if (!is.null(size)) weightNote(p),
      smoothObs = smoothCurve(tt, y), smoothFit = if (useFit) smoothCurve(tt, f$mean))
  }
  list(panels = shareXlim(unname(panels)), haveFit = haveFit, time = time, group = group, groups = groups)
}

# The time variable of a plot (E19): the `time` argument of the plot (checked: a
# numeric column of the model's data), else the model's time variable (`default`,
# from jsm(time = ) or detected when the model was built; see timeVariable()),
# else the detection rule of detectTimeVariable().
plotTimeVariable <- function(eqs, prep, time, default = NULL) {
  if (!is.null(time)) return(checkTimeColumn(prep, time))
  if (length(default) == 1L && !is.na(default) &&
      any(vapply(prep, function(df) default %in% names(df), NA))) return(as.character(default))
  cand <- timeCandidates(eqs, prep)
  if (!length(cand))
    stop("No numeric covariate to put on the horizontal axis; give `time`.", call. = FALSE)
  cand[[1L]]
}

# `time` must name one variable of the model's data that is numeric wherever it occurs.
checkTimeColumn <- function(prep, time, arg = "`time`")
  checkColumnName(time, arg, prep, numeric = TRUE)

# `value` must be the name of one column of (some of) the data frames `frames` (`allowNA`: NA
# passes, e.g. `group = NA` = no subject lines; `numeric`: the column must be numeric wherever it
# occurs, B10); `what` names the frames in the message. Returns the name.
checkColumnName <- function(value, arg, frames, what = "the model's data", allowNA = FALSE,
                            numeric = FALSE) {
  if (allowNA && length(value) == 1L && is.na(value)) return(NA_character_)
  if (!isString(value)) stop(arg, " must be the name of one variable.", call. = FALSE)
  has <- vapply(frames, function(df) value %in% names(df), NA)
  if (!any(has))
    stop(arg, " = '", value, "' is not a ", if (identical(what, "the predictions")) "column" else "variable",
         " of ", what, ".", call. = FALSE)
  if (numeric && !all(vapply(frames[has], function(df) is.numeric(df[[value]]), NA)))
    stop(arg, " = '", value, "' must be a numeric variable (it is ",
         paste(unique(vapply(frames[has], function(df) class(df[[value]])[1L], "")), collapse = " / "),
         ").", call. = FALSE)
  value
}

# The candidate time variables, best first: the numeric covariates (incl. link
# variables) of the equations, ordered by the number of equations using them (ties:
# first appearance in the order of `eqs`).
timeCandidates <- function(eqs, prep) {
  # link variables include the GROUPS of the linked long formulae: not time candidates (B27)
  covariates <- unique(unlist(lapply(eqs, function(p) p$variables$covariate), use.names = FALSE))
  cand <- lapply(eqs, function(p) {
    df <- prep[[p$response]]
    v <- intersect(unique(c(p$variables$covariate, p$linkVariables)), covariates)
    v[vapply(v, function(z) z %in% names(df) && is.numeric(df[[z]]), NA)]
  })
  all <- unlist(cand, use.names = FALSE)
  if (!length(all)) return(character())
  first <- unique(all)
  cnt <- as.vector(table(all)[first])
  first[order(-cnt, seq_along(first))]
}

# Default grouping variable: the grouping term of most random-effect blocks.
plotGroupVariable <- function(eqs, prep, group) {
  if (!is.null(group)) return(checkColumnName(group, "`group`", prep, allowNA = TRUE))
  g <- unlist(lapply(eqs, reGroups), use.names = FALSE)
  if (!length(g)) {
    g <- unlist(lapply(eqs, function(p) p$variables$group), use.names = FALSE)
  }
  if (!length(g)) return(NA_character_)
  first <- unique(g)
  first[which.max(table(g)[first])]
}
