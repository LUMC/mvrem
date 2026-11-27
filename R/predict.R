#' Predictions for new subjects or time points
#'
#' `predict()` gives the posterior of the fitted mean (or linear predictor) of
#' a response at new rows: new time points of subjects in the data (their
#' random effects are used) or new subjects (random effects from their
#' distribution, so wider intervals). Two methods, chosen automatically:
#'
#' * `method = "auto"` (default): `"sample"` when the model was fitted with
#'   INLA's configurations (`control.compute = list(config = TRUE)`, detected
#'   from the fit object), else `"refit"`. If sampling is not available for the
#'   model (see below), `"auto"` refits with a message. The method used is the
#'   attribute `method` of the result and is shown when it is printed.
#' * `method = "refit"`: INLA predicts a missing response, so the new
#'   rows are appended to the model's data with a missing response and the model
#'   is fitted again (one INLA run per `predict()` call: exact, but as slow as
#'   the fit). Outcome rows get their zero equations, so the linked longitudinal
#'   predictors are evaluated at the new rows as well.
#' * `method = "sample"`: no new INLA run. `nsamples` joint draws of the latent
#'   field and the hyperparameters are taken from the fitted posterior (INLA's
#'   Gaussian approximation at each integration point of the hyperparameters,
#'   points chosen by their posterior weight: the distribution of
#'   `INLA::inla.posterior.sample(skew.corr = FALSE)`, drawn in R, which is
#'   much faster). Per draw, the linear predictor of a new row is its fixed
#'   effects (design of the new rows with the fit's factor levels) plus its
#'   random effects: an existing subject has its sampled effects, a new subject
#'   gets effects drawn from \eqn{N(0, Q^{-1})}, \eqn{Q} built from the sampled
#'   precisions (and correlations of a correlated block; one draw per subject,
#'   shared by all its rows and responses). An outcome adds, for each `scaled()`
#'   link, the sampled \eqn{\gamma} times the linked longitudinal predictor at
#'   the row's values of that formula's groups and covariates (the matching of
#'   the zero equations). `type = "mean"` applies the inverse link per draw.
#'   The summaries are the mean, SD and quantiles of the draws; `mode` is the
#'   maximum of a kernel density estimate ([stats::density()]) of the draws.
#'   The random-effect blocks must have at most five correlated effects and
#'   estimated precisions, and the links no parameters; otherwise use
#'   `method = "refit"`. An explicit `method = "sample"` on a fit without the
#'   configurations is an error (with these hints).
#'
#' **Getting a fit with the configurations.** Fit with
#' `jsm(..., inlaOptions = list(control.compute = list(config = TRUE)))`, or
#' add them to a fitted model with
#' `m <- refit(m, inlaOptions = list(control.compute = list(config = TRUE)))`
#' (keeps the other options only if they are repeated: `inlaOptions` replaces
#' the model's) or `m <- update(m, inlaOptions = ...)` (repeats the whole call).
#'
#' **Cost, speed and accuracy.** The fit itself is not slower, but the fit
#' object keeps one sparse precision matrix per integration point and is
#' several times larger (the "multi" test model: 5.3 MB without, 25.2 MB with the
#' configurations). Prediction by sampling is then 5-10 times faster than a
#' refit (seconds instead of a full INLA run per `predict()` call). The two
#' methods agree within Monte-Carlo error, validated in the package tests
#' against the refit for new and existing subjects, longitudinal and outcome
#' responses, both `type`s: means within a few Monte-Carlo standard errors, SDs
#' and 95% bounds within a few percent of the refit's SD. Caveat: sampling uses
#' INLA's discrete integration points of the hyperparameters (no skewness
#' correction), so when a hyperparameter posterior is heavy-tailed (e.g. a
#' random-effect SD near zero, few groups) the tails of the predictions can
#' differ more from the refit's; check with `method = "refit"` if the intervals
#' matter, and increase `nsamples` for stable tail quantiles.
#'
#' `plot()` of the result draws the predicted trajectories: one panel per
#' response, the posterior mean per subject (solid) with its 95% interval
#' (dotted), over the observed data of the model (grey) when `data = TRUE`.
#'
#' @param object A fitted model (class `mvremModel`).
#' @param newdata A data frame with the variables of the equations to predict
#'   (covariates, grouping variables, `trials()` / `exposure()` inputs and, for
#'   outcomes, the variables of the linked longitudinal equations); the response
#'   itself is not needed (and ignored). Or a list of such data frames named by
#'   response. The variables of `weights()` are optional: a weight only enters the
#'   likelihood of an observed response, so it does not change a prediction; a
#'   variable used only by `weights()` and missing in `newdata` is set to 1 in the
#'   new rows.
#' @param response Responses to predict when `newdata` is one data frame
#'   (default: every response whose variables are all in `newdata`). The
#'   result lists the responses in the order the formulae were passed to
#'   [jsm()], whatever the order of `response` or of a `newdata` list.
#' @param type `"mean"` (default): the fitted mean of the response (proportion /
#'   rate scale for `trials()` / `exposure()` responses); `"linear"`: the linear
#'   predictor.
#' @param ... `predict()`: `verbose`, `threads`, `inlaOptions` of the new INLA
#'   run of `method = "refit"` (default: those of the model; with
#'   `method = "auto"` they select the refit). `print()`, `plot()`: ignored /
#'   passed to the data frame's print method.
#' @param method `"auto"` (default), `"refit"` or `"sample"`, see Details.
#' @param nsamples Number of posterior draws of `method = "sample"`.
#' @param seed Seed of the random draws of `method = "sample"` (`NULL`: the
#'   current state of R's random number generator); with a seed, the state of
#'   the generator is restored afterwards.
#' @return `predict()`: a data frame of class `mvremPredictions` with `response`,
#'   the columns of `newdata` and the posterior `mean`, `sd`, `q0.025`, `q0.5`,
#'   `q0.975`, `mode`; rows grouped by response in the order of the model;
#'   attributes `type`, `method` (the method used: `"refit"` or `"sample"`),
#'   `methodChosen` (`"auto"` if chosen automatically, else the method given),
#'   `nsamples` (sampling only) and `model`: the refitted model (`"refit"`) or
#'   `object` itself (`"sample"`).
#' @examples
#' \dontrun{
#' # d: a data frame with id, t, y, x1
#' new <- data.frame(id = rep(c("new", "S001"), each = 5), t = rep(0:4 / 2, 2))
#'
#' # a fit without INLA's configurations: predict() refits (one INLA run)
#' m <- jsm(form = list(y ~ t + scaled(x1), x1 ~ t + (1 + t | id)), data = d)
#' p <- predict(m, new)        # a new subject and an existing one, both responses
#' attr(p, "method")           # "refit"
#' plot(p)
#' predict(m, new, response = "x1", type = "linear")
#'
#' # a fit with the configurations: predict() samples (no INLA run, 5-10x faster)
#' mc <- jsm(form = list(y ~ t + scaled(x1), x1 ~ t + (1 + t | id)), data = d,
#'           inlaOptions = list(control.compute = list(config = TRUE)))
#' # or: mc <- update(m, inlaOptions = list(control.compute = list(config = TRUE)))
#' ps <- predict(mc, new)      # method "sample" chosen automatically
#' ps                          # the last printed line names the method
#' predict(mc, new, nsamples = 4000, seed = 1)   # more draws, reproducible
#' pr <- predict(mc, new, method = "refit")      # force the INLA run (e.g. to compare)
#' }
#' @name predict.mvremModel
#' @export
predict.mvremModel <- function(object, newdata, response = NULL, type = c("mean", "linear"), ...,
                               method = c("auto", "refit", "sample"), nsamples = 1000L, seed = NULL) {
  needFit(object)
  type <- match.arg(type)
  chosen <- match.arg(method)
  opts <- list(...)
  bad <- setdiff(names(opts), c("verbose", "threads", "inlaOptions"))
  if (length(bad)) stop("Unknown argument(s) to predict(): ", paste(bad, collapse = ", "), call. = FALSE)
  # "auto": sample when the fit keeps INLA's configurations, unless options of an INLA run
  # are given
  method <- if (chosen == "auto") { if (hasConfig(object) && !length(opts)) "sample" else "refit" } else chosen
  if (method == "sample") checkSampleArgs(opts, nsamples, seed)
  pr <- predictionModel(object, newdata, response)
  # E38: sampling does not follow a latent chain; an explicit "sample" says so first
  if (chosen == "sample" && length(latentParts(object$forms)))
    stop("`method = \"sample\"` does not support models with latent equations yet; use ",
         "`method = \"refit\"`.", call. = FALSE)
  if (method == "sample") needConfig(object)

  preds <- NULL
  x <- object
  if (method == "sample") {
    preds <- tryCatch(
      predictSampled(object, pr$model$inla, pr$newdata, pr$response, type, as.integer(nsamples), seed),
      mvremSampleUnsupported = function(e) {
        if (chosen != "auto") stop(e)
        message("predict(): sampling is not available for this model (", sub(";.*$", "", conditionMessage(e)),
                "); refitting with INLA instead.")
        NULL
      })
    if (is.null(preds)) method <- "refit"
  }
  if (method == "refit") {
    x <- pr$model
    x$options <- utils::modifyList(object$options, opts)
    x <- refit(x)
    if (is.null(x$fit))
      stop("INLA failed on the data with the new rows: ", firstLine(x$inlaError), call. = FALSE)
    preds <- refitPredictions(x, object, pr$newdata, pr$response, type)
  }

  out <- lapply(pr$response, function(r) {
    pred <- preds[[r]]
    rownames(pred) <- NULL
    cbind(data.frame(response = r, stringsAsFactors = FALSE), pr$newdata[[r]], pred)
  })
  # the summaries stay last: a column union in the order of the first part would put the
  # extra newdata columns of a later response after them
  res <- rbindFill(out, last = predictionStats)
  rownames(res) <- NULL
  structure(res, class = c("mvremPredictions", "data.frame"), model = x, type = type, method = method,
            methodChosen = chosen, nsamples = if (method == "sample") as.integer(nsamples))
}

# The summary columns of a prediction (as fittedValues()).
predictionStats <- c("mean", "sd", "q0.025", "q0.5", "q0.975", "mode")

# Arguments of predict(method = "sample"): no INLA-run options, a valid nsamples and seed.
checkSampleArgs <- function(opts, nsamples, seed) {
  if (length(opts))
    stop("`", paste(names(opts), collapse = "`, `"), "` only apply to `method = \"refit\"`.", call. = FALSE)
  if (!isNumber(nsamples) || nsamples < 2 ||
      nsamples != round(nsamples))
    stop("`nsamples` must be a single integer >= 2.", call. = FALSE)
  if (!is.null(seed) && !isNumber(seed))
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  invisible()
}

# The predictions of the refitted model `x` per response: predictionModel() APPENDS the new
# rows to each response's prepared data, so they are the last rows of fittedValues(x, r);
# checked against the number of rows of the original model.
refitPredictions <- function(x, object, newdata, response, type) {
  lapply(stats::setNames(response, response), function(r) {
    fv <- fittedValues(x, r, type = type)
    n <- nrow(newdata[[r]])
    old <- nrow(bothParts(object$data)[[r]])
    if (nrow(fv) != old + n)
      stop("Internal error: the refitted model of '", r, "' has ", nrow(fv), " rows, expected ",
           old, " + ", n, ".", call. = FALSE)
    fv[seq.int(old + 1L, old + n), predictionStats]
  })
}

#' @rdname predict.mvremModel
#' @export
print.mvremPredictions <- function(x, ...) {
  print.data.frame(x, ...)
  m <- attr(x, "method")
  if (!is.null(m)) {
    how <- if (identical(m, "sample")) paste0("by sampling (", attr(x, "nsamples"), " posterior draws)") else
      "by refitting with INLA"
    why <- if (identical(attr(x, "methodChosen"), "auto"))
      paste0("; chosen automatically", if (identical(m, "sample")) ": the fit has INLA's configurations" else
        if (!hasConfig(attr(x, "model"))) ": the fit has no INLA configurations") else ""
    cat("Predictions ", how, why, ".\n", sep = "")
  }
  invisible(x)
}

# The checked responses and new data of predict() and the model rebuilt on its
# data with the new rows appended (response missing; `fit` still the original's,
# replaced by refit()):
# list(model, newdata = list by response, response).
predictionModel <- function(object, newdata, response) {
  eqs <- bothParts(object$forms)
  prep <- bothParts(object$data)
  # variables used only by weights(): optional in newdata (set to 1 in the new rows)
  weightOnly <- function(p) if (is.null(p$aux$weights)) character() else
    setdiff(all.vars(p$aux$weights), c(p$variables$covariate, p$variables$group, p$linkVariables,
                                       unlist(lapply(p$aux[names(p$aux) != "weights"], all.vars))))
  needed <- function(p) setdiff(formulaVariables(p, response = FALSE), weightOnly(p))

  surv <- names(object$inla$survival)[object$inla$survival]
  if (is.data.frame(newdata)) {
    newdata <- plainDataFrame(newdata)
    if (is.null(response)) {
      response <- names(eqs)[vapply(eqs, function(p) all(needed(p) %in% names(newdata)), NA)]
      if (length(response) && all(response %in% surv))
        stop("`newdata` has only the variables of survival response(s) (",
             paste(response, collapse = ", "), "), which cannot be predicted.", call. = FALSE)
      response <- setdiff(response, surv)
      if (!length(response))
        stop("`newdata` has the variables of no equation; needed: ",
             paste(vapply(eqs, function(p) paste0(p$response, ": ", paste(needed(p), collapse = ", ")), ""),
                   collapse = "; "), ".", call. = FALSE)
    }
    response <- selectResponses(object, response, observed = TRUE)
    newdata <- stats::setNames(rep(list(newdata), length(response)), response)
  } else if (is.list(newdata) && !is.null(names(newdata)) && all(vapply(newdata, is.data.frame, NA))) {
    response <- selectResponses(object, names(newdata), observed = TRUE)
    newdata <- lapply(newdata, plainDataFrame)
  } else {
    stop("`newdata` must be a data frame or a list of data frames named by response.", call. = FALSE)
  }
  if (any(response %in% surv))
    stop("Prediction is not available for survival responses (", paste(intersect(response, surv), collapse = ", "),
         ").", call. = FALSE)
  for (r in response) {
    miss <- setdiff(needed(eqs[[r]]), names(newdata[[r]]))
    if (length(miss))
      stop("`newdata` for '", r, "' lacks the variable(s) ", paste0("`", miss, "`", collapse = ", "), ".",
           call. = FALSE)
    if (!nrow(newdata[[r]])) stop("`newdata` for '", r, "' has no rows.", call. = FALSE)
  }

  # grouping variables / terms may take new values (new subjects); any other factor or
  # character column must keep to the values seen in the fit, else the rebuilt model gets
  # a new fixed-effect column that only the new rows use (estimated from the prior alone)
  groupCols <- unique(unlist(lapply(eqs, function(p) c(p$variables$group, reGroups(p)))))
  seenValues <- function(v) unique(unlist(lapply(prep, function(df)
    if (v %in% names(df)) { if (is.factor(df[[v]])) levels(df[[v]]) else unique(as.character(df[[v]])) })))
  # the model's data per response with the new rows appended (response missing)
  full <- lapply(names(eqs), function(r) {
    base <- prep[[r]][setdiff(names(prep[[r]]), ".key")]
    if (!r %in% response) return(base)
    new <- newdata[[r]]
    for (v in setdiff(weightOnly(eqs[[r]]), names(new))) new[[v]] <- rep_len(1, nrow(new))
    for (v in setdiff(names(base), names(new))) new[[v]] <- NA
    new[[eqs[[r]]$variables$response]] <- NA
    for (v in setdiff(names(base), groupCols)) {
      if (!is.factor(base[[v]]) && !is.character(base[[v]])) next
      unseen <- setdiff(unique(as.character(new[[v]])), c(seenValues(v), NA))
      if (length(unseen))
        stop("Unknown level(s) of `", v, "` in `newdata` for response ", r, ": ",
             paste(unseen, collapse = ", "), " (the fit has: ", paste(seenValues(v), collapse = ", "),
             ").", call. = FALSE)
    }
    # factors become character so that rbind() keeps new levels (new subjects); prepareData()
    # re-creates the factors
    for (v in names(base)) if (is.factor(base[[v]]) && !is.factor(new[[v]])) new[[v]] <- as.character(new[[v]])
    rbind(base, new[names(base)])
  })
  names(full) <- names(eqs)
  dataOut <- full[names(object$forms$out)]
  dataLong <- full[names(object$forms$long)]
  x <- object
  x$data <- prepareData(object$forms, dataOut = dataOut, dataLong = dataLong, na.action = stats::na.pass)
  x$dataSummary <- summarizeData(object$forms, dataOut = dataOut, dataLong = dataLong)
  x$inla <- buildInlaFormula(object$forms, x$data, object$inla$naming,
                             modelZeroLogPrec(object))
  list(model = x, newdata = newdata, response = response)
}

#' @rdname predict.mvremModel
#' @param x The result of `predict()`.
#' @param response Responses to plot (default: all predicted ones), in the order of the model.
#' @param time,group Names of the time and grouping (subject) variables; default:
#'   those of [plot.mvremModel()].
#' @param ids Subjects (levels of `group`) to plot; default: all predicted subjects, and at most
#'   `n` observed ones.
#' @param data If `TRUE` (default), the observed data of the model are drawn in
#'   grey behind the predictions (at most `n` subjects).
#' @param n Maximum number of observed subjects drawn (default: 30, and 300 for
#'   the interactive widget).
#' @param ncol Number of panel columns (default: automatic).
#' @param interactive `NULL` (default): an HTML widget with mouse-over
#'   identification if the htmlwidgets package is installed, else base graphics;
#'   `TRUE` / `FALSE` force one of them (see [htmlPlot()]).
#' @param colorBy,shapeBy Name of a categorical variable of the data or the predictions (2 to 12
#'   levels for colour, 2 to 6 for shape) colouring / shaping the observed and predicted
#'   subjects, as in [plot.mvremModel()]; the widget offers all such variables in dropdowns
#'   whenever there are any, and clicking a level in their legends hides / shows that subset.
#' @export
plot.mvremPredictions <- function(x, response = NULL, time = NULL, group = NULL, data = TRUE, n = NULL,
                                  ids = NULL, ncol = NULL, interactive = NULL, colorBy = NULL,
                                  shapeBy = NULL, ...) {
  interactive <- useInteractive(interactive)
  n <- defaultSubjectCount(n, interactive)
  m <- attr(x, "model")
  eqs <- bothParts(m$forms)
  prep <- plotData(m)
  resp <- intersect(selectResponses(m, response, observed = TRUE), x$response)
  if (!length(resp)) stop("No predictions for the response(s) given.", call. = FALSE)
  time <- plotTimeVariable(eqs, prep, time, timeVariable(m))
  group <- plotGroupVariable(eqs, prep, group)
  checkColumnName(time, "`time`", list(x), what = "the predictions")
  linear <- identical(attr(x, "type"), "linear")
  # observed data: the model's rows without the appended prediction rows; not drawn for
  # predictions of the linear predictor, which is on another scale than the data
  obs <- if (isTRUE(data) && !linear) plotLinesData(m, resp, time, group, n, ids, fitted = FALSE)$panels
  # colorBy / shapeBy: the categorical variables of the data and of the new rows (a new clinic is a level)
  groups <- groupSpec(c(prep, list(x[setdiff(names(x), "response")])), colorBy, shapeBy)

  panels <- lapply(seq_along(resp), function(i) {
    r <- resp[i]
    d <- x[x$response == r, , drop = FALSE]
    gg <- subjectOfRows(d, group, fill = "prediction")
    if (!is.null(ids)) {                      # only these subjects' predictions
      d <- d[gg %in% as.character(ids), , drop = FALSE]
      gg <- gg[gg %in% as.character(ids)]
      if (!nrow(d)) return(NULL)
    }
    subjects <- subjectSeries(unique(gg), gg, TRUE, d[[time]],
                              list(mean = d$mean, lo = d$q0.025, hi = d$q0.975),
                              cv = groupColumns(d, groups))
    o <- if (!is.null(obs)) obs[[i]]
    # (subjects without an observed value of this response have empty series: dropped)
    observed <- if (!is.null(o) && !isTRUE(o$missing)) Filter(function(s) length(s$x) > 0L, o$subjects)
    ylab <- if (linear) paste0("linear predictor of ", eqs[[r]]$lhs) else if (!is.null(o) && !isTRUE(o$missing)) o$ylab else eqs[[r]]$lhs
    list(response = r, title = equationTitle(eqs[[r]]), note = if (length(observed)) o$note,
         xlab = time, ylab = ylab, subjects = subjects, observed = observed,
         xlim = range(c(d[[time]], unlist(lapply(observed, `[[`, "x"))), finite = TRUE),
         ylim = range(c(d$q0.025, d$q0.975, unlist(lapply(observed, `[[`, "y"))), finite = TRUE))
  })
  panels <- Filter(Negate(is.null), panels)
  if (!length(panels)) stop("No predictions for the subject(s) in `ids`.", call. = FALSE)
  panels <- shareXlim(panels)
  ncol <- autoNcol(length(panels), ncol)
  legend <- predictionsLegend(!is.null(obs))
  if (isTRUE(interactive)) return(htmlPredictions(panels, ncol, legend, group, groups))

  colData <- plotColors$base$subjectFaint
  colPred <- plotColors$accent
  withPanels(length(panels), ncol, function(i) {
    pn <- panels[[i]]
    panelFrame(pn)
    for (s in pn$observed) {
      sc <- seriesStyle(s, groups, colData, alpha = predictionFade)
      graphics::lines(s$x, s$y, col = sc$line, lwd = 0.8)
      graphics::points(s$x, s$y, col = sc$points, pch = sc$pch, cex = 0.5 * if (is.null(s$size)) 1 else s$size)
    }
    for (s in pn$subjects) {
      sc <- seriesStyle(s, groups, colPred, alpha = 1)
      graphics::lines(s$x, s$lo, col = sc$line, lwd = 1, lty = 3)
      graphics::lines(s$x, s$hi, col = sc$line, lwd = 1, lty = 3)
      graphics::lines(s$x, s$mean, col = sc$line, lwd = 2)
      graphics::points(s$x, s$mean, col = sc$points, pch = sc$pch, cex = 0.6)
    }
    if (i == 1L) { baseLegend(legend); groupLegends(groups) }
  })
  invisible(x)
}

# Legend of plot.mvremPredictions() in both renderers (see legendEntry()).
predictionsLegend <- function(haveObs) {
  w <- plotColors$widget
  c(if (haveObs) list(legendEntry("obs", "observations (subjects)", w$subject, "", plotColors$base$subjectFaint, 1, 1)),
    list(legendEntry("pred", "predicted mean", plotColors$accent, "", plotColors$accent, 2, 1),
         legendEntry("predInt", "95% interval", w$accentLight, "2 3", plotColors$accent, 1, 3)))
}

# Opacity of the coloured observations behind predictions (colorBy): faint, so that the predicted
# mean and bounds (opaque) stand out; used by the widget (`fade`) and the base plot.
predictionFade <- 0.2

# widget of plot.mvremPredictions(): observed lines (grey), predicted mean (solid blue) and
# 95% bounds (dotted light blue, without points) per subject, linked by subject id. The
# prediction lines carry their colours, as the widget's defaults are those of observed /
# fitted lines.
htmlPredictions <- function(panels, ncol, legend, group, groups = NULL) {
  ps <- lapply(panels, function(pn) {
    lines <- list()
    for (s in pn$observed)
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = "obs", x = r5(s$x), y = r5(s$y)),
                                       if (!is.null(s$size)) list(size = r5(s$size)),
                                       if (length(s$cv)) list(cv = s$cv, fade = predictionFade))   # colorBy / shapeBy
    for (s in pn$subjects) {
      cv <- if (length(s$cv)) list(cv = s$cv, fade = 1)   # predictions: opaque in their colour
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = "predInt", color = plotColors$widget$accentLight, dash = "2 3",
                                            points = FALSE, x = r5(s$x), y = r5(s$lo)), cv)
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = "predInt", color = plotColors$widget$accentLight, dash = "2 3",
                                            points = FALSE, x = r5(s$x), y = r5(s$hi)), cv)
      lines[[length(lines) + 1L]] <- c(list(id = s$id, kind = "pred", color = plotColors$accent, width = 2,
                                            x = r5(s$x), y = r5(s$mean)), cv)
    }
    widgetPanel(pn, lines)
  })
  subjectWidget(ps, ncol, legend, group, groups)
}
