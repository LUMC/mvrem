# Survival curves of survival outcomes (E31): the Kaplan-Meier estimate of the observed data
# and the model's marginal survival function, as a data frame (survivalCurves()) and as the
# panel of a survival response in plot() (survivalPanel(); both renderers).

#' Survival curves of a survival outcome
#'
#' The Kaplan-Meier estimate of the observed survival times and the survival
#' function of the fitted model, averaged over the subjects:
#' \deqn{\bar S(t) = \frac{1}{n} \sum_i \exp(-\lambda_i t^\alpha),}{S(t) = mean_i exp(-lambda_i t^alpha),}
#' with \eqn{\log \lambda_i} the linear predictor of subject \eqn{i} (incl. its
#' random effects and the `scaled()` associations; \eqn{\alpha = 1} for
#' [exponentialSurvFamily()]). Comparing the two curves checks the fit of the
#' baseline hazard and of the covariates together. `plot()` draws them in the
#' panel of a survival response.
#'
#' The model curve is computed at the posterior means of the linear predictors
#' and of \eqn{\alpha} (plug-in). When the fit keeps INLA's integration points
#' (`jsm(..., inlaOptions = list(control.compute = list(config = TRUE)))`, see
#' [predict.mvremModel()]) the curve is instead the posterior mean of
#' \eqn{\bar S(t)} over `nsamples` joint posterior draws, with a pointwise 95%
#' credible band.
#'
#' The association uses the marker at the covariate values of the subject's
#' survival row (see [weibullSurvFamily()]), so the curve shows the model as it
#' was fitted, not a current-value joint model.
#'
#' @param x A fitted `jsm` model with at least one survival outcome.
#' @param response Survival responses (default: all of them).
#' @param colorBy `NULL` or the name of a categorical variable of the survival data
#'   (2 to 12 levels, one value per subject; as in [plot.mvremModel()]): one Kaplan-Meier
#'   and one model curve per level (a missing value is the level `"(missing)"`).
#' @param nTimes Number of time points of the model curve (from 0 to the largest
#'   observed time).
#' @param nsamples Number of posterior draws for the band (fits with
#'   `config = TRUE` only).
#' @param seed Random seed of the draws (`NULL`, the default: the current state), as in
#'   [predict.mvremModel()].
#' @return A data frame of class `mvremSurvival`: `response`, `level` (the level
#'   of `colorBy`, `"all"` without it), `curve` (`"km"` or `"model"`), `time`,
#'   `surv`, `lower`, `upper` (the model's 95% band; for `"km"` the
#'   Greenwood interval on the log(-log) scale), `n.risk`, `n.event`,
#'   `n.censor` (Kaplan-Meier rows: one per distinct observed time, plus time 0).
#'   Attributes `band` (TRUE if the model curve has a band), `nsamples` and `colorBy`.
#' @seealso [plot.mvremModel()], [weibullSurvFamily()], [fittedValues()]
#' @export
survivalCurves <- function(x, response = NULL, colorBy = NULL, nTimes = 101L, nsamples = 200L,
                           seed = NULL) {
  needFit(x)
  checkCount(nTimes, "nTimes", 2L)
  checkCount(nsamples, "nsamples", 2L)
  surv <- survivalResponses(x$inla)
  if (!length(surv)) stop("The model has no survival outcome.", call. = FALSE)
  response <- if (is.null(response)) intersect(modelResponses(x), surv) else selectResponses(x, response)
  bad <- setdiff(response, surv)
  if (length(bad))
    stop("Not a survival response: ", paste(bad, collapse = ", "), ".", call. = FALSE)
  band <- hasConfig(x)
  dr <- NULL
  if (band) {
    if (!is.null(seed)) {
      old <- saveSeed()
      on.exit(restoreSeed(old), add = TRUE)
      set.seed(seed)
    }
    dr <- posteriorDraws(x$fit, nsamples)
  }
  out <- lapply(response, function(r) {
    sdat <- survivalData(x, r, colorBy)
    grid <- seq(0, max(sdat$time), length.out = nTimes)
    do.call(rbind, lapply(levels(sdat$level), function(l) {
      k <- sdat$level == l
      km <- kaplanMeier(sdat$time[k], sdat$event[k])
      mc <- modelSurvival(x, r, k, grid, dr, nsamples)
      rbind(data.frame(response = r, level = l, curve = "km", km, stringsAsFactors = FALSE),
            data.frame(response = r, level = l, curve = "model", time = grid, surv = mc$surv,
                       lower = mc$lower, upper = mc$upper, n.risk = NA_real_, n.event = NA_real_,
                       n.censor = NA_real_, stringsAsFactors = FALSE))
    }))
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  structure(res, class = c("mvremSurvival", "data.frame"), band = band,
            nsamples = if (band) nsamples, colorBy = colorBy)
}

#' @export
print.mvremSurvival <- function(x, ...) {
  print(as.data.frame(unclass(x)[names(x)], stringsAsFactors = FALSE), ...)
  cat(survivalNote(attr(x, "band"), attr(x, "nsamples")), "\n", sep = "")
  invisible(x)
}

# One line on how the model curve was computed (print, panel note).
survivalNote <- function(band, nsamples) {
  if (isTRUE(band)) sprintf("model: mean over subjects, 95%% band from %d posterior draws", nsamples)
  else "model: mean over subjects at posterior means (band: fit with config = TRUE)"
}

# The survival data of response r: time, event (0 / 1), subject id, level of `colorBy` (a
# factor with the levels of the data; a missing value is the level "(missing)", B25), in the
# rows of the prepared data (= the observation stack).
survivalData <- function(x, r, colorBy = NULL, group = NA_character_) {
  inl <- x$inla
  st <- inl$stacks[[r]]$main
  df <- plotData(x)[[r]]
  lev <- NULL
  if (!is.null(colorBy)) {
    cand <- groupCandidates(list(df))
    if (!isString(colorBy) || !colorBy %in% cand)
      stop("`colorBy` = '", paste(colorBy, collapse = ", "), "' is not a categorical variable (2 to ",
           maxColorLevels, " levels) of the data of '", r, "'",
           if (length(cand)) paste0("; candidates: ", paste(cand, collapse = ", ")), ".", call. = FALSE)
    v <- df[[colorBy]]
    lev <- if (is.factor(v)) droplevels(v) else factor(v)
    if (anyNA(lev))
      lev <- factor(ifelse(is.na(lev), "(missing)", as.character(lev)), c(levels(lev), "(missing)"))
  }
  data.frame(time = st[[inl$responses[[r]]]], event = st[[inl$inputs[[r]][["event"]]]],
             id = subjectOfRows(df, group, fill = NA_character_),
             level = if (is.null(lev)) factor(rep("all", nrow(st))) else lev,
             stringsAsFactors = FALSE)
}

# Kaplan-Meier estimate: one row per distinct observed time (and time 0) with the survival
# just after it, Greenwood's 95% interval on the log(-log) scale and the counts.
kaplanMeier <- function(time, event) {
  ok <- is.finite(time) & !is.na(event)
  time <- time[ok]; event <- event[ok]
  tt <- sort(unique(time))
  nr <- vapply(tt, function(u) sum(time >= u), 0)
  ne <- vapply(tt, function(u) sum(time == u & event == 1), 0)
  nc <- vapply(tt, function(u) sum(time == u & event == 0), 0)
  s <- cumprod(1 - ne / nr)
  v <- cumsum(ifelse(nr > ne, ne / (nr * (nr - ne)), 0))       # Greenwood: var(log S)
  se <- sqrt(v) / abs(log(s))
  z <- stats::qnorm(0.975)
  lo <- ifelse(s > 0 & s < 1, s^exp(z * se), NA_real_)
  hi <- ifelse(s > 0 & s < 1, s^exp(-z * se), NA_real_)
  data.frame(time = c(0, tt), surv = c(1, s), lower = c(1, lo), upper = c(1, hi),
             n.risk = c(length(time), nr), n.event = c(0, ne), n.censor = c(0, nc))
}

# The model's survival curve averaged over the rows `k` of response r at the times `grid`:
# with posterior draws `dr` the mean and 95% band of the average over draws, else the plug-in
# at the posterior means of the linear predictor and of alpha.
modelSurvival <- function(x, r, k, grid, dr = NULL, nsamples = 0L) {
  inl <- x$inla
  # the Weibull shape: INLA's row when estimated, else its fixed value (B26); 1 = exponential
  alphaName <- NULL
  alphaFixed <- 1
  if (identical(x$forms$out[[r]]$family$kind, "weibull")) {
    ht <- hyperTable(x)
    alphaName <- ht$.inla[ht$response == r & ht$type %in% "family" & ht$short %in% "alpha"]
    if (!length(alphaName)) {
      alphaName <- NULL
      alphaFixed <- fixedFamilyHyper(x, r, "alpha")
      if (is.null(alphaFixed)) stop("Internal: no Weibull shape found for '", r, "'.", call. = FALSE)
    }
  }
  pt <- function(eta, alpha) {                                  # rows x times
    lam <- exp(eta)
    colMeans(exp(-outer(lam, grid^alpha)))
  }
  if (is.null(dr)) {
    rows <- observationRows(as.character(stackRows(inl)$.block), r)
    eta <- x$fit$summary.linear.predictor$mean[rows][k]
    alpha <- if (is.null(alphaName)) alphaFixed else x$fit$summary.hyperpar[alphaName, "mean"]
    return(list(surv = pt(eta, alpha), lower = NA_real_, upper = NA_real_))
  }
  smp <- etaSampler(x, inl, dr, nsamples)
  eta <- smp$etaRows(r, inl$stacks[[r]]$main[k, , drop = FALSE])
  alpha <- if (is.null(alphaName)) rep(alphaFixed, nsamples) else smp$hyperRow(alphaName, "Weibull shape")
  S <- vapply(seq_len(nsamples), function(s) pt(eta[, s], alpha[s]), numeric(length(grid)))
  q <- apply(S, 1L, stats::quantile, c(0.025, 0.975), names = FALSE)
  list(surv = rowMeans(S), lower = q[1L, ], upper = q[2L, ])
}

# The panel of survival response r in plot(): Kaplan-Meier step curves (per level of the
# colour variable `colorBy` when it is a variable of the survival data), the censored times as
# open circles on them (one series per subject: linked to the subject's lines in the other
# panels), and with a fit the model curves and band. Drawn by drawSurvivalPanel() / htmlSurvival().
survivalPanel <- function(x, r, groups, group, haveFit) {
  p <- bothParts(x$forms)[[r]]
  df <- plotData(x)[[r]]
  colorBy <- groups$color
  if (!is.null(colorBy) && is.null(df[[colorBy]])) colorBy <- NULL
  sd <- survivalData(x, r, colorBy, group)
  cur <- if (haveFit) tryCatch(survivalCurves(x, r, colorBy), error = function(e) {
    warning("No model survival curve for '", r, "': ", conditionMessage(e), call. = FALSE)
    NULL
  })
  if (is.null(cur)) cur <- do.call(rbind, lapply(levels(sd$level), function(l) {
    k <- sd$level == l
    data.frame(response = r, level = l, curve = "km", kaplanMeier(sd$time[k], sd$event[k]),
               stringsAsFactors = FALSE)
  }))
  cv <- groupColumns(df, groups)
  # censored subjects: the open circle at the Kaplan-Meier value of its level just after its time
  cens <- which(sd$event == 0)
  ids <- ifelse(is.na(sd$id), paste0("row ", seq_len(nrow(sd))), sd$id)
  censored <- lapply(cens, function(i) {
    km <- cur[cur$curve == "km" & cur$level == as.character(sd$level[i]), ]
    list(id = ids[i], x = sd$time[i], y = km$surv[match(sd$time[i], km$time)],
         cv = if (length(cv)) lapply(cv, function(v) v[i]))
  })
  list(response = r, missing = FALSE, survival = TRUE, title = equationTitle(p),
       xlab = p$lhs, ylab = "survival probability", xlim = c(0, max(sd$time, na.rm = TRUE)),
       ylim = c(0, 1), curves = cur, by = colorBy, censored = censored,
       note = if (haveFit && any(cur$curve == "model"))
         survivalNote(attr(cur, "band"), attr(cur, "nsamples")))
}

# Step-function coordinates of a Kaplan-Meier curve (time, surv), continued to `tmax`.
kmSteps <- function(time, surv, tmax = max(time)) {
  n <- length(time)
  list(x = c(rep(time, each = 2L)[-1L], tmax), y = c(rep(surv, each = 2L)[-2L * n], surv[n]))
}

# Colour of the curves of level `l` (base graphics): the colour of the level of `by`, black
# for the Kaplan-Meier curve without it, the accent for the model.
survivalColor <- function(groups, by, l, curve) {
  if (!is.null(by)) return(colorOf(l, groups, plotColors$base$mean))
  if (curve == "km") plotColors$base$mean else plotColors$accent
}

# Legend entries of the survival panels (both renderers).
survivalLegend <- function(haveModel, band) {
  c(list(legendEntry("km", "Kaplan-Meier", plotColors$widget$mean, "", plotColors$base$mean, 2, 1),
         legendEntry("censored", "censored", plotColors$widget$mean, "", plotColors$base$mean, 1, 0, pch = 1)),
    if (haveModel) list(legendEntry("surv", "model survival", plotColors$accent, "6 3", plotColors$accent, 2, 2)),
    if (haveModel && band) list(legendEntry("survInt", "95% band", plotColors$widget$accentLight, "2 3",
                                            plotColors$accent, 1, 3)))
}

# The survival legend entries of the panels of plotLinesData() (none without a survival panel).
survivalEntries <- function(d) {
  sv <- Filter(function(pn) isTRUE(pn$survival), d$panels)
  if (!length(sv)) return(list())
  cur <- do.call(rbind, lapply(sv, `[[`, "curves"))
  survivalLegend(any(cur$curve == "model"), any(is.finite(cur$lower[cur$curve == "model"])))
}

# Base graphics: the survival panel `pn` of survivalPanel(), with its own legend (bottom left).
drawSurvivalPanel <- function(pn, groups) {
  panelFrame(pn)
  cur <- pn$curves
  for (l in unique(cur$level)) {
    km <- cur[cur$curve == "km" & cur$level == l, ]
    st <- kmSteps(km$time, km$surv, pn$xlim[2L])
    graphics::lines(st$x, st$y, col = survivalColor(groups, pn$by, l, "km"), lwd = 2)
    md <- cur[cur$curve == "model" & cur$level == l, ]
    if (nrow(md)) {
      col <- survivalColor(groups, pn$by, l, "model")
      graphics::lines(md$time, md$surv, col = col, lwd = 2, lty = 2)
      if (any(is.finite(md$lower))) {
        graphics::lines(md$time, md$lower, col = col, lwd = 1, lty = 3)
        graphics::lines(md$time, md$upper, col = col, lwd = 1, lty = 3)
      }
    }
  }
  for (s in pn$censored) {
    l <- if (is.null(pn$by)) "all" else s$cv[[pn$by]]
    graphics::points(s$x, s$y, pch = 1, cex = 0.7, col = survivalColor(groups, pn$by, l, "km"))
  }
  e <- survivalLegend(any(cur$curve == "model"), any(is.finite(cur$lower[cur$curve == "model"])))
  graphics::legend("bottomleft", bty = "n", cex = 0.75, seg.len = 1.5,
                   legend = vapply(e, `[[`, "", "label"), col = vapply(e, `[[`, "", "col"),
                   lwd = vapply(e, `[[`, 0, "lwd"), lty = vapply(e, `[[`, 0, "lty"),
                   pch = vapply(e, `[[`, 0, "pch"))
}

# Widget: the survival panel `pn` as a widget panel: the curves are hoverable elements (`id`)
# carrying the level of `by` (`cv`: coloured and hidden with it by the colour dropdown), the
# censored times open circles linked to their subjects.
htmlSurvival <- function(pn) {
  cur <- pn$curves
  curves <- list()
  for (l in unique(cur$level)) {
    lab <- if (is.null(pn$by)) "" else paste0(", ", pn$by, " = ", l)
    cv <- if (!is.null(pn$by)) list(cv = stats::setNames(list(l), pn$by), fade = 1)
    km <- cur[cur$curve == "km" & cur$level == l, ]
    st <- kmSteps(km$time, km$surv, pn$xlim[2L])
    curves[[length(curves) + 1L]] <- c(list(id = paste0("Kaplan-Meier", lab), kind = "km",
                                            label = paste0("Kaplan-Meier", lab), color = plotColors$widget$mean,
                                            dash = "", width = 2, x = r5(st$x), y = r5(st$y)), cv)
    md <- cur[cur$curve == "model" & cur$level == l, ]
    if (nrow(md)) {
      curves[[length(curves) + 1L]] <- c(list(id = paste0("model", lab), kind = "surv",
                                              label = paste0("model survival", lab), color = plotColors$accent,
                                              dash = "6 3", width = 2, x = r5(md$time), y = r5(md$surv)), cv)
      if (any(is.finite(md$lower)))
        for (b in c("lower", "upper"))
          curves[[length(curves) + 1L]] <- c(list(id = paste0("model ", b, " 95%", lab), kind = "survInt",
                                                  label = paste0("model survival, ", b, " 95% bound", lab),
                                                  color = plotColors$widget$accentLight, dash = "2 3", width = 1.3,
                                                  x = r5(md$time), y = r5(md[[b]])), cv)
    }
  }
  lines <- lapply(pn$censored, function(s)
    c(list(id = s$id, kind = "censored", x = r5(s$x), y = r5(s$y), color = plotColors$widget$mean,
           hollow = TRUE), if (length(s$cv)) list(cv = s$cv, tint = TRUE)))
  widgetPanel(pn, lines, curves, connect = FALSE)
}
