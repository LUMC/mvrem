#' Fitted values and fit criteria of a model
#'
#' `fittedValues()` returns the posterior summaries of the fitted means (or
#' linear predictors) of each observation, together with the data of the
#' observation, `residuals()` adds the residuals (observed value minus fitted
#' mean) and `modelFit()` gives the criteria for model comparison.
#'
#' A latent equation (see [latentSpec()]) has no observations; `fittedValues(x, "eta")`
#' gives the posterior of its latent field at its evaluation points (one row per
#' combination of its variables in the data of the outcomes linked to it; no `observed`
#' column; `type` makes no difference).
#'
#' Fitted means are on the scale of the response's mean: for `trials()`
#' responses a proportion and for `exposure()` responses a rate (the observed
#' value in column `observed` is on the same scale). Rows with a missing
#' response get a prediction. For a survival response the fitted value is
#' INLA's \eqn{\lambda = \exp(\eta)} (hazard scale, see [families]), not a
#' mean time, and `observed` is `NA`; survival responses have no residuals.
#'
#' The criteria need INLA to compute them: fit with `jsm(..., criteria = TRUE)`
#' (or `inlaOptions = list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))`).
#'
#' @inheritParams randomEffects
#' @param type `"mean"` (default): fitted mean of the response; `"linear"`: the
#'   linear predictor.
#' @return
#' * `fittedValues()`: a data frame with `response`, the variables of the
#'   equation's data, `observed` (for `type = "mean"`) and the posterior
#'   `mean`, `sd`, `q0.025`, `q0.5`, `q0.975`, `mode`; one row per observation.
#' * `residuals()`: the same (`type = "mean"`) with a column `residual` =
#'   `observed - mean` (`NA` for a missing response); see [plotResiduals()].
#' * `modelFit()`: a data frame with one row per observed response, an
#'   `observations` row (the observation rows of all responses together) and an
#'   `overall` row (INLA's totals over all rows of the model, which include the
#'   zero-equation and latent-equation rows: pseudo-observations whose number grows
#'   with the links of the model, so compare models with the `observations` row or the
#'   rows of the responses they share): `dic`, `waic`, `logScore` (minus the sum of log
#'   CPO, smaller is better), `cpoFailures` (observations whose CPO INLA flags as
#'   unreliable) and, for the overall row, the log marginal likelihood `mlik`.
#' @examples
#' \dontrun{
#' m <- jsm(list(y ~ t + scaled(x1)), list(x1 ~ t + (1 + t | id)), data = d, criteria = TRUE)
#' head(fittedValues(m, "y"))
#' head(residuals(m, "y"))
#' modelFit(m)
#' }
#' @name fittedValues
NULL

#' @rdname fittedValues
#' @export
fittedValues <- function(x, ...) UseMethod("fittedValues")
#' @rdname fittedValues
#' @export
fittedValues.mvremModel <- function(x, response = NULL, type = c("mean", "linear"), ...) {
  needFit(x)
  type <- match.arg(type)
  # latent equations (E38) only when named: they have no observations
  resp <- selectResponses(x, response, observed = is.null(response))
  latent <- intersect(resp, names(latentParts(x$forms)))
  resp <- setdiff(resp, latent)
  tab <- statTable(if (type == "mean") x$fit$summary.fitted.values else
    x$fit$summary.linear.predictor)
  st <- stackRows(x$inla)
  block <- as.character(st$.block)
  prep <- bothParts(x$data)
  out <- list()
  for (r in resp) {
    rows <- observationRows(block, r)
    df <- prep[[r]]
    df <- df[setdiff(names(df), ".key")]
    # one stack row per prepared data row (buildInlaFormula()); anything else is a bug
    if (length(rows) != nrow(df))
      stop("Internal error: ", length(rows), " stack rows for ", nrow(df), " data rows of '", r, "'.",
           call. = FALSE)
    part <- cbind(data.frame(response = r, stringsAsFactors = FALSE), df)
    if (type == "mean") {
      y <- st[[x$inla$responses[[r]]]][rows]
      aux <- x$inla$inputs[[r]]
      if ("Ntrials" %in% names(aux)) y <- y / st[[aux[["Ntrials"]]]][rows]
      if ("E" %in% names(aux)) y <- y / st[[aux[["E"]]]][rows]
      if (isTRUE(x$inla$survival[r])) y <- NA_real_   # a survival time is not on the scale of lambda
      part$observed <- y
    }
    out[[r]] <- cbind(part, tab[rows, , drop = FALSE])
  }
  for (r in latent) out[[r]] <- latentValues(x, r)
  out <- out[orderResponses(x, names(out))]
  out <- rbindFill(out)
  rownames(out) <- NULL
  out
}

# The responses with residuals: `response` (checked; NULL = all) without the
# survival responses (a time to event has no residual). Survival responses named
# explicitly give a warning, or an error when nothing is left.
residualResponses <- function(x, response) {
  given <- !is.null(response)
  response <- selectResponses(x, response, observed = TRUE)
  surv <- intersect(response, survivalResponses(x$inla))
  response <- setdiff(response, surv)
  if (!length(response))
    stop("No residuals: the selected responses are survival responses.", call. = FALSE)
  if (given && length(surv))
    warning("No residuals for survival response(s) ", paste(surv, collapse = ", "), "; skipped.", call. = FALSE)
  response
}

#' @rdname fittedValues
#' @param object A fitted model (class `mvremModel`), for `residuals()`.
#' @export
residuals.mvremModel <- function(object, response = NULL, ...) {
  response <- residualResponses(object, response)
  fv <- fittedValues(object, response, type = "mean")
  fv$residual <- fv$observed - fv$mean
  fv
}

#' @rdname fittedValues
#' @export
modelFit <- function(x, ...) UseMethod("modelFit")
#' @export
modelFit.mvremModel <- function(x, ...) {
  needFit(x)
  f <- x$fit
  if (is.null(f$dic) || is.null(f$waic) || is.null(f$cpo))
    stop("The model was fitted without model criteria. Fit it with `jsm(..., criteria = TRUE)` ",
         "or `refit(m, inlaOptions = list(control.compute = list(dic = TRUE, waic = TRUE, ",
         "cpo = TRUE)))`.", call. = FALSE)
  block <- as.character(stackRows(x$inla)$.block)
  one <- function(rows) {
    cpo <- f$cpo$cpo[rows]
    ok <- is.finite(cpo) & cpo > 0
    data.frame(dic = sum(f$dic$local.dic[rows], na.rm = TRUE),
               waic = sum(f$waic$local.waic[rows], na.rm = TRUE),
               logScore = -sum(log(cpo[ok])),
               cpoFailures = sum(f$cpo$failure[rows] > 0, na.rm = TRUE),
               mlik = NA_real_)
  }
  resp <- selectResponses(x, NULL, observed = TRUE)
  obsRows <- lapply(resp, function(r) observationRows(block, r))
  tab <- do.call(rbind, lapply(obsRows, one))
  # Q5 (E38): the observation rows of all responses together, and INLA's totals over all rows
  # (also the zero-equation and latent rows, pseudo-observations with a fixed precision)
  obs <- one(unlist(obsRows))
  all <- one(seq_along(block))
  all$dic <- f$dic$dic
  all$waic <- f$waic$waic
  all$mlik <- unname(f$mlik[1L])
  out <- cbind(response = c(resp, "observations", "overall"), rbind(tab, obs, all),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
