#' Options of a latent equation
#'
#' A latent equation, written `latent(eta) ~ <fixed> + <random effects> + scaled(x1) + ...`,
#' is an unobserved process: its linear predictor combines fixed and random effects and scaled
#' copies of the longitudinal responses named in its `scaled()` terms, and outcomes are linked
#' to it with `scaled(eta)`, like to a longitudinal response. It has no data and no family; its
#' options are given in the `family` list of [jsm()] / [collectFormulas()]:
#' `family = list(eta = latentSpec(...))`.
#'
#' The latent process is evaluated at the values of its variables (its own groups and
#' covariates and those of its longitudinal targets) found in the data of the outcomes linked
#' to it, through one zero equation per value combination, `0 = <latent predictor> - u`, with a
#' latent field `u` shared by all these outcomes (each copies it with its own coefficient). Its
#' longitudinal targets have one zero-equation field each towards it.
#'
#' The location of the latent process is not identified when the outcomes have intercepts;
#' `center = TRUE` constrains `u` to sum to zero over its evaluation points. Its scale and sign
#' are identified only by the priors (`eta` times c with the outcome coefficients divided by c
#' fits equally well), unless a coefficient is fixed with `scaled(x, coef = 1)`.
#'
#' @param center Constrain the latent field to sum to zero over its evaluation points
#'   (default `TRUE`).
#' @param centerTargets Constrain the zero-equation field of each longitudinal target of the
#'   latent equation (its linear predictor at the latent's evaluation points) to sum to zero
#'   (default `FALSE`).
#' @param logPrec Fixed log precision of the latent field and of the zero-equation fields of
#'   its targets towards it: a small value makes them diffuse, so that they follow the zero
#'   equations. `NULL` (default): the model's value for all zero-equation fields (argument
#'   `zeroLogPrec` of [jsm()], default -6).
#' @return An object of class `mvremLatentSpec`.
#' @examples
#' latentSpec()
#' latentSpec(center = TRUE, centerTargets = TRUE, logPrec = -15)
#' @export
latentSpec <- function(center = TRUE, centerTargets = FALSE, logPrec = NULL) {
  checkFlag(center, "center")
  checkFlag(centerTargets, "centerTargets")
  if (!is.null(logPrec) && !(isNumber(logPrec) && is.finite(logPrec)))
    stop("`logPrec` must be NULL or one finite number (a log precision).", call. = FALSE)
  structure(list(center = center, centerTargets = centerTargets, logPrec = logPrec),
            class = "mvremLatentSpec")
}

#' @export
format.mvremLatentSpec <- function(x, ...) {
  paste0("latentSpec(center = ", x$center, ", centerTargets = ", x$centerTargets,
         ", logPrec = ", if (is.null(x$logPrec)) "NULL" else format(x$logPrec), ")")
}

#' @export
print.mvremLatentSpec <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# Default fixed log precision of every zero-equation field (diffuse; see zeroFieldHyper()).
defaultZeroLogPrec <- -6

# The zero-field log precision a model was built with (models built before E38: the default).
modelZeroLogPrec <- function(x)
  if (is.null(x$inla$zeroLogPrec)) defaultZeroLogPrec else x$inla$zeroLogPrec

# The zero-field log precision of a latent equation's fields: its own `logPrec`, else the
# model's value.
latentLogPrec <- function(p, zeroLogPrec)
  if (is.null(p$spec$logPrec)) zeroLogPrec else p$spec$logPrec

# The latent equations named in an outcome's scaled() terms.
latentTargets <- function(o, latent) intersect(scaledTargets(o), names(latent))

# The variables a latent equation is evaluated at: its groups and covariates and those of its
# longitudinal targets (the order of linkedVariables()).
latentVariables <- function(p, long) linkedVariables(p, long)

# Evaluation points of each latent equation: one row per unique combination of its variables
# (latentVariables()) in the data of the outcomes linked to it, in the order of first appearance
# (outcomes in model order). Numeric columns stay numeric (covariates); others become
# character. A list named by latent response (empty if none).
latentPoints <- function(forms, data) {
  latent <- latentParts(forms)
  res <- list()
  for (L in names(latent)) {
    vars <- latentVariables(latent[[L]], forms$long)
    users <- names(forms$out)[vapply(forms$out, function(o) L %in% scaledTargets(o), NA)]
    cols <- lapply(vars, function(v) {
      vals <- lapply(users, function(o) data$out[[o]][[v]])
      if (all(vapply(vals, is.numeric, NA))) unlist(vals, use.names = FALSE)
      else unlist(lapply(vals, as.character), use.names = FALSE)
    })
    pts <- as.data.frame(stats::setNames(cols, vars), stringsAsFactors = FALSE, optional = TRUE)
    key <- do.call(paste, c(lapply(pts, as.character), sep = "\r"))
    pts <- pts[!duplicated(key), , drop = FALSE]
    rownames(pts) <- NULL
    res[[L]] <- pts
  }
  res
}

# The posterior of a latent equation's field at its evaluation points (fittedValues()):
# `response`, the point's variables, the summaries of the field matched by its labels.
latentValues <- function(x, r) {
  lf <- x$inla$latent[[r]]
  pts <- x$inla$latentPoints[[r]]
  sr <- x$fit$summary.random[[lf$index]]
  rows <- match(latentLabels(pts, r, lf$vars, x$inla$naming), sr$ID)
  if (anyNA(rows))
    stop("Internal error: latent field of '", r, "' lacks evaluation points.", call. = FALSE)
  tab <- statTable(sr[rows, setdiff(names(sr), "ID"), drop = FALSE])
  cbind(data.frame(response = rep(r, nrow(pts)), stringsAsFactors = FALSE), pts, tab)
}

# Is the scale of latent equation L fixed by a fixed coefficient of one of its links
# (scaled(L, coef = ) in an outcome, or scaled(x, coef = ) in the latent equation)?
latentScaleFixed <- function(forms, L) {
  hasCoef <- function(s) !is.null(s$coef)
  inLatent <- any(vapply(latentParts(forms)[[L]]$scaled, hasCoef, NA))
  inOutcomes <- any(unlist(lapply(forms$out, function(o)
    vapply(Filter(function(s) s$target == L, o$scaled), hasCoef, NA))))
  inLatent || inOutcomes
}
