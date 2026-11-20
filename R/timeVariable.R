# The time variable of a model (E19).
#
# The model itself needs no time variable: zero equations match outcome rows to the
# longitudinal predictor by the values of the long formula's groups and covariates,
# and a survival outcome evaluates the marker at the covariate values of its row.
# The time variable is used by what looks at the data over time: the plots
# (plot(), plotResiduals(), plot() of predictions) and, later, a time-varying
# survival association (E17). jsm(time = ) declares it; without it the numeric
# covariate used by most formulae is chosen when the model is built.

# list(variable, source = "given" / "detected" / "none", candidates) for jsm().
resolveTimeVariable <- function(forms, dat, time = NULL) {
  eqs <- bothParts(forms)
  prep <- bothParts(dat)
  cand <- timeCandidates(eqs, prep)
  if (!is.null(time)) {
    checkTimeColumn(prep, time)
    if (!time %in% cand)
      stop("`time` = '", time, "' is not a numeric covariate of any formula",
           if (length(cand)) paste0(" (candidates: ", paste(cand, collapse = ", "), ")"), ".",
           call. = FALSE)
    return(list(variable = time, source = "given", candidates = cand))
  }
  if (!length(cand)) return(list(variable = NA_character_, source = "none", candidates = cand))
  if (length(cand) > 1L)
    message("Time variable: '", cand[1L], "' (the numeric covariate of most formulae; other ",
            "candidates: ", paste(cand[-1L], collapse = ", "), "). Set it with jsm(time = ).")
  list(variable = cand[1L], source = "detected", candidates = cand)
}

#' Time variable of a model
#'
#' The variable on the horizontal axis of [plot()][plot.mvremModel],
#' [plotResiduals()] and the plot of [predictions][predict.mvremModel]: the one
#' given by `jsm(time = )`, else the numeric covariate (including the variables
#' an outcome needs for its links) used by most formulae, chosen when the model
#' is built (ties: the first in the order of the formulae; a message names the
#' other candidates). The `time` argument of each plot function still overrides
#' it for that plot. The model itself does not use it: a longitudinal process is
#' linked to an outcome at the values of all its covariates and groups.
#'
#' @param x An `mvremModel` (e.g. from [jsm()]).
#' @param ... Unused.
#' @return The name of the time variable, `NA` if the model has no numeric
#'   covariate. Attribute `source`: `"given"`, `"detected"` or `"none"`.
#' @examples
#' \dontrun{
#' m <- jsm(list(y ~ t + scaled(x), x ~ t + age + (1 | id)), data = d, time = "t", fit = FALSE)
#' timeVariable(m)
#' }
#' @export
timeVariable <- function(x, ...) UseMethod("timeVariable")

#' @export
timeVariable.mvremModel <- function(x, ...) {
  tv <- x$time
  if (is.null(tv))          # objects built before E19
    tv <- suppressMessages(resolveTimeVariable(x$forms, x$data))
  structure(tv$variable, source = tv$source)
}
