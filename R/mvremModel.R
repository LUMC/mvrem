#' Accessing a fitted (or only constructed) model
#'
#' Generic functions to inspect objects of class `mvremModel`, such as the
#' result of [jsm()]. They are generics so that other model types of the package
#' can provide the same interface. No knowledge of INLA is needed to use them.
#'
#' An `mvremModel` is a list with the components `call`, `forms` (result of
#' [collectFormulas()]), `data` (result of [prepareData()]), `display` (categorical
#' columns of the data that no formula uses, kept for colouring / shaping the plots only),
#' `dataSummary`
#' ([summarizeData()]), `inla` (result of [buildInlaFormula()]), `time` (the time variable:
#' `variable`, `source`, `candidates`; see [timeVariable()]), `options` (how
#' INLA is run), `fit` (the INLA result, `NULL` if not fitted) and `inlaError`.
#' Prefer the functions below to `$`.
#'
#' Every result (tables, plots, predictions, data summary) lists the responses
#' in the order in which the formulae were passed to [jsm()]: with one list
#' (`form`) the order of the list, outcomes and longitudinal responses
#' interleaved as given; with the pair `formOut` / `formLong` the outcomes, then
#' the longitudinal responses. A `response` argument of a result function only
#' selects responses; they still come in this order. `modelResponses()` gives it.
#'
#' @param x An `mvremModel`, e.g. from [jsm()].
#' @param n Number of rows of each part of the stacked data (default: all).
#' @param ... Further arguments of the methods.
#' @return
#' * `modelFormulas()`: the parsed formulae ([collectFormulas()] result).
#' * `modelData()`: the data of each equation ([prepareData()] result).
#' * `modelResponses()`: the response names in the order the formulae were
#'   passed (also for a [collectFormulas()] result).
#' * `dataSummary()`: rows supplied and complete, groups, missing values and
#'   zero equations ([summarizeData()] result).
#' * `inlaFormula()`: the INLA formula (an R formula); `inlaSpec()`: everything
#'   built for INLA ([buildInlaFormula()] result: families, `control.family`,
#'   terms, stacked data, naming).
#' * `inlaData()`: the stacked data as one data frame with a column `.block`
#'   telling the equation each row belongs to.
#' * `inlaArguments()`: the exact arguments `INLA::inla()` is called with, for
#'   users who want to run INLA themselves.
#' * `inlaFit()`: the result of `INLA::inla()` (`NULL` if the model is not fitted).
#' * `fixedEffects()`: posterior summaries of the fixed effects, one row per
#'   response and term.
#' * `hyperparameters()`: posterior summaries of the hyperparameters (residual
#'   and family parameters, precisions and correlations of random effects,
#'   association coefficients): `response`, `parameter` (a readable name), `type` (`"family"`,
#'   `"random"`, `"association"`, or `NA` for a row mvrem does not recognise) and the posterior
#'   `mean`, `sd`, `q0.025`, `q0.5`, `q0.975`, `mode`; rows grouped by response in the order the
#'   formulae were passed. INLA's parametrisation is used: precisions, not standard deviations.
#' * `associations()`: the association coefficients \eqn{\gamma} of the
#'   `scaled()` terms, one row per outcome and longitudinal response.
#' * `refit()`: the model fitted again (see [jsm()] for the options).
#' @name mvremModel
NULL

#' @rdname mvremModel
#' @export
modelFormulas <- function(x, ...) UseMethod("modelFormulas")
#' @export
modelFormulas.mvremModel <- function(x, ...) x$forms

#' @rdname mvremModel
#' @export
modelResponses <- function(x, ...) {
  forms <- if (inherits(x, "mvremModel")) x$forms else x
  if (!is.null(forms$order)) return(forms$order)
  # model objects made before the order was stored: outcomes, then longitudinal
  c(names(forms$out), names(forms$long))
}

# The responses `r` in the order of the model (a selection does not reorder).
orderResponses <- function(x, r) {
  all <- modelResponses(x)
  all[all %in% r]
}

# The responses selected by a `response` argument (NULL = all), checked, in the order of
# the model (the order in which the formulae were passed), whatever the order given.
# `observed = TRUE`: only responses with data (a latent equation, E38, has none): NULL selects
# the observed responses, a latent one named explicitly is an error.
selectResponses <- function(x, response, observed = FALSE) {
  allResp <- modelResponses(x)
  latent <- names(latentParts(if (inherits(x, "mvremModel")) x$forms else x))
  if (observed) allResp <- setdiff(allResp, latent)
  if (is.null(response)) return(allResp)
  if (observed && any(response %in% latent))
    stop("'", intersect(response, latent)[1L], "' is a latent equation: it has no observations.",
         call. = FALSE)
  bad <- setdiff(response, allResp)
  if (length(bad))
    stop("Unknown response(s): ", paste(bad, collapse = ", "), ". Responses of the model: ",
         paste(allResp, collapse = ", "), ".", call. = FALSE)
  orderResponses(x, response)
}

# The rows of a result table grouped by response in the order of the model (stable: the
# order within a response is kept; rows of no response last).
orderByResponse <- function(x, df, col = "response", first = NULL, then = NULL) {
  if (!NROW(df)) return(df)
  out <- df[responseOrder(x, df, col, first, then), , drop = FALSE]
  rownames(out) <- NULL
  out
}
# The permutation behind orderByResponse(): a stable sort by `first` (a key sorted before the
# response, e.g. residual rows last), the response in the model's order (unknown responses last)
# and `then` (a key within a response).
responseOrder <- function(x, df, col = "response", first = NULL, then = NULL) {
  keys <- list(first, match(df[[col]], modelResponses(x)), then, seq_len(nrow(df)))
  do.call(order, c(Filter(Negate(is.null), keys), list(na.last = TRUE)))
}

#' @rdname mvremModel
#' @export
modelData <- function(x, ...) UseMethod("modelData")
#' @export
modelData.mvremModel <- function(x, ...) x$data

#' @rdname mvremModel
#' @export
dataSummary <- function(x, ...) UseMethod("dataSummary")
#' @export
dataSummary.mvremModel <- function(x, ...) x$dataSummary

#' @rdname mvremModel
#' @export
inlaFormula <- function(x, ...) UseMethod("inlaFormula")
#' @export
inlaFormula.mvremModel <- function(x, ...) x$inla$formula

#' @rdname mvremModel
#' @export
inlaSpec <- function(x, ...) UseMethod("inlaSpec")
#' @export
inlaSpec.mvremModel <- function(x, ...) x$inla

#' @rdname mvremModel
#' @export
inlaData <- function(x, ...) UseMethod("inlaData")
#' @rdname mvremModel
#' @export
inlaData.mvremModel <- function(x, n = Inf, ...) inlaData.default(x$inla, n)
# a buildInlaFormula() result
#' @export
inlaData.default <- function(x, n = Inf, ...) {
  if (!isNumber(n) || n < 0) stop("`n` must be a single number >= 0.", call. = FALSE)
  if (!is.list(x) || is.null(x$stacks))
    stop("`x` must be a model or the result of buildInlaFormula().", call. = FALSE)
  stackRows(x, n)
}

#' @rdname mvremModel
#' @export
inlaArguments <- function(x, ...) UseMethod("inlaArguments")
#' @export
inlaArguments.mvremModel <- function(x, ...) inlaArgumentList(x$inla, x$options)

#' @rdname mvremModel
#' @export
inlaFit <- function(x, ...) UseMethod("inlaFit")
#' @export
inlaFit.mvremModel <- function(x, ...) x$fit

#' @rdname mvremModel
#' @export
fixedEffects <- function(x, ...) UseMethod("fixedEffects")
#' @export
fixedEffects.mvremModel <- function(x, ...) fixedEffectsTable(x)

#' @rdname mvremModel
#' @export
hyperparameters <- function(x, ...) UseMethod("hyperparameters")
#' @export
hyperparameters.mvremModel <- function(x, ...) {
  h <- orderByResponse(x, hyperTable(x))
  h[setdiff(names(h), hyperInternalColumns)]
}

#' @rdname mvremModel
#' @export
associations <- function(x, ...) UseMethod("associations")
#' @export
associations.mvremModel <- function(x, ...) {
  h <- associationRows(x)
  out <- cbind(outcome = h$response, longitudinal = h$link,
               h[setdiff(names(h), c("response", "parameter", "type", hyperInternalColumns))])
  # E38: with latent equations, which kind of link each row is
  latent <- names(latentParts(x$forms))
  if (length(latent) && nrow(out))
    out <- cbind(out[1:2], level = ifelse(out$outcome %in% latent, "marker -> latent",
                                          ifelse(out$longitudinal %in% latent, "latent -> outcome",
                                                 "marker -> outcome")),
                 out[-(1:2)])
  rownames(out) <- NULL
  out
}

# The association rows of hyperTable(), outcomes in the order of the model (each outcome's
# scaled() terms in the order written).
associationRows <- function(x) {
  h <- hyperTable(x)
  orderByResponse(x, h[!is.na(h$type) & h$type == "association", , drop = FALSE])
}

#' @rdname mvremModel
#' @param verbose,threads,inlaOptions Options of INLA's run for `refit()` (see
#'   [jsm()]); default: those of the model.
#' @export
refit <- function(x, ...) UseMethod("refit")
#' @rdname mvremModel
#' @export
refit.mvremModel <- function(x, ..., verbose = x$options$verbose, threads = x$options$threads,
                             inlaOptions = x$options$inlaOptions) {
  dots <- list(...)
  if (length(dots)) stop("Unknown argument(s) to refit(): ", paste(names(dots), collapse = ", "),
                         call. = FALSE)
  checkFlag(verbose, "verbose")
  checkInlaOptions(inlaOptions, threads)
  x$options <- list(verbose = verbose, threads = threads, inlaOptions = inlaOptions)
  res <- runInla(x)
  if (inherits(res, "error")) {
    x["fit"] <- list(NULL)          # x$fit <- NULL would remove the component from the list
    x$inlaError <- res
    warning("INLA failed: ", firstLine(res),
            "\nThe model was built but not fitted; change the model or the options and use ",
            "refit().", call. = FALSE)
  } else {
    x$fit <- res
    x["inlaError"] <- list(NULL)
  }
  x
}

#' @rdname mvremModel
#' @param object An `mvremModel`.
#' @param evaluate If `TRUE` (default) the changed call is evaluated, otherwise
#'   it is returned.
#' @details `update()` repeats the call that created the model with the given
#'   arguments replaced (e.g. `update(m, family = list(y = gammaFamily()))`).
#'   Giving `form` or `formOut` replaces all formulae: the stored formula
#'   arguments (`form`, `formOut`, `formLong`) are dropped first, so
#'   `update(m, form = list(...))` switches a pair-style model to one list of
#'   formulae and `update(m, formOut = ..., formLong = ...)` the other way;
#'   `formLong` alone replaces only the longitudinal formulae of a pair.
#' @export
update.mvremModel <- function(object, ..., evaluate = TRUE) {
  call <- object$call
  if (is.null(call)) stop("The model has no stored call.", call. = FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  if (any(c("form", "formOut") %in% names(extras)))
    for (a in intersect(c("form", "formOut", "formLong"), names(call))) call[[a]] <- NULL
  for (a in names(extras)) call[[a]] <- extras[[a]]
  if (evaluate) eval(call, parent.frame()) else call
}

#' @export
coef.mvremModel <- function(object, ...) {
  fe <- fixedEffects(object)
  stats::setNames(fe$mean, paste0(fe$response, ": ", fe$term))
}

# Argument checks: one non-missing string / number.
isString <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
isNumber <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
isFlag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
# rbind() of a list of data frames (rows collected in a loop), or `template` (the empty data
# frame with the typed columns) when there are none; row names dropped.
bindRows <- function(rows, template) {
  out <- if (length(rows)) do.call(rbind, rows) else template
  rownames(out) <- NULL
  out
}
# A list whose elements all have non-empty names (an empty list passes); `unique`: no repeats.
isNamedList <- function(x, unique = FALSE) {
  if (!is.list(x)) return(FALSE)
  nms <- names(x)
  if (!length(x)) return(TRUE)
  !is.null(nms) && !any(nms == "") && (!unique || !anyDuplicated(nms))
}

# `arg` must be one whole number >= min (a count such as `n`, `ncol`, `nsamples`); NULL
# passes when `nullOk`. Returns the integer value invisibly.
checkCount <- function(x, arg, min = 0L, nullOk = FALSE) {
  if (is.null(x) && nullOk) return(invisible(NULL))
  if (!isNumber(x) || x != round(x) || x < min)
    stop("`", arg, "` must be a single whole number >= ", min, ".", call. = FALSE)
  invisible(as.integer(x))
}
checkFlag <- function(x, arg) {
  if (!isFlag(x)) stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  invisible(x)
}

# The outcome and longitudinal parts of a two-part list (`forms`, `data`, a model's
# `x$forms` / `x$data`) as one list named by response: outcomes first, then the markers.
# This is the internal (structural) order; results follow modelResponses() instead.
bothParts <- function(z) c(z$out, z$long)

# ---------------------------------------------------------------------------
# printing

#' @export
print.mvremModel <- function(x, ...) {
  inl <- x$inla
  cat("Multivariate joint model (mvrem)\n\n")
  famText <- function(p) paste0(p$response, " (", p$family$family, ")")
  cat("Outcomes:      ", paste(vapply(x$forms$out, famText, ""), collapse = ", "), "\n", sep = "")
  cat("Longitudinal:  ", paste(vapply(x$forms$long, famText, ""), collapse = ", "), "\n", sep = "")
  latent <- latentParts(x$forms)
  if (length(latent))
    cat("Latent:        ", paste(names(latent), collapse = ", "), " (no observations)\n", sep = "")
  links <- vapply(c(x$forms$out, latent), function(o) {
    t <- scaledTargets(o)
    if (length(t)) paste0(o$response, " <- ", paste(t, collapse = ", ")) else NA_character_
  }, "")
  cat("Associations:  ", paste(links[!is.na(links)], collapse = "; "), "\n", sep = "")
  # Q3 (E38): the scale and sign of a latent process are fixed only by the priors unless a
  # coefficient of its links is fixed
  free <- names(latent)[vapply(names(latent), function(L) !latentScaleFixed(x$forms, L), NA)]
  if (length(free))
    cat("               (scale and sign of ", paste(free, collapse = ", "), " set only by the priors; ",
        "fix one with scaled(", free[1L], ", coef = 1))\n", sep = "")
  tv <- timeVariable(x)
  if (!is.na(tv))
    cat("Time variable: ", tv, if (identical(attr(tv, "source"), "detected"))
      " (detected; set with jsm(time = ))", "\n", sep = "")
  eq <- orderByResponse(x, x$dataSummary$equations)   # (older objects: summary long first)
  cat("Observations:  ", paste0(eq$response, ": ", eq$complete, collapse = ", "), "\n", sep = "")
  used <- vapply(eq$response, function(r) nrow(bothParts(x$data)[[r]]), 0L)
  inc <- eq$rows - eq$complete                    # rows with a missing value
  note <- function(v, what) if (any(v > 0))
    cat("               (", what, ": ", paste0(eq$response[v > 0], ": ", v[v > 0], collapse = ", "),
        ")\n", sep = "")
  note(inc * (used > eq$complete), "rows with missing values, kept")
  note(eq$rows - used, "rows dropped for missing values")
  if (!is.null(x$fit)) {
    cat("\nFitted with INLA. Use summary(x), fixedEffects(x), hyperparameters(x), ",
        "associations(x).\n", sep = "")
  } else if (!is.null(x$inlaError)) {
    cat("\nNot fitted: INLA failed (", firstLine(x$inlaError), "). Use refit(x).\n", sep = "")
  } else {
    cat("\nNot fitted (use refit(x)). Available: modelFormulas(x), modelData(x), ",
        "inlaFormula(x), inlaData(x), latexEquation(x).\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.mvremModel <- function(object, ...) {
  needFit(object)
  structure(list(call = object$call, model = object,
                 fixed = fixedEffects(object), hyper = hyperparameters(object),
                 associations = associations(object)),
            class = "summary.mvremModel")
}

#' @export
print.summary.mvremModel <- function(x, digits = 3, ...) {
  show <- function(df) {
    num <- vapply(df, is.numeric, NA)
    df[num] <- lapply(df[num], function(v)
      format(signif(v, digits), scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
    print(df, row.names = FALSE, right = TRUE)
  }
  cat("Call:\n", paste(deparse(x$call, width.cutoff = 100), collapse = "\n"), "\n\n", sep = "")
  cat("Fixed effects (posterior):\n")
  show(x$fixed[setdiff(names(x$fixed), "mode")])
  cat("\nAssociations (gamma, posterior):\n")
  show(x$associations[setdiff(names(x$associations), "mode")])
  other <- x$hyper[is.na(x$hyper$type) | x$hyper$type != "association", ]
  cat("\nHyperparameters (posterior; precisions):\n")
  show(other[setdiff(names(other), c("type", "mode"))])
  invisible(x)
}
