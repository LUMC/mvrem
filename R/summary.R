#' Summarize the data of a joint model
#'
#' Summarizes, for each equation of the model, how much of the supplied data
#' can be used: the number of rows supplied, the number of complete rows (no
#' missing value in the variables used by the formula; these are the rows kept
#' by [prepareData()] with `na.action = na.omit`), the number of groups per
#' grouping term among the complete rows, the missing values per variable and
#' the number of zero equations per (longitudinal, outcome) pair.
#'
#' @param forms Result of [collectFormulas()].
#' @param data,dataOut,dataLong Data as for [prepareData()]: `data` for all
#'   formulae (a data frame or a list of data frames named by response),
#'   `dataOut` / `dataLong` to override it for one part.
#' @return An object of class `mvremDataSummary`: a list with data frames
#'   * `equations`: one row per equation (in the order the formulae were passed)
#'     with `part`, `response`, `family`,
#'     `rows` (supplied), `complete`, `dropped` and `groups` (number of levels
#'     per grouping term among the complete rows, e.g. `"id: 50"`),
#'   * `missing`: one row per equation and variable used by it, with `part`,
#'     `response`, `variable` and `missing` (number of missing values among the
#'     supplied rows),
#'   * `zero`: one row per (longitudinal, outcome) pair with `long`, `out`, `vars`
#'     (matching variables) and `rows` (zero equations = unique combinations of
#'     the matching variables in the complete outcome rows).
#' @export
summarizeData <- function(forms, data = NULL, dataOut = NULL, dataLong = NULL) {
  raw <- list(out = resolveData(forms, "out", data, dataOut),
              long = resolveData(forms, "long", data, dataLong))
  # complete rows = what na.omit keeps; jsm() has prepared the data with its own na.action
  # (default na.pass), so the data are prepared (and validated) once more here
  prep <- prepareParts(forms, data, dataOut, dataLong, stats::na.omit, allowEmpty = TRUE)

  eqRows <- list()
  missRows <- list()
  # equations in the order the formulae were passed (modelResponses())
  for (r in selectResponses(forms, NULL, observed = TRUE)) {
    p <- bothParts(forms)[[r]]
    part <- p$part
    df <- raw[[part]][[r]]
    pd <- prep[[part]][[r]]
    groupTerms <- unique(reGroups(p))
    groups <- if (length(groupTerms)) {
      paste(paste0(groupTerms, ": ",
                   vapply(groupTerms, function(g) length(unique(pd[[g]])), 0L)),
            collapse = ", ")
    } else ""
    eqRows[[length(eqRows) + 1L]] <- data.frame(
      part = part, response = r, family = formulaFamily(p)$family,
      rows = nrow(df), complete = nrow(pd), dropped = nrow(df) - nrow(pd), groups = groups
    )
    vars <- formulaVariables(p)
    missRows[[length(missRows) + 1L]] <- data.frame(
      part = rep(part, length(vars)), response = rep(r, length(vars)), variable = vars,
      missing = vapply(vars, function(v) sum(is.na(df[[v]])), 0L, USE.NAMES = FALSE)
    )
  }

  zeroRows <- list()
  for (o in forms$out) for (s in o$scaled) {
    # a latent equation (E38) is evaluated at its variables and those of its targets
    l <- c(forms$long, latentParts(forms))[[s$target]]
    vars <- linkedVariables(l, forms$long)
    pd <- prep$out[[o$response]]
    combos <- unique(pd[, vars, drop = FALSE])
    zeroRows[[length(zeroRows) + 1L]] <- data.frame(
      long = l$response, out = o$response, vars = paste(vars, collapse = ", "),
      rows = nrow(combos)
    )
  }

  structure(list(
    equations = bindRows(eqRows, data.frame(part = character(), response = character(), family = character(),
                                            rows = integer(), complete = integer(), dropped = integer(),
                                            groups = character())),
    missing = bindRows(missRows, data.frame(part = character(), response = character(),
                                            variable = character(), missing = integer())),
    zero = bindRows(zeroRows, data.frame(long = character(), out = character(), vars = character(),
                                         rows = integer()))
  ), class = "mvremDataSummary")
}

#' @export
print.mvremDataSummary <- function(x, ...) {
  cat("Equations (rows supplied / complete / with missing values):\n")
  print(x$equations, row.names = FALSE)
  if (nrow(x$zero)) {
    cat("\nZero equations (unique combinations of matching variables in outcome rows):\n")
    print(x$zero, row.names = FALSE)
  }
  miss <- x$missing[x$missing$missing > 0, , drop = FALSE]
  cat("\nMissing values per variable:\n")
  if (nrow(miss)) print(miss, row.names = FALSE) else cat("  none\n")
  invisible(x)
}
