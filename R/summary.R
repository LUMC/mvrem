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
#' @param dataOut,dataLong Data as for [prepareData()]: one data frame per part
#'   or a list of data frames named by response.
#' @return An object of class `mvremDataSummary`: a list with data frames
#'   * `equations`: one row per equation with `part`, `response`, `family`,
#'     `rows` (supplied), `complete`, `dropped` and `groups` (number of levels
#'     per grouping term among the complete rows, e.g. `"id: 50"`),
#'   * `missing`: one row per equation and variable used by it, with `part`,
#'     `response`, `variable` and `missing` (number of missing values among the
#'     supplied rows),
#'   * `zero`: one row per (longitudinal, outcome) pair with `long`, `out`, `vars`
#'     (matching variables) and `rows` (zero equations = unique combinations of
#'     the matching variables in the complete outcome rows).
#' @export
summarizeData <- function(forms, dataOut, dataLong) {
  raw <- list(out = matchData(dataOut, forms$out, "dataOut"),
              long = matchData(dataLong, forms$long, "dataLong"))
  prep <- prepareData(dataOut, dataLong, forms, na.action = stats::na.omit)

  eqRows <- list()
  missRows <- list()
  for (part in c("long", "out")) for (p in forms[[part]]) {
    r <- p$response
    df <- raw[[part]][[r]]
    pd <- prep[[part]][[r]]
    groupTerms <- vapply(p$re, `[[`, "", "group")
    groupTerms <- unique(groupTerms)
    groups <- if (length(groupTerms)) {
      paste(paste0(groupTerms, ": ",
                   vapply(groupTerms, function(g) length(unique(pd[[g]])), 0L)),
            collapse = ", ")
    } else ""
    eqRows[[length(eqRows) + 1L]] <- data.frame(
      part = part, response = r, family = if (is.null(p$family)) defaultFamily()$family else p$family$family,
      rows = nrow(df), complete = nrow(pd), dropped = nrow(df) - nrow(pd), groups = groups
    )
    vars <- unique(c(unlist(p$variables, use.names = FALSE), p$linkVariables))
    missRows[[length(missRows) + 1L]] <- data.frame(
      part = rep(part, length(vars)), response = rep(r, length(vars)), variable = vars,
      missing = vapply(vars, function(v) sum(is.na(df[[v]])), 0L, USE.NAMES = FALSE)
    )
  }

  zeroRows <- list()
  for (o in forms$out) for (s in o$scaled) {
    l <- forms$long[[s$target]]
    vars <- c(l$variables$group, l$variables$covariate)
    pd <- prep$out[[o$response]]
    combos <- unique(pd[, vars, drop = FALSE])
    zeroRows[[length(zeroRows) + 1L]] <- data.frame(
      long = l$response, out = o$response, vars = paste(vars, collapse = ", "),
      rows = nrow(combos)
    )
  }

  bind <- function(x, cols) {
    if (length(x)) return(do.call(rbind, x))
    as.data.frame(stats::setNames(replicate(length(cols), character(), simplify = FALSE), cols))
  }
  structure(list(
    equations = bind(eqRows, c("part", "response", "family", "rows", "complete", "dropped",
                               "groups")),
    missing = bind(missRows, c("part", "response", "variable", "missing")),
    zero = bind(zeroRows, c("long", "out", "vars", "rows"))
  ), class = "mvremDataSummary")
}

#' @export
print.mvremDataSummary <- function(x, ...) {
  cat("Equations (rows supplied / complete / dropped because of missing values):\n")
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
