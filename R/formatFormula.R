#' Show the INLA formula annotated with its data
#'
#' Formats the INLA formula built by [buildInlaFormula()] one term per line,
#' grouped by response. Each group is preceded by a comment with the user's
#' formula and the first rows of the corresponding stacked data (observation
#' rows and, for longitudinal responses, zero-equation rows).
#'
#' @param inl Result of [buildInlaFormula()].
#' @param nrows Number of rows shown per part of the stacked data (`0` hides
#'   the data).
#' @param digits Number of significant digits of numbers in the data.
#' @return A character vector of class `mvremFormulaText` (one element per
#'   line), printed with `cat()`.
#' @export
formatInlaFormula <- function(inl, nrows = 4L, digits = 3L) {
  stopifnot(length(nrows) == 1L, nrows >= 0)
  out <- paste(inl$naming$responseObject, "~ -1 +")
  resp <- names(inl$terms)
  for (i in seq_along(resp)) {
    r <- resp[i]
    out <- c(out, "   #",
             paste0("   # ----- ", r, ": ", deparse1(inl$formulas[[r]]), " -----"),
             "   #")
    if (nrows > 0) {
      parts <- lapply(inl$stacks[[r]], utils::head, n = nrows)
      tab <- do.call(rbind, unname(parts))
      shown <- utils::capture.output(print(tab, row.names = FALSE, digits = digits))
      out <- c(out, paste0("   # ", shown), "   #")
    }
    terms <- vapply(inl$terms[[r]], deparse1, "", width.cutoff = 500L)
    last <- i == length(resp)
    out <- c(out, paste0(terms, ifelse(last & seq_along(terms) == length(terms), "", " +")))
  }
  structure(out, class = "mvremFormulaText")
}

#' @export
print.mvremFormulaText <- function(x, ...) {
  cat(x, sep = "\n")
  invisible(x)
}
