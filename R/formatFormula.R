#' Show the INLA formula annotated with its data
#'
#' Formats the INLA formula built by [buildInlaFormula()] one term per line,
#' grouped by response. Each group is preceded by a comment with the user's
#' formula, the likelihood family and `control.family` (also of the zero
#' equations of a longitudinal response), the family-input columns
#' (`# inputs: Ntrials = k_Ntrials`, `scale = x_scale` for `weights()`; not with
#' `nrows = 0`) and the first rows of the corresponding stacked data (observation
#' rows and, for longitudinal responses, zero-equation rows).
#'
#' @param inl Result of [buildInlaFormula()], or a model (class `mvremModel`,
#'   e.g. from [jsm()]).
#' @param ... Further arguments of the methods.
#' @param nrows Number of rows shown per part of the stacked data (`0` hides
#'   the data).
#' @param digits Number of significant digits of numbers in the data.
#' @return A character vector of class `mvremFormulaText` (one element per
#'   line), printed with `cat()`.
#' @export
formatInlaFormula <- function(inl, ...) UseMethod("formatInlaFormula")

#' @rdname formatInlaFormula
#' @export
formatInlaFormula.default <- function(inl, nrows = 4L, digits = 3L, ...) {
  stopifnot(length(nrows) == 1L, nrows >= 0)
  out <- paste(inl$naming$responseObject, "~ -1 +")
  resp <- names(inl$terms)
  for (i in seq_along(resp)) {
    r <- resp[i]
    isLatent <- r %in% names(inl$latent)       # E38: zero rows only
    out <- c(out, "   #",
             paste0("   # ----- ", r, ": ", deparse1(inl$formulas[[r]]), " -----"),
             if (!isLatent) familyComment("", inl, inl$responses[[r]]),
             if (nrows > 0) inputComment(inl$inputs[[r]]))
    if (r %in% names(inl$zeroResponses))
      out <- c(out, familyComment(if (isLatent) "latent equations: " else "zero equations: ",
                                  inl, inl$zeroResponses[[r]]))
    out <- c(out, "   #")
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

#' @rdname formatInlaFormula
#' @export
formatInlaFormula.mvremModel <- function(inl, ...) formatInlaFormula(inl$inla, ...)

familyComment <- function(prefix, inl, col) {
  ctrl <- inl$controlFamily[[col]]
  paste0("   # ", prefix, "family = \"", inl$families[[col]], "\"",
         if (length(ctrl)) paste0(", control.family = ", deparse1(ctrl, width.cutoff = 500L)))
}

# The family inputs of a response (INLA argument = stack column), e.g.
# `# inputs: scale = x_scale` for weights(); none -> no line. The `event` column is
# part of the inla.surv response, not an inla() argument.
inputComment <- function(inputs) {
  if (!length(inputs)) return(character())
  arg <- ifelse(names(inputs) == "event", "event (inla.surv)", names(inputs))
  note <- if (all(names(inputs) == "event")) "" else " (rows of other likelihoods: 1)"
  paste0("   # inputs: ", paste0(arg, " = ", inputs, collapse = ", "), note)
}

#' @export
print.mvremFormulaText <- function(x, ...) {
  cat(x, sep = "\n")
  invisible(x)
}
