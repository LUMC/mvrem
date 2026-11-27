#' All estimates of a fitted model in one table
#'
#' `estimatesTable()` joins the fixed effects ([fixedEffects()]), the association
#' coefficients ([associations()]), the standard deviations and correlations of
#' the random effects and the residual standard deviations
#' ([varianceComponents()]) and the other hyperparameters of the likelihoods
#' ([hyperparameters()]) into one data frame, one row per parameter, labelled
#' with the names of the model. It is the table of a paper's results section;
#' `format()` writes it as Markdown or LaTeX.
#'
#' @param x A fitted `mvremModel`, e.g. from [jsm()].
#' @param response Names of the responses to include (default: all).
#' @param style For `format()`: `"text"` (default; aligned columns, as printed),
#'   `"markdown"` (a pipe table), `"latex"` (a `tabular` environment for a LaTeX
#'   document) or `"array"` (the same table as a LaTeX `array`, to be placed in
#'   display math (`$$ ... $$`) of an HTML document: KaTeX and MathJax render it,
#'   `tabular` they do not).
#' @param digits Significant digits in `format()` / `print()`.
#' @param ... Ignored.
#' @return A data frame of class `mvremEstimates` with `response`, `type`
#'   (`"fixed"`, `"association"`, `"sd"`, `"correlation"`, `"residual"`,
#'   `"family"`, or `"other"` for hyperparameters mvrem does not recognise),
#'   `parameter` (e.g. `(Intercept)`, `t`, `gamma (x1 -> y)`, `SD t | id`, `corr (Intercept), t | id`,
#'   `residual SD`, `phi`), and the posterior `mean`, `sd`, `q0.025`, `q0.5`,
#'   `q0.975`; rows grouped by response (in the order the formulae were passed to [jsm()]). `format()` returns a
#'   character string (lines separated by newlines) with `mean (sd)` and the
#'   95% interval per row.
#' @examples
#' \dontrun{
#' estimatesTable(m)
#' cat(format(estimatesTable(m), "markdown"))
#' cat(format(estimatesTable(m, "y"), "latex", digits = 2))
#' cat("$$", format(estimatesTable(m), "array"), "$$", sep = "\n")   # in an R Markdown chunk
#' }
#' @export
estimatesTable <- function(x, ...) UseMethod("estimatesTable")

#' @rdname estimatesTable
#' @export
estimatesTable.mvremModel <- function(x, response = NULL, ...) {
  needFit(x)
  resp <- selectResponses(x, response)
  stats <- c("mean", "sd", "q0.025", "q0.5", "q0.975")
  rows <- function(r, type, parameter, df) {
    df <- as.data.frame(df)
    cbind(data.frame(response = as.character(r), type = rep(type, nrow(df)),
                     parameter = rep_len(as.character(parameter), nrow(df)), stringsAsFactors = FALSE),
          df[stats])
  }
  effLabel <- function(e) modelEffectLabel(x, e)
  out <- list()
  fe <- fixedEffects(x)
  out$fixed <- rows(fe$response, "fixed", effLabel(fe$term), fe)
  a <- associations(x)
  out$assoc <- rows(a$outcome, "association", paste0("gamma (", linkLabel(a$longitudinal, a$outcome), ")"), a)
  vc <- varianceComponentsOrNull(x, "the estimates table")$table
  if (!is.null(vc) && nrow(vc)) {
    sdr <- vc[vc$parameter == "sd" & vc$block != "residual", ]
    out$sd <- rows(sdr$response, "sd", paste0("SD ", effLabel(sdr$effect), " | ", sdr$group), sdr)
    cr <- vc[vc$parameter == "correlation", ]
    out$cor <- rows(cr$response, "correlation",
                    paste0("corr ", effLabel(cr$effect), ", ", effLabel(cr$effect2), " | ", cr$group), cr)
    res <- vc[vc$block == "residual", ]
    out$res <- rows(res$response, "residual", "residual SD", res)
  }
  h <- hyperTable(x)
  fam <- h[h$type %in% "family" & !h$short %in% "prec", ]
  out$fam <- rows(fam$response, "family", fam$parameter, fam)
  oth <- h[is.na(h$type), ]
  out$oth <- rows(rep("", nrow(oth)), "other", oth$parameter, oth)
  tab <- do.call(rbind, unname(out))
  tab <- tab[tab$response %in% c(resp, ""), , drop = FALSE]
  tab <- orderByResponse(x, tab, then = match(tab$type, c("fixed", "association", "sd", "correlation",
                                                           "residual", "family", "other")))
  structure(tab, class = c("mvremEstimates", "data.frame"))
}

#' @rdname estimatesTable
#' @export
format.mvremEstimates <- function(x, style = c("text", "markdown", "latex", "array"), digits = 3, ...) {
  style <- match.arg(style)
  num <- function(v) numText(v, digits)
  # a zero-row table gives zero rows (paste0() of empty vectors with constants has length 1)
  rowPaste <- function(...) if (nrow(x)) paste0(...) else character()
  est <- rowPaste(num(x$mean), " (", num(x$sd), ")")
  ci <- rowPaste("[", num(x$q0.025), ", ", num(x$q0.975), "]")
  resp <- ifelse(duplicated(x$response), "", x$response)
  if (style == "text") {
    if (!nrow(x)) return(" response type parameter mean (sd) 95% interval\n (no estimates)")
    df <- data.frame(response = resp, type = x$type, parameter = x$parameter,
                     `mean (sd)` = est, `95% interval` = ci, check.names = FALSE)
    return(paste(utils::capture.output(print(df, row.names = FALSE, right = FALSE)), collapse = "\n"))
  }
  if (style == "markdown") {
    md <- function(v) gsub("|", "\\|", v, fixed = TRUE)
    lines <- c("| response | parameter | mean (sd) | 95% interval |", "|---|---|---|---|",
               rowPaste("| ", md(resp), " | ", md(x$parameter), " | ", est, " | ", ci, " |"))
    return(paste(lines, collapse = "\n"))
  }
  tex <- function(v) gsub("->", "$\\\\to$", gsub("([_&%$#{}])", "\\\\\\1", v))
  if (style == "array") {       # the tabular as display math: KaTeX / MathJax render it
    txt <- function(v) ifelse(nzchar(v), paste0("\\text{", tex(v), "}"), "")
    lines <- c("\\begin{array}{llrr}", "\\hline",
               "\\text{response} & \\text{parameter} & \\text{mean (sd)} & \\text{95\\% interval} \\\\",
               "\\hline",
               rowPaste(txt(resp), " & ", txt(x$parameter), " & ", gsub(" ", "\\\\ ", est), " & ",
                        gsub(" ", "\\\\ ", ci), " \\\\"),
               "\\hline", "\\end{array}")
    return(paste(lines, collapse = "\n"))
  }
  lines <- c("\\begin{tabular}{llrr}", "\\hline",
             "response & parameter & mean (sd) & 95\\% interval \\\\", "\\hline",
             rowPaste(tex(resp), " & ", tex(x$parameter), " & ", est, " & ", ci, " \\\\"),
             "\\hline", "\\end{tabular}")
  paste(lines, collapse = "\n")
}

#' @rdname estimatesTable
#' @export
print.mvremEstimates <- function(x, digits = 3, ...) {
  cat(format(x, "text", digits = digits), "\n", sep = "")
  invisible(x)
}
