#' LaTeX equations of a joint model
#'
#' Produces the LaTeX equations of formulae processed by [collectFormulas()]
#' (or [parseJsmFormula()]), optionally using the data prepared by
#' [prepareData()] to expand factor covariates into indicator terms.
#'
#' Observations are indexed by \eqn{i}; random effects by the level of their
#' grouping factor at observation \eqn{i}, e.g. `id[i]`. Fixed effects are
#' \eqn{\beta^{(r)}_k} (\eqn{k = 0} for the intercept), random effects
#' \eqn{b^{(r)}_k} share the index of the fixed effect of the same column (other
#' columns continue the numbering). Longitudinal responses get a residual
#' \eqn{\varepsilon^{(r)}_i} (Gaussian family). Each family (see [families])
#' writes its link of the mean, e.g. \eqn{\mathrm{logit}(\mu^{(r)}_i)}, and a
#' distribution line, e.g. \eqn{y_i \sim \mathrm{Beta}(\mu_i \phi, (1 - \mu_i)\phi)},
#' followed by the priors of its hyperparameters given by prior constructors;
#' [inlaFamily()] and a standalone outcome formula without family use a generic
#' link \eqn{g_r(\mathrm{E}[\cdot])}. In a whole model every outcome \eqn{y} is
#' associated with every longitudinal response \eqn{x} (with or without a
#' `scaled(x)` term); the association term is \eqn{\gamma_{y,x}\,u^{(x \to y)}_{\mathrm{id}[i],t_i}}:
#' the outcome is matched to the longitudinal submodel by the values of the
#' grouping variables and covariates of \eqn{x} (here `id`, `t`), not by row.
#' \eqn{u^{(x \to y)}_{\mathrm{id},t}} is defined by the zero equation
#' \eqn{0 = \eta^{(x)}(\mathrm{id}, t) - u^{(x \to y)}_{\mathrm{id},t}} for each
#' combination \eqn{(\mathrm{id}, t) \in \mathcal{D}_y} present in the outcome
#' data: the linear predictor of \eqn{x} (fixed and random effects). Zero
#' equations are only produced for a whole model, ordered as: longitudinal
#' equations, zero equations, outcome equations, then the distributions
#' (family, residuals, random effects; longitudinal before outcome) and
#' finally those of \eqn{u}. Normal distributions are parametrised by
#' precision, as in INLA: \eqn{\mathcal{N}(0, \tau^{-1})} and
#' \eqn{\mathcal{N}(\mathbf{0}, \mathbf{Q}^{-1})}.
#' A `prior` of a random-effect block (see [priors]) adds a line after the
#' block's distribution, e.g. \eqn{(\tau)^{-1/2} \sim \mathrm{Exp}(\lambda)} or
#' \eqn{\mathbf{Q} \sim \mathcal{W}_k(df, \mathbf{V})}; blocks without
#' `prior` (INLA's default) get none.
#'
#' @param x A parsed formula (element of `collectFormulas()$out` / `$long`) or
#'   the result of [collectFormulas()].
#' @param data `NULL`, or the data of step 2: a data frame (element of
#'   `prepareData()$out` / `$long`) for a single formula, the result of
#'   [prepareData()] for a whole model.
#' @return A character string of class `mvremLatex`: the body of a LaTeX
#'   `aligned` environment (lines separated by `\\`). Use [format()] with
#'   `delim = "$$"` to wrap it for display.
#' @export
latexEquation <- function(x, data = NULL) {
  if (all(c("response", "part") %in% names(x))) {
    if (!is.null(data) && !is.data.frame(data))
      stop("For a single formula, `data` must be a data frame.", call. = FALSE)
    eq <- latexFormula(x, data)
    lines <- c(eq$eq, eq$dist)
  } else if (all(c("out", "long") %in% names(x))) {
    if (!is.null(data) && !all(c("out", "long") %in% names(data)))
      stop("For a whole model, `data` must be the result of prepareData().", call. = FALSE)
    parsed <- c(x$out, x$long)
    dfs <- if (is.null(data)) list() else c(data$out, data$long)
    eq <- lapply(parsed, function(p) latexFormula(p, dfs[[p$response]], longs = x$long))
    zero <- character()
    uDist <- character()
    for (o in x$out) for (s in o$scaled) {
      l <- x$long[[s$target]]
      u <- uSymbol(l, o$response, obs = NULL)
      eta <- latexFormula(l, dfs[[l$response]], obs = NULL)$eta
      vars <- matchVarsLatex(l, obs = NULL)
      zero <- c(zero, paste0("0 &= ", paste(c(eta, paste0("- ", u)), collapse = " + "),
                             ", \\quad \\left(", paste(vars, collapse = ", "),
                             "\\right) \\in \\mathcal{D}_{", latexName(o$response), "}"))
      uDist <- c(uDist, paste0(u, " &\\sim \\mathcal{N}\\left(0, \\left(\\tau_{",
                               uTag(s$target, o$response), "}\\right)^{-1}\\right)"))
    }
    zero <- gsub("+ - ", "- ", zero, fixed = TRUE)
    eqLong <- eq[names(x$long)]
    eqOut <- eq[names(x$out)]
    lines <- c(unlist(lapply(eqLong, `[[`, "eq"), use.names = FALSE), zero,
               unlist(lapply(eqOut, `[[`, "eq"), use.names = FALSE),
               unlist(lapply(c(eqLong, eqOut), `[[`, "dist"), use.names = FALSE), uDist)
  } else {
    stop("`x` must be a parsed formula or the result of collectFormulas().", call. = FALSE)
  }
  structure(paste0("\\begin{aligned}\n", paste(lines, collapse = " \\\\\n"),
                   "\n\\end{aligned}"),
            class = "mvremLatex")
}

#' @export
format.mvremLatex <- function(x, delim = "", ...) {
  paste0(delim, unclass(x), delim)
}

#' @export
print.mvremLatex <- function(x, ...) {
  cat(format(x, ...), "\n", sep = "")
  invisible(x)
}

uTag <- function(long, out) paste0(latexName(long), " \\to ", latexName(out))
# Variables matching the outcome data to the zero equation of long formula l
# (its grouping variables and covariates), at observation obs (NULL: no index).
matchVarsLatex <- function(l, obs = "i") {
  vapply(c(l$variables$group, l$variables$covariate), function(v) {
    if (v %in% l$variables$group) {
      paste0("\\mathrm{", escapeLatex(v), "}", if (!is.null(obs)) paste0("[", obs, "]"))
    } else {
      latexName(v, obs)
    }
  }, "", USE.NAMES = FALSE)
}

# u of a (long, out) pair, indexed by the values of the matching variables;
# l = parsed long formula, or its name when unknown (then indexed by obs).
uSymbol <- function(l, out, obs = "i") {
  if (is.character(l)) return(paste0("u^{(", uTag(l, out), ")}_{", obs, "}"))
  paste0("u^{(", uTag(l$response, out), ")}_{", paste(matchVarsLatex(l, obs), collapse = ","), "}")
}

# One formula: `eq` (equation line), `dist` (family / residual / RE
# distribution lines) and `eta` (fixed and random terms of the linear predictor).
# obs: observation index ("i"), or NULL to write variables without index (used
# for zero equations, which are indexed by covariate values); longs: parsed
# long formulae (to index u in association terms by the matching variables).
latexFormula <- function(p, df, obs = "i", longs = NULL) {
  r <- p$response
  tag <- paste0("(", latexName(r), ")")
  fixedCols <- formulaColumns(p$fixedFormula, df)
  coefIdx <- seq_along(fixedCols) - as.integer("(Intercept)" %in% fixedCols)
  names(coefIdx) <- fixedCols
  nextIdx <- if (length(coefIdx)) max(coefIdx) + 1L else 0L

  terms <- vapply(seq_along(fixedCols), function(k)
    coefTerm(paste0("\\beta^{", tag, "}_{", coefIdx[k], "}"), fixedCols[k], df, obs), "")

  dist <- character()
  for (blk in p$re) {
    lhs <- stats::as.formula(call("~", blk$term[[2L]]), env = environment(p$formula))
    effects <- formulaColumns(lhs, df)
    idx <- integer(length(effects))
    for (e in seq_along(effects)) {
      if (effects[e] %in% names(coefIdx)) {
        idx[e] <- coefIdx[[effects[e]]]
      } else {
        idx[e] <- nextIdx
        coefIdx[[effects[e]]] <- nextIdx
        nextIdx <- nextIdx + 1L
      }
    }
    g <- paste0("\\mathrm{", escapeLatex(blk$group), "}")
    terms <- c(terms, vapply(seq_along(effects), function(e)
      coefTerm(paste0("b^{", tag, "}_{", idx[e], ",", g,
                      if (!is.null(obs)) paste0("[", obs, "]"), "}"), effects[e], df, obs), ""))
    bl <- paste0("b^{", tag, "}_{", idx, ",l}")
    if (blk$correlated && length(effects) > 1L) {
      dist <- c(dist, paste0("\\left(", paste(bl, collapse = ", "), "\\right)^\\top",
                             " &\\sim \\mathcal{N}\\left(\\mathbf{0}, \\left(\\mathbf{Q}^{", tag, "}_{", g,
                             "}\\right)^{-1}\\right), \\quad l = 1, \\ldots, n_{", g, "}"))
    } else {
      dist <- c(dist, paste0(bl, " &\\sim \\mathcal{N}\\left(0, \\left(\\tau^{", tag, "}_{",
                             idx, ",", g, "}\\right)^{-1}\\right), \\quad l = 1, \\ldots, n_{",
                             g, "}"))
    }
    if (!is.null(df)) checkBlockPrior(blk$prior, length(effects), blk$correlated, blk$term)
    tauIdx <- if (length(idx) == 1L) idx else "k"
    dist <- c(dist, priorLatex(blk$prior, length(effects),
                               tau = paste0("\\tau^{", tag, "}_{", tauIdx, ",", g, "}"),
                               Q = paste0("\\mathbf{Q}^{", tag, "}_{", g, "}"), idx = idx))
  }

  eta <- terms
  for (s in p$scaled) {
    l <- if (!is.null(longs[[s$target]])) longs[[s$target]] else s$target
    terms <- c(terms, paste0("\\gamma_{", latexName(r), ",", latexName(s$target), "}\\,",
                             uSymbol(l, r, obs)))
  }

  # family (see familyLatex()); unknown (standalone outcome) -> generic link
  fam <- p$family
  if (is.null(fam) && p$part == "long") fam <- defaultFamily()
  y <- exprLatex(p$responseExpr)
  if (is.null(fam)) {
    lhs <- paste0("g_{", latexName(r), "}\\left(\\mathrm{E}\\left[", y, "\\right]\\right)")
  } else {
    fl <- familyLatex(fam, y, tag, latexName(r), p$aux)
    lhs <- fl$lhs
    terms <- c(terms, fl$residual)
    dist <- c(fl$dist, dist)
  }
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "0"
  list(eq = paste0(lhs, " &= ", rhs), dist = dist, eta = eta)
}

# Columns of the design matrix of the right-hand side of f: from the data if
# given (factor levels expanded), otherwise the term labels.
formulaColumns <- function(f, df) {
  tt <- stats::delete.response(stats::terms(f))
  if (!is.null(df)) return(colnames(stats::model.matrix(tt, df)))
  c(if (attr(tt, "intercept") == 1L) "(Intercept)", attr(tt, "term.labels"))
}

# coefficient times covariate column (intercept: coefficient only)
coefTerm <- function(coef, col, df, obs = "i") {
  if (col == "(Intercept)") return(coef)
  paste0(coef, "\\, ", columnLatex(col, df, obs))
}

# A design-matrix column, e.g. "I(t^2)", "gv" (factor g, level v), "gv:t".
columnLatex <- function(col, df, obs = "i") {
  pieces <- strsplit(col, ":", fixed = TRUE)[[1L]]
  out <- vapply(pieces, function(piece) {
    if (!is.null(df)) for (v in names(df)) {
      if ((is.factor(df[[v]]) || is.character(df[[v]])) && startsWith(piece, v)) {
        lev <- substring(piece, nchar(v) + 1L)
        if (lev %in% as.character(unique(df[[v]])))
          return(paste0("\\mathbb{I}\\left[", latexName(v, obs), " = \\text{",
                        escapeLatex(lev), "}\\right]"))
      }
    }
    e <- tryCatch(str2lang(piece), error = function(err) NULL)
    if (is.null(e)) paste0("\\mathrm{", escapeLatex(piece), "}") else exprLatex(e, obs)
  }, "")
  paste(out, collapse = " \\cdot ")
}

escapeLatex <- function(x) gsub("([_&%$#{}])", "\\\\\\1", x)

# Variable name: "t" -> t_{i}, "x1" -> x_{1,i}, "age" -> \mathit{age}_{i}.
latexName <- function(name, idx = NULL) {
  if (grepl("^[A-Za-z]$", name)) {
    base <- name
    sub <- idx
  } else if (grepl("^[A-Za-z][0-9]+$", name)) {
    base <- substr(name, 1L, 1L)
    sub <- paste(c(substring(name, 2L), idx), collapse = ",")
  } else {
    base <- paste0("\\mathit{", escapeLatex(name), "}")
    sub <- idx
  }
  if (length(sub)) paste0(base, "_{", sub, "}") else base
}

# R expression -> LaTeX; variables get the observation index obs (NULL: none).
exprLatex <- function(e, obs = "i") {
  if (is.name(e)) return(latexName(as.character(e), obs))
  if (!is.call(e)) return(format(e))
  fn <- as.character(e[[1L]])
  a <- as.list(e)[-1L]
  arg <- function(k) exprLatex(a[[k]], obs)
  paren <- function(x) paste0("\\left(", x, "\\right)")
  switch(fn,
    "(" = paren(arg(1L)),
    "I" = arg(1L),
    "+" = if (length(a) == 1L) arg(1L) else paste(arg(1L), "+", arg(2L)),
    "-" = if (length(a) == 1L) paste0("-", arg(1L)) else paste(arg(1L), "-", arg(2L)),
    "*" = , ":" = paste(arg(1L), "\\cdot", arg(2L)),
    "/" = paste0("\\frac{", arg(1L), "}{", arg(2L), "}"),
    "^" = paste0("{", arg(1L), "}^{", arg(2L), "}"),
    "sqrt" = paste0("\\sqrt{", arg(1L), "}"),
    "abs" = paste0("\\left|", arg(1L), "\\right|"),
    "log" = , "exp" = , "sin" = , "cos" = , "tan" = paste0("\\", fn, paren(arg(1L))),
    paste0("\\operatorname{", escapeLatex(fn), "}",
           paren(paste(vapply(seq_along(a), arg, ""), collapse = ", ")))
  )
}
