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
#' link \eqn{g_r(\mathrm{E}[\cdot])}. In a whole model an outcome \eqn{y} is
#' associated with each longitudinal response \eqn{x} named in a `scaled(x)` term
#' of its formula; the association term is \eqn{\gamma_{y,x}\,u^{(x \to y)}_{\mathrm{id}[i],t_i}}:
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
#' @param x A parsed formula (element of `collectFormulas()$out` / `$long`),
#'   the result of [collectFormulas()], or a model (class `mvremModel`, e.g. from
#'   [jsm()]; its formulae and data are used).
#' @param data `NULL`, or the data of step 2: a data frame (element of
#'   `prepareData()$out` / `$long`) for a single formula, the result of
#'   [prepareData()] for a whole model. Not used for a model object.
#' @param annotate For a fitted model: `"none"` (default) writes the symbols
#'   only; `"tooltip"` annotates every symbol with its posterior estimate and
#'   further information for an HTML page rendered with KaTeX (hovering a symbol
#'   shows them; see Details); `"inline"` writes the posterior mean and, on a
#'   second line, the sd under each estimated parameter (`\\underset{\\substack{...}}`,
#'   valid in LaTeX documents with amsmath too).
#' @param digits Significant digits of the estimates.
#' @param ... Further arguments of the methods.
#' @details With `annotate = "tooltip"` the symbols are wrapped in KaTeX's
#'   `\htmlData{...}{}` (data attributes on the rendered span), which KaTeX
#'   only accepts with the render option `trust: true`; the package's KaTeX
#'   header (`vignettes/katex/katex-header.html` in the source repository)
#'   sets it and adds the tooltip script `vignettes/katex/mvrem-tips.js`. The tooltip of
#'   a parameter gives its posterior mean, sd and 95% interval (fixed effects,
#'   associations, family hyperparameters, precisions; a covariance matrix
#'   \eqn{\mathbf{Q}} its SDs and correlations); random effects and residuals
#'   give their SD; the response, covariates, fitted means, numbers of levels,
#'   the zero-equation field \eqn{u} and its fixed precision get descriptive
#'   notes. Such output is not valid LaTeX outside KaTeX.
#' @return A character string of class `mvremLatex`: the body of a LaTeX
#'   `aligned` environment (lines separated by `\\`). Use [format()] with
#'   `delim = "$$"` to wrap it for display.
#' @export
latexEquation <- function(x, ...) UseMethod("latexEquation")

#' @rdname latexEquation
#' @export
latexEquation.default <- function(x, data = NULL, ...) {
  if (all(c("response", "part") %in% names(x))) {
    if (!is.null(data) && !is.data.frame(data))
      stop("For a single formula, `data` must be a data frame.", call. = FALSE)
    eq <- latexFormula(x, data)
    lines <- c(eq$eq, eq$dist)
  } else if (all(c("out", "long") %in% names(x))) {
    if (!is.null(data) && !all(c("out", "long") %in% names(data)))
      stop("For a whole model, `data` must be the result of prepareData().", call. = FALSE)
    latent <- latentParts(x)                   # E38: latent equations
    targets <- c(x$long, latentLatexTargets(latent, x$long))
    parsed <- c(bothParts(x), latent)
    dfs <- if (is.null(data)) list() else bothParts(data)
    eq <- lapply(parsed, function(p) latexFormula(p, dfs[[p$response]], longs = targets,
                                                  obs = if (p$part == "latent") NULL else "i"))
    zero <- character()
    uDist <- character()
    for (o in c(x$out, latent)) for (s in o$scaled) {
      if (s$target %in% names(latent)) next      # the latent's own line defines its field
      l <- x$long[[s$target]]
      u <- uSymbol(l, o$response, obs = NULL)
      eta <- latexFormula(l, dfs[[l$response]], obs = NULL)$eta
      vars <- matchVarsLatex(l, obs = NULL)
      zero <- c(zero, paste0("0 &= ", paste(c(eta, paste0("- ", u)), collapse = " + "),
                             ", \\quad \\left(", paste(vars, collapse = ", "),
                             "\\right) \\in ", annot(paste0("\\mathcal{D}_{", latexName(o$response), "}"),
                                                     annotKey("D", o$response, s$target))))
      uDist <- c(uDist, paste0(u, " &\\sim ", normalSymbol(), "\\left(0, \\left(",
                               annot(paste0("\\tau_{", uTag(s$target, o$response), "}"),
                                     annotKey("utau", s$target, o$response)),
                               "\\right)^{-1}\\right)"))
    }
    for (L in names(latent)) {
      u <- uSymbol(targets[[L]], L, obs = NULL)
      uDist <- c(uDist, paste0(u, " &\\sim ", normalSymbol(), "\\left(0, \\left(\\tau_{",
                               latexName(L), "}\\right)^{-1}\\right)",
                               if (isTRUE(latent[[L]]$spec$center))
                                 paste0(", \\quad \\textstyle\\sum ", "u^{(", latexName(L), ")} = 0")))
    }
    zero <- gsub("+ - ", "- ", zero, fixed = TRUE)
    eqLong <- eq[names(x$long)]
    eqLatent <- eq[names(latent)]
    eqOut <- eq[names(x$out)]
    lines <- c(unlist(lapply(eqLong, `[[`, "eq"), use.names = FALSE), zero,
               unlist(lapply(eqLatent, `[[`, "eq"), use.names = FALSE),
               unlist(lapply(eqOut, `[[`, "eq"), use.names = FALSE),
               unlist(lapply(c(eqLong, eqLatent, eqOut), `[[`, "dist"), use.names = FALSE), uDist)
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

# The latent equations as targets of links in LaTeX (E38): the parsed formula with its
# variables extended by those of its own targets (its field is indexed by all of them) and
# `latentField = TRUE`, so that uSymbol() writes u^{(eta)} without an arrow.
latentLatexTargets <- function(latent, long) {
  lapply(latent, function(p) {
    tg <- long[scaledTargets(p)]
    p$variables$group <- unique(c(p$variables$group, unlist(lapply(tg, function(l) l$variables$group))))
    p$variables$covariate <- unique(c(p$variables$covariate,
                                      unlist(lapply(tg, function(l) l$variables$covariate))))
    p$latentField <- TRUE
    p
  })
}
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
  idx <- paste(matchVarsLatex(l, obs), collapse = ",")
  # the field of a latent equation (E38) is shared by its outcomes: u^{(eta)}
  if (isTRUE(l$latentField)) return(paste0("u^{(", latexName(l$response), ")}_{", idx, "}"))
  annot(paste0("u^{(", uTag(l$response, out), ")}_{", idx, "}"), annotKey("u", l$response, out))
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

  terms <- vapply(seq_along(fixedCols), function(k)
    coefTerm(annot(paste0("\\beta^{", tag, "}_{", coefIdx[k], "}"), annotKey("fixed", r, fixedCols[k])),
             fixedCols[k], df, obs), "")

  dist <- character()
  for (blk in p$re) {
    effects <- blockEffects(blk, p, df)
    coefIdx <- extendCoefIndex(coefIdx, effects)
    b <- reBlockLatex(blk, p, df, effects, unname(coefIdx[effects]), tag, obs)
    terms <- c(terms, b$terms)
    dist <- c(dist, b$dist)
  }

  eta <- terms
  for (s in p$scaled) {
    l <- if (!is.null(longs[[s$target]])) longs[[s$target]] else s$target
    gamma <- annot(paste0("\\gamma_{", latexName(r), ",", latexName(s$target), "}"),
                   annotKey("assoc", r, s$target))
    terms <- c(terms, paste0(gamma, "\\,", uSymbol(l, r, obs)))
    dist <- c(dist, coefPriorLatex(s$prior, gamma))
  }

  # a latent equation (E38) defines its field at its evaluation points (no family)
  if (p$part == "latent") {
    self <- if (!is.null(longs[[r]])) longs[[r]] else latentLatexTargets(list(p), list())[[1L]]
    rhs <- if (length(terms)) paste(terms, collapse = " + ") else "0"
    vars <- matchVarsLatex(self, obs = NULL)
    return(list(eq = paste0(uSymbol(self, r, obs = NULL), " &= ", rhs, ", \\quad \\left(",
                            paste(vars, collapse = ", "), "\\right) \\in \\mathcal{D}_{",
                            latexName(r), "}"),
                dist = dist, eta = eta))
  }

  # family (see familyLatex()); unknown (standalone outcome) -> generic link
  fam <- p$family
  if (is.null(fam) && p$part == "long") fam <- defaultFamily()
  y <- exprLatex(p$responseExpr)
  if (is.null(fam)) {
    lhs <- paste0("g_{", latexName(r), "}\\left(\\mathrm{E}\\left[", y, "\\right]\\right)")
  } else {
    fl <- familyLatex(fam, y, tag, latexName(r), p$aux, resp = r)
    lhs <- fl$lhs
    terms <- c(terms, fl$residual)
    dist <- c(fl$dist, dist)
  }
  rhs <- if (length(terms)) paste(terms, collapse = " + ") else "0"
  list(eq = paste0(lhs, " &= ", rhs), dist = dist, eta = eta)
}

# Left-hand side of a random-effect block as a one-sided formula (`~ 1 + t`), and the
# block's design columns (effects) on the data df (term labels without data).
blockLhs <- function(blk, p) stats::as.formula(call("~", blk$term[[2L]]), env = environment(p$formula))
blockEffects <- function(blk, p, df) formulaColumns(blockLhs(blk, p), df)

# Coefficient indices (named by design column) extended by the effects of a random-effect
# block: an effect with a fixed-effect column shares its index k, the others get the next
# free ones.
extendCoefIndex <- function(coefIdx, effects) {
  for (e in setdiff(effects, names(coefIdx)))
    coefIdx[[e]] <- if (length(coefIdx)) max(coefIdx) + 1L else 0L
  coefIdx
}

# One random-effect block of formula p: `terms` (b_{k,g[i]} x covariate, for the linear
# predictor) and `dist` (its distribution line, then its prior lines). idx = the effects'
# coefficient indices, tag = "(r)".
reBlockLatex <- function(blk, p, df, effects, idx, tag, obs) {
  r <- p$response
  g <- paste0("\\mathrm{", escapeLatex(blk$group), "}")
  bKey <- annotKey("b", r, blk$group, effects)
  terms <- vapply(seq_along(effects), function(e)
    coefTerm(annot(paste0("b^{", tag, "}_{", idx[e], ",", g,
                          if (!is.null(obs)) paste0("[", obs, "]"), "}"), bKey[e]), effects[e], df, obs), "")
  sym <- list(
    tau = function(k) paste0("\\tau^{", tag, "}_{", k, ",", g, "}"),   # precision of effect(s) k
    precKey = annotKey("prec", r, blk$group, effects),
    Q = annot(paste0("\\mathbf{Q}^{", tag, "}_{", g, "}"), blockKey(r, blk$group, effects)))
  bl <- annot(paste0("b^{", tag, "}_{", idx, ",l}"), bKey)
  nLev <- annot(paste0("n_{", g, "}"), annotKey("nlev", r, blk$group))
  correlated <- blk$correlated && length(effects) > 1L
  dist <- if (correlated) {
    paste0("\\left(", paste(bl, collapse = ", "), "\\right)^\\top",
           " &\\sim ", normalSymbol(), "\\left(\\mathbf{0}, \\left(", sym$Q,
           "\\right)^{-1}\\right), \\quad l = 1, \\ldots, ", nLev)
  } else {
    paste0(bl, " &\\sim ", normalSymbol(), "\\left(0, \\left(",
           annot(sym$tau(idx), sym$precKey), "\\right)^{-1}\\right), \\quad l = 1, \\ldots, ", nLev)
  }
  list(terms = terms, dist = c(dist, reBlockPriorLatex(blk, p, df, effects, idx, correlated, sym)))
}

# Prior lines of a random-effect block (none without a prior). With data the prior is
# checked against the block and its data-dependent parts are resolved; a data-scaled scalar
# prior of an uncorrelated block then gives one line per effect with its own numbers.
reBlockPriorLatex <- function(blk, p, df, effects, idx, correlated, sym) {
  prior <- blk$prior
  if (!is.null(df)) {
    checkBlockPrior(prior, length(effects), blk$correlated, blk$term)
    if (!correlated && isScalarDataPrior(prior)) {
      priors <- resolveEffectPriors(prior, p, df, effects, blk)
      return(unlist(lapply(seq_along(effects), function(e)
        priorLatex(priors[[e]], 1L, tau = annot(sym$tau(idx[e]), sym$precKey[e]), Q = NULL, idx = idx[e])),
        use.names = FALSE))
    }
    prior <- resolveBlockPrior(prior, p, df, effects, blk)
  }
  tau <- sym$tau(if (length(idx) == 1L) idx else "k")
  if (length(effects) == 1L) tau <- annot(tau, sym$precKey)
  priorLatex(prior, length(effects), tau = tau, Q = sym$Q, idx = idx)
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
          return(paste0("\\mathbb{I}\\left[", annot(latexName(v, obs), annotKey("cov", v)), " = \\text{",
                        escapeLatex(lev, "text"), "}\\right]"))
      }
    }
    e <- tryCatch(str2lang(piece), error = function(err) NULL)
    if (is.null(e)) paste0("\\mathrm{", escapeLatex(piece), "}") else exprLatex(e, obs)
  }, "")
  paste(out, collapse = " \\cdot ")
}

# Escapes every LaTeX special character of a raw string (never already-escaped text).
# `mode` is the mode of the enclosing command: "text" for \text{} (factor levels, inline
# items), "math" for \mathrm{} / \mathit{} / \operatorname{} (names), because KaTeX accepts
# \textasciicircum / \textasciitilde / \textbackslash only in text mode and \backslash only in math mode.
escapeLatex <- function(x, mode = c("math", "text")) {
  mode <- match.arg(mode)
  map <- c("_" = "\\_", "&" = "\\&", "%" = "\\%", "$" = "\\$", "#" = "\\#", "{" = "\\{", "}" = "\\}",
           if (mode == "text") c("^" = "\\textasciicircum{}", "~" = "\\textasciitilde{}", "\\" = "\\textbackslash{}")
           else c("^" = "\\hat{}", "~" = "\\tilde{}", "\\" = "\\backslash{}"))
  # Character by character (not nested gsub): the replacements contain { } and \, which a later
  # substitution would escape again.
  vapply(x, function(s) {
    if (is.na(s)) return(s)
    ch <- strsplit(s, "", fixed = TRUE)[[1L]]
    paste(ifelse(ch %in% names(map), map[ch], ch), collapse = "")
  }, "", USE.NAMES = FALSE)
}

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
  if (is.name(e)) return(annot(latexName(as.character(e), obs), annotKey("cov", as.character(e))))
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

#' @rdname latexEquation
#' @export
latexEquation.mvremModel <- function(x, annotate = c("none", "tooltip", "inline"), digits = 3, ...) {
  annotate <- match.arg(annotate)
  if (annotate == "none") return(latexEquation(x$forms, x$data))
  needFit(x)
  old <- latexState$est
  latexState$est <- list(mode = annotate, digits = digits, keys = estimateLookup(x))
  on.exit(latexState$est <- old, add = TRUE)
  latexEquation(x$forms, x$data)
}

# ---------------------------------------------------------------------------
# Estimates in the equations: while latexEquation(model, annotate = ...) runs,
# latexState$est holds the mode and a lookup (key -> what to show); the symbol
# builders call annot(symbol, key), which wraps the symbol when the key is known.
#
# Why a package-level environment rather than an argument: the symbols are built
# deep in a chain of helpers (latexFormula(), reBlockLatex(), familyLatex(),
# priorLatex(), exprLatex(), columnLatex(), ...) that also serve the plain
# latexEquation() of formulae without a model. Threading the lookup through all of
# them would add an argument to every helper for one optional feature; the state is
# set for one call by latexEquation.mvremModel() and restored on exit (also on error),
# and without it annot() returns the symbol unchanged.

latexState <- new.env(parent = emptyenv())

# Keys of the lookup, "<kind>|<part>|...": built by annotKey() on both sides, by
# estimateLookup() (producer) and by the symbol builders (consumers). Kinds and parts:
#   fixed|r|col  assoc|out|long  family|r|short  prec|r|g|effect  block|r|g|eff1|eff2..
#   b|r|g|effect  eps|r  resp|r  mean|r  cov|v  nlev|r|g  u|long|out  utau|long|out
#   D|out|long  dist|<name>
# Vectorised over the parts like paste(); a NULL part counts as "" (a standalone
# formula has no response name).
annotKey <- function(kind, ...) {
  parts <- lapply(list(...), function(x) if (is.null(x)) "" else x)
  do.call(paste, c(list(kind), parts, sep = "|"))
}
distKey <- function(name) annotKey("dist", name)
blockKey <- function(r, g, effects) annotKey("block", r, g, paste(effects, collapse = "|"))

# The band of an entry with a posterior credible interval [lo, hi]; other entries
# (observed / fitted values) have band "range". Only credible intervals are shown
# inline and their band is implied in tooltips.
credibleBand <- "95%"
isCredible <- function(e) identical(e$band, credibleBand)

# An entry of the lookup: par (title), item / mean / sd / lo / hi (parallel vectors,
# one line each in the tooltip), band (label of [lo, hi]: credibleBand or "range"),
# note (free text).
estEntry <- function(par, mean = NULL, sd = NULL, lo = NULL, hi = NULL, item = NULL,
                     band = credibleBand, note = NULL)
  list(par = par, mean = mean, sd = sd, lo = lo, hi = hi, item = item, band = band, note = note)

# A symbol annotated with the estimate stored under `key` in latexState$est (see
# estimateLookup()): vectorised over sym / key; without an active lookup, or an unknown key,
# the plain symbol; mode "tooltip" wraps it in \htmlData{}, "inline" (credible-band entries
# only) in \underset{}.
annot <- function(sym, key) {
  est <- latexState$est
  if (is.null(est)) return(sym)
  if (length(sym) > 1L) return(unname(mapply(annot, sym, key)))
  if (is.null(e <- est$keys[[key]])) return(sym)
  if (est$mode == "inline") {
    if (is.null(e$mean) || !isCredible(e)) return(sym)
    return(paste0("\\underset{\\substack{", inlineEstimate(e, est$digits), "}}{", sym, "}"))
  }
  paste0("\\htmlData{", htmlDataString(e, est$digits), "}{", sym, "}")
}

# Lines of a \substack under the symbol: the mean and "(sd)" on two lines ("SD 0.485" /
# "(0.027)" for a random effect); one line per item ("SD t 0.218 (0.0167)") for a covariance matrix.
inlineEstimate <- function(e, digits = 3) {
  sc <- function(x) paste0("\\scriptstyle ", x)
  num <- function(v) numLatex(v, digits)
  lines <- if (length(e$mean) == 1L) {
    c(sc(paste0(if (!is.null(e$item)) paste0("\\text{", escapeLatex(e$item, "text"), "}\\ "), num(e$mean))),
      if (!is.null(e$sd)) sc(paste0("(", num(e$sd), ")")))
  } else {
    vapply(seq_along(e$mean), function(i) sc(paste0("\\text{", escapeLatex(e$item[i], "text"), "}\\ ",
                                                   num(e$mean[i]),
                                                   if (!is.null(e$sd)) paste0("\\,(", num(e$sd[i]), ")"))), "")
  }
  paste(lines, collapse = " \\\\ ")
}

# key=value pairs for \htmlData: values may not contain commas, braces or "="
# (KaTeX splits on them); vectors are joined with ";" (split again by mvrem-tips.js).
htmlDataString <- function(e, digits) {
  txt <- function(v) {
    v <- gsub("[^][A-Za-z0-9 ():.^+*/|'<>-]", " ", v, perl = TRUE)
    trimws(gsub(" +", " ", v))
  }
  num <- function(v) numText(v, digits)
  pairs <- c(par = txt(e$par))
  if (!is.null(e$item)) pairs["item"] <- paste(txt(e$item), collapse = ";")
  for (f in c("mean", "sd", "lo", "hi")) if (!is.null(e[[f]])) pairs[f] <- paste(num(e[[f]]), collapse = ";")
  if (!is.null(e$lo) && !isCredible(e)) pairs["band"] <- txt(e$band)
  if (!is.null(e$note)) pairs["note"] <- txt(e$note)
  paste0(names(pairs), "=", pairs, collapse = ",")
}

# One line about the prior of a hyperparameter (tooltip text: no commas / "=" / ";").
priorNote <- function(prior) {
  if (is.null(prior)) return("prior: INLA default")
  if (!inherits(prior, "mvremPrior")) {
    if (is.list(prior)) return(paste0("prior: INLA ", if (!is.null(prior$prior)) prior$prior else "hyper list",
                                      if (!is.null(prior$param)) paste0(" with param ", paste(numText(prior$param), collapse = " ")),
                                      if (isTRUE(prior$fixed)) " (fixed)"))
    return(NULL)
  }
  a <- prior$args
  if (isScalarDataPrior(prior)) {
    r <- dataPriorRules
    return(switch(prior$dist,
      gamma = paste0("prior on the precision: Gamma with shape ", numText(a$shape),
                     " and mean 1/s^2 (s: the lme4 SD of the effect)"),
      pc = paste0("PC prior on the SD: P(SD > u) is ", numText(a$alpha), " with u ", numText(r$pcU),
                  " times the lme4 SD of the effect"),
      halfNormal = paste0("half-normal prior on the SD with scale ", numText(r$halfNormalScale),
                          " times the lme4 SD of the effect"),
      halfCauchy = paste0("half-Cauchy prior on the SD with scale ", numText(r$halfCauchyScale),
                          " times the lme4 SD of the effect")))
  }
  fromData <- if (identical(a$source, "data") && !is.null(a$sdData))
    paste0(" - from the data: the lme4 SD of the effect is ", numText(a$sdData))
  v <- priorParams(prior)      # block size unknown: a Wishart's default df is written "k + 2"
  out <- switch(prior$dist,
    gamma = paste0("prior on the precision: Gamma with shape ", numText(v$shape), " and rate ", numText(v$rate),
                   " (INLA loggamma)"),
    pc = paste0("PC prior on the SD: P(SD > ", numText(v$u), ") is ", numText(v$alpha), " (SD is Exp(",
                numText(v$rate), "))"),
    halfNormal = paste0("half-normal prior on the SD with scale ", numText(v$scale)),
    halfCauchy = paste0("half-Cauchy prior on the SD with scale ", numText(v$scale)),
    wishart = paste0("Wishart prior on the precision matrix with ", if (is.null(v$df)) "k + 2" else numText(v$df),
                     " degrees of freedom (E[Q] is df times the scale matrix",
                     if (is.character(a$scale) || identical(a$source, "data"))
                       " - here the precision matrix of an lme4 fit of the formula divided by df",
                     ")"),
    normal = paste0("normal prior on the coefficient with mean ", numText(v$mean), " and SD ", numText(v$sd)),
    paste0("prior: ", prior$dist))
  paste0(out, fromData)
}

# Explanations of the distributions and priors written in the equations (keys distKey(<name>)).
distributionNotes <- function() {
  n <- function(par, note) estEntry(par, note = note)
  notes <- list(
    normal = n("normal distribution", "parametrised by its precision tau (variance 1/tau) as in INLA"),
    t = n("Student t distribution", "nu degrees of freedom and precision tau (of the standardised t)"),
    gamma = n("gamma distribution", "shape phi and rate phi/mu: mean mu and variance mu^2/phi (phi is the precision parameter)"),
    beta = n("beta distribution", "mean mu and precision phi: variance mu(1 - mu)/(1 + phi)"),
    poisson = n("Poisson distribution", "mean E lambda: lambda is the rate per unit of the exposure E"),
    negbin = n("negative binomial distribution", "mean E mu and size kappa: variance mu + mu^2/kappa (overdispersed Poisson)"),
    binomial = n("binomial distribution", "n trials with success probability p"),
    bernoulli = n("Bernoulli distribution", "one trial with success probability p"),
    betabinomial = n("beta-binomial distribution", "n trials; the success probability varies around p with overdispersion rho"),
    zero = n("point mass at zero", "a structural zero with probability pi (zero inflation)"),
    truncated = n("zero-truncated distribution", "the count distribution restricted to values above zero (type 0 zero inflation)"),
    inla = n("INLA likelihood", "see inla.doc() for its parametrisation"),
    weibull = n("Weibull distribution", "shape alpha and scale lambda (INLA weibullsurv variant 0): hazard alpha lambda t^(alpha - 1); right-censored survival times"),
    exp = n("exponential distribution", "rate lambda: constant hazard lambda (INLA exponentialsurv); right-censored survival times"),
    event = n("event indicator", "1: the event was observed at the time; 0: right-censored at the time"),
    priorGamma = n("gamma prior", "on the precision: shape a and rate b (mean a/b; INLA loggamma)"),
    priorPc = n("penalised-complexity prior", "the SD is exponential with the rate shown: P(SD > u) is alpha"),
    priorHalfNormal = n("half-normal prior", "on the SD: normal with the precision shown truncated to positive values"),
    priorHalfCauchy = n("half-Cauchy prior", "on the SD with the scale shown"),
    priorWishart = n("Wishart prior", "on the precision matrix Q with the degrees of freedom and scale matrix V shown: E[Q] is df V"),
    priorNormal = n("normal prior", "on the association coefficient: mean and precision (1/SD^2) shown"))
  stats::setNames(notes, distKey(names(notes)))
}
# A number with `digits` significant digits as text (NA as `na`); every display of an estimate
# goes through it (estimatesTable(), the LaTeX tooltips, numLatex()).
# The annotated normal-distribution symbol (random effects, u fields, Gaussian priors).
normalSymbol <- function() annot("\\mathcal{N}", distKey("normal"))
numText <- function(v, digits = 3, na = "NA")
  ifelse(is.na(v), na, trimws(formatC(signif(v, digits), digits = digits, format = "g")))

# The prior of the random-effect block of parsed formula p that contains `effect` for group g.
blockPrior <- function(p, g, effect, df, random = NULL) {
  # the prior as resolved by buildInlaFormula() (data-scaled priors: per effect)
  for (b in random) if (identical(b$group, g) && effect %in% b$effects) return(b$prior)
  for (blk in p$re) {
    if (blk$group != g) next
    if (effect %in% blockEffects(blk, p, df)) return(blk$prior)
  }
  NULL
}

# The lookup for a fitted model (keys: annotKey(), as used by the symbol builders):
# the entries of each source, combined in this order.
estimateLookup <- function(x) {
  ctx <- lookupContext(x)
  mergeEntries(distributionNotes(), fixedEntries(x, ctx), associationEntries(x, ctx),
               familyEntries(ctx), precisionEntries(x, ctx),
               concatEntries(lapply(names(ctx$eqs), function(r) responseEntries(x, ctx, r))),
               covariateEntries(x, ctx), zeroEntries(x))
}

# Combine named lists of entries; a later entry replaces an earlier one of the same key.
mergeEntries <- function(...) {
  keys <- list()
  for (part in list(...)) if (length(part)) keys[names(part)] <- part
  keys
}

# Concatenate a list of named lists of entries (unnamed: c() would prefix the keys).
concatEntries <- function(parts) do.call(c, unname(parts))

# One entry as a named list (key -> estEntry(...)).
keyedEntry <- function(key, ...) stats::setNames(list(estEntry(...)), key)

# Entries over the elements of idx: key(i) -> entry(i).
entriesOver <- function(idx, key, entry) {
  out <- lapply(idx, entry)
  names(out) <- vapply(idx, key, "", USE.NAMES = FALSE)
  out
}

# What the builders share: the parsed formulae and prepared data (named by response), the
# hyperparameter and variance-component tables, and how to label / look up an SD.
lookupContext <- function(x) {
  vc <- varianceComponentsOrNull(x, "the equation annotations")
  vt <- if (is.null(vc)) NULL else vc$table
  list(
    prep = bothParts(x$data), eqs = bothParts(x$forms), h = hyperTable(x), vt = vt,
    effLabel = function(e) effectLabel(e, x$inla$naming),
    # the SD row of `effect` | g of response r (g = NA: the residual), or NULL
    sdOf = function(r, g, effect) {
      if (is.null(vt)) return(NULL)
      k <- which(vt$parameter == "sd" & vt$response == r & vt$group %in% g & vt$effect == effect)
      if (length(k) == 1L) vt[k, ] else NULL
    },
    sdNote = function(row) if (!is.null(row)) paste0("SD ", numText(row$mean), " (", numText(row$sd), ")"))
}

# Lookup entries of the fixed effects (`fixed|r|col`), from fixedEffects().
fixedEntries <- function(x, ctx) {
  fe <- fixedEffects(x)
  entriesOver(seq_len(nrow(fe)), function(i) annotKey("fixed", fe$response[i], fe$term[i]), function(i)
    estEntry(paste0("fixed effect of ", ctx$effLabel(fe$term[i]), " in ", fe$response[i]),
             fe$mean[i], fe$sd[i], fe$q0.025[i], fe$q0.975[i]))
}

# Lookup entries of the associations (`assoc|out|long`), from associations(), noting the prior.
associationEntries <- function(x, ctx) {
  a <- associations(x)
  entriesOver(seq_len(nrow(a)), function(i) annotKey("assoc", a$outcome[i], a$longitudinal[i]), function(i) {
    sc <- Filter(function(s) s$target == a$longitudinal[i], ctx$eqs[[a$outcome[i]]]$scaled)
    pr <- if (length(sc)) sc[[1L]]$prior
    estEntry(paste0("association gamma (", linkLabel(a$longitudinal[i], a$outcome[i]), ")"),
             a$mean[i], a$sd[i], a$q0.025[i], a$q0.975[i],
             note = if (is.null(pr)) paste0(priorNote(defaultCoefPrior()), " (default)") else priorNote(pr))
  })
}

# Likelihood hyperparameters, keyed by the family's short name (precisions: "prec").
familyEntries <- function(ctx) {
  h <- ctx$h
  short <- function(i) h$short[i]
  entriesOver(which(h$type %in% "family"), function(i) annotKey("family", h$response[i], short(i)), function(i) {
    res <- if (short(i) == "prec") ctx$sdOf(h$response[i], NA, "residual")
    fam <- ctx$eqs[[h$response[i]]]$family
    estEntry(paste0(h$parameter[i], " of ", h$response[i]), h$mean[i], h$sd[i], h$q0.025[i], h$q0.975[i],
             note = paste(c(ctx$sdNote(res), priorNote(fam$hyper[[short(i)]])), collapse = ". "))
  })
}

# Precisions of the random effects (not the correlations).
precisionEntries <- function(x, ctx) {
  h <- ctx$h
  rows <- which(h$type %in% "random" & !is.na(h$effect) & is.na(h$effect2))
  entriesOver(rows, function(i) annotKey("prec", h$response[i], h$group[i], h$effect[i]), function(i) {
    r <- h$response[i]
    eff <- h$effect[i]
    estEntry(paste0("precision of random effect ", ctx$effLabel(eff), " | ", h$group[i], " in ", r),
             h$mean[i], h$sd[i], h$q0.025[i], h$q0.975[i],
             note = paste(c(ctx$sdNote(ctx$sdOf(r, h$group[i], eff)),
                            priorNote(blockPrior(ctx$eqs[[r]], h$group[i], eff, ctx$prep[[r]], x$inla$random[[r]]))),
                          collapse = ". "))
  })
}

# Per response: its random-effect blocks, residual, response and fitted mean.
responseEntries <- function(x, ctx, r) {
  p <- ctx$eqs[[r]]
  res <- ctx$sdOf(r, NA, "residual")
  # response and fitted mean: observed / posterior values over the observations
  y <- x$inla$stacks[[r]]$main[[x$inla$responses[[r]]]]
  fv <- fittedValues(x, r)
  c(concatEntries(lapply(seq_along(p$re), function(i) blockEntries(ctx, r, p, ctx$prep[[r]], p$re[[i]], i))),
    if (!is.null(res)) keyedEntry(annotKey("eps", r), paste0("residual of ", r), res$mean, res$sd, res$q0.025,
                                  res$q0.975, item = "SD"),
    keyedEntry(annotKey("resp", r),
               paste0("response ", p$lhs, " (", if (p$part == "out") "outcome" else "longitudinal", ")"),
               mean(y, na.rm = TRUE), lo = min(y, na.rm = TRUE), hi = max(y, na.rm = TRUE), item = "observed",
               band = "range",
               note = paste0(length(y), " observations", if (anyNA(y)) paste0(" (", sum(is.na(y)), " missing)"))),
    keyedEntry(annotKey("mean", r),
               paste0(if (isTRUE(x$inla$survival[r])) "fitted hazard scale lambda of " else "fitted mean of ",
                      r, " (one per observation: fittedValues())"),
               mean(fv$mean), lo = min(fv$mean), hi = max(fv$mean), item = "posterior means", band = "range"))
}

# One random-effect block of response r: each effect's SD, the number of levels and, for
# the block, its SDs and correlations (the matrix Q).
blockEntries <- function(ctx, r, p, df, blk, i) {
  g <- blk$group
  effects <- blockEffects(blk, p, df)
  b <- entriesOver(effects, function(e) annotKey("b", r, g, e), function(e) {
    row <- ctx$sdOf(r, g, e)
    estEntry(paste0("random effect ", ctx$effLabel(e), " | ", g, " of ", r),
             row$mean, row$sd, row$q0.025, row$q0.975, item = if (!is.null(row)) "SD",
             note = "one value per level: randomEffects()")
  })
  nl <- if (!is.null(df[[g]])) length(unique(df[[g]]))
  nlev <- keyedEntry(annotKey("nlev", r, g), paste0("levels of ", g, " in the data of ", r),
                     note = if (!is.null(nl)) paste0(nl, " levels"))
  vt <- ctx$vt
  k <- if (!is.null(vt)) which(vt$response == r & vt$blockId %in% i)
  block <- if (length(k))
    keyedEntry(blockKey(r, g, effects), paste0("covariance of the random effects | ", g, " of ", r),
               vt$mean[k], vt$sd[k], vt$q0.025[k], vt$q0.975[k], note = priorNote(blk$prior),
               item = ifelse(vt$parameter[k] == "sd", paste0("SD ", ctx$effLabel(vt$effect[k])),
                             paste0("corr ", ctx$effLabel(vt$effect[k]), ":", ctx$effLabel(vt$effect2[k]))))
  c(b, nlev, block)
}

# Covariates and family inputs: summary over the data frames that contain them.
covariateEntries <- function(x, ctx) {
  vars <- x$forms$variables
  concatEntries(lapply(unique(vars$variable[vars$role %in% c("covariate", "aux")]), function(v) {
    vals <- unlist(lapply(ctx$prep, function(df) if (!is.null(df[[v]]) && is.numeric(df[[v]])) df[[v]]))
    if (length(vals))
      return(keyedEntry(annotKey("cov", v), paste0("variable ", v), mean(vals, na.rm = TRUE),
                        lo = min(vals, na.rm = TRUE), hi = max(vals, na.rm = TRUE), item = "values", band = "range"))
    lev <- unique(unlist(lapply(ctx$prep, function(df) if (!is.null(df[[v]])) as.character(df[[v]]))))
    if (length(lev))
      keyedEntry(annotKey("cov", v), paste0("variable ", v),
                 note = paste0(length(lev), " levels: ", paste(utils::head(sort(lev), 6), collapse = " "),
                               if (length(lev) > 6) " ..."))
  }))
}

# Zero equations: u, its fixed precision, the matched combinations.
zeroEntries <- function(x) {
  concatEntries(lapply(x$inla$zero, function(z) {
    n <- length(x$inla$levels[[z$levels]])
    vars <- paste(z$vars, collapse = " ")
    c(keyedEntry(annotKey("u", z$long, z$out), paste0("zero-equation field u (", linkLabel(z$long, z$out), ")"),
                 note = paste0("linear predictor of ", z$long, " at the ", n, " combinations of (", vars,
                               ") in the data of ", z$out, ". Its precision is fixed")),
      keyedEntry(annotKey("utau", z$long, z$out), paste0("precision of u (", linkLabel(z$long, z$out), ")"),
                 note = "fixed at exp(-6) (not estimated)"),
      keyedEntry(annotKey("D", z$out, z$long), paste0("combinations of (", vars, ") in the data of ", z$out),
                 note = paste0(n, " unique combinations")))
  }))
}
