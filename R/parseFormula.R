#' Collect information about all formulae of a joint model
#'
#' Parses every outcome and longitudinal formula and returns, per formula, the
#' original formula, its lme4-compatible simplification, and the extra
#' parameters given to `re()` and `scaled()` terms.
#'
#' @param formOut List of outcome formulae. May contain `re()` and `scaled()` terms.
#' @param formLong List of longitudinal formulae. May contain `re()` terms only.
#' @param familyOut,familyLong Likelihood families: named character vectors
#'   (names = response variables of `formOut` / `formLong`, values = INLA
#'   likelihood names, see `names(INLA::inla.models()$likelihood)`), e.g.
#'   `c(y = "beta")`. Responses not listed get `"gaussian"`. Values are
#'   validated when INLA is installed.
#' @return A list with elements `out` and `long`, each a list named by response
#'   variable, whose elements are the result of [parseJsmFormula()] with an
#'   added `family` (outcome elements additionally have `linkVariables`:
#'   grouping variables and covariates of their `scaled()` targets not already
#'   used by the outcome formula), and
#'   `variables`, a data frame with one row per variable and role, with columns
#'   `variable`, `role` (`"response"`, `"covariate"` or `"group"`), `part` and
#'   `formula` (response of the formula it occurs in).
#' @export
collectFormulas <- function(formOut, formLong, familyOut = NULL, familyLong = NULL) {
  if (inherits(formOut, "formula")) formOut <- list(formOut)
  if (inherits(formLong, "formula")) formLong <- list(formLong)

  out <- lapply(formOut, parseJsmFormula, part = "out")
  long <- lapply(formLong, parseJsmFormula, part = "long")
  names(out) <- vapply(out, `[[`, "", "response")
  names(long) <- vapply(long, `[[`, "", "response")

  allResp <- c(names(out), names(long))
  dup <- unique(allResp[duplicated(allResp)])
  if (length(dup))
    stop("Response(s) defined more than once: ", paste(dup, collapse = ", "), call. = FALSE)

  for (o in out) {
    bad <- setdiff(vapply(o$scaled, `[[`, "", "target"), names(long))
    if (length(bad))
      stop("scaled() in formula for '", o$response, "' refers to unknown longitudinal ",
           "response(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # variables of the scaled() targets (other than their response) must also be
  # available in the outcome data: the zero equations are evaluated at outcome rows
  for (r in names(out)) {
    targets <- vapply(out[[r]]$scaled, `[[`, "", "target")
    out[[r]]$linkVariables <- setdiff(
      unique(unlist(lapply(long[targets], function(l) l$variables[c("group", "covariate")]),
                    use.names = FALSE)),
      unlist(out[[r]]$variables, use.names = FALSE))
    if (is.null(out[[r]]$linkVariables)) out[[r]]$linkVariables <- character()
  }

  famOut <- resolveFamilies(familyOut, names(out), "familyOut", names(long))
  famLong <- resolveFamilies(familyLong, names(long), "familyLong", names(out))
  for (r in names(out)) out[[r]]$family <- famOut[[r]]
  for (r in names(long)) long[[r]]$family <- famLong[[r]]

  list(out = out, long = long, variables = variableTable(c(out, long)))
}

#' Likelihood families available in INLA
#'
#' @return Character vector of INLA likelihood names, or `NULL` if INLA is not
#'   installed.
#' @export
inlaFamilies <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) return(NULL)
  names(INLA::inla.models()$likelihood)
}

defaultFamily <- "gaussian"

# Validate a named family vector and complete it with the default family.
resolveFamilies <- function(family, resp, argName, otherResp) {
  res <- stats::setNames(rep(defaultFamily, length(resp)), resp)
  if (is.null(family)) return(res)
  if (!is.character(family) || anyNA(family))
    stop("`", argName, "` must be a named character vector.", call. = FALSE)
  nms <- names(family)
  if (length(family) && (is.null(nms) || any(nms == "") || anyDuplicated(nms)))
    stop("`", argName, "` must have unique, non-empty names (the response variables).",
         call. = FALSE)
  unknown <- setdiff(nms, resp)
  if (length(unknown)) {
    other <- intersect(unknown, otherResp)
    stop("`", argName, "` names do not match any response of this part: ",
         paste(unknown, collapse = ", "),
         if (length(other)) paste0(" (", paste(other, collapse = ", "),
                                   " belongs to the other part)"),
         call. = FALSE)
  }
  known <- inlaFamilies()
  if (is.null(known)) {
    warning("INLA is not installed: likelihood families in `", argName,
            "` are not validated.", call. = FALSE)
  } else {
    bad <- setdiff(family, known)
    if (length(bad))
      stop("Unknown INLA likelihood family in `", argName, "`: ",
           paste(bad, collapse = ", "), ". See inlaFamilies().", call. = FALSE)
  }
  res[nms] <- family
  res
}

# One row per (variable, role, formula) across all parsed formulae.
variableTable <- function(parsed) {
  rows <- lapply(parsed, function(p) {
    v <- p$variables
    n <- lengths(v)
    data.frame(
      variable = unlist(v, use.names = FALSE),
      role = rep(names(v), n),
      part = rep(p$part, sum(n)),
      formula = rep(p$response, sum(n)),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, c(list(data.frame(variable = character(), role = character(),
                                          part = character(), formula = character())),
                          rows))
  rownames(res) <- NULL
  res
}

# Data variables used by a formula, by role. all.vars() skips function names,
# so log(y), I(t^2), poly(t, 2), ns(t, df = 3) etc. yield only the variables.
collectVariables <- function(lhs, fixed, re) {
  allVars <- function(exprs) unique(unlist(lapply(exprs, all.vars)))
  list(
    response = all.vars(lhs),
    covariate = allVars(c(fixed, lapply(re, function(b) b$term[[2L]]))),
    group = allVars(lapply(re, function(b) b$term[[3L]]))
  )
}

#' Parse a single joint-model formula
#'
#' @param formula A two-sided formula. Random effects may be given in lme4
#'   notation, `(1 + t | id)`, or as `re(1 + t | id, ...)` where `...` holds extra
#'   (named) parameters, e.g. `prior = pcSd(u = 1)` (see [priors]). Outcome
#'   formulae may contain `scaled(x, ...)` terms.
#' @param part Either `"out"` (outcome) or `"long"` (longitudinal).
#' @return A list with elements
#'   * `response`: name of the response variable (e.g. `y` for `log(y) ~ ...`);
#'     used to name formulae, data frames and `scaled()` targets,
#'   * `lhs`: the left-hand side as written (e.g. `log(y)`),
#'   * `part`: `"out"` or `"long"`,
#'   * `formula`: the formula as given by the user,
#'   * `variables`: list with character vectors `response` (variables in the
#'     left-hand side, e.g. `y` for `log(y)`), `covariate` (variables in fixed-effect
#'     terms and on the left of `|` in random-effect terms, e.g. `t` for `I(t^2)`)
#'     and `group` (variables on the right of `|`, e.g. `id`, or `a` and `b` for `a:b`),
#'   * `fixedFormula`: the fixed-effects part only (no random effects, no `scaled()`),
#'   * `lme4Formula`: formula with `re()` replaced by lme4 bar terms and `scaled()` removed,
#'   * `re`: list of random-effect blocks, each with `term` (the lme4 bar expression),
#'     `group` (grouping factor name), `correlated` (`FALSE` for `||`) and
#'     `params` (named list of extra arguments, evaluated in the formula environment)
#'     and `prior` (the `prior` argument, see [priors], or `NULL`),
#'   * `scaled`: (outcome formulae only) list of associations, each with `target`
#'     (longitudinal response name) and `params`.
#' @export
parseJsmFormula <- function(formula, part = c("out", "long")) {
  part <- match.arg(part)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("Each model formula must be a two-sided formula.", call. = FALSE)

  env <- environment(formula)
  lhs <- deparse1(formula[[2L]])
  response <- all.vars(formula[[2L]])
  if (length(response) != 1L)
    stop("The left-hand side of a formula must contain exactly one variable (found '",
         lhs, "').", call. = FALSE)
  terms <- flattenTerms(formula[[3L]])

  keep <- list()
  fixed <- list()
  fixedTerms <- list()
  re <- list()
  scaled <- list()

  for (tm in terms) {
    e <- tm$expr
    if (is.call(e) && identical(e[[1L]], quote(re))) {
      blk <- parseRe(e, env, response)
      re[[length(re) + 1L]] <- blk
      keep[[length(keep) + 1L]] <- list(expr = call("(", blk$term), sign = tm$sign)
    } else if (is.call(e) && identical(e[[1L]], quote(scaled))) {
      if (part == "long")
        stop("scaled() is only allowed in outcome formulae (found in formula for '",
             response, "').", call. = FALSE)
      scaled[[length(scaled) + 1L]] <- parseScaled(e, env, response)
    } else if (is.call(e) && identical(e[[1L]], quote(`(`)) && isBar(e[[2L]])) {
      re[[length(re) + 1L]] <- makeReBlock(e[[2L]], list())
      keep[[length(keep) + 1L]] <- tm
    } else {
      if (any(c("re", "scaled") %in% all.names(e)))
        stop("re() and scaled() must be top-level terms (formula for '", response, "').",
             call. = FALSE)
      keep[[length(keep) + 1L]] <- tm
      fixed[[length(fixed) + 1L]] <- e
      fixedTerms[[length(fixedTerms) + 1L]] <- tm
    }
  }

  lme4Formula <- buildFormula(formula[[2L]], keep, env)

  res <- list(
    response = response,
    lhs = lhs,
    part = part,
    formula = formula,
    lme4Formula = lme4Formula,
    fixedFormula = buildFormula(formula[[2L]], fixedTerms, env),
    variables = collectVariables(formula[[2L]], fixed, re),
    re = re
  )
  if (part == "out") res$scaled <- scaled
  res
}

# Split a RHS expression into additive terms, keeping track of the sign.
flattenTerms <- function(e, sign = 1L) {
  if (is.call(e) && identical(e[[1L]], quote(`+`))) {
    if (length(e) == 2L) return(flattenTerms(e[[2L]], sign))
    return(c(flattenTerms(e[[2L]], sign), flattenTerms(e[[3L]], sign)))
  }
  if (is.call(e) && identical(e[[1L]], quote(`-`))) {
    if (length(e) == 2L) return(flattenTerms(e[[2L]], -sign))
    return(c(flattenTerms(e[[2L]], sign), flattenTerms(e[[3L]], -sign)))
  }
  list(list(expr = e, sign = sign))
}

isBar <- function(e) {
  is.call(e) && (identical(e[[1L]], quote(`|`)) || identical(e[[1L]], quote(`||`)))
}

makeReBlock <- function(bar, params, prior = NULL) {
  list(
    term = bar,
    group = deparse1(bar[[3L]]),
    correlated = identical(bar[[1L]], quote(`|`)),
    params = params,
    prior = prior
  )
}

# Evaluate the extra (named) arguments of an re()/scaled() call.
extraParams <- function(args, env, what, response) {
  if (!length(args)) return(list())
  nms <- names(args)
  if (is.null(nms) || any(nms == ""))
    stop("All extra arguments of ", what, "() must be named (formula for '",
         response, "').", call. = FALSE)
  lapply(args, eval, envir = env)
}

parseRe <- function(e, env, response) {
  args <- as.list(e)[-1L]
  if (!length(args) || !isBar(args[[1L]]))
    stop("The first argument of re() must be a random-effects term like `1 + t | id` ",
         "(formula for '", response, "').", call. = FALSE)
  if (!is.null(names(args)) && names(args)[1L] != "")
    stop("The first argument of re() must not be named (formula for '", response, "').",
         call. = FALSE)
  params <- extraParams(args[-1L], env, "re", response)
  prior <- params[["prior"]]
  params <- params[names(params) != "prior"]
  if (!length(params)) params <- list()
  if (!is.null(prior) && !inherits(prior, "mvremPrior"))
    stop("`prior` in re() must be created by a prior constructor (gammaPrec(), pcSd(), ",
         "halfNormalSd(), halfCauchySd(), wishartPrec()); raw INLA settings go in `hyper` ",
         "(formula for '", response, "').", call. = FALSE)
  makeReBlock(args[[1L]], params, prior)
}

parseScaled <- function(e, env, response) {
  args <- as.list(e)[-1L]
  x <- if (length(args)) args[[1L]] else NULL
  if (!is.name(x))
    stop("The first argument of scaled() must be the name of a longitudinal ",
         "response, e.g. scaled(x1) (formula for '", response, "').", call. = FALSE)
  if (!is.null(names(args)) && names(args)[1L] != "")
    stop("The first argument of scaled() must not be named (formula for '", response,
         "').", call. = FALSE)
  params <- extraParams(args[-1L], env, "scaled", response)
  if ("prior" %in% names(params))
    stop("`prior` is not supported in scaled() yet (formula for '", response, "').",
         call. = FALSE)
  list(target = as.character(x), params = params)
}

buildFormula <- function(lhs, terms, env) {
  rhs <- NULL
  for (tm in terms) {
    if (is.null(rhs)) {
      rhs <- if (tm$sign > 0) tm$expr else call("-", tm$expr)
    } else {
      rhs <- call(if (tm$sign > 0) "+" else "-", rhs, tm$expr)
    }
  }
  if (is.null(rhs)) rhs <- 1
  f <- call("~", lhs, rhs)
  stats::as.formula(f, env = env)
}
