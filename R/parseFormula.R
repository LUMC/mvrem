#' Collect information about all formulae of a joint model
#'
#' Parses every outcome and longitudinal formula and returns, per formula, the
#' original formula, its lme4-compatible simplification, and the extra
#' parameters given to `re()` and `scaled()` terms.
#'
#' Every outcome is associated with every longitudinal response (a scaled copy
#' of its linear predictor, see [buildInlaFormula()]). A `scaled(x, ...)` term
#' is only needed to give non-default parameters of the association with `x`.
#'
#' @param formOut List of outcome formulae. May contain `re()` and `scaled()` terms.
#' @param formLong List of longitudinal formulae. May contain `re()` terms only.
#' @param familyOut,familyLong Likelihood families: named lists (names =
#'   response variables of `formOut` / `formLong`) of family objects (see
#'   [families]), e.g. `list(y = betaFamily(phi = gammaPrec(1, 0.1)))`. An
#'   element may also be an INLA likelihood name (`"beta"`; a named character
#'   vector is accepted too), meaning the family's constructor with its
#'   defaults ([inlaFamily()] for likelihoods without a constructor; names are
#'   validated when INLA is installed). Responses not listed get
#'   [gaussianFamily()].
#' @return A list with elements `out` and `long`, each a list named by response
#'   variable, whose elements are the result of [parseJsmFormula()] with an
#'   added `family` (an `mvremFamily` object; `trials()` / `exposure()` on the
#'   left-hand side are checked against it). In outcome elements, `scaled` holds
#'   one association per longitudinal response (in the order of `formLong`;
#'   `params` from the formula's `scaled()` term, or empty), and
#'   `linkVariables` the grouping variables and covariates of all longitudinal
#'   formulae not already used by the outcome formula (the outcome data must
#'   contain them). Also returned:
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
    targets <- vapply(o$scaled, `[[`, "", "target")
    bad <- setdiff(targets, names(long))
    if (length(bad))
      stop("scaled() in formula for '", o$response, "' refers to unknown longitudinal ",
           "response(s): ", paste(bad, collapse = ", "), call. = FALSE)
    dupT <- unique(targets[duplicated(targets)])
    if (length(dupT))
      stop("scaled() in formula for '", o$response, "' is given more than once for: ",
           paste(dupT, collapse = ", "), call. = FALSE)
  }

  # every outcome is associated with every longitudinal response; scaled() only
  # supplies the parameters. The variables of the long formulae (other than
  # their response) must also be available in the outcome data: the zero
  # equations are evaluated at outcome rows
  longVars <- unique(unlist(lapply(long, function(l) l$variables[c("group", "covariate")]),
                            use.names = FALSE))
  for (r in names(out)) {
    given <- out[[r]]$scaled
    names(given) <- vapply(given, `[[`, "", "target")
    out[[r]]$scaled <- unname(lapply(names(long), function(l)
      list(target = l, params = if (is.null(given[[l]])) list() else given[[l]]$params)))
    out[[r]]$linkVariables <- setdiff(longVars, unlist(out[[r]]$variables, use.names = FALSE))
    if (is.null(out[[r]]$linkVariables)) out[[r]]$linkVariables <- character()
  }

  famOut <- resolveFamilies(familyOut, names(out), "familyOut", names(long))
  famLong <- resolveFamilies(familyLong, names(long), "familyLong", names(out))
  for (r in names(out)) out[[r]]$family <- famOut[[r]]
  for (r in names(long)) long[[r]]$family <- famLong[[r]]
  for (p in c(out, long)) checkAux(p)

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

# Validate a named list / character vector of families and complete it with
# the default family; elements become mvremFamily objects.
resolveFamilies <- function(family, resp, argName, otherResp) {
  res <- stats::setNames(lapply(resp, function(r) defaultFamily()), resp)
  if (is.null(family)) return(res)
  if (inherits(family, "mvremFamily"))
    stop("`", argName, "` must be a named list of families, e.g. list(y = ",
         family$kind, "Family()).", call. = FALSE)
  isName <- function(f) is.character(f) && length(f) == 1L && !is.na(f)
  ok <- if (is.character(family)) !anyNA(family) else
    is.list(family) && all(vapply(family, function(f) inherits(f, "mvremFamily") || isName(f), NA))
  if (!ok)
    stop("`", argName, "` must be a named list of family objects (see ?families) or ",
         "INLA likelihood names.", call. = FALSE)
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
  for (r in nms) {
    f <- family[[r]]
    res[[r]] <- if (inherits(f, "mvremFamily")) f else familyFromName(f, paste0("`", argName, "`"))
  }
  res
}

# trials() / exposure() on the left-hand side must be used by the family.
checkAux <- function(p) {
  bad <- setdiff(names(p$aux), p$family$aux)
  if (length(bad))
    stop(paste0(bad, "()", collapse = ", "), " on the left-hand side of the formula for '",
         p$response, "' is not used by family '", p$family$family, "'",
         if (length(p$family$aux)) paste0(" (it takes ", paste0(p$family$aux, "()"), ")"),
         ".", call. = FALSE)
  invisible()
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
collectVariables <- function(lhs, fixed, re, aux) {
  allVars <- function(exprs) unique(unlist(lapply(exprs, all.vars)))
  vars <- list(
    response = all.vars(lhs),
    covariate = allVars(c(fixed, lapply(re, function(b) b$term[[2L]]))),
    group = allVars(lapply(re, function(b) b$term[[3L]]))
  )
  vars$aux <- if (length(aux)) allVars(aux) else character()
  vars
}

#' Parse a single joint-model formula
#'
#' @param formula A two-sided formula. The left-hand side is the response
#'   (an expression of one variable, e.g. `y` or `log(y)`), optionally
#'   followed by `|` and family inputs (see [families]):
#'   `y | trials(n) ~ ...`, `y | exposure(pt) ~ ...`. Random effects may be given in lme4
#'   notation, `(1 + t | id)`, or as `re(1 + t | id, ...)` where `...` holds extra
#'   (named) parameters, e.g. `prior = pcSd(u = 1)` (see [priors]). Outcome
#'   formulae may contain `scaled(x, ...)` terms.
#' @param part Either `"out"` (outcome) or `"long"` (longitudinal).
#' @return A list with elements
#'   * `response`: name of the response variable (e.g. `y` for `log(y) ~ ...`);
#'     used to name formulae, data frames and `scaled()` targets,
#'   * `lhs`: the response expression as written (e.g. `log(y)`, without the
#'     `| trials(n)` part),
#'   * `responseExpr`: the response expression as a call / name,
#'   * `aux`: named list of the family inputs of the left-hand side
#'     (`trials`, `exposure`; expressions), empty if none,
#'   * `part`: `"out"` or `"long"`,
#'   * `formula`: the formula as given by the user,
#'   * `variables`: list with character vectors `response` (variables in the
#'     left-hand side, e.g. `y` for `log(y)`), `covariate` (variables in fixed-effect
#'     terms and on the left of `|` in random-effect terms, e.g. `t` for `I(t^2)`)
#'     `group` (variables on the right of `|`, e.g. `id`, or `a` and `b` for `a:b`)
#'     and `aux` (variables of `trials()` / `exposure()`),
#'   * `fixedFormula`: the fixed-effects part only (no random effects, no `scaled()`),
#'   * `lme4Formula`: formula with `re()` replaced by lme4 bar terms and `scaled()` removed,
#'   * `re`: list of random-effect blocks, each with `term` (the lme4 bar expression),
#'     `group` (grouping factor name), `correlated` (`FALSE` for `||`) and
#'     `params` (named list of extra arguments, evaluated in the formula environment)
#'     and `prior` (the `prior` argument, see [priors], or `NULL`),
#'   * `scaled`: (outcome formulae only) list of the formula's `scaled()` terms,
#'     each with `target` (longitudinal response name) and `params`
#'     ([collectFormulas()] completes it to all longitudinal responses).
#' @export
parseJsmFormula <- function(formula, part = c("out", "long")) {
  part <- match.arg(part)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("Each model formula must be a two-sided formula.", call. = FALSE)

  env <- environment(formula)
  lhsParts <- parseLhs(formula[[2L]])
  responseExpr <- lhsParts$response
  lhs <- deparse1(responseExpr)
  response <- all.vars(responseExpr)
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

  lme4Formula <- buildFormula(responseExpr, keep, env)

  res <- list(
    response = response,
    lhs = lhs,
    responseExpr = responseExpr,
    aux = lhsParts$aux,
    part = part,
    formula = formula,
    lme4Formula = lme4Formula,
    fixedFormula = buildFormula(responseExpr, fixedTerms, env),
    variables = collectVariables(responseExpr, fixed, re, lhsParts$aux),
    re = re
  )
  if (part == "out") res$scaled <- scaled
  res
}

# Family inputs allowed on the left-hand side: y | trials(n) + exposure(E).
auxTerms <- c("trials", "exposure")

# Split the left-hand side into the response expression and the family inputs.
parseLhs <- function(e) {
  if (!(is.call(e) && identical(e[[1L]], quote(`|`))))
    return(list(response = e, aux = list()))
  aux <- list()
  for (tm in flattenTerms(e[[3L]])) {
    a <- tm$expr
    fn <- if (is.call(a) && is.name(a[[1L]])) as.character(a[[1L]]) else ""
    if (tm$sign < 0 || !fn %in% auxTerms)
      stop("Unsupported term on the left-hand side after `|`: ", deparse1(a), "; use ",
           paste0(auxTerms, "()", collapse = " or "), ".", call. = FALSE)
    if (length(a) != 2L || (!is.null(names(a)) && names(a)[2L] != ""))
      stop(fn, "() takes exactly one unnamed argument (found ", deparse1(a), ").",
           call. = FALSE)
    if (!is.null(aux[[fn]]))
      stop(fn, "() is given more than once on the left-hand side.", call. = FALSE)
    aux[[fn]] <- a[[2L]]
  }
  list(response = e[[2L]], aux = aux)
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
