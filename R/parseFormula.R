#' Collect information about all formulae of a joint model
#'
#' Parses every outcome and longitudinal formula and returns, per formula, the
#' original formula, its lme4-compatible simplification, and the extra
#' parameters given to `re()` and `scaled()` terms.
#'
#' An outcome is associated with a longitudinal response `x` (a scaled copy of
#' its linear predictor, see [buildInlaFormula()]) only through a `scaled(x, ...)`
#' term in its formula; at least one such term is required in the model
#' (otherwise an error). Further arguments of `scaled(x, ...)` set non-default
#' parameters of the association. Every longitudinal response must be named in
#' some `scaled()` term (otherwise an error): all equations are linked.
#'
#' The formulae are given either as ONE list for both parts (`form` alone, the
#' primary style): a formula with a `scaled()` term is an outcome, a formula
#' without one is longitudinal and must be the target of some `scaled()` term
#' (see [splitFormulas()]); or as the pair of an outcome and a longitudinal list
#' (`form` + `formLong`, or `formOut` + `formLong`), which is the only way to
#' give an outcome without a `scaled()` term. Within each part the formulae keep
#' the order given; `order` records the order of all formulae as passed, which
#' every result of a model follows (see [modelResponses()]).
#'
#' @param form One list of formulae for both parts (also a single formula):
#'   outcomes are the formulae with `scaled()` terms. With `formLong`, the list
#'   of outcome formulae (the pair style, as `formOut`).
#' @param formLong List of longitudinal formulae (pair style). May contain `re()`
#'   terms only.
#' @param formOut List of outcome formulae (pair style, with `formLong`; an
#'   alternative name of `form`, giving both is an error). May contain `re()`
#'   and `scaled()` terms.
#' @param family Likelihood families of all responses (outcome and
#'   longitudinal): a named list (names = response variables) of family objects
#'   (see [families]) or INLA likelihood names, see below.
#' @param familyOut,familyLong Likelihood families of the outcome /
#'   longitudinal responses only, in the same form as `family` (names = response
#'   variables of `formOut` / `formLong`); entries override those of `family`
#'   for the same response. Named lists of family objects (see
#'   [families]), e.g. `list(y = betaFamily(phi = gammaPrec(1, 0.1)))`. An
#'   element may also be an INLA likelihood name (`"beta"`; a named character
#'   vector is accepted too), meaning the family's constructor with its
#'   defaults ([inlaFamily()] for likelihoods without a constructor; names are
#'   validated when INLA is installed). Responses not listed get
#'   [gaussianFamily()].
#' @return A list with elements `out` and `long`, each a list named by response
#'   variable, whose elements are the result of [parseJsmFormula()] with an
#'   added `family` (an `mvremFamily` object; `trials()` / `exposure()` / `event()` / `weights()` on the
#'   left-hand side are checked against it). In outcome elements, `scaled` holds
#'   one association per `scaled()` term of the formula (in the order written;
#'   `target`, `params`, `prior`), and `linkVariables` the grouping variables and
#'   covariates of the linked longitudinal formulae not already used by the
#'   outcome formula (the outcome data must contain them). Also returned:
#'   `variables`, a data frame with one row per variable and role, with columns
#'   `variable`, `role` (`"response"`, `"covariate"` or `"group"`), `part` and
#'   `formula` (response of the formula it occurs in), and `order`, the
#'   response names in the order the formulae were passed (one list: as
#'   listed; pair style: the outcomes, then the longitudinal responses).
#' @export
collectFormulas <- function(form = NULL, formLong = NULL, family = NULL, familyOut = NULL,
                            familyLong = NULL, formOut = NULL) {
  parts <- formulaParts(form, formLong, formOut)
  formOut <- parts$out
  formLong <- parts$long

  out <- lapply(formOut, parseJsmFormula, part = "out")
  long <- lapply(formLong, parseJsmFormula, part = "long")
  names(out) <- vapply(out, `[[`, "", "response")
  names(long) <- vapply(long, `[[`, "", "response")

  allResp <- c(names(out), names(long))
  dup <- unique(allResp[duplicated(allResp)])
  if (length(dup))
    stop("Response(s) defined more than once: ", paste(dup, collapse = ", "), call. = FALSE)

  checkScaledTargets(out, long)
  out <- lapply(out, addLinkVariables, long = long)
  fams <- attachFamilies(out, long, family, familyOut, familyLong)
  out <- fams$out
  long <- fams$long

  # the order in which the formulae were passed (one list: as listed; pair: outcomes, then
  # longitudinal): every result of a model follows it (modelResponses())
  order <- character(length(parts$isOut))
  order[parts$isOut] <- names(out)
  order[!parts$isOut] <- names(long)

  list(out = out, long = long, variables = variableTable(c(out, long)), order = order)
}

# The longitudinal responses an outcome is linked to: the targets of its scaled() terms,
# in the order written (empty for a longitudinal formula or an outcome without scaled()).
scaledTargets <- function(p) vapply(p$scaled, `[[`, "", "target")

# collectFormulas(): every scaled() target is a longitudinal response, named once per
# outcome; the model has at least one scaled() term and every longitudinal response is
# the target of one (all equations are linked).
checkScaledTargets <- function(out, long) {
  for (o in out) {
    targets <- scaledTargets(o)
    bad <- setdiff(targets, names(long))
    if (length(bad))
      stop("scaled() in formula for '", o$response, "' refers to unknown longitudinal ",
           "response(s): ", paste(bad, collapse = ", "), call. = FALSE)
    dupT <- unique(targets[duplicated(targets)])
    if (length(dupT))
      stop("scaled() in formula for '", o$response, "' is given more than once for: ",
           paste(dupT, collapse = ", "), call. = FALSE)
  }
  if (!any(vapply(out, function(o) length(o$scaled) > 0L, NA)))
    stop("No scaled() term in the outcome formulae: at least one is required to link an ",
         "outcome to a longitudinal response, e.g. y ~ t + scaled(x1) with x1 ~ t + (1 | id).",
         call. = FALSE)
  unlinked <- setdiff(names(long), unlist(lapply(out, scaledTargets)))
  if (length(unlinked))
    stop("Longitudinal response(s) not named in any scaled() term: ",
         paste(unlinked, collapse = ", "), ". Every longitudinal equation must be linked to an ",
         "outcome, e.g. y ~ t + scaled(", unlinked[1L], ").", call. = FALSE)
  invisible()
}

# An outcome's `linkVariables`: the groups and covariates of its linked longitudinal
# formulae that the outcome formula does not use. The outcome data must contain them,
# because the zero equations of a link are evaluated at the outcome rows.
addLinkVariables <- function(o, long) {
  linkVars <- unique(unlist(lapply(long[scaledTargets(o)], function(l)
    l$variables[c("group", "covariate")]), use.names = FALSE))
  o$linkVariables <- setdiff(linkVars, unlist(o$variables, use.names = FALSE))
  if (is.null(o$linkVariables)) o$linkVariables <- character()
  o
}

# Each response's family (`family`, overridden per part by `familyOut` / `familyLong`),
# checked against the left-hand-side inputs; survival families only for outcomes.
attachFamilies <- function(out, long, family, familyOut, familyLong) {
  famOut <- resolveFamilies(familyOut, names(out), "familyOut", names(long), family)
  famLong <- resolveFamilies(familyLong, names(long), "familyLong", names(out), family)
  for (r in names(out)) out[[r]]$family <- famOut[[r]]
  for (r in names(long)) long[[r]]$family <- famLong[[r]]
  for (p in c(out, long)) checkAux(p)
  survLong <- names(long)[vapply(long, function(p) "event" %in% p$family$aux, NA)]
  if (length(survLong))
    stop("Survival family for longitudinal response(s) ", paste0("'", survLong, "'", collapse = ", "),
         ": survival families are only for outcomes (`formOut`).", call. = FALSE)
  list(out = out, long = long)
}

# The outcome and longitudinal formula lists from the formula arguments of
# collectFormulas() / jsm(): `form` alone = one list for both parts (split by
# splitFormulas()); `form` or `formOut` with `formLong` = the pair style. `isOut`: per
# formula in the order passed, whether it is an outcome.
formulaParts <- function(form = NULL, formLong = NULL, formOut = NULL) {
  if (!is.null(form) && !is.null(formOut))
    stop("Give the formulae either as `form` (one list, or the outcomes with `formLong`) ",
         "or as `formOut` with `formLong`, not both.", call. = FALSE)
  pairName <- if (is.null(form)) "formOut" else "form"
  if (is.null(form)) form <- formOut
  if (is.null(form))
    stop("No model formulae: give `form` (one list of outcome and longitudinal formulae) ",
         "or `formOut` with `formLong`.", call. = FALSE)
  if (inherits(form, "formula")) form <- list(form)
  if (inherits(formLong, "formula")) formLong <- list(formLong)
  if (!is.null(formLong))
    return(list(out = form, long = formLong,
                isOut = rep(c(TRUE, FALSE), c(length(form), length(formLong)))))
  if (pairName == "formOut")
    stop("`formOut` needs `formLong`; to give all formulae in one list use `form`.", call. = FALSE)
  res <- splitFormulas(form)
  res$isOut <- vapply(unname(form), function(f) hasCall(f[[length(f)]], "scaled"), NA)
  res
}

#' Split one list of formulae into outcome and longitudinal formulae
#'
#' A formula with a `scaled()` term is an outcome, a formula without one is
#' longitudinal. Used by [collectFormulas()] and [jsm()] when the formulae are
#' given as one list (`form` without `formLong`); the parts are then validated
#' as in the pair style.
#'
#' Errors: a formula without `scaled()` that no `scaled()` term names (it would
#' be neither linked nor an outcome; an outcome without `scaled()` needs the
#' pair `formOut` / `formLong`), and a `scaled()` target that is itself an
#' outcome (has `scaled()` terms). Other problems (unknown targets, duplicated
#' responses, syntax) are reported by [collectFormulas()].
#'
#' @param form A list of formulae (or a single formula).
#' @return `list(out, long)`: the outcome and the longitudinal formulae, each in
#'   the order given.
#' @examples
#' splitFormulas(list(x1 ~ t + (1 | id), y ~ t + scaled(x1)))
#' @export
splitFormulas <- function(form) {
  if (inherits(form, "formula")) form <- list(form)
  if (!is.list(form) || !length(form) ||
      !all(vapply(form, function(f) inherits(f, "formula"), NA)))
    stop("`form` must be a formula or a non-empty list of formulae.", call. = FALSE)
  form <- unname(form)
  isOut <- vapply(form, function(f) hasCall(f[[length(f)]], "scaled"), NA)
  parsed <- lapply(seq_along(form), function(i)
    parseJsmFormula(form[[i]], if (isOut[i]) "out" else "long"))
  resp <- vapply(parsed, `[[`, "", "response")
  res <- list(out = form[isOut], long = form[!isOut])
  # duplicated responses: reported by collectFormulas()
  if (anyDuplicated(resp)) return(res)
  showForm <- function(i) paste0("`", deparse1(form[[i]]), "`")

  for (i in which(isOut)) {
    targets <- scaledTargets(parsed[[i]])
    bad <- intersect(targets, resp[isOut])
    if (length(bad))
      stop(paste0("scaled(", bad, ") in ", showForm(i), " targets ",
                  vapply(match(bad, resp), showForm, ""),
                  ", which is an outcome (it has a scaled() term); a scaled() target must be a ",
                  "longitudinal formula", collapse = "\n"), ".", call. = FALSE)
  }
  targets <- unlist(lapply(parsed[isOut], scaledTargets))
  # unknown targets (e.g. a misspelt response): reported by collectFormulas()
  if (!all(targets %in% resp)) return(res)
  unlinked <- which(!isOut & !resp %in% targets)
  one <- length(unlinked) == 1L
  if (length(unlinked))
    stop(paste0(vapply(unlinked, showForm, ""), collapse = ", "),
         if (one) " has" else " have", " no scaled() term and ", if (one) "is" else "are",
         " not the target of one: link ", if (one) "it" else "each", " to an outcome with ",
         "scaled(", resp[unlinked[1L]], ") in the outcome's formula, or give an outcome ",
         "without scaled() in the pair `formOut` / `formLong`.", call. = FALSE)
  res
}

#' Likelihood families available in INLA
#'
#' @return Character vector of INLA likelihood names, or `NULL` if INLA is not
#'   installed.
#' @export
inlaFamilies <- function() {
  if (!hasInla()) return(NULL)
  names(INLA::inla.models()$likelihood)
}

# Validate a named list / character vector of families (argument `argName`)
# and return it as a list; NULL stays NULL.
checkFamilyList <- function(family, argName) {
  if (is.null(family)) return(NULL)
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
  as.list(family)
}

# Families of the responses `resp` of one part: the default family, replaced by
# the entries of `shared` (the `family` argument, for responses of either
# part) and then by those of `family` (the part's own argument `argName`).
# Elements become mvremFamily objects.
resolveFamilies <- function(family, resp, argName, otherResp, shared = NULL) {
  res <- stats::setNames(lapply(resp, function(r) defaultFamily()), resp)
  toFamily <- function(f, arg)
    if (inherits(f, "mvremFamily")) f else familyFromName(f, paste0("`", arg, "`"))

  shared <- checkFamilyList(shared, "family")
  unknownShared <- setdiff(names(shared), c(resp, otherResp))
  if (length(unknownShared))
    stop("`family` names do not match any response: ",
         paste(unknownShared, collapse = ", "), call. = FALSE)
  for (r in intersect(names(shared), resp)) res[[r]] <- toFamily(shared[[r]], "family")

  family <- checkFamilyList(family, argName)
  unknown <- setdiff(names(family), resp)
  if (length(unknown)) {
    other <- intersect(unknown, otherResp)
    stop("`", argName, "` names do not match any response of this part: ",
         paste(unknown, collapse = ", "),
         if (length(other)) paste0(" (", paste(other, collapse = ", "),
                                   " belongs to the other part)"),
         call. = FALSE)
  }
  for (r in names(family)) res[[r]] <- toFamily(family[[r]], argName)
  res
}

# trials() / exposure() / event() / weights() on the left-hand side must be used by the family.
checkAux <- function(p) {
  accepted <- c(p$family$aux, if (isTRUE(p$family$weights)) "weights")
  bad <- setdiff(names(p$aux), accepted)
  # B4: one "(it takes ...)" listing the accepted inputs, not one per input; weights()
  # (taken by many families) only when it is itself rejected
  takes <- if ("weights" %in% bad) accepted else p$family$aux
  if (length(bad))
    stop(paste0(bad, "()", collapse = ", "), " on the left-hand side of the formula for '",
         p$response, "' is not used by family '", p$family$family, "'",
         if (length(takes)) paste0(" (it takes ", paste0(takes, "()", collapse = ", "), ")"),
         ".", call. = FALSE)
  if ("event" %in% p$family$aux && is.null(p$aux$event))
    stop("Survival family '", p$family$family, "' of '", p$response, "' needs the event ",
         "indicator on the left-hand side: `", p$lhs, " | event(status) ~ ...`.", call. = FALSE)
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
#'   `y | trials(n) ~ ...`, `y | exposure(pt) ~ ...`, `time | event(status) ~ ...` (survival
#'   families: `status` 1 = event at `time`, 0 = right-censored), `y | weights(w) ~ ...`
#'   (observation weights: the precision of row i is multiplied by \eqn{w_i}; gaussian,
#'   Student-t, gamma, beta and [inlaFamily()]). Random effects may be given in lme4
#'   notation, `(1 + t | id)`, or as `re(1 + t | id, ...)` where `...` holds extra
#'   (named) parameters, e.g. `prior = pcSd(u = 1)` (see [priors]). Outcome
#'   formulae may contain `scaled(x, ...)` terms, with `prior = normalCoef(...)`
#'   for the association coefficient.
#' @param part Either `"out"` (outcome) or `"long"` (longitudinal).
#' @return A list with elements
#'   * `response`: name of the response variable (e.g. `y` for `log(y) ~ ...`);
#'     used to name formulae, data frames and `scaled()` targets,
#'   * `lhs`: the response expression as written (e.g. `log(y)`, without the
#'     `| trials(n)` part),
#'   * `responseExpr`: the response expression as a call / name,
#'   * `aux`: named list of the family inputs of the left-hand side
#'     (`trials`, `exposure`, `event`, `weights`; expressions), empty if none,
#'   * `part`: `"out"` or `"long"`,
#'   * `formula`: the formula as given by the user,
#'   * `variables`: list with character vectors `response` (variables in the
#'     left-hand side, e.g. `y` for `log(y)`), `covariate` (variables in fixed-effect
#'     terms and on the left of `|` in random-effect terms, e.g. `t` for `I(t^2)`)
#'     `group` (variables on the right of `|`, e.g. `id`, or `a` and `b` for `a:b`)
#'     and `aux` (variables of `trials()` / `exposure()` / `event()` / `weights()`),
#'   * `fixedFormula`: the fixed-effects part only (no random effects, no `scaled()`),
#'   * `lme4Formula`: formula with `re()` replaced by lme4 bar terms and `scaled()` removed,
#'   * `re`: list of random-effect blocks, each with `term` (the lme4 bar expression),
#'     `group` (grouping factor name), `correlated` (`FALSE` for `||`) and
#'     `params` (named list of extra arguments, evaluated in the formula environment)
#'     and `prior` (the `prior` argument, see [priors], or `NULL`),
#'   * `scaled`: (outcome formulae only) list of the formula's `scaled()` terms,
#'     each with `target` (longitudinal response name), `params` and `prior`
#'     (a `normalCoef()` prior of the association coefficient, or `NULL`).
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
    # B3: any number of brackets around a bar term is a random-effect block (as in lme4)
    inner <- stripParens(e)
    if (is.call(e) && identical(e[[1L]], quote(re))) {
      blk <- parseRe(e, env, response)
      re[[length(re) + 1L]] <- blk
      keep[[length(keep) + 1L]] <- list(expr = call("(", lme4Bar(blk$term)), sign = tm$sign)
    } else if (is.call(e) && identical(e[[1L]], quote(scaled))) {
      if (part == "long")
        stop("scaled() is only allowed in outcome formulae (found in formula for '",
             response, "').", call. = FALSE)
      scaled[[length(scaled) + 1L]] <- parseScaled(e, env, response)
    } else if (isBar(inner)) {
      re[[length(re) + 1L]] <- makeReBlock(inner, list())
      keep[[length(keep) + 1L]] <- list(expr = call("(", lme4Bar(inner)), sign = tm$sign)
    } else {
      # B2: only nested CALLS re(...) / scaled(...) are errors; a variable named `re` is a covariate
      if (hasCall(e, c("re", "scaled")))
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

# Family inputs allowed on the left-hand side: y | trials(n) + exposure(E),
# time | event(status) for survival families, y | weights(w) (observation
# weights: INLA `scale`, families that honour it).
auxTerms <- c("trials", "exposure", "event", "weights")

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

# Remove any number of enclosing parentheses: ((1 | id)) -> 1 | id.
stripParens <- function(e) {
  while (is.call(e) && identical(e[[1L]], quote(`(`)) && length(e) == 2L) e <- e[[2L]]
  e
}

# TRUE if the expression contains a call to one of `fns` (a variable of that name is not a call).
hasCall <- function(e, fns) {
  if (!is.call(e)) return(FALSE)
  if (is.name(e[[1L]]) && as.character(e[[1L]]) %in% fns) return(TRUE)
  any(vapply(as.list(e), hasCall, NA, fns = fns))
}

# The bar term for lme4: an intercept-only `||` block becomes `|` (the same model), because
# lme4's double-bar expansion of (1 || id) yields the unparsable `(1|id)+(0+|id)` (B1).
lme4Bar <- function(bar) {
  if (!identical(bar[[1L]], quote(`||`))) return(bar)
  lhsTerms <- stats::terms(stats::as.formula(call("~", bar[[2L]])))
  if (length(attr(lhsTerms, "term.labels"))) return(bar)
  bar[[1L]] <- quote(`|`)
  bar
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
  if (length(args)) args[[1L]] <- stripParens(args[[1L]])  # re((1 | id)) as in lme4 (B3)
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
  if (!is.null(prior) && !(inherits(prior, "mvremPrior") && prior$target != "coef"))
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
  prior <- params[["prior"]]
  params <- params[names(params) != "prior"]
  if (!length(params)) params <- list()
  if (!is.null(prior) && !(inherits(prior, "mvremPrior") && prior$target == "coef"))
    stop("`prior` in scaled() must be normalCoef(); raw INLA settings go in `hyper` ",
         "(formula for '", response, "').", call. = FALSE)
  list(target = as.character(x), params = params, prior = prior)
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
