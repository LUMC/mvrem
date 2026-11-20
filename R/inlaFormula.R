#' Naming scheme for INLA formula variables
#'
#' Templates used to name the columns, latent fields and index values of the
#' INLA formula. Placeholders in braces are filled in by [buildInlaFormula()]:
#' `{resp}` response, `{term}` design-matrix term (intercept shown as
#' `intercept`), `{terms}` terms of a random-effect block joined by `sep`,
#' `{group}` grouping term, `{level}` grouping level, `{long}` / `{out}`
#' longitudinal / outcome response, `{vars}` grouping variables and covariates of
#' the longitudinal formula joined by `sep`, `{values}` their values joined by
#' `valueSep`. Names are made syntactic by replacing other characters by `.`.
#' A placeholder a template cannot use is an error. `interceptLabel` is not a
#' template: it is how the intercept is written in result tables and plots
#' (`"(Intercept)"` as base R; the INLA names use `intercept`, default `1`).
#'
#' The default reproduces names like `y_1`, `y_t_id`, `x1_TO_y_id_t`, `y_Ntrials`,
#' `y_FROM_x1_id_t` and index values like `"y|1|S001"`, `"x1|TO|y|S001|0.5"`.
#' The scheme in use can be set globally with `options(mvrem.naming = namingScheme(...))`.
#'
#' @param ... Templates overriding the defaults, e.g. `zeroWeight = "w_{long}"`.
#' @return A named list of templates.
#' @export
namingScheme <- function(...) {
  defaults <- list(
    responseObject = "MX",           # response matrix of the INLA model
    levelsObject = "LEVELS",         # list holding the `values` of the f() terms
    response = "{resp}_{resp}",      # response column
    zeroResponse = "{resp}_{resp}_0eq",
    trials = "{resp}_Ntrials",       # trials() column (INLA `Ntrials`)
    exposure = "{resp}_E",           # exposure() column (INLA `E`)
    event = "{resp}_event",          # event() column (status of an `inla.surv` response)
    weights = "{resp}_scale",        # weights() column (INLA `scale`)
    intercept = "1",                 # the intercept in INLA names (`y_1`)
    interceptLabel = "(Intercept)",  # the intercept in result tables and plots (not a template)
    fixed = "{resp}_{term}",         # fixed-effect covariate column (also RE weight)
    reIndex = "{resp}_{term}_{group}",
    reLevels = "{resp}_{terms}_{group}",
    reValue = "{resp}|{term}|{level}",
    zeroIndex = "{long}_TO_{out}_{vars}",
    zeroWeight = "{long}_WEIGHT",
    zeroValue = "{long}|TO|{out}|{values}",
    scaledIndex = "{out}_FROM_{long}_{vars}",
    sep = "_",
    valueSep = "|"
  )
  user <- list(...)
  bad <- setdiff(names(user), names(defaults))
  if (length(bad))
    stop("Unknown naming template(s): ", paste(bad, collapse = ", "), call. = FALSE)
  res <- utils::modifyList(defaults, user)
  for (tpl in names(res)) {
    if (!is.character(res[[tpl]]) || length(res[[tpl]]) != 1L || is.na(res[[tpl]]))
      stop("Naming template `", tpl, "` must be a single string.", call. = FALSE)
    allowed <- templatePlaceholders[[tpl]]
    found <- regmatches(res[[tpl]], gregexpr("\\{[^{}]*\\}", res[[tpl]]))[[1L]]
    unknown <- setdiff(found, paste0("{", allowed, "}"))
    if (length(unknown))
      stop("Unknown placeholder(s) ", paste(unknown, collapse = ", "), " in naming template `", tpl,
           "`", if (length(allowed)) paste0(" (allowed: ", paste0("{", allowed, "}", collapse = ", "), ")"),
           ".", call. = FALSE)
  }
  res
}

# The placeholders each template may use (what fillName() / fillValue() supply for it).
templatePlaceholders <- list(
  response = "resp", zeroResponse = "resp", trials = "resp", exposure = "resp", event = "resp",
  weights = "resp",
  fixed = c("resp", "term"), reIndex = c("resp", "term", "group"), reLevels = c("resp", "terms", "group"),
  reValue = c("resp", "term", "level"), zeroIndex = c("long", "out", "vars"), zeroWeight = "long",
  zeroValue = c("long", "out", "values"), scaledIndex = c("out", "long", "vars")
)

# The label of a random / fixed effect in result tables and plots: the intercept as
# `naming$interceptLabel` (default `(Intercept)`, as base R), other effects as their column.
effectLabel <- function(e, naming) ifelse(e == "(Intercept)", naming$interceptLabel, e)

#' Build the INLA formula of a joint model
#'
#' Translates the parsed formulae ([collectFormulas()]) and the prepared data
#' ([prepareData()]) into a single INLA formula, using the names given by
#' `naming` (see [namingScheme()]). Terms, per response:
#' * **fixed**: one covariate column per design-matrix column (`y_1`, `y_t`);
#' * **random**: per random-effect block; a single effect or uncorrelated
#'   effects (`||`) give one `iid` term per effect, correlated effects (`|`)
#'   give an `iid<k>d` term (`iidkd` for k > 5) for the first effect and `copy`
#'   terms for the others. Slopes use the fixed covariate column as weight.
#'   Indices are labels (`"y|t|S001"`) listed in `LEVELS`;
#' * **zero-equation** (longitudinal responses, one per linked outcome; an outcome
#'   is linked to the longitudinal responses of its `scaled()` terms): latent field `x1_TO_y_id_t` with weight `x1_WEIGHT`
#'   (equation `0 = eta_x1 - u`), indexed by the grouping variables and
#'   covariates of the longitudinal formula at the outcome rows; the precision
#'   of `u` is fixed and small (`initial = -6`, diffuse), so `u` follows the
#'   zero equations;
#' * **scaled** (outcome responses, one per `scaled()` term):
#'   `y_FROM_x1_id_t`, a scaled copy of the corresponding zero-equation field.
#'   The prior of the scaling coefficient is set by `scaled(x1, prior =
#'   normalCoef(...))` (default: `normalCoef(0, 10)`); other parameters by the
#'   further arguments of `scaled(x1, ...)`.
#'
#' A `prior` given in `re(...)` / `scaled(...)` (see [priors]) sets the `hyper`
#' argument of the term. Other parameters given in `re(...)` / `scaled(...)` override
#' the defaults (nested lists such as `hyper` are merged; a warning is given when
#' `hyper` replaces the distribution set by `prior`).
#'
#' The likelihood of each response, its link and hyperparameter priors come
#' from its family (see [families]), giving `families` and `controlFamily`.
#' `trials()` / `exposure()` / `event()` / `weights()` on the left-hand side add
#' a column to the response's stacked data, right after the response.
#'
#' @param forms Result of [collectFormulas()].
#' @param data Result of [prepareData()].
#' @param naming Naming scheme, see [namingScheme()].
#' @return A list with elements
#'   * `responses`: response column per response, `zeroResponses` per
#'     longitudinal response,
#'   * `families`: INLA likelihood family per response column (named by the
#'     columns of `responses`, then `zeroResponses`; zero equations are
#'     `"gaussian"`), in the order of the likelihoods of the INLA model,
#'   * `controlFamily`: INLA `control.family`, one entry per element of
#'     `families` (see [controlFamily()]); zero equations have a fixed
#'     precision (placeholder: log precision 10),
#'   * `inputs`: per response, the stack columns of its family inputs, named
#'     by INLA argument (`c(Ntrials = "y_Ntrials")`, `c(E = "z_E")`,
#'     `c(scale = "x1_scale")` for `weights()`; empty if
#'     none). Rows of other likelihoods (and zero equations) hold `NA` there,
#'     to be replaced by 1 when the INLA arguments are assembled,
#'   * `fixed`: per response, `formula` and `columns` (INLA name -> design-matrix column),
#'   * `random`: per response, a list of blocks with `group`, `lhs` (formula of
#'     the block's effects), `effects`
#'     (design-matrix columns), `index`, `weight` (`NA` for the intercept),
#'     `levels` (name in `LEVELS`) and `terms` (`f()` calls),
#'   * `zero`: per (long, out) pair, named `"<long>_<out>"`: `long`, `out`,
#'     `vars`, `index`, `weight`, `levels` and `term`,
#'   * `scaled`: per pair, named `"<out>_<long>"`: `out`, `long`, `index`, `term`,
#'   * `levels`: the `LEVELS` list (index values of all labelled terms),
#'   * `naming`: the naming scheme used,
#'   * `terms`: per response, the terms of its equation(s) in the formula,
#'   * `formulas`: the user's formulae,
#'   * `stacks`: per response, the stacked data with the formula's column names:
#'     `main` (observation rows) and, for longitudinal responses, one element
#'     per outcome (zero-equation rows: one per unique
#'     combination of the long formula's grouping variables and covariates in
#'     the outcome data; the outcome is matched to them by these values),
#'   * `formula`: the INLA formula; its environment holds `LEVELS`.
#' @export
buildInlaFormula <- function(forms, data, naming = getOption("mvrem.naming", namingScheme())) {
  nm <- naming
  groups <- groupLevels(data)
  responses <- vapply(bothParts(forms), function(p)
    fillName(nm$response, resp = p$response), "")
  # zero-equation likelihoods only for the long responses linked to an outcome
  linked <- unique(unlist(lapply(forms$out, scaledTargets), use.names = FALSE))
  zeroResponses <- vapply(forms$long[intersect(names(forms$long), linked)], function(p)
    fillName(nm$zeroResponse, resp = p$response), "")

  links <- linkTerms(forms, data, nm)
  eqs <- equationTerms(forms, data, groups, nm, links)
  inputs <- lapply(bothParts(forms), function(p) {
    a <- names(p$aux)
    stats::setNames(vapply(a, function(x) fillName(nm[[x]], resp = p$response), ""),
                    inputArguments[a])
  })
  checkInlaNames(nm, responses, zeroResponses, eqs$fixed, inputs, eqs$random, links$zero,
                 links$scaled)
  lik <- likelihoodSpec(forms, responses, zeroResponses)

  # The f() terms refer to their `values` as LEVELS$<index> (levelsRef()) instead of
  # embedding the level vectors in the formula: the formula stays printable and short,
  # and INLA evaluates LEVELS in the formula's environment, which holds the list.
  levels <- c(links$levels, eqs$levels)
  env <- new.env(parent = globalenv())
  assign(nm$levelsObject, levels, envir = env)
  rhs <- Reduce(function(a, b) call("+", a, b), unlist(eqs$terms, use.names = FALSE), quote(-1))

  res <- list(
    responses = responses,
    zeroResponses = zeroResponses,
    families = lik$families,
    controlFamily = lik$control,
    survival = lik$survival,
    inputs = inputs,
    fixed = eqs$fixed,
    random = eqs$random,
    zero = links$zero,
    scaled = links$scaled,
    levels = levels,
    groups = groups,
    naming = nm,
    terms = eqs$terms,
    formulas = lapply(bothParts(forms), `[[`, "formula"),
    formula = stats::as.formula(call("~", as.name(nm$responseObject), rhs), env = env)
  )
  res$stacks <- buildStacks(forms, data, res)
  res
}

# The terms of the links, one zero-equation field and one scaled copy per scaled() term
# (= per (outcome, long) pair): list(zero, scaled, levels), `levels` = the zero fields'
# LEVELS entries (unique combinations of the long formula's groups + covariates in the
# outcome data).
linkTerms <- function(forms, data, nm) {
  zero <- list()
  scaled <- list()
  levels <- list()
  for (o in forms$out) for (s in o$scaled) {
    l <- forms$long[[s$target]]
    vars <- c(l$variables$group, l$variables$covariate)
    varsTxt <- paste(syntactic(vars), collapse = nm$sep)
    index <- fillName(nm$zeroIndex, long = l$response, out = o$response, vars = varsTxt)
    weight <- fillName(nm$zeroWeight, long = l$response)
    levels[[index]] <- unique(zeroLabels(data$out[[o$response]], l$response, o$response,
                                         vars, nm))
    zero[[paste0(l$response, "_", o$response)]] <- list(
      long = l$response, out = o$response, vars = vars, index = index, weight = weight,
      levels = index,
      term = fTerm(index, weight,
                   list(model = "iid", values = levelsRef(nm, index), hyper = zeroFieldHyper),
                   list())
    )
    sIndex <- fillName(nm$scaledIndex, long = l$response, out = o$response, vars = varsTxt)
    coefPrior <- if (is.null(s$prior)) defaultCoefPrior() else s$prior
    if (!is.null(s$prior))
      warnHyperOverride(coefHyper(coefPrior), s$params[["hyper"]],
                        paste0("scaled(", l$response, ") in the formula for ", o$response))
    scaled[[paste0(o$response, "_", l$response)]] <- list(
      out = o$response, long = l$response, index = sIndex, prior = s$prior,
      term = fTerm(sIndex, NULL, list(copy = index, hyper = coefHyper(coefPrior)), s$params)
    )
  }
  list(zero = zero, scaled = scaled, levels = levels)
}

# Per response the fixed and random terms, and the formula terms in order: fixed columns,
# random f() terms, then the outcome's scaled copies / the marker's zero fields.
# list(fixed, random, terms, levels), `levels` = the random fields' LEVELS entries.
equationTerms <- function(forms, data, groups, nm, links) {
  fixed <- list()
  random <- list()
  terms <- list()
  levels <- list()
  for (p in bothParts(forms)) {
    r <- p$response
    df <- bothParts(data)[[r]]
    fixed[[r]] <- fixedTerms(p, df, nm)
    rnd <- randomTerms(p, df, groups, nm)
    for (b in rnd) levels[[b$levels]] <- b$values
    random[[r]] <- lapply(rnd, function(b) b[names(b) != "values"])
    linkCalls <- if (p$part == "out") Filter(function(s) s$out == r, links$scaled) else
      Filter(function(z) z$long == r, links$zero)
    terms[[r]] <- unname(c(lapply(names(fixed[[r]]$columns), as.name),
                           unlist(lapply(random[[r]], `[[`, "terms"), recursive = FALSE),
                           lapply(linkCalls, `[[`, "term")))
  }
  list(fixed = fixed, random = random, terms = terms, levels = levels)
}

# Every column / object name the formula and the stacks use must be unique: a collision
# with a response, zero-response or weight column silently overwrote data (B7).
checkInlaNames <- function(nm, responses, zeroResponses, fixed, inputs, random, zero, scaled) {
  zeroWeights <- vapply(zero, `[[`, "", "weight")[!duplicated(vapply(zero, `[[`, "", "long"))]
  allIdx <- c(nm$responseObject, nm$levelsObject, responses, zeroResponses,
              unlist(lapply(fixed, function(f) names(f$columns))),
              unlist(inputs, use.names = FALSE),
              unlist(lapply(random, function(r) lapply(r, `[[`, "index"))),
              vapply(zero, `[[`, "", "index"), unname(zeroWeights),
              vapply(scaled, `[[`, "", "index"))
  dup <- unique(allIdx[duplicated(allIdx)])
  if (length(dup))
    stop("Naming scheme produces duplicated INLA names: ", paste(dup, collapse = ", "),
         call. = FALSE)
  invisible()
}

# The INLA likelihoods, named by response column: the responses' families, then one
# Gaussian likelihood per zero-response column (fixed high precision, zeroControlFamily).
likelihoodSpec <- function(forms, responses, zeroResponses) {
  fams <- lapply(bothParts(forms), function(p)
    if (is.null(p$family)) defaultFamily() else p$family)
  list(
    families = c(stats::setNames(vapply(fams, `[[`, "", "family"), responses),
                 stats::setNames(rep("gaussian", length(zeroResponses)), zeroResponses)),
    control = c(stats::setNames(lapply(fams, controlFamily), responses),
                stats::setNames(rep(list(zeroControlFamily), length(zeroResponses)), zeroResponses)),
    survival = stats::setNames(vapply(fams, function(f) "event" %in% f$aux, NA), names(responses))
  )
}

# Prior of the zero-equation field u: fixed, small precision (diffuse), so that u
# is determined by the zero equations (an estimated precision shrinks u, and with
# it the long predictor, to 0).
zeroFieldHyper <- list(prec = list(initial = -6, fixed = TRUE))

# Labels of the zero-equation field at the rows of df.
zeroLabels <- function(df, long, out, vars, nm) {
  valueTxt <- do.call(paste, c(lapply(vars, function(v) as.character(df[[v]])),
                               sep = nm$valueSep))
  fillValue(nm$zeroValue, long = long, out = out, values = valueTxt)
}

# Values of a grouping term (e.g. "id", "site:id") at the rows of df.
groupValues <- function(df, group) {
  if (!is.null(df[[group]])) return(as.character(df[[group]]))
  vars <- all.vars(str2lang(group))
  do.call(paste, c(lapply(vars, function(v) as.character(df[[v]])), sep = ":"))
}

# Design matrix of the right-hand side of f at the rows of df, with factor
# levels taken from ref (the data the formula belongs to).
designMatrix <- function(f, df, ref = df) {
  tt <- stats::delete.response(stats::terms(f))
  xlev <- stats::.getXlevels(tt, stats::model.frame(tt, ref))
  stats::model.matrix(tt, stats::model.frame(tt, df, xlev = xlev, na.action = stats::na.pass))
}

# Columns of one response's equation at the rows of df: fixed covariates,
# random-effect indices and their weights.
equationColumns <- function(p, inl, df, ref) {
  nm <- inl$naming
  r <- p$response
  X <- designMatrix(p$fixedFormula, df, ref)
  cols <- as.data.frame(X[, inl$fixed[[r]]$columns, drop = FALSE], optional = TRUE)
  names(cols) <- names(inl$fixed[[r]]$columns)
  for (b in inl$random[[r]]) {
    Z <- designMatrix(b$lhs, df, ref)
    g <- groupValues(df, b$group)
    for (e in seq_along(b$effects)) {
      label <- termLabel(b$effects[e], nm)
      cols[[b$index[e]]] <- fillValue(nm$reValue, resp = r, term = label, level = g)
      if (!is.na(b$weight[e])) cols[[b$weight[e]]] <- Z[, b$effects[e]]
    }
  }
  cols
}

# Family input columns (trials(), exposure(), event(), weights()) of a response at the rows of df.
inputColumns <- function(p, inl, df) {
  cols <- inl$inputs[[p$response]]
  res <- data.frame(row.names = seq_len(nrow(df)))
  for (a in names(cols)) {
    aux <- p$aux[[names(inputArguments)[inputArguments == a]]]
    res[[cols[[a]]]] <- as.numeric(eval(aux, df, environment(p$formula)))
  }
  res
}

# Left-hand-side inputs -> names of the stack columns in `inl$inputs` (the
# INLA argument for trials / exposure / weights; `event` is part of the inla.surv response).
inputArguments <- c(trials = "Ntrials", exposure = "E", event = "event", weights = "scale")

# Stacked data of each response's equations, with the column names used in
# the INLA formula. Outcome stacks: response, equation columns and scaled
# indices. Longitudinal stacks: the observation rows followed by the
# zero-equation rows (one set per outcome).
buildStacks <- function(forms, data, inl) {
  nm <- inl$naming
  stacks <- list()
  for (p in forms$out) {
    r <- p$response
    df <- data$out[[r]]
    st <- data.frame(eval(p$responseExpr, df, environment(p$formula)))
    names(st) <- inl$responses[[r]]
    st <- cbind(st, inputColumns(p, inl, df), equationColumns(p, inl, df, df))
    for (sc in Filter(function(s) s$out == r, inl$scaled)) {
      z <- inl$zero[[paste0(sc$long, "_", r)]]
      st[[sc$index]] <- zeroLabels(df, sc$long, r, z$vars, nm)
    }
    stacks[[r]] <- list(main = st)
  }
  for (p in forms$long) {
    r <- p$response
    ref <- data$long[[r]]
    zs <- Filter(function(z) z$long == r, inl$zero)
    st <- data.frame(eval(p$responseExpr, ref, environment(p$formula)))
    names(st) <- inl$responses[[r]]
    if (length(zs)) st[[inl$zeroResponses[[r]]]] <- NA_real_
    st <- cbind(st, inputColumns(p, inl, ref), equationColumns(p, inl, ref, ref))
    if (length(zs)) st[[zs[[1L]]$weight]] <- NA_real_
    for (z in zs) st[[z$index]] <- NA_character_
    parts <- list(main = st)
    for (z in zs) {
      # one zero equation per unique combination of the matching variables
      # (groups + covariates of the long formula) present in the outcome data
      od <- data$out[[z$out]]
      od <- od[!duplicated(zeroLabels(od, r, z$out, z$vars, nm)), , drop = FALSE]
      zr <- data.frame(NA_real_, 0, row.names = NULL)[rep(1L, nrow(od)), , drop = FALSE]
      names(zr) <- c(inl$responses[[r]], inl$zeroResponses[[r]])
      for (col in inl$inputs[[r]]) zr[[col]] <- NA_real_
      zr <- cbind(zr, equationColumns(p, inl, od, ref))
      zr[[z$weight]] <- -1
      for (z2 in zs) zr[[z2$index]] <- NA_character_
      zr[[z$index]] <- zeroLabels(od, r, z$out, z$vars, nm)
      rownames(zr) <- NULL
      parts[[z$out]] <- zr
    }
    stacks[[r]] <- parts
  }
  stacks
}

# Union of the levels of each grouping factor over all prepared data frames.
groupLevels <- function(data) {
  res <- list()
  for (df in bothParts(data)) for (g in names(df)[vapply(df, is.factor, TRUE)]) {
    res[[g]] <- union(res[[g]], levels(df[[g]]))
  }
  lapply(res, sort)
}

syntactic <- function(x) gsub("[^A-Za-z0-9_.]+", ".", x)

# Fill {placeholders}; vector arguments give vector results.
fillValue <- function(template, ...) {
  args <- list(...)
  n <- max(lengths(args), 1L)
  res <- rep_len(template, n)
  for (a in names(args)) {
    v <- rep_len(as.character(args[[a]]), n)
    res <- vapply(seq_len(n), function(i)
      gsub(paste0("{", a, "}"), v[i], res[i], fixed = TRUE), "")
  }
  res
}

fillName <- function(template, ...) syntactic(fillValue(template, ...))

termLabel <- function(col, nm) ifelse(col == "(Intercept)", nm$intercept, col)

levelsRef <- function(nm, name) call("$", as.name(nm$levelsObject), as.name(name))

fixedTerms <- function(p, df, nm) {
  f <- p$fixedFormula
  X <- stats::model.matrix(stats::delete.response(stats::terms(f)), df)
  cols <- colnames(X)
  list(formula = f,
       columns = stats::setNames(cols, fillName(nm$fixed, resp = p$response,
                                                term = termLabel(cols, nm))))
}

randomTerms <- function(p, df, groups, nm) {
  r <- p$response
  blocks <- lapply(p$re, function(blk) {
    lhs <- stats::as.formula(call("~", blk$term[[2L]]), env = environment(p$formula))
    effects <- colnames(stats::model.matrix(lhs, df))
    labels <- termLabel(effects, nm)
    k <- length(effects)
    lev <- groups[[blk$group]]
    index <- fillName(nm$reIndex, resp = r, term = labels, group = blk$group)
    weight <- ifelse(effects == "(Intercept)", NA_character_,
                     fillName(nm$fixed, resp = r, term = labels))
    valuesOf <- function(e) fillValue(nm$reValue, resp = r, term = labels[e], level = lev)
    checkBlockPrior(blk$prior, k, blk$correlated, blk$term)
    what <- paste0("block `", deparse1(blk$term), "`")

    if (k == 1L || !blk$correlated) {
      # one prior per effect: a data-scaled prior differs between the effects of `||`
      priors <- resolveEffectPriors(blk$prior, p, df, effects, blk)
      hypers <- lapply(priors, function(pr) if (!is.null(pr)) list(hyper = priorHyper(pr, 1L)))
      warnHyperOverride(hypers[[1L]]$hyper, blk$params[["hyper"]], what)
      levName <- vapply(seq_len(k), function(e)
        fillName(nm$reLevels, resp = r, terms = labels[e], group = blk$group), "")
      terms <- lapply(seq_len(k), function(e)
        fTerm(index[e], weight[e],
              c(list(model = "iid", values = levelsRef(nm, levName[e]),
                     n = call("length", levelsRef(nm, levName[e]))), hypers[[e]]),
              blk$params))
      values <- lapply(seq_len(k), valuesOf)
      out <- lapply(seq_len(k), function(e)
        list(group = blk$group, lhs = lhs, effects = effects[e], index = index[e],
             weight = weight[e],
             levels = levName[e], terms = terms[e], values = values[[e]], prior = priors[[e]]))
      return(out)
    }

    prior <- resolveBlockPrior(blk$prior, p, df, effects, blk)
    hyper <- if (!is.null(prior)) list(hyper = priorHyper(prior, k))
    warnHyperOverride(hyper$hyper, blk$params[["hyper"]], what)

    levName <- fillName(nm$reLevels, resp = r, terms = paste(labels, collapse = nm$sep),
                        group = blk$group)
    model <- if (k <= 5L) paste0("iid", k, "d") else "iidkd"
    defaults <- list(model = model)
    if (model == "iidkd") defaults$order <- k
    defaults <- c(defaults, list(values = levelsRef(nm, levName),
                                 n = call("length", levelsRef(nm, levName))), hyper)
    first <- fTerm(index[1L], weight[1L], defaults, blk$params)
    copies <- lapply(seq_len(k)[-1L], function(e)
      fTerm(index[e], weight[e], list(copy = index[1L]), list()))
    list(list(group = blk$group, lhs = lhs, effects = effects, index = index, weight = weight,
              levels = levName, terms = c(list(first), copies),
              values = unlist(lapply(seq_len(k), valuesOf)), prior = prior))
  })
  unlist(blocks, recursive = FALSE)
}

# Warn when a user `hyper` replaces the distribution set by `prior`; what =
# description of the block / family.
warnHyperOverride <- function(fromPrior, user, what) {
  if (is.null(fromPrior) || !is.list(user)) return(invisible())
  over <- names(fromPrior)[vapply(names(fromPrior), function(h)
    is.list(user[[h]]) && any(c("prior", "param") %in% names(user[[h]])), NA)]
  if (length(over))
    warning("`hyper` overrides the `prior` of ", what, " (",
            paste(over, collapse = ", "), ").", call. = FALSE)
  invisible()
}

# f(index, [weight], ...) call; user parameters override defaults by name.
# Values (lists, vectors) are embedded as the calls that build them: terms()
# deparses embedded objects without their names (`list(list("pc.prec", ...))`),
# so INLA would silently ignore e.g. `hyper`.
fTerm <- function(index, weight, defaults, params) {
  args <- utils::modifyList(defaults, params, keep.null = TRUE)
  args <- lapply(args, function(a)
    if (is.language(a) || (is.atomic(a) && length(a) == 1L && is.null(attributes(a)))) a
    else str2lang(deparse1(a)))
  as.call(c(list(quote(f), as.name(index)),
            if (!is.null(weight) && !is.na(weight)) list(as.name(weight)),
            args))
}
