#' Likelihood families
#'
#' Constructors of the likelihood family of a response, given in `family`
#' (or `familyOut` / `familyLong`) of [collectFormulas()], e.g.
#' `family = list(y = betaFamily(phi = gammaPrec(1, 0.1)))`. Each sets the
#' INLA likelihood, its link and the priors of its hyperparameters, i.e. the
#' response's entry of INLA's `control.family`.
#'
#' | Constructor | INLA likelihood | Response | Hyperparameters | LHS input |
#' |---|---|---|---|---|
#' | `gaussianFamily()` | `gaussian` | real | `prec` (residual precision \eqn{\tau_\varepsilon}) | `weights()` |
#' | `studentTFamily()` | `t` | real | `prec`, `dof` (\eqn{\nu}) | `weights()` |
#' | `lognormalFamily()` | `lognormal` | positive | `prec` (of \eqn{\log y}) | |
#' | `gammaFamily()` | `gamma` | positive | `prec` (\eqn{\phi}, Var = \eqn{\mu^2/\phi}) | `weights()` |
#' | `betaFamily()` | `beta` | in (0, 1) | `phi` (\eqn{\phi}, Var = \eqn{\mu(1-\mu)/(1+\phi)}) | `weights()` |
#' | `poissonFamily()` | `poisson` | count | | `exposure()` |
#' | `negBinomialFamily()` | `nbinomial` | count | `size` (\eqn{\kappa}, Var = \eqn{\mu + \mu^2/\kappa}) | `exposure()` |
#' | `binomialFamily()` | `binomial` | 0, ..., n | | `trials()` (default 1) |
#' | `betaBinomialFamily()` | `betabinomial` | 0, ..., n | `rho` (overdispersion) | `trials()` (default 1) |
#' | `zeroInflatedPoissonFamily()` | `zeroinflatedpoisson<type>` | count | `prob` (\eqn{\pi}) | `exposure()` |
#' | `zeroInflatedNegBinomialFamily()` | `zeroinflatednbinomial<type>` | count | `size`, `prob` | `exposure()` |
#' | `weibullSurvFamily()` | `weibullsurv` | time >= 0 (0 only censored) | `alpha` (shape \eqn{\alpha}) | `event()` (required) |
#' | `exponentialSurvFamily()` | `exponentialsurv` | time >= 0 (0 only censored) | | `event()` (required) |
#' | `inlaFamily()` | any (non-survival) | | any | `trials()`, `exposure()`, `weights()` |
#'
#' **Survival outcomes.** `weibullSurvFamily()` and `exponentialSurvFamily()`
#' model a time to event with right censoring, written as
#' `time | event(status) ~ ...` (`status` 1 = the event was observed at `time`,
#' 0 = censored at `time`), one row per subject. INLA's `weibullsurv` (variant 0)
#' has hazard \eqn{h(t) = \alpha \lambda t^{\alpha - 1}} and survival
#' \eqn{S(t) = \exp(-\lambda t^\alpha)} with \eqn{\log \lambda} the linear
#' predictor (the exponential is \eqn{\alpha = 1}), so an association
#' coefficient is a log hazard ratio. The marker enters the hazard at the values
#' of its covariates given in the outcome row (e.g. `t = time`). Survival
#' responses are not plotted, have no residuals and cannot be predicted; their
#' fitted value ([fittedValues()]) is \eqn{\lambda} (the hazard scale), not a
#' mean time. A time of 0 is accepted for a censored row only.
#'
#' **Hyperparameter priors.** A hyperparameter argument is `NULL` (INLA's
#' default prior), a prior constructor (see [priors]) or a raw INLA `hyper`
#' entry, e.g. `list(prior = "loggamma", param = c(1, 0.01))` or
#' `list(initial = log(100), fixed = TRUE)`. The constructors are accepted for
#' precision-like parameters (`prec`, `phi`: [gammaPrec()], [pcSd()],
#' [halfNormalSd()], [halfCauchySd()], placing the prior on the parameter
#' \eqn{\tau} or on \eqn{\tau^{-1/2}}) and, for `size`, [gammaPrec()]
#' (\eqn{\kappa \sim \mathrm{Gamma}}); other parameters take raw INLA entries.
#'
#' **Inputs on the left-hand side.** Binomial-type families take the number of
#' trials and count families an exposure (INLA's `Ntrials` / `E`), given in the
#' formula as in brms: `y | trials(n) ~ ...`, `y | exposure(pt) ~ ...` (any
#' expression of the data in the parentheses). Without `trials()` the binomial
#' is a Bernoulli (\eqn{n = 1}); without `exposure()`, \eqn{E = 1}.
#'
#' **Observation weights.** `y | weights(w) ~ ...` (INLA's `scale`; positive
#' and finite) multiplies the precision of observation \eqn{i} by \eqn{w_i}:
#' gaussian and Student-t \eqn{\varepsilon_i \sim N(0, (w_i \tau_\varepsilon)^{-1})}
#' (e.g. \eqn{w_i} = the inverse of a known relative measurement variance), gamma
#' and beta \eqn{\phi_i = w_i \phi}. These are the likelihoods whose INLA
#' documentation lists `scale` (`INLA::inla.doc("gaussian")`, `"T"`, `"gamma"`,
#' `"beta"`); `lognormal` does not document it, so `lognormalFamily()` does not
#' take `weights()`. `inlaFamily()` passes `weights()` on unchecked: INLA ignores
#' `scale` for likelihoods that do not use it. Without `weights()`, \eqn{w_i = 1}.
#' These are not likelihood weights (INLA's `weights`, which multiply the
#' log-likelihood). They can be combined with another input of the same family
#' (`inlaFamily()`: `y | trials(n) + weights(w) ~ ...`).
#'
#' Responses are checked against the support of the family by [prepareData()].
#'
#' @param prec,phi,size,rho,prob,dof,alpha Prior of the hyperparameter (see Details;
#'   `alpha`: raw INLA entry only, INLA's default is the `pc.alphaw` prior).
#' @param link Link function (INLA `control.link` model), e.g. `"logit"`,
#'   `"probit"`, `"cloglog"`. Validated against INLA's list when INLA is installed.
#' @param type Zero-inflation type: `1` (\eqn{y = 0} with probability \eqn{\pi},
#'   otherwise from the count distribution) or `0` (\eqn{P(y = 0) = \pi},
#'   positive values from the zero-truncated count distribution).
#' @param family For `inlaFamily()`: name of an INLA likelihood (see
#'   [inlaFamilies()]); survival likelihoods only through `weibullSurvFamily()`
#'   / `exponentialSurvFamily()`.
#' @param hyper For `inlaFamily()`: named list of hyperparameter priors, keyed
#'   by INLA's names (e.g. `prec`, or `theta1`).
#' @param control Further `control.family` entries (named list), merged into
#'   the ones set by the other arguments, e.g. `list(beta.censor.value = 1e-4)`.
#' @return An object of class `mvremFamily`.
#' @examples
#' betaFamily(phi = gammaPrec(1, 0.1))
#' binomialFamily(link = "probit")
#' gaussianFamily(prec = pcSd(u = 1))
#' @name families
NULL

#' @rdname families
#' @export
gaussianFamily <- function(prec = NULL, link = "identity", control = list()) {
  newFamily("gaussian", "gaussian", link, "identity",
            list(prec = prec), list(prec = hyperSpec("\\tau^{#r#}_{\\varepsilon}")),
            control, aux = character(), support = "real", weights = TRUE)
}

#' @rdname families
#' @export
studentTFamily <- function(prec = NULL, dof = NULL, control = list()) {
  newFamily("t", "t", "identity", "identity",
            list(prec = prec, dof = dof),
            list(prec = hyperSpec("\\tau^{#r#}_{\\varepsilon}"),
                 dof = hyperSpec("\\nu^{#r#}", character())),
            control, aux = character(), support = "real", weights = TRUE)
}

#' @rdname families
#' @export
lognormalFamily <- function(prec = NULL, control = list()) {
  newFamily("lognormal", "lognormal", "identity", "identity",
            list(prec = prec), list(prec = hyperSpec("\\tau^{#r#}_{\\varepsilon}")),
            control, aux = character(), support = "positive")
}

#' @rdname families
#' @export
gammaFamily <- function(prec = NULL, link = "log", control = list()) {
  newFamily("gamma", "gamma", link, "log",
            list(prec = prec), list(prec = hyperSpec("\\phi^{#r#}")),
            control, aux = character(), support = "positive", weights = TRUE)
}

#' @rdname families
#' @export
betaFamily <- function(phi = NULL, link = "logit", control = list()) {
  newFamily("beta", "beta", link, "logit",
            list(phi = phi), list(phi = hyperSpec("\\phi^{#r#}")),
            control, aux = character(), support = "unit", weights = TRUE)
}

#' @rdname families
#' @export
poissonFamily <- function(link = "log", control = list()) {
  newFamily("poisson", "poisson", link, "log", list(), list(),
            control, aux = "exposure", support = "count")
}

#' @rdname families
#' @export
negBinomialFamily <- function(size = NULL, link = "log", control = list()) {
  newFamily("nbinomial", "nbinomial", link, "log",
            list(size = size), list(size = hyperSpec("\\kappa^{#r#}", "gamma")),
            control, aux = "exposure", support = "count")
}

#' @rdname families
#' @export
binomialFamily <- function(link = "logit", control = list()) {
  newFamily("binomial", "binomial", link, "logit", list(), list(),
            control, aux = "trials", support = "binomial")
}

#' @rdname families
#' @export
betaBinomialFamily <- function(rho = NULL, link = "logit", control = list()) {
  newFamily("betabinomial", "betabinomial", link, "logit",
            list(rho = rho), list(rho = hyperSpec("\\rho^{#r#}", character())),
            control, aux = "trials", support = "binomial")
}

#' @rdname families
#' @export
zeroInflatedPoissonFamily <- function(type = 1, prob = NULL, control = list()) {
  type <- zeroInflationType(type)
  fam <- newFamily(paste0("zeroinflatedpoisson", type), "zip", "log", "log",
                   list(prob = prob), list(prob = hyperSpec("\\pi^{#r#}", character())),
                   control, aux = "exposure", support = "count")
  fam$type <- type
  fam
}

#' @rdname families
#' @export
zeroInflatedNegBinomialFamily <- function(type = 1, size = NULL, prob = NULL,
                                          control = list()) {
  type <- zeroInflationType(type)
  fam <- newFamily(paste0("zeroinflatednbinomial", type), "zinb", "log", "log",
                   list(size = size, prob = prob),
                   list(size = hyperSpec("\\kappa^{#r#}", "gamma"),
                        prob = hyperSpec("\\pi^{#r#}", character())),
                   control, aux = "exposure", support = "count")
  fam$type <- type
  fam
}

#' @rdname families
#' @export
weibullSurvFamily <- function(alpha = NULL, link = "log", control = list()) {
  newFamily("weibullsurv", "weibull", link, "log",
            list(alpha = alpha), list(alpha = hyperSpec("\\alpha^{#r#}", character())),
            control, aux = "event", support = "survival")
}

#' @rdname families
#' @export
exponentialSurvFamily <- function(link = "log", control = list()) {
  newFamily("exponentialsurv", "exponential", link, "log", list(), list(),
            control, aux = "event", support = "survival")
}

#' @rdname families
#' @export
inlaFamily <- function(family, link = NULL, hyper = list(), control = list()) {
  if (!isString(family))
    stop("`family` must be a single INLA likelihood name.", call. = FALSE)
  checkInlaFamily(family, "inlaFamily()")
  if (!is.list(hyper) || (length(hyper) && (is.null(names(hyper)) || any(names(hyper) == ""))))
    stop("`hyper` must be a named list.", call. = FALSE)
  specs <- lapply(names(hyper), function(h)
    hyperSpec(paste0("\\theta^{#r#}_{\\mathrm{", escapeLatex(h), "}}")))
  newFamily(family, "inla", link, NULL, hyper, stats::setNames(specs, names(hyper)),
            control, aux = c("trials", "exposure"), support = "any", weights = TRUE)
}

zeroInflationType <- function(type) {
  if (!(identical(type, 0) || identical(type, 1) || identical(type, 0L) || identical(type, 1L)))
    stop("`type` must be 0 or 1.", call. = FALSE)
  as.integer(type)
}

# Scalar prior constructors accepted for precision-like hyperparameters.
precDists <- c("gamma", "pc", "halfNormal", "halfCauchy")

# LaTeX symbol (with tagPlaceholder for the response tag) and accepted constructors.
hyperSpec <- function(symbol, dists = precDists) list(symbol = symbol, dists = dists)

# Placeholder of the response tag "(r)" in the symbols of hyperSpec(), filled in by
# hyperSymbol(). It cannot occur otherwise in a symbol: "#" is TeX's macro-parameter
# character, and escapeLatex() writes a "#" in a user's name (inlaFamily() hyper names)
# as "\#", so "r#" never appears in escaped text.
tagPlaceholder <- "#r#"
hyperSymbol <- function(fam, h, tag) gsub(tagPlaceholder, tag, fam$specs[[h]]$symbol, fixed = TRUE)

hyperFields <- c("prior", "param", "initial", "fixed", "to.theta", "from.theta")

# aux: LHS inputs the family takes (trials / exposure / event); weights: whether it takes
# weights() (INLA `scale`; kept apart from `aux` as it is a precision factor, not an input
# of the distribution).
newFamily <- function(family, kind, link, linkDefault, hyper, specs, control, aux, support,
                      weights = FALSE) {
  if (!is.null(link)) {
    if (!isString(link))
      stop("`link` must be a single character string.", call. = FALSE)
    known <- familyLinks(family)
    if (!is.null(known) && !link %in% known)
      stop("Link '", link, "' is not available for INLA likelihood '", family, "'; ",
           "available: ", paste(setdiff(known, "default"), collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(control) || (length(control) && (is.null(names(control)) || any(names(control) == ""))))
    stop("`control` must be a named list of control.family entries.", call. = FALSE)
  hyper <- hyper[!vapply(hyper, is.null, NA)]
  for (h in names(hyper)) checkHyper(hyper[[h]], h, specs[[h]]$dists, family)
  structure(list(family = family, kind = kind, link = link, linkDefault = linkDefault,
                 hyper = hyper, specs = specs, control = control, aux = aux,
                 support = support, weights = weights),
            class = "mvremFamily")
}

checkHyper <- function(x, name, dists, family) {
  what <- paste0("Hyperparameter `", name, "` of family '", family, "'")
  if (inherits(x, "mvremPrior")) {
    if (isScalarDataPrior(x))
      stop(what, ": data-scaled priors (\"data\") are for random-effect blocks only; ",
           "give a number.", call. = FALSE)
    if (!x$dist %in% dists)
      stop(what, " does not accept this prior constructor; ",
           if (length(dists)) paste0("use ", paste0(priorConstructor(dists), "()", collapse = ", "),
                                     " or ") else "use ",
           "a raw INLA entry, e.g. list(prior = ..., param = ...).", call. = FALSE)
  } else if (is.list(x) && length(x) && !is.null(names(x)) && all(names(x) %in% hyperFields)) {
    return(invisible())
  } else {
    stop(what, " must be NULL, a prior constructor (see ?priors) or a list with ",
         "elements among ", paste(hyperFields, collapse = ", "), ".", call. = FALSE)
  }
  invisible()
}

priorConstructor <- function(dist) {
  c(gamma = "gammaPrec", pc = "pcSd", halfNormal = "halfNormalSd",
    halfCauchy = "halfCauchySd", wishart = "wishartPrec")[dist]
}

# Is INLA installed? The single check of the package (tests mock it to cover the
# no-INLA paths: validation warnings, `fit = FALSE` only).
hasInla <- function() requireNamespace("INLA", quietly = TRUE)

# Links INLA knows for a likelihood (NULL when INLA is not installed).
familyLinks <- function(family) {
  if (!hasInla()) return(NULL)
  INLA::inla.models()$likelihood[[family]]$link
}

checkInlaFamily <- function(family, where) {
  if (!hasInla()) {
    warning("INLA is not installed: likelihood family '", family, "' in ", where,
            " is not validated.", call. = FALSE)
    return(invisible())
  }
  lik <- INLA::inla.models()$likelihood
  if (!family %in% names(lik))
    stop("Unknown INLA likelihood family in ", where, ": ", family, ". See inlaFamilies().",
         call. = FALSE)
  if (isTRUE(lik[[family]]$survival))
    stop("Survival likelihood '", family, "' (", where, ") is not supported; use ",
         "weibullSurvFamily() or exponentialSurvFamily() with `time | event(status) ~ ...`.",
         call. = FALSE)
  invisible()
}

# Family object of an INLA likelihood name (the constructor with its defaults).
familyFromName <- function(name, where) {
  switch(name,
    gaussian = gaussianFamily(),
    t = studentTFamily(),
    lognormal = lognormalFamily(),
    gamma = gammaFamily(),
    beta = betaFamily(),
    poisson = poissonFamily(),
    nbinomial = negBinomialFamily(),
    binomial = binomialFamily(),
    betabinomial = betaBinomialFamily(),
    zeroinflatedpoisson0 = zeroInflatedPoissonFamily(0),
    zeroinflatedpoisson1 = zeroInflatedPoissonFamily(1),
    zeroinflatednbinomial0 = zeroInflatedNegBinomialFamily(0),
    zeroinflatednbinomial1 = zeroInflatedNegBinomialFamily(1),
    weibullsurv = weibullSurvFamily(),
    exponentialsurv = exponentialSurvFamily(),
    {
      checkInlaFamily(name, where)
      suppressWarnings(inlaFamily(name))
    }
  )
}

defaultFamily <- function() gaussianFamily()

#' INLA `control.family` entry of a family
#'
#' @param x An `mvremFamily` object (see [families]).
#' @return A named list: `control.link` (when the link differs from INLA's
#'   default), `hyper` (hyperparameters with a prior) and the `control`
#'   entries of the family.
#' @examples
#' controlFamily(betaFamily(phi = gammaPrec(1, 0.1), link = "probit"))
#' @export
controlFamily <- function(x) {
  if (!inherits(x, "mvremFamily"))
    stop("`x` must be a family object, see ?families.", call. = FALSE)
  ctrl <- list()
  if (!is.null(x$link) && !identical(x$link, x$linkDefault))
    ctrl$control.link <- list(model = x$link)
  if (length(x$hyper))
    ctrl$hyper <- lapply(x$hyper, function(h)
      if (inherits(h, "mvremPrior")) priorHyper(h, 1L)[[1L]] else h)
  warnHyperOverride(ctrl$hyper, x$control$hyper, paste0("family '", x$family, "'"))
  utils::modifyList(ctrl, x$control)
}

# control.family of the zero equations: Gaussian with fixed, large precision
# (placeholder value, log precision 10).
zeroControlFamily <- list(hyper = list(prec = list(initial = 10, fixed = TRUE)))

#' @export
format.mvremFamily <- function(x, ...) {
  paste0("<mvrem family> ", x$family,
         if (!is.null(x$link)) paste0(", link = ", x$link),
         "; control.family = ", deparse1(controlFamily(x)))
}

#' @export
print.mvremFamily <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# Check the response (and trials / exposure / event / weights) values of a formula against the
# support of its family; df = prepared data of the formula.
checkResponse <- function(p, df) {
  fam <- p$family
  if (is.null(fam) || !nrow(df)) return(invisible())
  env <- environment(p$formula)
  isInt <- function(v) abs(v - round(v)) < 1e-8
  check <- function(what, v, ok, msg) {
    if (!is.numeric(v) && !is.logical(v)) stop(what, " must be numeric.", call. = FALSE)
    ok[is.na(ok)] <- TRUE
    # Inf / NaN give NA in the integer test and must not pass as missing values (B8)
    ok[is.infinite(v) | is.nan(v)] <- FALSE
    if (!all(ok)) stop(what, " ", msg, " (", sum(!ok), " value(s) not).", call. = FALSE)
  }
  n <- if (!is.null(p$aux$trials)) eval(p$aux$trials, df, env)
  if (!is.null(n)) check(paste0("trials() of '", p$lhs, "'"), n, n >= 1 & isInt(n),
                         "must be positive integers")
  if (!is.null(p$aux$exposure)) {
    E <- eval(p$aux$exposure, df, env)
    check(paste0("exposure() of '", p$lhs, "'"), E, E > 0 & is.finite(E), "must be positive")
  }
  if (!is.null(p$aux$weights)) {
    w <- eval(p$aux$weights, df, env)
    check(paste0("weights() of '", p$lhs, "'"), w, w > 0 & is.finite(w), "must be positive and finite")
  }
  ev <- NULL
  if (!is.null(p$aux$event)) {
    ev <- eval(p$aux$event, df, env)
    check(paste0("event() of '", p$lhs, "'"), ev, ev == 0 | ev == 1,
          "must be 0 (censored) or 1 (event)")
  }
  if (fam$support == "any") return(invisible())
  y <- eval(p$responseExpr, df, env)
  what <- paste0("Response '", p$lhs, "' (family '", fam$family, "')")
  switch(fam$support,
    real = check(what, y, rep(TRUE, length(y)), "must be finite"),
    positive = check(what, y, y > 0, "must be positive"),
    survival = check(what, y, y > 0 | (y == 0 & ev == 0),
                     "must be positive (0 only for a censored time)"),
    unit = check(what, y, y > 0 & y < 1, "must be in (0, 1)"),
    count = check(what, y, y >= 0 & isInt(y), "must be non-negative integers"),
    binomial = if (is.null(n)) {
      check(what, y, y == 0 | y == 1, "must be 0 or 1 (use trials() for numbers of successes)")  # NA stays NA (missing responses pass)
    } else {
      check(what, y, y >= 0 & isInt(y) & y <= n, "must be integers between 0 and trials()")
    }
  )
  invisible()
}

# LaTeX of a family: `lhs` of the equation, `residual` (term added to the
# linear predictor, or NULL), `dist` (distribution and prior lines).
# y = LaTeX of the response, tag = "(r)", aux = parsed trials() / exposure() / event() /
# weights() (weights multiply the precision-like parameter: w_i \tau, w_i \phi).
# resp: the response name (keys of the estimate annotations, see annot()).
#
# Three shapes: residual families (gaussian with identity link, lognormal, t) write
# `y = eta + eps` and the distribution of eps; mean families write `link(m) = eta` and
# `y ~ <distribution of the mean m>`, with m and the distribution from familyMeans and
# familyDistributions(); inlaFamily() writes a generic link and the INLA name.
familyLatex <- function(fam, y, tag, rName, aux, obs = "i", resp = NULL) {
  sym <- function(h) annot(hyperSymbol(fam, h, tag), annotKey("family", resp, h))
  y <- annot(y, annotKey("resp", resp))
  auxLatex <- function(a) if (!is.null(aux[[a]])) exprLatex(aux[[a]], obs)
  paren <- function(...) paste0("\\left(", paste0(...), "\\right)")
  # precision-like parameter of row `obs`, times the observation weight of weights()
  w <- auxLatex("weights")
  wsym <- function(h) if (is.null(w)) sym(h) else paste0(w, "\\, ", sym(h))
  normal <- function(m) paste0(annot("\\mathcal{N}", distKey("normal")), paren(m, ", \\left(", wsym("prec"), "\\right)^{-1}"))
  residual <- NULL
  if (fam$kind %in% c("lognormal", "t") || (fam$kind == "gaussian" && identical(fam$link, "identity"))) {
    residual <- annot(paste0("\\varepsilon^{", tag, "}_{", obs, "}"), annotKey("eps", resp))
    lhs <- if (fam$kind == "lognormal") paste0("\\log", paren(y)) else y
    d <- paste0(residual, " &\\sim ", if (fam$kind == "t") {
      paste0(annot("t", distKey("t")), "_{", sym("dof"), "}", paren("0, \\left(", wsym("prec"), "\\right)^{-1}"))
    } else {
      normal("0")
    })
  } else if (fam$kind == "inla") {
    lhs <- paste0("g_{", rName, "}\\left(\\mathrm{E}\\left[", y, "\\right]\\right)")
    d <- paste0(y, " &\\sim ", annot(paste0("\\mathrm{", escapeLatex(fam$family), "}"), distKey("inla")),
                if (!is.null(w)) paste0(", \\quad ", "\\text{scale}", " ", w))
  } else {
    m <- annot(paste0(familyMeans[[fam$kind]], "^{", tag, "}_{", obs, "}"), annotKey("mean", resp))
    lhs <- linkLatex(fam$link, m)
    dists <- familyDistributions(fam, sym, wsym, auxLatex, paren, normal)
    d <- paste0(y, " &\\sim ", dists[[fam$kind]](m))
  }
  priors <- unlist(lapply(names(fam$hyper), function(h) {
    pr <- fam$hyper[[h]]
    if (inherits(pr, "mvremPrior")) priorLatex(pr, 1L, tau = sym(h), Q = NULL, idx = 0L)
  }), use.names = FALSE)
  list(lhs = lhs, residual = residual, dist = c(d, priors))
}

# Symbol of the mean of a mean family (kind -> symbol; superscript tag and index added).
familyMeans <- c(gaussian = "\\mu", gamma = "\\mu", beta = "\\mu", poisson = "\\lambda", nbinomial = "\\mu",
                 binomial = "p", betabinomial = "p", zip = "\\lambda", zinb = "\\mu", weibull = "\\lambda",
                 exponential = "\\lambda")

# Distributions of the mean families (kind -> function of the mean symbol m giving the
# right-hand side of `y ~ ...`); the helpers are familyLatex()'s: sym / wsym = annotated
# hyperparameter symbol (wsym: times the weight), auxLatex = LHS input, paren, normal.
familyDistributions <- function(fam, sym, wsym, auxLatex, paren, normal) {
  distName <- function(name) annot(paste0("\\mathrm{", name, "}"), distKey(tolower(name)))
  # the mean times the exposure of exposure()
  scaledMean <- function(m) {
    E <- auxLatex("exposure")
    if (is.null(E)) m else paste0(E, "\\, ", m)
  }
  # count distribution, zero-truncated ("^{+}") for zero-inflation type 0
  count <- function(name, args, truncated = FALSE)
    paste0(distName(name), if (truncated) annot("^{+}", distKey("truncated")), paren(args))
  zeroInflated <- function(name, args)
    paste0(sym("prob"), "\\,", annot("\\delta_{0}", distKey("zero")), " + \\left(1 - ", sym("prob"), "\\right) ",
           count(name, args, fam$type == 0L))
  nbArgs <- function(m) paste0(scaledMean(m), ", ", sym("size"))
  trials <- auxLatex("trials")
  ev <- auxLatex("event")
  survivalStatus <- if (is.null(ev)) "" else paste0(", \\quad ", annot(ev, distKey("event")), " \\in \\{0, 1\\}")
  list(
    gaussian = function(m) normal(m),
    gamma = function(m) paste0(distName("Gamma"), paren(wsym("prec"), ", ", wsym("prec"), " / ", m)),
    beta = function(m) paste0(distName("Beta"), paren(m, "\\, ", wsym("phi"), ", \\left(1 - ", m, "\\right) ", wsym("phi"))),
    poisson = function(m) paste0(distName("Poisson"), paren(scaledMean(m))),
    nbinomial = function(m) count("NegBin", nbArgs(m)),
    binomial = function(m) if (is.null(trials)) count("Bernoulli", m) else count("Binomial", paste0(trials, ", ", m)),
    betabinomial = function(m) paste0(distName("BetaBinomial"), paren(if (is.null(trials)) "1" else trials, ", ", m, ", ", sym("rho"))),
    zip = function(m) zeroInflated("Poisson", scaledMean(m)),
    zinb = function(m) zeroInflated("NegBin", nbArgs(m)),
    weibull = function(m) paste0(distName("Weibull"), paren(sym("alpha"), ", ", m), ", \\quad S(t) = \\exp",
                                 paren("-", m, "\\, t^{", sym("alpha"), "}"), survivalStatus),
    exponential = function(m) paste0(distName("Exp"), paren(m), ", \\quad S(t) = \\exp", paren("-", m, "\\, t"),
                                     survivalStatus))
}

linkLatex <- function(link, x) {
  if (is.null(link)) link <- "default"
  arg <- paste0("\\left(", x, "\\right)")
  switch(link,
    identity = x,
    log = paste0("\\log", arg),
    probit = paste0("\\Phi^{-1}", arg),
    paste0("\\operatorname{", escapeLatex(link), "}", arg)
  )
}
