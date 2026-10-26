#' Likelihood families
#'
#' Constructors of the likelihood family of a response, given in `familyOut` /
#' `familyLong` of [collectFormulas()], e.g.
#' `familyOut = list(y = betaFamily(phi = gammaPrec(1, 0.1)))`. Each sets the
#' INLA likelihood, its link and the priors of its hyperparameters, i.e. the
#' response's entry of INLA's `control.family`.
#'
#' | Constructor | INLA likelihood | Response | Hyperparameters | LHS input |
#' |---|---|---|---|---|
#' | `gaussianFamily()` | `gaussian` | real | `prec` (residual precision \eqn{\tau_\varepsilon}) | |
#' | `studentTFamily()` | `t` | real | `prec`, `dof` (\eqn{\nu}) | |
#' | `lognormalFamily()` | `lognormal` | positive | `prec` (of \eqn{\log y}) | |
#' | `gammaFamily()` | `gamma` | positive | `prec` (\eqn{\phi}, Var = \eqn{\mu^2/\phi}) | |
#' | `betaFamily()` | `beta` | in (0, 1) | `phi` (\eqn{\phi}, Var = \eqn{\mu(1-\mu)/(1+\phi)}) | |
#' | `poissonFamily()` | `poisson` | count | | `exposure()` |
#' | `negBinomialFamily()` | `nbinomial` | count | `size` (\eqn{\kappa}, Var = \eqn{\mu + \mu^2/\kappa}) | `exposure()` |
#' | `binomialFamily()` | `binomial` | 0, ..., n | | `trials()` (default 1) |
#' | `betaBinomialFamily()` | `betabinomial` | 0, ..., n | `rho` (overdispersion) | `trials()` (default 1) |
#' | `zeroInflatedPoissonFamily()` | `zeroinflatedpoisson<type>` | count | `prob` (\eqn{\pi}) | `exposure()` |
#' | `zeroInflatedNegBinomialFamily()` | `zeroinflatednbinomial<type>` | count | `size`, `prob` | `exposure()` |
#' | `inlaFamily()` | any (non-survival) | | any | `trials()`, `exposure()` |
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
#' Responses are checked against the support of the family by [prepareData()].
#'
#' @param prec,phi,size,rho,prob,dof Prior of the hyperparameter (see Details).
#' @param link Link function (INLA `control.link` model), e.g. `"logit"`,
#'   `"probit"`, `"cloglog"`. Validated against INLA's list when INLA is installed.
#' @param type Zero-inflation type: `1` (\eqn{y = 0} with probability \eqn{\pi},
#'   otherwise from the count distribution) or `0` (\eqn{P(y = 0) = \pi},
#'   positive values from the zero-truncated count distribution).
#' @param family For `inlaFamily()`: name of an INLA likelihood (see
#'   [inlaFamilies()]); survival likelihoods are not supported yet.
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
            list(prec = prec), list(prec = hyperSpec("\\tau^{T}_{\\varepsilon}")),
            control, aux = character(), support = "real")
}

#' @rdname families
#' @export
studentTFamily <- function(prec = NULL, dof = NULL, control = list()) {
  newFamily("t", "t", "identity", "identity",
            list(prec = prec, dof = dof),
            list(prec = hyperSpec("\\tau^{T}_{\\varepsilon}"),
                 dof = hyperSpec("\\nu^{T}", character())),
            control, aux = character(), support = "real")
}

#' @rdname families
#' @export
lognormalFamily <- function(prec = NULL, control = list()) {
  newFamily("lognormal", "lognormal", "identity", "identity",
            list(prec = prec), list(prec = hyperSpec("\\tau^{T}_{\\varepsilon}")),
            control, aux = character(), support = "positive")
}

#' @rdname families
#' @export
gammaFamily <- function(prec = NULL, link = "log", control = list()) {
  newFamily("gamma", "gamma", link, "log",
            list(prec = prec), list(prec = hyperSpec("\\phi^{T}")),
            control, aux = character(), support = "positive")
}

#' @rdname families
#' @export
betaFamily <- function(phi = NULL, link = "logit", control = list()) {
  newFamily("beta", "beta", link, "logit",
            list(phi = phi), list(phi = hyperSpec("\\phi^{T}")),
            control, aux = character(), support = "unit")
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
            list(size = size), list(size = hyperSpec("\\kappa^{T}", "gamma")),
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
            list(rho = rho), list(rho = hyperSpec("\\rho^{T}", character())),
            control, aux = "trials", support = "binomial")
}

#' @rdname families
#' @export
zeroInflatedPoissonFamily <- function(type = 1, prob = NULL, control = list()) {
  type <- zeroInflationType(type)
  fam <- newFamily(paste0("zeroinflatedpoisson", type), "zip", "log", "log",
                   list(prob = prob), list(prob = hyperSpec("\\pi^{T}", character())),
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
                   list(size = hyperSpec("\\kappa^{T}", "gamma"),
                        prob = hyperSpec("\\pi^{T}", character())),
                   control, aux = "exposure", support = "count")
  fam$type <- type
  fam
}

#' @rdname families
#' @export
inlaFamily <- function(family, link = NULL, hyper = list(), control = list()) {
  if (!is.character(family) || length(family) != 1L || is.na(family))
    stop("`family` must be a single INLA likelihood name.", call. = FALSE)
  checkInlaFamily(family, "inlaFamily()")
  if (!is.list(hyper) || (length(hyper) && (is.null(names(hyper)) || any(names(hyper) == ""))))
    stop("`hyper` must be a named list.", call. = FALSE)
  specs <- lapply(names(hyper), function(h)
    hyperSpec(paste0("\\theta^{T}_{\\mathrm{", escapeLatex(h), "}}")))
  newFamily(family, "inla", link, NULL, hyper, stats::setNames(specs, names(hyper)),
            control, aux = c("trials", "exposure"), support = "any")
}

zeroInflationType <- function(type) {
  if (!(identical(type, 0) || identical(type, 1) || identical(type, 0L) || identical(type, 1L)))
    stop("`type` must be 0 or 1.", call. = FALSE)
  as.integer(type)
}

# Scalar prior constructors accepted for precision-like hyperparameters.
precDists <- c("gamma", "pc", "halfNormal", "halfCauchy")

# LaTeX symbol ("T" is replaced by the response tag) and accepted constructors.
hyperSpec <- function(symbol, dists = precDists) list(symbol = symbol, dists = dists)

hyperFields <- c("prior", "param", "initial", "fixed", "to.theta", "from.theta")

newFamily <- function(family, kind, link, linkDefault, hyper, specs, control, aux, support) {
  if (!is.null(link)) {
    if (!is.character(link) || length(link) != 1L || is.na(link))
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
                 support = support),
            class = "mvremFamily")
}

checkHyper <- function(x, name, dists, family) {
  what <- paste0("Hyperparameter `", name, "` of family '", family, "'")
  if (inherits(x, "mvremPrior")) {
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

# Links INLA knows for a likelihood (NULL when INLA is not installed).
familyLinks <- function(family) {
  if (!requireNamespace("INLA", quietly = TRUE)) return(NULL)
  INLA::inla.models()$likelihood[[family]]$link
}

checkInlaFamily <- function(family, where) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    warning("INLA is not installed: likelihood family '", family, "' in ", where,
            " is not validated.", call. = FALSE)
    return(invisible())
  }
  lik <- INLA::inla.models()$likelihood
  if (!family %in% names(lik))
    stop("Unknown INLA likelihood family in ", where, ": ", family, ". See inlaFamilies().",
         call. = FALSE)
  if (isTRUE(lik[[family]]$survival))
    stop("Survival likelihood '", family, "' (", where, ") is not supported yet.",
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

# Check the response (and trials / exposure) values of a formula against the
# support of its family; df = prepared data of the formula.
checkResponse <- function(p, df) {
  fam <- p$family
  if (is.null(fam) || !nrow(df)) return(invisible())
  env <- environment(p$formula)
  isInt <- function(v) abs(v - round(v)) < 1e-8
  check <- function(what, v, ok, msg) {
    if (!is.numeric(v) && !is.logical(v)) stop(what, " must be numeric.", call. = FALSE)
    ok[is.na(ok)] <- TRUE
    if (!all(ok)) stop(what, " ", msg, " (", sum(!ok), " value(s) not).", call. = FALSE)
  }
  n <- if (!is.null(p$aux$trials)) eval(p$aux$trials, df, env)
  if (!is.null(n)) check(paste0("trials() of '", p$lhs, "'"), n, n >= 1 & isInt(n),
                         "must be positive integers")
  if (!is.null(p$aux$exposure)) {
    E <- eval(p$aux$exposure, df, env)
    check(paste0("exposure() of '", p$lhs, "'"), E, E > 0 & is.finite(E), "must be positive")
  }
  if (fam$support == "any") return(invisible())
  y <- eval(p$responseExpr, df, env)
  what <- paste0("Response '", p$lhs, "' (family '", fam$family, "')")
  switch(fam$support,
    real = check(what, y, rep(TRUE, length(y)), ""),
    positive = check(what, y, y > 0, "must be positive"),
    unit = check(what, y, y > 0 & y < 1, "must be in (0, 1)"),
    count = check(what, y, y >= 0 & isInt(y), "must be non-negative integers"),
    binomial = if (is.null(n)) {
      check(what, y, y %in% c(0, 1), "must be 0 or 1 (use trials() for numbers of successes)")
    } else {
      check(what, y, y >= 0 & isInt(y) & y <= n, "must be integers between 0 and trials()")
    }
  )
  invisible()
}

# LaTeX of a family: `lhs` of the equation, `residual` (term added to the
# linear predictor, or NULL), `dist` (distribution and prior lines).
# y = LaTeX of the response, tag = "(r)", aux = parsed trials()/exposure().
familyLatex <- function(fam, y, tag, rName, aux, obs = "i") {
  sym <- function(h) gsub("T", tag, fam$specs[[h]]$symbol, fixed = TRUE)
  mean <- function(s) paste0(s, "^{", tag, "}_{", obs, "}")
  link <- function(x) linkLatex(fam$link, x)
  auxLatex <- function(a) if (!is.null(aux[[a]])) exprLatex(aux[[a]], obs)
  scaledMean <- function(m) {
    E <- auxLatex("exposure")
    if (is.null(E)) m else paste0(E, "\\, ", m)
  }
  paren <- function(...) paste0("\\left(", paste0(...), "\\right)")
  eps <- paste0("\\varepsilon^{", tag, "}_{", obs, "}")
  normal <- function(m) paste0("\\mathcal{N}", paren(m, ", \\left(", sym("prec"), "\\right)^{-1}"))
  # count distribution, zero-truncated ("^{+}") for zero-inflation type 0
  count <- function(name, args, truncated = FALSE)
    paste0("\\mathrm{", name, "}", if (truncated) "^{+}", paren(args))
  zeroInflated <- function(name, args)
    paste0(sym("prob"), "\\,\\delta_{0} + \\left(1 - ", sym("prob"), "\\right) ",
           count(name, args, fam$type == 0L))
  nbArgs <- function(m) paste0(scaledMean(m), ", ", sym("size"))
  residual <- NULL
  lhs <- NULL
  d <- switch(fam$kind,
    gaussian = if (identical(fam$link, "identity")) {
      lhs <- y
      residual <- eps
      paste0(eps, " &\\sim ", normal("0"))
    } else {
      m <- mean("\\mu")
      lhs <- link(m)
      paste0(y, " &\\sim ", normal(m))
    },
    lognormal = {
      lhs <- paste0("\\log", paren(y))
      residual <- eps
      paste0(eps, " &\\sim ", normal("0"))
    },
    t = {
      lhs <- y
      residual <- eps
      paste0(eps, " &\\sim t_{", sym("dof"), "}", paren("0, \\left(", sym("prec"), "\\right)^{-1}"))
    },
    gamma = {
      m <- mean("\\mu")
      lhs <- link(m)
      paste0(y, " &\\sim \\mathrm{Gamma}", paren(sym("prec"), ", ", sym("prec"), " / ", m))
    },
    beta = {
      m <- mean("\\mu")
      lhs <- link(m)
      paste0(y, " &\\sim \\mathrm{Beta}", paren(m, "\\, ", sym("phi"), ", \\left(1 - ", m,
                                                  "\\right) ", sym("phi")))
    },
    poisson = {
      m <- mean("\\lambda")
      lhs <- link(m)
      paste0(y, " &\\sim \\mathrm{Poisson}", paren(scaledMean(m)))
    },
    nbinomial = {
      m <- mean("\\mu")
      lhs <- link(m)
      paste0(y, " &\\sim ", count("NegBin", nbArgs(m)))
    },
    binomial = {
      m <- mean("p")
      lhs <- link(m)
      n <- auxLatex("trials")
      paste0(y, " &\\sim ", if (is.null(n)) count("Bernoulli", m) else
        count("Binomial", paste0(n, ", ", m)))
    },
    betabinomial = {
      m <- mean("p")
      lhs <- link(m)
      n <- auxLatex("trials")
      paste0(y, " &\\sim \\mathrm{BetaBinomial}",
             paren(if (is.null(n)) "1" else n, ", ", m, ", ", sym("rho")))
    },
    zip = {
      m <- mean("\\lambda")
      lhs <- link(m)
      paste0(y, " &\\sim ", zeroInflated("Poisson", scaledMean(m)))
    },
    zinb = {
      m <- mean("\\mu")
      lhs <- link(m)
      paste0(y, " &\\sim ", zeroInflated("NegBin", nbArgs(m)))
    },
    inla = {
      lhs <- paste0("g_{", rName, "}\\left(\\mathrm{E}\\left[", y, "\\right]\\right)")
      paste0(y, " &\\sim \\mathrm{", escapeLatex(fam$family), "}")
    }
  )
  priors <- unlist(lapply(names(fam$hyper), function(h) {
    pr <- fam$hyper[[h]]
    if (inherits(pr, "mvremPrior")) priorLatex(pr, 1L, tau = sym(h), Q = NULL, idx = 0L)
  }), use.names = FALSE)
  list(lhs = lhs, residual = residual, dist = c(d, priors))
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
