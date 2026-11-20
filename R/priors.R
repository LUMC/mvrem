#' Prior distributions for random-effect blocks and associations
#'
#' Constructors of prior distributions for the hyperparameters of a
#' random-effect block, given as `re(..., prior = <constructor>)`, and for the
#' association coefficient \eqn{\gamma} of a `scaled()` term, given as
#' `scaled(x, prior = normalCoef(...))`. Each name
#' tells the distribution and the quantity it is placed on: the precision
#' \eqn{\tau}, the standard deviation \eqn{\sigma = \tau^{-1/2}}, the
#' precision matrix \eqn{\mathbf{Q}} of a correlated block, or a coefficient.
#'
#' * `gammaPrec()`: \eqn{\tau \sim \mathrm{Gamma}(shape, rate)} (INLA `loggamma`;
#'   the defaults are INLA's default for `iid` effects).
#' * `pcSd()`: penalised-complexity prior with \eqn{P(\sigma > u) = \alpha},
#'   i.e. \eqn{\sigma \sim \mathrm{Exp}(-\log(\alpha) / u)} (INLA `pc.prec`).
#' * `halfNormalSd()`: \eqn{\sigma \sim \mathcal{N}^{+}(0, scale^2)} (INLA
#'   `logtnormal` with mean 0 and precision \eqn{1/scale^2}).
#' * `halfCauchySd()`: \eqn{\sigma \sim \mathrm{C}^{+}(0, scale)}, the
#'   half-Cauchy prior of Gelman (2006) (INLA `expression:` prior on
#'   \eqn{\log \tau}).
#' * `wishartPrec()`: \eqn{\mathbf{Q} \sim \mathcal{W}_k(df, \mathbf{V})} with
#'   \eqn{\mathrm{E}[\mathbf{Q}] = df\,\mathbf{V}} (INLA `wishart<k>d` /
#'   `wishartkd`, whose parameter \eqn{\mathbf{R} = \mathbf{V}^{-1}}). With
#'   `scale = "data"` the scale matrix is taken from the data: the formula is
#'   fitted with lme4 (`lmer()` for Gaussian and lognormal responses, `glmer()`
#'   for binomial and Poisson ones, with `trials()` / `exposure()`; other
#'   families are an error) and \eqn{\mathbf{V} = \hat{\mathbf{Q}} / df} with
#'   \eqn{\hat{\mathbf{Q}}} the precision matrix of the block in that fit, so
#'   that \eqn{\mathrm{E}[\mathbf{Q}] = \hat{\mathbf{Q}}}. INLA's default
#'   (identity scale) is informative in the data's units and shrinks small
#'   random-effect SDs towards \eqn{1/\sqrt{df}}, which also shrinks the
#'   correlations; a data-scaled prior avoids this. A singular lme4 fit is
#'   regularised with a warning (SDs floored at 1% of the largest, correlations
#'   clamped to \eqn{\pm 0.95}).
#' * `normalCoef()`: \eqn{\gamma \sim \mathcal{N}(mean, sd^2)} for the
#'   association coefficient of a `scaled()` term (INLA `normal` prior on the
#'   `beta` of the copied field, with precision \eqn{1/sd^2}); the defaults are
#'   mvrem's default for `scaled()`. Only `scaled()` accepts it.
#'
#' Scalar priors are for single effects and uncorrelated (`||`) blocks; the same
#' prior (for a data-scaled one: the same rule) is used for every effect of the block (split the block into several
#' `re()` terms for different priors).
#'
#' **Data-scaled scalar priors.** `pcSd(u = "data")`, `halfNormalSd(scale =
#' "data")`, `halfCauchySd(scale = "data")` and `gammaPrec(rate = "data")` take
#' their scale from the data, separately for every effect of the block: the
#' formula is fitted with lme4 (as for `wishartPrec(scale = "data")`; Gaussian,
#' lognormal, binomial and Poisson responses, other families are an error) and
#' \eqn{s} is the SD of the effect in that fit. Rules:
#' * `pcSd(u = "data", alpha = 0.05)`: \eqn{u = 6 s}, i.e. \eqn{P(\sigma > 6 s) = 0.05}
#'   (\eqn{\sigma \sim \mathrm{Exp}(0.5 / s)}, prior median \eqn{1.4 s});
#'   `alpha` can be changed.
#' * `halfNormalSd(scale = "data")`: \eqn{scale = 2 s} (prior median \eqn{1.35 s}).
#' * `halfCauchySd(scale = "data")`: \eqn{scale = 2 s} (prior median \eqn{2 s}).
#' * `gammaPrec(shape, rate = "data")`: \eqn{rate = shape \cdot s^2}, so that
#'   \eqn{\mathrm{E}[\tau] = 1 / s^2}; `shape` (default 1, an exponential
#'   precision) sets the dispersion (coefficient of variation
#'   \eqn{1/\sqrt{shape}}).
#'
#' An SD below the floor \eqn{0.5\,\hat\sigma / \mathrm{rms}(w)} (\eqn{\hat\sigma}
#' the residual SD of the lme4 fit, 1 on the link scale of a binomial or Poisson
#' response; \eqn{\mathrm{rms}(w)} the root mean square of the effect's covariate,
#' 1 for an intercept) is raised to it with a warning: a singular lme4 fit
#' (SD 0) says only that the effect is not distinguishable from the residual
#' noise, and the prior must still allow SDs of that order. INLA's default
#' (`loggamma(1, 5e-5)` on the precision) puts most of its mass on SDs near 0
#' and collapses small SDs; fixed priors such as `pcSd(1)` work only when their
#' scale suits the units of the data. The rules were chosen by simulation
#' (`tools/simDataPriors.R`: 40 subjects x 4 visits, SDs 0.2 / 0.1; slope SD
#' bias -72% with INLA's default, between -10% and +3% with the data-scaled priors).
#'
#' A correlated block with \eqn{k > 1}
#' effects accepts only `wishartPrec()`. The dimension \eqn{k} is known once the
#' data is (factor levels expand), so `wishartPrec()` is checked against the
#' block by [buildInlaFormula()].
#'
#' @param shape,rate Shape and rate of the gamma distribution; `rate = "data"`:
#'   from the data (see above).
#' @param mean,sd Mean and standard deviation of the normal prior of a
#'   coefficient.
#' @param u,alpha Upper bound and tail probability: \eqn{P(\sigma > u) = \alpha};
#'   `u = "data"`: from the data (see above; `alpha` then defaults to 0.05).
#' @param scale For `halfNormalSd()` / `halfCauchySd()`: scale of the
#'   distribution of \eqn{\sigma}, or `"data"` (see above). For `wishartPrec()`: the scale matrix
#'   \eqn{\mathbf{V}}, a symmetric positive-definite \eqn{k \times k} matrix,
#'   a positive number \eqn{s} meaning \eqn{s\,\mathbf{I}_k}, or `"data"` (see
#'   above).
#' @param df Degrees of freedom of the Wishart distribution, \eqn{df > k + 1};
#'   `NULL` means \eqn{k + 2}.
#' @return An object of class `mvremPrior`.
#' @examples
#' pcSd(u = 1)
#' wishartPrec(df = 4, scale = diag(c(1, 0.1)))
#' wishartPrec(scale = "data")
#' pcSd(u = "data")
#' gammaPrec(rate = "data")
#' normalCoef(mean = 0, sd = 1)
#' @name priors
NULL

#' @rdname priors
#' @export
gammaPrec <- function(shape = 1, rate = 5e-5) {
  checkPositive(shape, "shape")
  checkPositiveOrData(rate, "rate")
  newPrior("gamma", "prec", list(shape = shape, rate = rate))
}

#' @rdname priors
#' @export
pcSd <- function(u, alpha = if (isData(u)) dataPriorRules$pcAlpha else 0.01) {
  checkPositiveOrData(u, "u")
  if (!isNumber(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  newPrior("pc", "sd", list(u = u, alpha = alpha))
}

#' @rdname priors
#' @export
halfNormalSd <- function(scale) {
  checkPositiveOrData(scale, "scale")
  newPrior("halfNormal", "sd", list(scale = scale))
}

#' @rdname priors
#' @export
halfCauchySd <- function(scale) {
  checkPositiveOrData(scale, "scale")
  newPrior("halfCauchy", "sd", list(scale = scale))
}

#' @rdname priors
#' @export
wishartPrec <- function(df = NULL, scale = 1) {
  if (!is.null(df)) checkPositive(df, "df")
  if (is.character(scale)) {
    if (length(scale) != 1L || is.na(scale) || scale != "data")
      stop("`scale` must be a matrix, a positive number or \"data\".", call. = FALSE)
  } else if (is.matrix(scale)) {
    # is.finite() also catches NA / NaN; Inf entries would make eigen() fail (B9)
    if (!is.numeric(scale) || !all(is.finite(scale)) || nrow(scale) != ncol(scale) ||
        !isSymmetric(unname(scale)) ||
        any(eigen(scale, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("`scale` must be a symmetric positive-definite matrix.", call. = FALSE)
    scale <- unname(scale)
  } else {
    checkPositive(scale, "scale")
  }
  newPrior("wishart", "precMatrix", list(df = df, scale = scale))
}

#' @rdname priors
#' @export
normalCoef <- function(mean = 0, sd = 10) {
  if (!isNumber(mean) || !is.finite(mean))
    stop("`mean` must be a single finite number.", call. = FALSE)
  checkPositive(sd, "sd")
  newPrior("normal", "coef", list(mean = mean, sd = sd))
}

newPrior <- function(dist, target, args) {
  structure(list(dist = dist, target = target, args = args), class = "mvremPrior")
}

# Is x the string "data" (a data-scaled prior parameter)?
isData <- function(x) is.character(x) && length(x) == 1L && identical(x, "data")

checkPositiveOrData <- function(x, name) {
  if (is.character(x) && !isData(x))
    stop("`", name, "` must be a single positive number or \"data\".", call. = FALSE)
  if (!isData(x)) checkPositive(x, name)
}

checkPositive <- function(x, name) {
  if (!isNumber(x) || !is.finite(x) || x <= 0)
    stop("`", name, "` must be a single positive number.", call. = FALSE)
}

#' @export
format.mvremPrior <- function(x, ...) {
  a <- x$args
  num <- function(v) vapply(v, function(z) format(signif(z, 4)), "")
  if (isScalarDataPrior(x)) return(paste("<mvrem prior>", dataPriorText(x, num)))
  v <- priorParams(x)          # block size unknown: a Wishart's default df is written "k + 2"
  txt <- switch(x$dist,
    gamma = paste0("Gamma prior on precision: tau ~ Gamma(shape = ", num(v$shape),
                   ", rate = ", num(v$rate), ")   [INLA: loggamma, param = c(",
                   num(v$shape), ", ", num(v$rate), ")]"),
    pc = paste0("PC prior on SD: P(sigma > ", num(v$u), ") = ", num(v$alpha),
                ", i.e. sigma ~ Exp(", num(v$rate), ")   [INLA: pc.prec, param = c(",
                num(v$u), ", ", num(v$alpha), ")]"),
    halfNormal = paste0("Half-normal prior on SD: sigma ~ N+(0, sd = ", num(v$scale),
                        ")   [INLA: logtnormal, param = c(0, ", num(v$prec), ")]"),
    halfCauchy = paste0("Half-Cauchy prior on SD: sigma ~ C+(0, ", num(v$scale),
                        ")   [INLA: expression: prior on log precision]"),
    wishart = paste0("Wishart prior on precision matrix: Q ~ W(df = ",
                     if (is.null(v$df)) "k + 2" else num(v$df), ", scale = ",
                     scaleText(a$scale, num), if (identical(a$source, "data")) " (from the data)",
                     ")   [INLA: wishart<k>d / wishartkd, R = scale^-1]"),
    normal = paste0("Normal prior on coefficient: gamma ~ N(mean = ", num(v$mean), ", sd = ",
                    num(v$sd), ")   [INLA: beta, normal, param = c(", num(v$mean), ", ",
                    num(v$prec), ")]")
  )
  if (identical(a$source, "data") && !is.null(a$sdData))
    txt <- paste0(txt, "   (from the data: lme4 SD ", num(a$sdData), ")")
  paste("<mvrem prior>", txt)
}

# format() of a scalar prior with a "data" parameter (not yet resolved).
dataPriorText <- function(x, num) {
  a <- x$args
  r <- dataPriorRules
  each <- " x the lme4 SD of each effect)"
  switch(x$dist,
    pc = paste0("PC prior on SD: P(sigma > u) = ", num(a$alpha), " with u = data (", num(r$pcU), each,
                "   [INLA: pc.prec]"),
    halfNormal = paste0("Half-normal prior on SD: sigma ~ N+(0, sd = data (", num(r$halfNormalScale), each,
                        ")   [INLA: logtnormal]"),
    halfCauchy = paste0("Half-Cauchy prior on SD: sigma ~ C+(0, scale = data (", num(r$halfCauchyScale), each,
                        ")   [INLA: expression: prior on log precision]"),
    gamma = paste0("Gamma prior on precision: tau ~ Gamma(shape = ", num(a$shape),
                   ", rate = data (shape x the lme4 variance of each effect: E[tau] = 1 / SD^2))",
                   "   [INLA: loggamma]"))
}

#' @export
print.mvremPrior <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

scaleText <- function(s, num) {
  if (is.character(s)) return("\"data\" (lme4 fit of the formula: E[Q] = its precision matrix)")
  if (!is.matrix(s)) return(paste(num(s), "* I"))
  if (all(s[row(s) != col(s)] == 0)) return(paste0("diag(", paste(num(diag(s)), collapse = ", "), ")"))
  paste0("rbind(", paste(apply(s, 1L, function(r) paste0("c(", paste(num(r), collapse = ", "), ")")),
                         collapse = ", "), ")")
}

# The numbers a numeric (resolved) prior implies, derived here once for every renderer:
# format(), priorHyper(), coefHyper(), priorLatex(), coefPriorLatex(), priorNote().
#   gamma:      shape, rate          tau ~ Gamma(shape, rate)
#   pc:         u, alpha, rate       sigma ~ Exp(rate), rate = -log(alpha) / u
#   halfNormal: scale, prec          sigma ~ N+(0, 1 / prec), prec = 1 / scale^2
#   halfCauchy: scale                sigma ~ C+(0, scale)
#   normal:     mean, sd, prec       gamma ~ N(mean, 1 / prec), prec = 1 / sd^2
#   wishart:    df, scale            see wishartArgs() (k = block size, NULL if unknown)
priorParams <- function(prior, k = NULL) {
  a <- prior$args
  switch(prior$dist,
    gamma = list(shape = a$shape, rate = a$rate),
    pc = list(u = a$u, alpha = a$alpha, rate = -log(a$alpha) / a$u),
    halfNormal = list(scale = a$scale, prec = 1 / a$scale^2),
    halfCauchy = list(scale = a$scale),
    normal = list(mean = a$mean, sd = a$sd, prec = 1 / a$sd^2),
    wishart = wishartArgs(prior, k))
}

# Wishart settings resolved for a block of k effects: df (default k + 2) and the k x k
# scale matrix V (NULL for scale = "data"). With k = NULL (block unknown, as in format())
# df is NULL unless given and V NULL unless given as a matrix.
wishartArgs <- function(prior, k = NULL) {
  a <- prior$args
  df <- if (!is.null(a$df)) a$df else if (!is.null(k)) k + 2
  V <- if (is.matrix(a$scale)) a$scale else if (is.character(a$scale) || is.null(k)) NULL else diag(a$scale, k)
  list(df = df, scale = V)
}

# Check that `prior` suits a block of k effects (correlated: `|` block).
checkBlockPrior <- function(prior, k, correlated, term) {
  if (is.null(prior)) return(invisible())
  blk <- paste0("`", deparse1(term), "`")
  if (prior$target == "coef")
    stop("normalCoef() is a prior for the coefficient of scaled(); block ", blk,
         " needs gammaPrec(), pcSd(), halfNormalSd(), halfCauchySd() or wishartPrec().",
         call. = FALSE)
  if (correlated && k > 1L) {
    if (prior$dist != "wishart")
      stop("Correlated block ", blk, " (k = ", k, ") needs wishartPrec(); use `||` or ",
           "separate re() blocks for scalar priors.", call. = FALSE)
    if (is.matrix(prior$args$scale) && nrow(prior$args$scale) != k)
      stop("wishartPrec() scale of block ", blk, " must be a ", k, " x ", k, " matrix.",
           call. = FALSE)
    df <- wishartArgs(prior, k)$df
    if (df <= k + 1)
      stop("wishartPrec() of block ", blk, ": df must be > k + 1 = ", k + 1L, ".",
           call. = FALSE)
  } else if (prior$dist == "wishart") {
    stop("wishartPrec() needs a correlated block with >= 2 effects (block ", blk,
         "); use gammaPrec(), pcSd(), halfNormalSd() or halfCauchySd().", call. = FALSE)
  }
  invisible()
}

# INLA `hyper` list of a prior for a block of k effects.
priorHyper <- function(prior, k) {
  if (isScalarDataPrior(prior))
    stop(priorCall(prior), " must be resolved with resolveEffectPriors() first.")
  if (prior$dist == "wishart" && is.character(prior$args$scale))
    stop("wishartPrec(scale = \"data\") must be resolved with resolveBlockPrior() first.")
  v <- priorParams(prior, k)
  prec <- function(p, param) list(prec = c(list(prior = p), if (!is.null(param)) list(param = param)))
  switch(prior$dist,
    gamma = prec("loggamma", c(v$shape, v$rate)),
    pc = prec("pc.prec", c(v$u, v$alpha)),
    halfNormal = prec("logtnormal", c(0, v$prec)),
    halfCauchy = prec(halfCauchyExpression(v$scale), NULL),
    wishart = {
      R <- solve(v$scale)          # INLA's parameter is R = V^-1
      list(theta1 = list(prior = if (k <= 5L) paste0("wishart", k, "d") else "wishartkd",
                         param = c(v$df, diag(R), R[lower.tri(R)])))
    }
  )
}

# Default prior of the coefficient of a scaled() term (the INLA `beta` of the copy).
defaultCoefPrior <- function() normalCoef(mean = 0, sd = 10)

# The prior of a random-effect block with its data-dependent parts filled in:
# wishartPrec(scale = "data") becomes a wishartPrec() whose scale matrix is
# V = Q / df with Q the precision matrix of the block in an lme4 fit of the
# formula (so E[Q_prior] = Q). p = the parsed formula, data = its prepared data,
# effects = the block's design columns, blk = the block entry of p$re. Other
# priors are returned unchanged.
resolveBlockPrior <- function(prior, p, data, effects, blk) {
  if (is.null(prior) || prior$dist != "wishart" || !is.character(prior$args$scale)) return(prior)
  S <- blockCovariance(p, data, effects, blk)
  df <- wishartArgs(prior, length(effects))$df
  V <- solve(S) / df
  V <- (V + t(V)) / 2
  dimnames(V) <- NULL
  newPrior("wishart", "precMatrix", list(df = prior$args$df, scale = V, source = "data"))
}

# Covariance matrix of a random-effect block estimated with lme4 (lmer /
# glmer) from the formula p on the data (see lme4BlockFit()).
blockCovariance <- function(p, data, effects, blk) {
  what <- dataPriorWhat(p, blk)
  fit <- lme4BlockFit(p, data, what, "wishartPrec(scale = \"data\")", "a numeric scale")
  vc <- lme4::VarCorr(fit)
  S <- NULL
  for (nm in names(vc)) {
    if (varCorrGroup(nm) != blk$group) next
    m <- as.matrix(vc[[nm]])
    if (nrow(m) == length(effects) && setequal(rownames(m), effects)) { S <- m[effects, effects, drop = FALSE]; break }
  }
  if (is.null(S)) stop("No random-effect block `", deparse1(blk$term), "` in the lme4 fit for ", what, ".",
                       call. = FALSE)
  regularizeCovariance(unname(S), what)
}

# Grouping factor of a VarCorr() element. lme4 names the elements by grouping factor and,
# when several terms share one (the effects of a `||` block are split into one term each;
# two re() blocks may have the same group), makes the names unique: "id", "id.1", "id.2".
# The suffix is stripped to match the block's group.
varCorrGroup <- function(nm) sub("\\.[0-9]+$", "", nm)

dataPriorWhat <- function(p, blk)
  paste0("the data-scaled prior of block `", deparse1(blk$term), "` (formula of '", p$response, "')")

# lme4 fit (lmer / glmer) of the formula p (its lme4Formula) on the data, for
# a data-scaled prior. Supported families: gaussian (also with a non-identity
# link, fitted with the identity), lognormal (log response), binomial (trials()
# -> cbind), poisson (exposure() -> offset). lme4's warnings are collected into
# one warning; what = description of the block, prior / alt = the prior
# constructor and the alternative named in the error for other families.
lme4BlockFit <- function(p, data, what, prior, alt) {
  fam <- p$family
  kind <- if (is.null(fam)) "gaussian" else fam$kind
  f <- p$lme4Formula
  lhs <- p$responseExpr
  gl <- NULL                       # NULL: lmer(), else the glmer() family
  switch(kind,
    gaussian = NULL,
    lognormal = lhs <- call("log", lhs),
    binomial = {
      if (!is.null(p$aux$trials)) lhs <- call("cbind", lhs, call("-", p$aux$trials, lhs))
      gl <- stats::binomial(link = if (is.null(fam$link)) "logit" else fam$link)
    },
    poisson = {
      if (!is.null(p$aux$exposure)) f[[3L]] <- call("+", f[[3L]], call("offset", call("log", p$aux$exposure)))
      gl <- stats::poisson()
    },
    stop(prior, " needs a Gaussian, lognormal, binomial or Poisson response ",
         "(estimated with lme4); give ", alt, " for ", what, ".", call. = FALSE))
  f[[2L]] <- lhs
  msgs <- character()
  fit <- withCallingHandlers(
    tryCatch({
      if (is.null(gl)) {
        lme4::lmer(f, data = data, na.action = stats::na.omit, control = lme4Control("lmer"))
      } else {
        lme4::glmer(f, data = data, family = gl, na.action = stats::na.omit, control = lme4Control("glmer"))
      }
    }, error = function(e) stop("lme4 fit for ", what, " failed: ", conditionMessage(e), call. = FALSE)),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  if (length(msgs))
    warning("lme4 fit for ", what, ": ", paste(unique(trimws(msgs)), collapse = "; "), call. = FALSE)
  fit
}

# A positive-definite version of an lme4 block covariance: SDs floored at 1% of
# the largest SD, correlations clamped to +-0.95 (a singular lme4 fit has zero
# variances or correlations of +-1, which no Wishart can be centred on).
regularizeCovariance <- function(S, what) {
  sds <- sqrt(pmax(diag(S), 0))
  if (max(sds) <= 0) stop("All random-effect variances of the lme4 fit for ", what, " are zero.", call. = FALSE)
  k <- nrow(S)
  R <- diag(k)
  for (i in seq_len(k)) for (j in seq_len(k)) if (i != j && sds[i] > 0 && sds[j] > 0)
    R[i, j] <- S[i, j] / (sds[i] * sds[j])
  R2 <- pmin(pmax(R, -0.95), 0.95); diag(R2) <- 1
  sds2 <- pmax(sds, 0.01 * max(sds))
  S2 <- R2 * tcrossprod(sds2)
  ev <- eigen(S2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8 * max(ev)) {   # k > 2 with clamped correlations may still be singular
    S2 <- S2 + diag(1e-4 * max(ev), k)
  }
  if (any(sds2 != sds) || any(R2 != R) || !isTRUE(all.equal(S2, R2 * tcrossprod(sds2))))
    warning("The lme4 fit for ", what, " is singular; its covariance matrix was regularised ",
            "(SDs >= 1% of the largest, |correlations| <= 0.95) before centring the prior.", call. = FALSE)
  S2
}

# Rules of the data-scaled scalar priors (E15, chosen by simulation: see
# tools/simDataPriors.R and CLAUDE.md). s = the lme4 SD of the effect, floored
# at `floor` times the reference scale (residual SD of lmer, 1 on the link scale
# of glmer).
dataPriorRules <- list(pcU = 6, pcAlpha = 0.05, halfNormalScale = 2, halfCauchyScale = 2,
                       floor = 0.5)

# Is the prior a scalar prior with a "data" parameter (not yet resolved)?
isScalarDataPrior <- function(prior)
  !is.null(prior) && prior$dist %in% c("gamma", "pc", "halfNormal", "halfCauchy") &&
    any(vapply(prior$args[setdiff(names(prior$args), "source")], isData, NA))

# The priors of the effects of an uncorrelated block (a single effect or a `||`
# block), one per effect: a scalar prior with a "data" parameter is resolved per
# effect from the SD of that effect in an lme4 fit of the formula (see
# dataScalarPrior()); other priors (also NULL) are repeated unchanged.
resolveEffectPriors <- function(prior, p, data, effects, blk) {
  if (!isScalarDataPrior(prior)) return(rep(list(prior), length(effects)))
  s <- effectSds(p, data, effects, blk, priorCall(prior))
  lapply(seq_along(effects), function(e) dataScalarPrior(prior, s[[e]]))
}

# The prior with the lme4 SD s filled in; args$source = "data", args$sdData = s.
dataScalarPrior <- function(prior, s) {
  a <- prior$args
  r <- dataPriorRules
  out <- switch(prior$dist,
    pc = pcSd(u = r$pcU * s, alpha = a$alpha),
    halfNormal = halfNormalSd(r$halfNormalScale * s),
    halfCauchy = halfCauchySd(r$halfCauchyScale * s),
    gamma = gammaPrec(shape = a$shape, rate = a$shape * s^2))
  out$args$source <- "data"
  out$args$sdData <- s
  out
}

# How the prior was written (for messages), e.g. pcSd(u = "data").
priorCall <- function(prior) {
  a <- prior$args
  switch(prior$dist,
    pc = "pcSd(u = \"data\")", halfNormal = "halfNormalSd(scale = \"data\")",
    halfCauchy = "halfCauchySd(scale = \"data\")", gamma = "gammaPrec(rate = \"data\")",
    paste0(priorConstructor(prior$dist), "()"))
}

# SDs of the effects of an uncorrelated block in an lme4 fit of the formula (one
# fit for all effects); SDs below `floor` times the reference scale (lmer:
# residual SD; glmer: 1 on the link scale) are raised to it with a warning.
effectSds <- function(p, data, effects, blk, prior) {
  what <- dataPriorWhat(p, blk)
  fit <- lme4BlockFit(p, data, what, prior, "a numeric parameter")
  vc <- lme4::VarCorr(fit)
  grp <- vapply(names(vc), varCorrGroup, "")
  s <- vapply(effects, function(e) {
    hit <- NA_real_
    for (exact in c(TRUE, FALSE)) for (j in which(grp == blk$group)) {
      m <- as.matrix(vc[[j]])
      if (!e %in% rownames(m) || (exact && nrow(m) != 1L)) next
      hit <- sqrt(max(m[e, e], 0))
      break
    }
    if (is.na(hit)) stop("No random effect `", e, "` | ", blk$group, " in the lme4 fit for ", what, ".",
                         call. = FALSE)
    hit
  }, 0)
  glmm <- inherits(fit, "glmerMod")
  ref <- if (glmm) 1 else stats::sigma(fit)
  # the floor of an effect with covariate w: f x ref / rms(w), so that its contribution
  # b w has an SD of f x ref (intercept: rms 1)
  X <- stats::model.matrix(blockLhs(blk, p), stats::model.frame(fit))
  rms <- vapply(effects, function(e) if (e %in% colnames(X)) sqrt(mean(X[, e]^2)) else 1, 0)
  rms[!is.finite(rms) | rms <= 0] <- 1
  floor <- unname(dataPriorRules$floor * ref / rms)
  low <- s < floor
  if (any(low))
    warning("The lme4 fit for ", what, " gives a small (or zero: singular fit) SD for ",
            paste0("`", effects[low], "`", collapse = ", "), " (", paste(format(signif(s[low], 3)), collapse = ", "),
            "); the prior uses the floor ", paste(format(signif(floor[low], 3)), collapse = ", "), " (",
            dataPriorRules$floor * 100, "% of ", if (glmm) "1 on the link scale" else "the residual SD",
            ", divided by the root mean square of the effect's covariate) instead.", call. = FALSE)
  unname(pmax(s, floor))
}

# INLA `hyper` list of the coefficient prior of a scaled() term.
coefHyper <- function(prior) {
  v <- priorParams(prior)
  list(beta = list(fixed = FALSE, prior = "normal", param = c(v$mean, v$prec)))
}

# LaTeX line of the prior of an association coefficient (gamma = its symbol).
coefPriorLatex <- function(prior, gamma) {
  if (is.null(prior)) return(character())
  v <- priorParams(prior)
  paste0(gamma, " &\\sim ", annot("\\mathcal{N}", distKey("priorNormal")), "\\left(", numLatex(v$mean),
         ", \\left(", numLatex(v$prec), "\\right)^{-1}\\right)")
}

# Density of sigma = exp(-theta / 2) ~ C+(0, s), on theta = log precision
# (Jacobian |d sigma / d theta| = sigma / 2).
halfCauchyExpression <- function(s) {
  paste0("expression: sigma = exp(-theta/2); s = ", format(s, digits = 15),
         "; logdens = log(2) - log(pi) - log(s) - log(1 + (sigma/s)^2);",
         " return(logdens + log(sigma) - log(2));")
}

# LaTeX lines of the prior of an RE block: tau = precision symbol(s) (scalar
# priors) or Q (wishart); idx = effect indices for the "k = ..." suffix.
priorLatex <- function(prior, k, tau, Q, idx) {
  if (is.null(prior)) return(character())
  suffix <- if (length(idx) > 1L) paste0(", \\quad k = ", paste(idx, collapse = ", ")) else ""
  sd <- paste0("\\left(", tau, "\\right)^{-1/2}")
  if (isScalarDataPrior(prior)) return(dataPriorLatex(prior, tau, sd, suffix, length(idx) > 1L))
  a <- prior$args
  if (prior$dist == "wishart" && is.matrix(a$scale)) k <- nrow(a$scale)
  v <- priorParams(prior, k)
  switch(prior$dist,
    gamma = paste0(tau, " &\\sim ", annot("\\mathrm{Gamma}", distKey("priorGamma")), "\\left(", numLatex(v$shape), ", ",
                   numLatex(v$rate), "\\right)", suffix),
    pc = paste0(sd, " &\\sim ", annot("\\mathrm{Exp}", distKey("priorPc")), "\\left(", numLatex(v$rate),
                "\\right), \\quad P\\left(", sd, " > ", numLatex(v$u), "\\right) = ",
                numLatex(v$alpha), suffix),
    halfNormal = paste0(sd, " &\\sim ", annot("\\mathcal{N}^{+}", distKey("priorHalfNormal")), "\\left(0, \\left(",
                        numLatex(v$prec), "\\right)^{-1}\\right)", suffix),
    halfCauchy = paste0(sd, " &\\sim ", annot("\\mathrm{C}^{+}", distKey("priorHalfCauchy")), "\\left(0, ",
                        numLatex(v$scale), "\\right)", suffix),
    wishart = {
      # the scale as written (s I_k for a number, \hat{V} for "data"), not wishartArgs()'s matrix
      V <- if (is.character(a$scale)) "\\hat{\\mathbf{V}}_{\\mathrm{data}}" else matrixLatex(a$scale, k)
      paste0(Q, " &\\sim ", annot(paste0("\\mathcal{W}_{", k, "}"), distKey("priorWishart")), "\\left(",
             numLatex(v$df), ", ", V, "\\right)")
    }
  )
}

# LaTeX line of a scalar prior with a "data" parameter (no data to resolve it):
# the lme4 SD is the symbol \hat{s}_{data} (\hat{s}_{k,data} for several effects).
dataPriorLatex <- function(prior, tau, sd, suffix, several) {
  a <- prior$args
  r <- dataPriorRules
  sHat <- paste0("\\hat{s}_{", if (several) "k,", "\\mathrm{data}}")
  times <- function(cc) if (cc == 1) sHat else paste0(numLatex(cc), "\\,", sHat)
  switch(prior$dist,
    gamma = paste0(tau, " &\\sim ", annot("\\mathrm{Gamma}", distKey("priorGamma")), "\\left(", numLatex(a$shape), ", ",
                   times(a$shape), "^{2}\\right)", suffix),
    # rate -log(alpha) / u with u = pcU s written as a fraction (cf. priorParams())
    pc = paste0(sd, " &\\sim ", annot("\\mathrm{Exp}", distKey("priorPc")), "\\left(\\frac{", numLatex(-log(a$alpha)), "}{",
                times(r$pcU), "}\\right), \\quad P\\left(", sd, " > ", times(r$pcU), "\\right) = ",
                numLatex(a$alpha), suffix),
    halfNormal = paste0(sd, " &\\sim ", annot("\\mathcal{N}^{+}", distKey("priorHalfNormal")), "\\left(0, \\left(\\left(",
                        times(r$halfNormalScale), "\\right)^{-2}\\right)^{-1}\\right)", suffix),
    halfCauchy = paste0(sd, " &\\sim ", annot("\\mathrm{C}^{+}", distKey("priorHalfCauchy")), "\\left(0, ",
                        times(r$halfCauchyScale), "\\right)", suffix))
}

# Number with `digits` significant digits; 5e-05 -> 5 \cdot 10^{-5}.
numLatex <- function(x, digits = 3) {
  s <- formatC(signif(x, digits), digits = digits, format = "g")
  s <- trimws(s)
  if (grepl("e", s)) {
    parts <- strsplit(s, "e", fixed = TRUE)[[1L]]
    s <- paste0(if (parts[1L] != "1") paste0(parts[1L], " \\cdot "), "10^{",
                as.integer(parts[2L]), "}")
  }
  s
}

matrixLatex <- function(s, k) {
  if (!is.matrix(s)) {
    I <- paste0("\\mathbf{I}_{", k, "}")
    return(if (s == 1) I else paste0(numLatex(s), "\\,", I))
  }
  rows <- apply(s, 1L, function(r) paste(vapply(r, numLatex, ""), collapse = " & "))
  paste0("\\begin{pmatrix}", paste(rows, collapse = " \\\\ "), "\\end{pmatrix}")
}
