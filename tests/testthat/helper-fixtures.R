# Shared fixtures for the tests: simulated data and models fitted at most once per
# test run (cached in `.fx`). The cache is read-only for the tests: never modify a
# fixture in place (copy it first).
#
#   fxData()              list(long, surv): simulated data frames (see below)
#   fxSmallData()         a tiny data frame (4 subjects x 3 times; id, site, t, x, y, g (factor),
#                         k, n, cnt, pt) for structural tests that need no fit
#   fxSpec(name)          the jsm() arguments of a reference model
#   fxForms(name)         collectFormulas() result of a reference model
#   fxModel(name, fit)    jsm() model, fitted with criteria = TRUE; fit = FALSE builds it
#                         without INLA (and without criteria); skips if INLA is missing or
#                         the fit fails (the failure is cached, the fit is not retried)
#   fxWarnings(name, fit) the warnings raised while building / fitting the fixture
#   fxPredict(name)       predict() of the fitted "multi" model on fxNewdata (cached)
#   fxSurvivalArm()       the "surv" model with a subject-level factor `arm` (A / B) in the
#                         survival data (not in any formula: a display column), fitted with
#                         config = TRUE (survival curves per level, with a band; E31)
#   withNullDevice(expr)  draw base graphics to a null device
#
# Reference models (name):
#   "beta"   y (beta) ~ t + scaled(x1);  x1 ~ t + (1 + t | id)            [one data frame]
#   "multi"  k | trials(n) (binomial, logit) ~ t + scaled(x1) + scaled(x2);
#            cnt | exposure(pt) (poisson) ~ grp + scaled(x2);
#            x1 ~ t + (1 + t | id);  x2 ~ t + (1 | id)  (x2 poisson)         [one data frame]
#   "surv"   time | event(status) (weibull) ~ 1 + scaled(x1);  x1 ~ t + (1 + t | id)
#                                                                          [dataOut / dataLong]
#   "priors" every prior kind, an outcome RE on `site`, two blocks on one group, a `||`
#            block, non-default family hyperparameters:
#            yg (gaussian, prec = pcSd) ~ t + re(1 | site, pcSd) + scaled(x1, normalCoef(0, 2))
#                                         + scaled(x2) + scaled(x3);
#            x1 ~ t + re(1 + t | id, wishartPrec(scale = "data"));
#            x2 (negbin, size = gammaPrec) ~ t + re(1 | id, halfCauchySd) + re(0 + t | id, halfNormalSd);
#            x3 (student t, prec = gammaPrec) ~ t + re(1 + t || id, gammaPrec)   [one data frame]

.fx <- new.env(parent = emptyenv())

fxData <- function() {
  if (!is.null(.fx$data)) return(.fx$data)
  set.seed(2026)
  n <- 40
  ids <- sprintf("S%02d", seq_len(n))
  long <- data.frame(id = rep(ids, each = 4), t = rep(c(0, 0.5, 1, 1.5), n), stringsAsFactors = FALSE)
  i <- match(long$id, ids)
  long$site <- factor(c("A", "B", "C", "D")[(i - 1L) %% 4L + 1L])
  b0 <- rnorm(n, sd = 0.5); b1 <- rnorm(n, sd = 0.2); c0 <- rnorm(n, sd = 0.3); c1 <- rnorm(n, sd = 0.15)
  s0 <- c(A = -0.2, B = 0.1, C = 0.3, D = -0.2)
  eta1 <- 1 + 0.5 * long$t + b0[i] + b1[i] * long$t
  eta2 <- 0.5 + 0.3 * long$t + c0[i]
  long$x1 <- eta1 + rnorm(nrow(long), sd = 0.2)
  long$x2 <- rpois(nrow(long), exp(eta2))
  long$x3 <- 0.5 - 0.2 * long$t + c0[i] + c1[i] * long$t + rnorm(nrow(long), sd = 0.3)
  mu <- plogis(-1 + 0.8 * eta1)
  long$y <- rbeta(nrow(long), 20 * mu, 20 * (1 - mu))
  long$yg <- 0.5 * eta1 + 0.3 * eta2 - 0.4 * long$x3 + s0[as.character(long$site)] + rnorm(nrow(long), sd = 0.3)
  long$n <- sample(5:10, nrow(long), replace = TRUE)
  long$k <- rbinom(nrow(long), long$n, plogis(-0.5 + 0.5 * eta1 + 0.3 * eta2))
  long$grp <- factor(sample(c("a", "b", "c"), nrow(long), replace = TRUE))
  long$pt <- runif(nrow(long), 0.5, 2)
  long$cnt <- rpois(nrow(long), long$pt * exp(-0.2 + 0.4 * eta2))
  rownames(long) <- NULL
  # survival: one row per subject, hazard driven by the baseline level of x1
  lam <- exp(-2.5 + 0.8 * (1 + b0))
  Tm <- (-log(runif(n)) / lam)^(1 / 1.5)
  Cn <- runif(n, 1, 4)
  surv <- data.frame(id = ids, t = 0, time = pmin(Tm, Cn), status = as.integer(Tm <= Cn),
                     stringsAsFactors = FALSE)
  .fx$data <- list(long = long, surv = surv)
  .fx$data
}

fxSmallData <- function() {
  if (!is.null(.fx$small)) return(.fx$small)
  set.seed(7)
  d <- data.frame(id = rep(c("S1", "S2", "S3", "S4"), each = 3), site = rep(c("A", "A", "B", "B"), each = 3),
                  t = rep(c(0, 1, 2), 4), stringsAsFactors = FALSE)
  d$x <- 1 + 0.5 * d$t + rnorm(12, sd = 0.3)
  d$y <- 0.5 * d$x + rnorm(12, sd = 0.3)
  d$g <- factor(c("u", "v", "w", "v", "w", "u", "w", "u", "v", "u", "u", "v"))  # not a function of t
  d$n <- 5L
  d$k <- rbinom(12, 5, 0.4)
  d$cnt <- rpois(12, 2)
  d$pt <- runif(12, 0.5, 2)
  .fx$small <- d
  d
}

fxNewdata <- data.frame(id = c("N1", "N2", "N3", "S01", "S01"), t = c(0, 1, 2, 2, 3),
                        n = 5L, pt = 1, grp = "a", stringsAsFactors = FALSE)

fxSpec <- function(name) {
  switch(name,
    beta = list(formOut = list(y ~ t + scaled(x1)),
                formLong = list(x1 ~ t + (1 + t | id)),
                data = fxData()$long, family = list(y = betaFamily())),
    multi = list(formOut = list(k | trials(n) ~ t + scaled(x1) + scaled(x2),
                                cnt | exposure(pt) ~ grp + scaled(x2)),
                 formLong = list(x1 ~ t + (1 + t | id), x2 ~ t + (1 | id)),
                 data = fxData()$long,
                 family = list(k = binomialFamily(), cnt = poissonFamily(), x2 = poissonFamily())),
    surv = list(formOut = list(time | event(status) ~ 1 + scaled(x1)),
                formLong = list(x1 ~ t + (1 + t | id)),
                dataOut = fxData()$surv, dataLong = fxData()$long[c("id", "t", "x1")],
                family = list(time = weibullSurvFamily())),
    priors = list(formOut = list(yg ~ t + re(1 | site, prior = pcSd(1)) + scaled(x1, prior = normalCoef(0, 2)) +
                                   scaled(x2) + scaled(x3)),
                  formLong = list(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")),
                                  x2 ~ t + re(1 | id, prior = halfCauchySd(1)) + re(0 + t | id, prior = halfNormalSd(0.5)),
                                  x3 ~ t + re(1 + t || id, prior = gammaPrec(1, 0.01))),
                  data = fxData()$long,
                  family = list(yg = gaussianFamily(prec = pcSd(1)),
                                x2 = negBinomialFamily(size = gammaPrec(1, 0.01)),
                                x3 = studentTFamily(prec = gammaPrec(1, 5e-5)))),
    stop("unknown fixture ", name))
}

fxForms <- function(name) {
  s <- fxSpec(name)
  collectFormulas(s$formOut, s$formLong, family = s$family)
}

fxModel <- function(name, fit = TRUE) {
  key <- paste0(name, if (fit) "_fit" else "_nofit")
  if (fit) testthat::skip_if_not_installed("INLA")
  if (is.null(.fx[[key]])) {
    s <- fxSpec(name)
    warns <- character()
    m <- withCallingHandlers(
      suppressMessages(do.call(jsm, c(s, list(fit = fit, criteria = fit)))),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    .fx[[paste0(key, "_warnings")]] <- warns
    .fx[[key]] <- if (fit && is.null(inlaFit(m))) "failed" else m
  }
  if (identical(.fx[[key]], "failed")) testthat::skip(paste("INLA fit of fixture", name, "failed"))
  .fx[[key]]
}

# An UNFITTED reference model with extra jsm() arguments (e.g. fxBuild("beta", threads = 2)).
fxBuild <- function(name, ...) suppressMessages(do.call(jsm, c(fxSpec(name), list(fit = FALSE), list(...))))

# A small balanced grid: n subjects (S01, ...) at the visit times `times`, with a seed; the
# responses are added by the test.
fxGridData <- function(n, times, seed = 1) {
  set.seed(seed)
  data.frame(id = rep(sprintf("S%02d", seq_len(n)), each = length(times)), t = rep(times, n),
             stringsAsFactors = FALSE)
}

fxSurvivalArm <- function() {
  testthat::skip_if_not_installed("INLA")
  if (is.null(.fx$survArm)) {
    s <- fxSpec("surv")
    s$dataOut$arm <- factor(rep(c("A", "B"), length.out = nrow(s$dataOut)))
    m <- suppressWarnings(suppressMessages(do.call(jsm, c(s, list(
      inlaOptions = list(control.compute = list(config = TRUE)))))))
    .fx$survArm <- if (is.null(inlaFit(m))) "failed" else m
  }
  if (identical(.fx$survArm, "failed")) testthat::skip("INLA fit of fixture surv + arm failed")
  .fx$survArm
}

fxWarnings <- function(name, fit = TRUE) {
  fxModel(name, fit)
  .fx[[paste0(name, if (fit) "_fit" else "_nofit", "_warnings")]]
}

fxPredict <- function(name = "multi") {
  key <- paste0(name, "_predict")
  if (is.null(.fx[[key]])) {
    m <- fxModel(name)
    .fx[[key]] <- suppressMessages(predict(m, fxNewdata))
  }
  .fx[[key]]
}

# Draw base graphics to a null device.
withNullDevice <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
