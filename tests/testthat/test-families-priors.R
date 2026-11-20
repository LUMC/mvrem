# Extended tests of R/families.R and R/priors.R: family constructors, control.family,
# families in collectFormulas() / prepareData() / buildInlaFormula(), response checks,
# prior constructors, INLA hyper lists, block-prior checks, the data-scaled Wishart prior,
# priors in the INLA formula and the fitted "priors" fixture.
# Helpers: helper-families.R (checkResp, halfCauchyLogDens, hyperOf, withWarnings,
# dataBlock, refFit, refBlock, expectBlockCov) and helper-fixtures.R.

# ---------------------------------------------------------------------------
# family constructors

test_that("family constructors: likelihood, kind, default link, aux and support", {
  spec <- list(
    list(gaussianFamily(), "gaussian", "gaussian", "identity", character(), "real"),
    list(studentTFamily(), "t", "t", "identity", character(), "real"),
    list(lognormalFamily(), "lognormal", "lognormal", "identity", character(), "positive"),
    list(gammaFamily(), "gamma", "gamma", "log", character(), "positive"),
    list(betaFamily(), "beta", "beta", "logit", character(), "unit"),
    list(poissonFamily(), "poisson", "poisson", "log", "exposure", "count"),
    list(negBinomialFamily(), "nbinomial", "nbinomial", "log", "exposure", "count"),
    list(binomialFamily(), "binomial", "binomial", "logit", "trials", "binomial"),
    list(betaBinomialFamily(), "betabinomial", "betabinomial", "logit", "trials", "binomial"),
    list(weibullSurvFamily(), "weibullsurv", "weibull", "log", "event", "survival"),
    list(exponentialSurvFamily(), "exponentialsurv", "exponential", "log", "event", "survival"),
    list(zeroInflatedPoissonFamily(), "zeroinflatedpoisson1", "zip", "log", "exposure", "count"),
    list(zeroInflatedNegBinomialFamily(), "zeroinflatednbinomial1", "zinb", "log", "exposure", "count"))
  for (s in spec) {
    f <- s[[1]]
    expect_s3_class(f, "mvremFamily")
    expect_equal(f$family, s[[2]])
    expect_equal(f$kind, s[[3]])
    expect_equal(f$link, s[[4]])
    expect_equal(f$linkDefault, s[[4]])
    expect_equal(f$aux, s[[5]])
    expect_equal(f$support, s[[6]])
    expect_length(f$hyper, 0)                        # NULL hyper arguments are dropped
    expect_equal(f$control, list())
    expect_equal(controlFamily(f), list())           # INLA defaults: nothing set
  }
})

test_that("family constructors: hyperparameter specs (names, accepted priors, LaTeX symbols)", {
  expect_named(studentTFamily()$specs, c("prec", "dof"))
  expect_equal(studentTFamily()$specs$dof$dists, character())
  expect_equal(gammaFamily()$specs$prec$symbol, "\\phi^{#r#}")
  expect_equal(betaFamily()$specs$phi$symbol, "\\phi^{#r#}")
  expect_equal(negBinomialFamily()$specs$size$dists, "gamma")
  expect_named(zeroInflatedNegBinomialFamily()$specs, c("size", "prob"))
  expect_equal(poissonFamily()$specs, list())
  expect_equal(binomialFamily()$specs, list())
  expect_equal(exponentialSurvFamily()$specs, list())
  expect_named(weibullSurvFamily()$specs, "alpha")
})

test_that("family constructors: precision-like hyperparameters accept the four scalar priors only", {
  scalar <- list(gammaPrec(1, 0.1), pcSd(1), halfNormalSd(2), halfCauchySd(3))
  ctors <- list(gaussian = function(h) gaussianFamily(prec = h),
                t = function(h) studentTFamily(prec = h),
                lognormal = function(h) lognormalFamily(prec = h),
                gamma = function(h) gammaFamily(prec = h),
                beta = function(h) betaFamily(phi = h))
  for (nm in names(ctors)) {
    for (pr in scalar) {
      f <- ctors[[nm]](pr)
      h <- names(f$hyper)
      expect_equal(controlFamily(f)$hyper[[h]], priorHyper(pr, 1L)$prec)
    }
    expect_error(ctors[[nm]](wishartPrec()), "does not accept this prior constructor; use gammaPrec\\(\\)")
    expect_error(ctors[[nm]](normalCoef()), "does not accept this prior constructor")
  }
})

test_that("negBinomialFamily() / zeroInflatedNegBinomialFamily(): `size` accepts gammaPrec() only", {
  expect_equal(controlFamily(negBinomialFamily(size = gammaPrec(2, 1)))$hyper$size,
               list(prior = "loggamma", param = c(2, 1)))
  expect_equal(controlFamily(zeroInflatedNegBinomialFamily(size = gammaPrec(3, 1)))$hyper$size,
               list(prior = "loggamma", param = c(3, 1)))
  for (pr in list(pcSd(1), halfNormalSd(1), halfCauchySd(1), wishartPrec(), normalCoef())) {
    expect_error(negBinomialFamily(size = pr), "use gammaPrec\\(\\) or a raw INLA entry")
    expect_error(zeroInflatedNegBinomialFamily(size = pr), "Hyperparameter `size` .* use gammaPrec\\(\\)")
  }
  f <- zeroInflatedNegBinomialFamily(0, size = gammaPrec(1, 1), prob = list(initial = -1))
  expect_equal(controlFamily(f), list(hyper = list(size = list(prior = "loggamma", param = c(1, 1)),
                                                   prob = list(initial = -1))))
})

test_that("family constructors: `rho`, `prob`, `dof` and `alpha` take raw INLA entries only", {
  raw <- list(prior = "gaussian", param = c(0, 1))
  for (pr in list(gammaPrec(), pcSd(1))) {
    expect_error(betaBinomialFamily(rho = pr), "`rho` .* use a raw INLA entry")
    expect_error(zeroInflatedPoissonFamily(prob = pr), "`prob` of family 'zeroinflatedpoisson1'")
    expect_error(zeroInflatedNegBinomialFamily(prob = pr), "`prob` of family 'zeroinflatednbinomial1'")
    expect_error(studentTFamily(dof = pr), "`dof` of family 't' does not accept")
    expect_error(weibullSurvFamily(alpha = pr), "`alpha` of family 'weibullsurv'")
  }
  expect_equal(controlFamily(betaBinomialFamily(rho = raw))$hyper$rho, raw)
  expect_equal(controlFamily(zeroInflatedPoissonFamily(prob = raw))$hyper$prob, raw)
  expect_equal(controlFamily(studentTFamily(dof = raw))$hyper$dof, raw)
  expect_equal(controlFamily(weibullSurvFamily(alpha = list(initial = 0.4, fixed = TRUE)))$hyper$alpha,
               list(initial = 0.4, fixed = TRUE))
})

test_that("family constructors: raw hyper lists accept all INLA fields and nothing else", {
  full <- list(prior = "loggamma", param = c(1, 1), initial = 2, fixed = FALSE,
               to.theta = function(x) log(x), from.theta = function(x) exp(x))
  expect_equal(controlFamily(gaussianFamily(prec = full))$hyper$prec, full)
  expect_equal(controlFamily(gaussianFamily(prec = list(initial = 4, fixed = TRUE)))$hyper$prec,
               list(initial = 4, fixed = TRUE))
  msg <- "`prec` of family 'gaussian' must be NULL, a prior constructor"
  expect_error(gaussianFamily(prec = list()), msg)
  expect_error(gaussianFamily(prec = list(1, 2)), msg)
  expect_error(gaussianFamily(prec = list(prior = "loggamma", foo = 1)), msg)
  expect_error(gaussianFamily(prec = list(prior = "loggamma", 1)), msg)
  expect_error(gaussianFamily(prec = "loggamma"), msg)
  expect_error(gaussianFamily(prec = 3), msg)
  expect_error(gaussianFamily(prec = c(prior = 1)), msg)
  expect_error(betaFamily(phi = NA), "`phi` of family 'beta' must be NULL")
  expect_error(weibullSurvFamily(alpha = list(shape = 1)), "`alpha` of family 'weibullsurv' must be NULL")
})

test_that("controlFamily(): control.link only for a non-default link", {
  expect_equal(controlFamily(gaussianFamily(link = "log")), list(control.link = list(model = "log")))
  expect_equal(controlFamily(gammaFamily(link = "log")), list())
  expect_equal(controlFamily(betaFamily(link = "logit")), list())
  expect_equal(controlFamily(binomialFamily(link = "cloglog")), list(control.link = list(model = "cloglog")))
  expect_equal(controlFamily(poissonFamily(link = "log")), list())
  expect_equal(controlFamily(betaBinomialFamily(link = "probit"))$control.link, list(model = "probit"))
  expect_equal(controlFamily(negBinomialFamily(link = "log")), list())
  expect_equal(controlFamily(weibullSurvFamily(link = "log")), list())
  expect_equal(controlFamily(exponentialSurvFamily()), list())
})

test_that("family constructors: `link` must be a single string", {
  expect_error(gaussianFamily(link = 1), "`link` must be a single character string")
  expect_error(gaussianFamily(link = c("identity", "log")), "single character string")
  expect_error(binomialFamily(link = NA_character_), "single character string")
})

test_that("family constructors: links are validated against INLA's list", {
  skip_if_not_installed("INLA")
  expect_error(poissonFamily(link = "logit"), "Link 'logit' is not available for INLA likelihood 'poisson'")
  expect_error(gammaFamily(link = "nonsense"), "not available for INLA likelihood 'gamma'")
  expect_error(betaFamily(link = "identity"), "not available for INLA likelihood 'beta'")
  expect_error(exponentialSurvFamily(link = "probit"), "not available for INLA likelihood 'exponentialsurv'")
  expect_error(weibullSurvFamily(link = "logit"), "not available for INLA likelihood 'weibullsurv'")
  # the error message lists the available links, without INLA's "default" entry
  err <- tryCatch(binomialFamily(link = "foo"), error = conditionMessage)
  expect_match(err, "available: ")
  expect_false(grepl("default", err))
  expect_match(err, "logit")
})

test_that("controlFamily(): `control` entries follow control.link and hyper", {
  f <- betaFamily(phi = gammaPrec(1, 0.1), link = "probit", control = list(beta.censor.value = 1e-4))
  expect_equal(controlFamily(f), list(control.link = list(model = "probit"),
                                      hyper = list(phi = list(prior = "loggamma", param = c(1, 0.1))),
                                      beta.censor.value = 1e-4))
  f <- betaBinomialFamily(rho = list(prior = "gaussian", param = c(0, 1)), control = list(variant = 1))
  expect_equal(controlFamily(f), list(hyper = list(rho = list(prior = "gaussian", param = c(0, 1))), variant = 1))
  expect_no_warning(cf <- controlFamily(gaussianFamily(prec = pcSd(1), control = list(variant = 0))))
  expect_equal(cf$variant, 0)
})

test_that("controlFamily(): a control hyper adding fields to a prior is merged without warning", {
  f <- gaussianFamily(prec = pcSd(1), control = list(hyper = list(prec = list(initial = 2))))
  expect_no_warning(cf <- controlFamily(f))
  expect_equal(cf$hyper$prec, list(prior = "pc.prec", param = c(1, 0.01), initial = 2))
})

test_that("controlFamily(): a control hyper of another parameter is added without warning", {
  f <- studentTFamily(prec = gammaPrec(), control = list(hyper = list(dof = list(prior = "normal", param = c(0, 1)))))
  expect_no_warning(cf <- controlFamily(f))
  expect_named(cf$hyper, c("prec", "dof"))
  expect_equal(cf$hyper$dof, list(prior = "normal", param = c(0, 1)))
})

test_that("controlFamily(): a control hyper overriding `prior` / `param` of a prior warns and wins", {
  f <- gaussianFamily(prec = pcSd(1), control = list(hyper = list(prec = list(prior = "loggamma"))))
  expect_warning(cf <- controlFamily(f), "`hyper` overrides the `prior` of family 'gaussian' \\(prec\\)")
  expect_equal(cf$hyper$prec$prior, "loggamma")
  expect_warning(controlFamily(gammaFamily(prec = gammaPrec(), control = list(hyper = list(prec = list(param = c(1, 1)))))),
                 "overrides the `prior` of family 'gamma' \\(prec\\)")
  # several overridden parameters are listed together
  f <- zeroInflatedNegBinomialFamily(size = gammaPrec(), prob = list(initial = 0),
                                     control = list(hyper = list(size = list(param = c(1, 1)),
                                                                 prob = list(prior = "gaussian"))))
  expect_warning(controlFamily(f), "family 'zeroinflatednbinomial1' \\(size, prob\\)")
})

test_that("controlFamily(): a control.link in `control` overrides the `link` argument", {
  f <- binomialFamily(link = "probit", control = list(control.link = list(model = "logit")))
  expect_equal(controlFamily(f)$control.link, list(model = "logit"))
})

test_that("family constructors: `control` must be a named list", {
  for (bad in list("a", 1, list(1), list(a = 1, 2), c(a = 1)))
    expect_error(gaussianFamily(control = bad), "`control` must be a named list of control.family entries")
  expect_error(weibullSurvFamily(control = list(1)), "`control` must be a named list")
  expect_error(inlaFamily("gaussian", control = list(1)), "`control` must be a named list")
})

test_that("zeroInflatedPoissonFamily() / zeroInflatedNegBinomialFamily(): type 0 / 1, invalid types", {
  for (type in list(0, 1, 0L, 1L)) {
    f <- zeroInflatedPoissonFamily(type)
    expect_equal(f$family, paste0("zeroinflatedpoisson", type))
    expect_identical(f$type, as.integer(type))
    g <- zeroInflatedNegBinomialFamily(type)
    expect_equal(g$family, paste0("zeroinflatednbinomial", type))
    expect_identical(g$type, as.integer(type))
  }
  for (bad in list(2, -1, 0.5, "1", TRUE, NA, c(0, 1), NULL, 1 + 0i)) {
    expect_error(zeroInflatedPoissonFamily(bad), "`type` must be 0 or 1")
    expect_error(zeroInflatedNegBinomialFamily(bad), "`type` must be 0 or 1")
  }
  expect_equal(zeroInflationType(0), 0L)
})

test_that("controlFamily(): rejects non-family objects", {
  expect_error(controlFamily(list(family = "gaussian")), "`x` must be a family object")
  expect_error(controlFamily("gaussian"), "`x` must be a family object")
  expect_error(controlFamily(gammaPrec()), "`x` must be a family object")
})

test_that("format() / print() of families", {
  expect_equal(format(gaussianFamily()), "<mvrem family> gaussian, link = identity; control.family = list()")
  txt <- format(betaFamily(phi = gammaPrec(1, 0.1), link = "probit"))
  expect_match(txt, "^<mvrem family> beta, link = probit; control.family = list\\(")
  expect_match(txt, 'control.link = list(model = "probit")', fixed = TRUE)
  expect_match(txt, 'hyper = list(phi = list(prior = "loggamma", param = c(1, 0.1)))', fixed = TRUE)
  expect_match(format(zeroInflatedPoissonFamily(0)), "^<mvrem family> zeroinflatedpoisson0, link = log;")
  f <- weibullSurvFamily()
  expect_output(res <- withVisible(print(f)), "^<mvrem family> weibullsurv, link = log; control.family = list\\(\\)$")
  expect_false(res$visible)
  expect_identical(res$value, f)
})

# ---------------------------------------------------------------------------
# inlaFamily() and the INLA family list

test_that("inlaFamily(): `family` must be a single name, `hyper` a named list", {
  for (bad in list(c("gaussian", "t"), NA_character_, 1, NULL, list("gaussian")))
    expect_error(inlaFamily(bad), "`family` must be a single INLA likelihood name")
  expect_error(inlaFamily("gaussian", hyper = "prec"), "`hyper` must be a named list")
  expect_error(inlaFamily("gaussian", hyper = list(1)), "`hyper` must be a named list")
  expect_error(inlaFamily("gaussian", hyper = list(prec = gammaPrec(), gammaPrec())), "`hyper` must be a named list")
})

test_that("inlaFamily(): generic kind, no link, both inputs, any support", {
  skip_if_not_installed("INLA")
  f <- inlaFamily("logistic")
  expect_null(f$link)
  expect_null(f$linkDefault)
  expect_equal(f$kind, "inla")
  expect_equal(f$support, "any")
  expect_equal(f$aux, c("trials", "exposure"))
  expect_equal(controlFamily(f), list())
  expect_equal(format(f), "<mvrem family> logistic; control.family = list()")   # no link shown
})

test_that("inlaFamily(): hyper entries of any name, scalar priors or raw lists; NULL dropped", {
  skip_if_not_installed("INLA")
  f <- inlaFamily("gaussian", hyper = list(theta_1 = list(initial = 1)))
  expect_equal(f$specs$theta_1$symbol, "\\theta^{#r#}_{\\mathrm{theta\\_1}}")   # generic symbol, escaped
  expect_equal(controlFamily(f)$hyper, list(theta_1 = list(initial = 1)))
  expect_equal(controlFamily(inlaFamily("gaussian", hyper = list(prec = halfNormalSd(1))))$hyper$prec,
               list(prior = "logtnormal", param = c(0, 1)))
  expect_error(inlaFamily("gaussian", hyper = list(prec = wishartPrec())), "does not accept this prior constructor")
  expect_error(inlaFamily("gaussian", hyper = list(prec = normalCoef())), "does not accept")
  expect_error(inlaFamily("gaussian", hyper = list(prec = list(foo = 1))), "must be NULL, a prior constructor")
  expect_length(inlaFamily("gaussian", hyper = list(prec = NULL))$hyper, 0)
})

test_that("inlaFamily(): a link is always passed on (INLA's default is unknown to mvrem)", {
  skip_if_not_installed("INLA")
  expect_equal(controlFamily(inlaFamily("gaussian", link = "identity"))$control.link, list(model = "identity"))
  expect_equal(controlFamily(inlaFamily("logistic", link = "identity", hyper = list(prec = gammaPrec(1, 1)))),
               list(control.link = list(model = "identity"),
                    hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))))
  expect_error(inlaFamily("gaussian", link = "nonsense"), "not available for INLA likelihood 'gaussian'")
})

test_that("inlaFamily(): unknown and survival likelihoods are rejected", {
  skip_if_not_installed("INLA")
  expect_error(inlaFamily("weibullsurv"), "Survival likelihood 'weibullsurv' \\(inlaFamily\\(\\)\\) is not supported")
  expect_error(inlaFamily("exponentialsurv"), "is not supported")
  expect_error(inlaFamily("coxph"), "use weibullSurvFamily\\(\\) or exponentialSurvFamily\\(\\)")
  expect_error(inlaFamily("nonsense"), "Unknown INLA likelihood family in inlaFamily\\(\\): nonsense")
})

# The branch taken when INLA is NOT installed (links not validated, a warning instead of
# the family-name check) cannot be exercised here: familyLinks() and checkInlaFamily() call
# requireNamespace("INLA") directly, without a wrapper that local_mocked_bindings() could
# replace (testthat cannot create a binding for a base function in the package namespace).

test_that("inlaFamilies(): lists INLA's likelihoods", {
  skip_if_not_installed("INLA")
  fams <- inlaFamilies()
  expect_type(fams, "character")
  expect_true(all(c("gaussian", "t", "lognormal", "gamma", "beta", "poisson", "nbinomial",
                    "binomial", "betabinomial", "zeroinflatedpoisson0", "zeroinflatedpoisson1",
                    "zeroinflatednbinomial0", "zeroinflatednbinomial1", "weibullsurv",
                    "exponentialsurv") %in% fams))
})

test_that("familyFromName(): mapped names give the constructors with their defaults", {
  mapped <- list(gaussian = gaussianFamily(), t = studentTFamily(), lognormal = lognormalFamily(),
                 gamma = gammaFamily(), beta = betaFamily(), poisson = poissonFamily(),
                 nbinomial = negBinomialFamily(), binomial = binomialFamily(),
                 betabinomial = betaBinomialFamily(),
                 zeroinflatedpoisson0 = zeroInflatedPoissonFamily(0),
                 zeroinflatedpoisson1 = zeroInflatedPoissonFamily(1),
                 zeroinflatednbinomial0 = zeroInflatedNegBinomialFamily(0),
                 zeroinflatednbinomial1 = zeroInflatedNegBinomialFamily(1),
                 weibullsurv = weibullSurvFamily(), exponentialsurv = exponentialSurvFamily())
  for (nm in names(mapped)) {
    f <- familyFromName(nm, "`family`")
    expect_equal(f, mapped[[nm]])
    expect_equal(f$family, nm)
  }
  expect_equal(familyFromName("exponentialsurv", "x")$kind, "exponential")
  expect_equal(defaultFamily(), gaussianFamily())
})

test_that("familyFromName(): other names give inlaFamily(); unknown / survival names are errors naming `where`", {
  skip_if_not_installed("INLA")
  expect_no_warning(f <- familyFromName("logistic", "`family`"))
  expect_equal(f, inlaFamily("logistic"))
  expect_error(familyFromName("nonsense", "`familyOut`"),
               "Unknown INLA likelihood family in `familyOut`: nonsense\\. See inlaFamilies\\(\\)")
  expect_error(familyFromName("coxph", "`family`"), "Survival likelihood 'coxph' \\(`family`\\)")
})

test_that("familyLinks() / checkInlaFamily(): use INLA's model list", {
  skip_if_not_installed("INLA")
  expect_true(all(c("logit", "probit", "cloglog") %in% familyLinks("binomial")))
  expect_true("log" %in% familyLinks("poisson"))
  expect_null(familyLinks("nonsense"))
  expect_invisible(checkInlaFamily("gaussian", "x"))
})

test_that("collectFormulas(): INLA names in `family` arguments are validated", {
  skip_if_not_installed("INLA")
  expect_error(collectFormulas(y ~ t + scaled(x), x ~ t, family = c(y = "betta")),
               "Unknown INLA likelihood family in `family`: betta")
  expect_error(collectFormulas(y ~ t + scaled(x), x ~ t, familyLong = c(x = "coxph")),
               "Survival likelihood 'coxph' \\(`familyLong`\\) is not supported")
  # a survival name maps to its constructor, which then needs event() on the left-hand side
  expect_error(collectFormulas(y ~ t + scaled(x), x ~ t, family = c(y = "weibullsurv")),
               "Survival family 'weibullsurv' of 'y' needs the event indicator")
})

# ---------------------------------------------------------------------------
# checkResponse()

test_that("checkResponse(): no family or no rows -> nothing checked", {
  p <- parseJsmFormula(y ~ t, "out")
  expect_null(checkResponse(p, data.frame(y = "a", t = 1)))
  expect_null(checkResp(betaFamily(), y ~ t, data.frame(y = numeric(), t = numeric())))
  expect_invisible(checkResp(gaussianFamily(), y ~ t, data.frame(y = 1, t = 1)))
})

test_that("checkResponse(): real support accepts any number, logical, NA; rejects non-numeric", {
  expect_silent(checkResp(gaussianFamily(), y ~ t, data.frame(y = c(-1e10, 0, 1e10, NA), t = 1)))
  expect_silent(checkResp(studentTFamily(), y ~ t, data.frame(y = c(TRUE, FALSE), t = 1)))
  expect_error(checkResp(gaussianFamily(), y ~ t, data.frame(y = c("a", "b"), t = 1)),
               "Response 'y' \\(family 'gaussian'\\) must be numeric")
  expect_error(checkResp(gaussianFamily(), y ~ t, data.frame(y = factor(c("a", "b")), t = 1)), "must be numeric")
  # the transformed response is checked, the message names the LHS as written
  expect_error(checkResp(gammaFamily(), log(y) ~ t, data.frame(y = c(1, 2), t = 1)),
               "Response 'log\\(y\\)' \\(family 'gamma'\\) must be positive \\(1 value\\(s\\) not\\)")
})

test_that("checkResponse(): positive support (gamma, lognormal)", {
  expect_silent(checkResp(gammaFamily(), y ~ t, data.frame(y = c(1e-10, 5, NA), t = 1)))
  expect_error(checkResp(gammaFamily(), y ~ t, data.frame(y = c(0, -1, 2), t = 1)),
               "must be positive \\(2 value\\(s\\) not\\)")
  expect_error(checkResp(lognormalFamily(), y ~ t, data.frame(y = 0, t = 1)), "family 'lognormal'\\) must be positive")
})

test_that("checkResponse(): unit support excludes 0 and 1", {
  expect_silent(checkResp(betaFamily(), y ~ t, data.frame(y = c(1e-6, 0.5, 1 - 1e-6, NA), t = 1)))
  expect_error(checkResp(betaFamily(), y ~ t, data.frame(y = 0, t = 1)), "must be in \\(0, 1\\)")
  expect_error(checkResp(betaFamily(), y ~ t, data.frame(y = 1, t = 1)), "must be in \\(0, 1\\)")
  expect_error(checkResp(betaFamily(), y ~ t, data.frame(y = c(TRUE, FALSE), t = 1)),
               "must be in \\(0, 1\\) \\(2 value\\(s\\) not\\)")
})

test_that("checkResponse(): count support", {
  for (fam in list(poissonFamily(), negBinomialFamily(), zeroInflatedPoissonFamily(0),
                   zeroInflatedNegBinomialFamily(1))) {
    expect_silent(checkResp(fam, y ~ t, data.frame(y = c(0, 1, 1000, 3 + 1e-10, NA), t = 1)))
    expect_error(checkResp(fam, y ~ t, data.frame(y = c(2.5, -1, 0), t = 1)),
                 "must be non-negative integers \\(2 value\\(s\\) not\\)")
  }
  expect_silent(checkResp(poissonFamily(), y ~ t, data.frame(y = c(TRUE, FALSE), t = 1)))
})

test_that("checkResponse(): an infinite count or trials value is rejected", {
  expect_error(checkResp(poissonFamily(), y ~ t, data.frame(y = c(1, Inf), t = 1)),
               "must be non-negative integers \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, data.frame(y = c(1, 2), n = c(3, Inf), t = 1)),
               "trials\\(\\) of 'y' must be positive integers")
  # -Inf and NaN are not missing values either; NA still is (B8)
  expect_error(checkResp(poissonFamily(), y ~ t, data.frame(y = c(-Inf, NaN, NA, 2), t = 1)),
               "must be non-negative integers \\(2 value\\(s\\) not\\)")
  expect_error(checkResp(poissonFamily(), y | exposure(E) ~ t, data.frame(y = 1:2, E = c(NaN, 1), t = 1)),
               "exposure\\(\\) of 'y' must be positive \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, data.frame(y = c(Inf, 2), n = c(3, 3), t = 1)),
               "must be integers between 0 and trials\\(\\) \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(gammaFamily(), y ~ t, data.frame(y = c(1, Inf), t = 1)),
               "must be positive \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(gaussianFamily(), y ~ t, data.frame(y = c(1, Inf, NaN, NA), t = 1)),
               "Response 'y' \\(family 'gaussian'\\) must be finite \\(2 value\\(s\\) not\\)")
  expect_silent(checkResp(gaussianFamily(), y ~ t, data.frame(y = c(TRUE, NA), t = 1)))
})

test_that("checkResponse(): binomial without trials() is 0 / 1", {
  for (fam in list(binomialFamily(), betaBinomialFamily())) {
    expect_silent(checkResp(fam, y ~ t, data.frame(y = c(0, 1, NA), t = 1)))
    expect_silent(checkResp(fam, y ~ t, data.frame(y = c(TRUE, FALSE), t = 1)))
    expect_error(checkResp(fam, y ~ t, data.frame(y = c(0.5, 2, 1), t = 1)),
                 "must be 0 or 1 \\(use trials\\(\\) for numbers of successes\\) \\(2 value\\(s\\) not\\)")
  }
})

test_that("checkResponse(): binomial with trials()", {
  ok <- data.frame(y = c(0, 3, 5, NA), n = c(1, 3, 10, 2), t = 1)
  expect_silent(checkResp(binomialFamily(), y | trials(n) ~ t, ok))           # y = trials accepted
  expect_silent(checkResp(betaBinomialFamily(), y | trials(2 * n) ~ t, ok))   # expression of the data
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, y = c(0, 4, 5, NA))),
               "must be integers between 0 and trials\\(\\) \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, y = c(-1, 0.5, 5, NA))),
               "between 0 and trials\\(\\) \\(2 value\\(s\\) not\\)")
  # trials: positive integers, numeric; NA passes
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, n = c(0, 3, 10, 2))),
               "trials\\(\\) of 'y' must be positive integers \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, n = c(1, 3.5, 10, 2))),
               "trials\\(\\) of 'y' must be positive integers")
  expect_error(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, n = c("a", "b", "c", "d"))),
               "trials\\(\\) of 'y' must be numeric")
  expect_silent(checkResp(binomialFamily(), y | trials(n) ~ t, transform(ok, n = c(NA, 3, 10, 2))))
})

test_that("checkResponse(): exposure() must be positive and finite", {
  ok <- data.frame(y = c(0, 1, 2), pt = c(0.1, 1, 100), t = 1)
  expect_silent(checkResp(poissonFamily(), y | exposure(pt) ~ t, ok))
  for (bad in list(c(0, 1, 2), c(-1, 1, 2), c(Inf, 1, 2)))
    expect_error(checkResp(negBinomialFamily(), y | exposure(pt) ~ t, transform(ok, pt = bad)),
                 "exposure\\(\\) of 'y' must be positive \\(1 value\\(s\\) not\\)")
  expect_error(checkResp(poissonFamily(), y | exposure(pt) ~ t, transform(ok, pt = c("a", "b", "c"))),
               "exposure\\(\\) of 'y' must be numeric")
  # the aux check runs before the response check
  expect_error(checkResp(poissonFamily(), y | exposure(pt) ~ t, data.frame(y = -1, pt = 0, t = 1)),
               "exposure\\(\\)")
})

test_that("checkResponse(): survival times and event status", {
  for (fam in list(weibullSurvFamily(), exponentialSurvFamily())) {
    ok <- data.frame(time = c(0.1, 2, 0, NA), st = c(1, 0, 0, 1), t = 1)
    expect_silent(checkResp(fam, time | event(st) ~ t, ok))                    # 0 censored accepted
    expect_silent(checkResp(fam, time | event(st) ~ t, transform(ok, st = c(TRUE, FALSE, FALSE, TRUE))))
    expect_error(checkResp(fam, time | event(st) ~ t, transform(ok, st = c(1, 0, 1, 1))),
                 "must be positive \\(0 only for a censored time\\) \\(1 value\\(s\\) not\\)")
    expect_error(checkResp(fam, time | event(st) ~ t, transform(ok, time = c(-1, -2, 0, 1))),
                 "0 only for a censored time\\) \\(2 value\\(s\\) not\\)")
    expect_error(checkResp(fam, time | event(st) ~ t, transform(ok, st = c(1, 0.5, 0, 2))),
                 "event\\(\\) of 'time' must be 0 \\(censored\\) or 1 \\(event\\) \\(2 value\\(s\\) not\\)")
    expect_error(checkResp(fam, time | event(st) ~ t, transform(ok, st = c("1", "0", "0", "1"))),
                 "event\\(\\) of 'time' must be numeric")
  }
})

test_that("checkResponse(): support 'any' skips the response but checks inputs", {
  skip_if_not_installed("INLA")
  expect_silent(checkResp(inlaFamily("logistic"), y ~ t, data.frame(y = c("a", "b"), t = 1)))
  expect_error(checkResp(inlaFamily("binomial"), y | trials(n) ~ t, data.frame(y = 1, n = 0, t = 1)),
               "trials\\(\\) of 'y' must be positive integers")
  expect_error(checkResp(inlaFamily("poisson"), y | exposure(e) ~ t, data.frame(y = 1, e = -1, t = 1)),
               "exposure\\(\\) of 'y' must be positive")
})

test_that("prepareData(): responses, trials and exposure are checked against the family", {
  d <- fxSmallData()
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x) + scaled(cnt)),
                       list(x ~ t + (1 | id), cnt | exposure(pt) ~ t),
                       family = list(k = binomialFamily(), cnt = poissonFamily()))
  expect_silent(prepareData(f, d))
  bad <- function(col, value) {
    dd <- d
    dd[[col]][1] <- value
    prepareData(f, dd)
  }
  expect_error(bad("k", 6), "Response 'k' \\(family 'binomial'\\) must be integers between 0 and trials\\(\\)")
  expect_error(bad("n", 0), "trials\\(\\) of 'k' must be positive integers")
  expect_error(bad("cnt", -1), "Response 'cnt' \\(family 'poisson'\\) must be non-negative integers")
  expect_error(bad("pt", -1), "exposure\\(\\) of 'cnt' must be positive")
  expect_silent(bad("k", NA))                      # a missing response is predicted, not rejected
  # a survival outcome: status 0 / 1, positive times (0 only if censored)
  ds <- data.frame(id = 1:4, time = c(1, 2, 0.5, 3), status = c(1, 0, 1, 1), x = 1:4, t = c(0, 1, 0, 1))
  fs <- collectFormulas(list(time | event(status) ~ 1 + scaled(x)), list(x ~ t + (1 | id)),
                        family = list(time = weibullSurvFamily()))
  expect_silent(prepareData(fs, ds))
  expect_error(prepareData(fs, transform(ds, status = c(1, 2, 1, 1))), "event\\(\\) of 'time' must be 0 \\(censored\\)")
  expect_error(prepareData(fs, transform(ds, time = c(0, 2, 0.5, 3))), "Response 'time'.*0 only for a censored time")
})

# ---------------------------------------------------------------------------
# families in the INLA formula

test_that("buildInlaFormula(): family hyper and link reach controlFamily; exposure is NA in zero-equation rows", {
  d <- fxSmallData()
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x) + scaled(cnt)),
                       list(x ~ t + (1 | id), cnt | exposure(pt) ~ t),
                       family = list(k = binomialFamily(link = "probit"),
                                     x = gaussianFamily(prec = pcSd(1)),
                                     cnt = zeroInflatedPoissonFamily()))
  inl <- buildInlaFormula(f, prepareData(f, d))
  expect_equal(inl$families, c(k_k = "binomial", x_x = "gaussian", cnt_cnt = "zeroinflatedpoisson1",
                               x_x_0eq = "gaussian", cnt_cnt_0eq = "gaussian"))
  expect_equal(inl$controlFamily$k_k, list(control.link = list(model = "probit")))
  expect_equal(inl$controlFamily$x_x, list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))))
  expect_equal(inl$controlFamily$cnt_cnt, list())
  expect_equal(inl$stacks$cnt$main$cnt_E, d$pt)
  expect_true(all(is.na(inl$stacks$cnt$k$cnt_E)))
  expect_named(inl$stacks$cnt$k, names(inl$stacks$cnt$main))
  # the family lines of the formatted formula
  txt <- formatInlaFormula(inl, nrows = 1)
  expect_true(any(grepl('# family = "binomial", control.family = list(control.link = list(model = "probit"))',
                        txt, fixed = TRUE)))
  expect_true(any(grepl('# family = "zeroinflatedpoisson1"$', txt)))
  # naming of the input columns
  inl2 <- buildInlaFormula(f, prepareData(f, d), naming = namingScheme(exposure = "E_{resp}"))
  expect_equal(inl2$inputs$cnt, c(E = "E_cnt"))
  expect_equal(names(inl2$stacks$cnt$main)[1:3], c("cnt_cnt", "cnt_cnt_0eq", "E_cnt"))
  expect_error(buildInlaFormula(f, prepareData(f, d), naming = namingScheme(exposure = "{resp}_1")),
               "duplicated INLA names: cnt_1")
})

# ---------------------------------------------------------------------------
# prior constructors

test_that("prior constructors: class, dist, target and arguments", {
  expect_equal(unclass(gammaPrec(2, 3)), list(dist = "gamma", target = "prec", args = list(shape = 2, rate = 3)))
  expect_equal(gammaPrec()$args, list(shape = 1, rate = 5e-5))
  expect_equal(unclass(pcSd(0.5, 0.05)), list(dist = "pc", target = "sd", args = list(u = 0.5, alpha = 0.05)))
  expect_equal(pcSd(2)$args$alpha, 0.01)
  expect_equal(unclass(halfNormalSd(2)), list(dist = "halfNormal", target = "sd", args = list(scale = 2)))
  expect_equal(unclass(halfCauchySd(25)), list(dist = "halfCauchy", target = "sd", args = list(scale = 25)))
  expect_equal(unclass(normalCoef(-1, 2)), list(dist = "normal", target = "coef", args = list(mean = -1, sd = 2)))
  expect_equal(normalCoef()$args, list(mean = 0, sd = 10))
  expect_equal(defaultCoefPrior(), normalCoef(0, 10))
  expect_equal(gammaPrec(1L, 2L)$args, list(shape = 1L, rate = 2L))   # integers are numeric
  expect_error(pcSd(), "argument \"u\" is missing")
})

test_that("prior constructors: positive arguments reject zero, negative, NA, Inf, non-numeric and vectors", {
  bads <- list(0, -1, NA, NA_real_, Inf, -Inf, NaN, "1", TRUE, c(1, 2), numeric(), NULL, list(1))
  for (bad in bads) {
    expect_error(gammaPrec(shape = bad), "`shape` must be a single positive number")
    expect_error(gammaPrec(rate = bad), "`rate` must be a single positive number")
    expect_error(pcSd(u = bad), "`u` must be a single positive number")
    expect_error(halfNormalSd(bad), "`scale` must be a single positive number")
    expect_error(halfCauchySd(bad), "`scale` must be a single positive number")
    expect_error(normalCoef(sd = bad), "`sd` must be a single positive number")
  }
  for (bad in bads[!vapply(bads, is.null, TRUE)])                 # df = NULL is the default (k + 2)
    expect_error(wishartPrec(df = bad), "`df` must be a single positive number")
})

test_that("pcSd(): alpha in (0, 1); normalCoef(): a finite mean", {
  for (bad in list(0, 1, -0.1, 1.5, NA, Inf, "0.5", c(0.1, 0.2), numeric()))
    expect_error(pcSd(1, alpha = bad), "`alpha` must be in \\(0, 1\\)")
  expect_equal(pcSd(1, alpha = 1e-10)$args$alpha, 1e-10)
  for (bad in list(Inf, -Inf, NA, NaN, "0", c(0, 1), numeric(), TRUE))
    expect_error(normalCoef(mean = bad), "`mean` must be a single finite number")
  expect_equal(normalCoef(mean = -3)$args$mean, -3)
})

test_that("wishartPrec(): df, scalar, matrix and \"data\" scales", {
  expect_equal(unclass(wishartPrec()), list(dist = "wishart", target = "precMatrix", args = list(df = NULL, scale = 1)))
  expect_equal(wishartPrec(df = 2.5)$args$df, 2.5)               # checked against k only per block
  expect_equal(wishartPrec(scale = 0.3)$args$scale, 0.3)
  expect_equal(wishartPrec(scale = "data")$args$scale, "data")
  # matrices: dimnames dropped, 1 x 1 accepted, integer matrices accepted
  V <- matrix(c(2, 0.5, 0.5, 1), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_identical(wishartPrec(scale = V)$args$scale, unname(V))
  expect_equal(wishartPrec(scale = matrix(2))$args$scale, matrix(2))
  expect_equal(wishartPrec(scale = diag(2L))$args$scale, diag(2L))
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2)))
    expect_error(wishartPrec(scale = bad), "`scale` must be a single positive number")
  for (bad in list("Data", "foo", c("data", "data"), NA_character_, "", matrix(c("1", "0", "0", "1"), 2)))
    expect_error(wishartPrec(scale = bad), "`scale` must be a matrix, a positive number or \"data\"")
  for (bad in list(matrix(1, 2, 3), matrix(c(1, NA, NA, 1), 2), matrix(c(TRUE, FALSE, FALSE, TRUE), 2),
                   matrix(c(1, 0, 0, 0), 2), matrix(c(1, 0, 0, -1), 2),         # singular, indefinite
                   matrix(c(1, 2, 2, 1), 2),                                     # not positive definite
                   matrix(c(1, 0.5, 0.4, 1), 2), matrix(c(1, 2, 0, 1), 2)))      # not symmetric
    expect_error(wishartPrec(scale = bad), "`scale` must be a symmetric positive-definite matrix")
})

test_that("wishartPrec(): a matrix with infinite entries is rejected with the package's message", {
  for (bad in list(matrix(c(1, Inf, Inf, 1), 2), matrix(c(-Inf, 0, 0, 1), 2), matrix(c(1, NaN, NaN, 1), 2)))
    expect_error(wishartPrec(scale = bad), "`scale` must be a symmetric positive-definite matrix")
})

test_that("format() / print() of priors", {
  txt <- format(gammaPrec())
  expect_match(txt, "^<mvrem prior> Gamma prior on precision: tau ~ Gamma\\(shape = 1, rate = 5e-05\\)")
  expect_match(txt, "[INLA: loggamma, param = c(1, 5e-05)]", fixed = TRUE)
  expect_match(format(gammaPrec(1.23456, 2)), "shape = 1.235, rate = 2", fixed = TRUE)   # 4 significant digits
  # sigma ~ Exp(lambda) with lambda = -log(alpha) / u = -log(0.05) / 2 = 1.498
  expect_match(format(pcSd(2, 0.05)), "P(sigma > 2) = 0.05, i.e. sigma ~ Exp(1.498)", fixed = TRUE)
  expect_match(format(pcSd(2, 0.05)), "[INLA: pc.prec, param = c(2, 0.05)]", fixed = TRUE)
  expect_match(format(halfCauchySd(5)), "Half-Cauchy prior on SD: sigma ~ C+(0, 5)", fixed = TRUE)
  expect_match(format(halfCauchySd(5)), "[INLA: expression: prior on log precision]", fixed = TRUE)
  # logtnormal precision = 1 / scale^2 = 1 / 4
  expect_match(format(halfNormalSd(2)), "sigma ~ N+(0, sd = 2)   [INLA: logtnormal, param = c(0, 0.25)]", fixed = TRUE)
  # beta precision = 1 / sd^2 = 4
  expect_match(format(normalCoef(1, 0.5)), "gamma ~ N(mean = 1, sd = 0.5)   [INLA: beta, normal, param = c(1, 4)]",
               fixed = TRUE)
  expect_match(format(wishartPrec()), "df = k + 2, scale = 1 * I", fixed = TRUE)
  expect_match(format(wishartPrec(df = 5, scale = 2)), "Q ~ W(df = 5, scale = 2 * I)", fixed = TRUE)
  expect_match(format(wishartPrec(df = 4, scale = diag(c(1, 0.1)))), "scale = diag(1, 0.1)", fixed = TRUE)
  expect_match(format(wishartPrec(scale = matrix(c(1, 0.5, 0.5, 2), 2))),
               "scale = rbind(c(1, 0.5), c(0.5, 2))", fixed = TRUE)
  expect_match(format(wishartPrec(scale = "data")), "df = k + 2, scale = \"data\" (lme4 fit", fixed = TRUE)
  p <- pcSd(1)
  expect_output(res <- withVisible(print(p)), "^<mvrem prior> PC prior on SD")
  expect_false(res$visible)
  expect_identical(res$value, p)
})

# ---------------------------------------------------------------------------
# INLA hyper lists

test_that("priorHyper(): scalar priors, independent of k", {
  for (k in c(1L, 3L)) {
    expect_equal(priorHyper(gammaPrec(2, 0.5), k), list(prec = list(prior = "loggamma", param = c(2, 0.5))))
    expect_equal(priorHyper(pcSd(3, 0.1), k), list(prec = list(prior = "pc.prec", param = c(3, 0.1))))
    # logtnormal precision = 1 / scale^2 = 4
    expect_equal(priorHyper(halfNormalSd(0.5), k), list(prec = list(prior = "logtnormal", param = c(0, 4))))
  }
  h <- priorHyper(halfCauchySd(1 / 3), 1L)
  expect_named(h, "prec")
  expect_named(h$prec, "prior")
  expect_equal(h$prec$prior, halfCauchyExpression(1 / 3))
  expect_match(h$prec$prior, "s = 0.333333333333333;", fixed = TRUE)    # 15 significant digits
})

test_that("halfCauchySd(): the INLA expression is the half-Cauchy density of sigma on log precision", {
  for (s in c(0.1, 1, 25)) {
    ld <- halfCauchyLogDens(s)
    # a proper density of theta = log precision
    dth <- function(th) vapply(th, function(x) exp(ld(x)), 0)
    th0 <- -2 * log(s)                                     # theta at sigma = s
    tot <- integrate(dth, th0 - 80, th0 + 80, rel.tol = 1e-10)$value
    expect_equal(tot, 1, tolerance = 1e-5)
    # transformed back to sigma = exp(-theta / 2) (|d theta / d sigma| = 2 / sigma): 2 * dcauchy
    sig <- c(0.01, 0.5, 1, 3, 100) * s
    dens <- vapply(-2 * log(sig), function(x) exp(ld(x)), 0) * 2 / sig
    expect_equal(dens, 2 * stats::dcauchy(sig, 0, s), tolerance = 1e-10)
    # P(sigma < s) = 1/2 for the half-Cauchy: P(theta > -2 log s)
    expect_equal(integrate(dth, th0, th0 + 80, rel.tol = 1e-10)$value, 0.5, tolerance = 1e-5)
  }
})

test_that("priorHyper(): Wishart parameter vector = df, diag(R), lower triangle of R (R = V^-1)", {
  # k = 2, non-diagonal
  V <- matrix(c(2, 0.5, 0.5, 1), 2)
  R <- solve(V)
  expect_equal(priorHyper(wishartPrec(df = 6, scale = V), 2L),
               list(theta1 = list(prior = "wishart2d", param = c(6, R[1, 1], R[2, 2], R[1, 2]))))
  # k = 2, scalar scale s: R = I / s, df = k + 2
  expect_equal(priorHyper(wishartPrec(scale = 4), 2L)$theta1$param, c(4, 0.25, 0.25, 0))
  # k = 3, non-diagonal, order R12, R13, R23
  V3 <- matrix(c(1, 0.2, 0.1, 0.2, 2, 0.3, 0.1, 0.3, 3), 3)
  R3 <- solve(V3)
  expect_equal(priorHyper(wishartPrec(scale = V3), 3L)$theta1,
               list(prior = "wishart3d", param = c(5, diag(R3), R3[1, 2], R3[1, 3], R3[2, 3])))
  # k = 4 and 5: wishart<k>d; scalar scales expand to s I
  expect_equal(priorHyper(wishartPrec(scale = 2), 4L)$theta1$prior, "wishart4d")
  h5 <- priorHyper(wishartPrec(scale = 2), 5L)$theta1
  expect_equal(h5$prior, "wishart5d")
  expect_equal(h5$param, c(7, rep(0.5, 5), rep(0, 10)))
  # k = 6, non-diagonal: wishartkd, 1 + 6 + 15 parameters in the order R12, ..., R16, R23, ...
  set.seed(3)
  A <- matrix(rnorm(36), 6)
  V6 <- crossprod(A) + diag(6)
  R6 <- solve(V6)
  h6 <- priorHyper(wishartPrec(df = 9, scale = V6), 6L)$theta1
  expect_equal(h6$prior, "wishartkd")
  upper <- unlist(lapply(1:5, function(i) R6[i, (i + 1):6]))
  expect_equal(h6$param, c(9, diag(R6), upper))
  expect_length(h6$param, 22)
  expect_equal(priorHyper(wishartPrec(scale = 1), 6L)$theta1$param[1], 8)   # df = k + 2
  # unresolved data scale
  expect_error(priorHyper(dataPrior(), 2L), "must be resolved with resolveBlockPrior\\(\\) first")
})

test_that("wishartArgs(): df default k + 2, scale as a k x k matrix", {
  expect_equal(wishartArgs(wishartPrec(), 3L), list(df = 5, scale = diag(3)))
  expect_equal(wishartArgs(wishartPrec(df = 7, scale = 0.5), 2L), list(df = 7, scale = diag(0.5, 2)))
  V <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(wishartArgs(wishartPrec(scale = V), 2L)$scale, V)
  expect_equal(wishartArgs(dataPrior(), 2L), list(df = 4, scale = NULL))
  # block size unknown (format(), priorNote()): df / scale only when given
  expect_equal(wishartArgs(wishartPrec(scale = 0.5)), list(df = NULL, scale = NULL))
  expect_equal(wishartArgs(wishartPrec(df = 7, scale = V)), list(df = 7, scale = V))
})

test_that("priorParams(): the derived numbers of every numeric prior", {
  expect_equal(priorParams(gammaPrec(2, 0.1)), list(shape = 2, rate = 0.1))
  expect_equal(priorParams(pcSd(0.5, 0.01)), list(u = 0.5, alpha = 0.01, rate = -log(0.01) / 0.5))
  expect_equal(priorParams(halfNormalSd(2)), list(scale = 2, prec = 0.25))
  expect_equal(priorParams(halfCauchySd(3)), list(scale = 3))
  expect_equal(priorParams(normalCoef(1, 2)), list(mean = 1, sd = 2, prec = 0.25))
  expect_equal(priorParams(wishartPrec(), 2L), list(df = 4, scale = diag(2)))
})

test_that("familyLatex(): the tag placeholder does not hit a capital T of an inlaFamily() hyper name", {
  skip_if_not_installed("INLA")
  f <- inlaFamily("gaussian", hyper = list(Theta = gammaPrec(1, 1)))
  expect_equal(hyperSymbol(f, "Theta", "(y)"), "\\theta^{(y)}_{\\mathrm{Theta}}")
  expect_match(familyLatex(f, "y_{i}", "(y)", "y", list())$dist[2L], "\\theta^{(y)}_{\\mathrm{Theta}} &", fixed = TRUE)
})

test_that("coefHyper(): normal prior on beta with precision 1 / sd^2", {
  expect_equal(coefHyper(normalCoef(1, 0.5)),
               list(beta = list(fixed = FALSE, prior = "normal", param = c(1, 4))))
  expect_equal(coefHyper(defaultCoefPrior()),
               list(beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.01))))
  expect_equal(coefHyper(normalCoef(-2, 100))$beta$param, c(-2, 1e-4))
})

# ---------------------------------------------------------------------------
# checkBlockPrior()

test_that("checkBlockPrior(): accepted combinations", {
  tm <- quote(1 + t | id)
  expect_null(checkBlockPrior(NULL, 2L, TRUE, tm))
  for (pr in list(gammaPrec(), pcSd(1), halfNormalSd(1), halfCauchySd(1))) {
    expect_invisible(checkBlockPrior(pr, 1L, TRUE, quote(1 | id)))          # single effect
    expect_invisible(checkBlockPrior(pr, 3L, FALSE, quote(1 + t + s || id))) # uncorrelated block
  }
  expect_invisible(checkBlockPrior(wishartPrec(), 2L, TRUE, tm))             # df = k + 2 = 4 > 3
  expect_invisible(checkBlockPrior(wishartPrec(df = 3.01), 2L, TRUE, tm))
  expect_invisible(checkBlockPrior(wishartPrec(scale = diag(2)), 2L, TRUE, tm))
  expect_invisible(checkBlockPrior(wishartPrec(scale = 7), 6L, TRUE, tm))
  expect_invisible(checkBlockPrior(dataPrior(), 2L, TRUE, tm))
})

test_that("checkBlockPrior(): rejected combinations", {
  tm <- quote(1 + t | id)
  expect_error(checkBlockPrior(normalCoef(), 2L, TRUE, tm),
               "normalCoef\\(\\) is a prior for the coefficient of scaled\\(\\); block `1 \\+ t \\| id`")
  expect_error(checkBlockPrior(normalCoef(), 1L, FALSE, quote(1 | id)), "normalCoef\\(\\) is a prior")
  for (pr in list(gammaPrec(), pcSd(1), halfNormalSd(1), halfCauchySd(1)))
    expect_error(checkBlockPrior(pr, 2L, TRUE, tm), "Correlated block `1 \\+ t \\| id` \\(k = 2\\) needs wishartPrec\\(\\)")
  expect_error(checkBlockPrior(wishartPrec(df = 3), 2L, TRUE, tm), "df must be > k \\+ 1 = 3")
  expect_error(checkBlockPrior(wishartPrec(df = 5), 4L, TRUE, tm), "df must be > k \\+ 1 = 5")
  expect_error(checkBlockPrior(wishartPrec(scale = diag(3)), 2L, TRUE, tm), "must be a 2 x 2 matrix")
  expect_error(checkBlockPrior(wishartPrec(scale = matrix(1)), 2L, TRUE, tm), "must be a 2 x 2 matrix")
  expect_error(checkBlockPrior(wishartPrec(), 1L, TRUE, quote(1 | id)),
               "wishartPrec\\(\\) needs a correlated block with >= 2 effects")
  expect_error(checkBlockPrior(wishartPrec(), 2L, FALSE, quote(1 + t || id)), "needs a correlated block")
  expect_error(checkBlockPrior(dataPrior(), 2L, FALSE, quote(1 + t || id)), "needs a correlated block")
})

# ---------------------------------------------------------------------------
# data-scaled Wishart prior

test_that("resolveBlockPrior(): NULL and non-data priors are returned unchanged", {
  b <- dataBlock(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), gaussianFamily())
  blk <- b$p$re[[1]]
  expect_null(resolveBlockPrior(NULL, b$p, b$df, eff2, blk))
  for (pr in list(pcSd(1), wishartPrec(scale = diag(2)), wishartPrec(df = 5, scale = 2)))
    expect_identical(resolveBlockPrior(pr, b$p, b$df, eff2, blk), pr)
})

test_that("resolveBlockPrior(): \"data\" -> V = S^-1 / df, so that E[Q] = df V = S^-1", {
  b <- dataBlock(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), gaussianFamily())
  blk <- b$p$re[[1]]
  S <- vcBlock(refFit("lmerX1"))
  r <- resolveBlockPrior(dataPrior(), b$p, b$df, eff2, blk)
  expect_s3_class(r, "mvremPrior")
  expect_null(r$args$df)                                 # df stays default (k + 2 = 4)
  expect_equal(r$args$scale, solve(S) / 4, tolerance = 1e-6)
  expect_true(isSymmetric(r$args$scale))
  expect_null(dimnames(r$args$scale))
  expect_equal(r$args$source, "data")
  # the INLA parameters: R = V^-1 = df S
  h <- priorHyper(r, 2L)$theta1
  expect_equal(h$param, c(4, 4 * S[1, 1], 4 * S[2, 2], 4 * S[1, 2]), tolerance = 1e-6)
  # a user df is kept and scales V
  r6 <- resolveBlockPrior(wishartPrec(df = 6, scale = "data"), b$p, b$df, eff2, blk)
  expect_equal(r6$args$df, 6)
  expect_equal(r6$args$scale, solve(S) / 6, tolerance = 1e-6)
  expect_invisible(checkBlockPrior(r6, 2L, TRUE, blk$term))
})

test_that("blockCovariance(): lmer for a gaussian response, whatever its link; no family = gaussian", {
  b <- dataBlock(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), gaussianFamily(link = "log"))
  ref <- refBlock(refFit("lmerX1"))
  expectBlockCov(b, ref)
  pn <- b$p
  pn$family <- NULL
  expectBlockCov(list(p = pn, df = b$df), ref)
})

test_that("blockCovariance(): the matrix follows the order of `effects`", {
  b <- dataBlock(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), gaussianFamily())
  expectBlockCov(b, refBlock(refFit("lmerX1"), rev(eff2)), effects = rev(eff2))
})

test_that("blockCovariance(): lmer on the log of a lognormal response", {
  b <- dataBlock(pt ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), lognormalFamily())
  expectBlockCov(b, refBlock(refFit("lmerLogPt")))
})

test_that("blockCovariance(): glmer for binomial with trials() and the family's link", {
  b <- dataBlock(k | trials(n) ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")),
                 binomialFamily(link = "probit"))
  S <- expectBlockCov(b, refBlock(refFit("glmerBinomProbit")))
  expect_true(all(eigen(S, only.values = TRUE)$values > 0))
})

test_that("blockCovariance(): glmer for poisson with exposure() as an offset", {
  b <- dataBlock(cnt | exposure(pt) ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), poissonFamily())
  S <- expectBlockCov(b, refBlock(refFit("glmerPois")))
  expect_true(all(eigen(S, only.values = TRUE)$values > 0))
})

test_that("blockCovariance(): the block is matched among several blocks of the same group", {
  b <- dataBlock(x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")) + (0 + x3 | id), gaussianFamily())
  S <- expectBlockCov(b, refBlock(refFit("lmerX1X3")))
  expect_equal(dim(S), c(2L, 2L))
})

test_that("blockCovariance(): errors for other families, missing blocks and failed fits", {
  b <- dataBlock(pt ~ t + re(1 + t | id, prior = wishartPrec(scale = "data")), lognormalFamily())
  blk <- b$p$re[[1]]
  for (fam in list(gammaFamily(), betaFamily(), negBinomialFamily(), studentTFamily(),
                   betaBinomialFamily(), zeroInflatedPoissonFamily())) {
    p <- b$p
    p$family <- fam
    expect_error(blockCovariance(p, b$df, eff2, blk),
                 "needs a Gaussian, lognormal, binomial or Poisson response .* block `1 \\+ t \\| id` \\(formula of 'pt'\\)")
  }
  expect_error(blockCovariance(b$p, b$df, c("(Intercept)", "x3"), blk),
               "No random-effect block `1 \\+ t \\| id` in the lme4 fit")
  expect_error(blockCovariance(b$p, b$df, "(Intercept)", blk), "No random-effect block")
  expect_error(blockCovariance(b$p, b$df[setdiff(names(b$df), "t")], eff2, blk),
               "lme4 fit for the data-scaled prior of block `1 \\+ t \\| id` \\(formula of 'pt'\\) failed: ")
})

test_that("regularizeCovariance(): positive-definite input unchanged, singular input regularised", {
  S <- matrix(c(0.25, 0.05, 0.05, 0.04), 2)
  expect_silent(res <- regularizeCovariance(S, "X"))
  expect_equal(res, S)
  S3 <- matrix(c(1, 0.3, 0.1, 0.3, 2, 0.2, 0.1, 0.2, 0.5), 3)
  expect_silent(expect_equal(regularizeCovariance(S3, "X"), S3))
  msg <- "The lme4 fit for X is singular; its covariance matrix was regularised"
  # correlation +1 / -1 -> +-0.95
  expect_warning(res <- regularizeCovariance(matrix(c(4, 2, 2, 1), 2), "X"), msg)
  expect_equal(stats::cov2cor(res)[1, 2], 0.95)
  expect_equal(diag(res), c(4, 1))
  expect_warning(res <- regularizeCovariance(matrix(c(1, -1, -1, 1), 2), "X"), msg)
  expect_equal(stats::cov2cor(res)[1, 2], -0.95)
  # a zero variance -> SD 1% of the largest (2 * 0.01 = 0.02), correlation 0
  expect_warning(res <- regularizeCovariance(diag(c(4, 0)), "X"), msg)
  expect_equal(res, diag(c(4, 0.02^2)))
  # a slightly negative variance is treated as zero
  expect_warning(res <- regularizeCovariance(diag(c(1, -1e-12)), "X"), msg)
  expect_equal(res, diag(c(1, 1e-4)))
  # k = 3 with clamped correlations still singular -> ridge, positive definite
  expect_warning(res <- regularizeCovariance(matrix(1, 3, 3), "X"), msg)
  ev <- eigen(res, only.values = TRUE)$values
  expect_true(min(ev) > 1e-8 * max(ev))
  expect_true(isSymmetric(res))
  expect_error(regularizeCovariance(matrix(0, 2, 2), "the fit Y"),
               "All random-effect variances of the lme4 fit for the fit Y are zero")
})

# ---------------------------------------------------------------------------
# priors in the INLA formula (the "priors" fixture, built without INLA)

test_that("buildInlaFormula(): every prior kind reaches the hyper of its block; copies carry none", {
  inl <- fxModel("priors", fit = FALSE)$inla
  expect_equal(hyperOf(inl$random$yg[[1]]$terms[[1]]), list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
  x1 <- inl$random$x1[[1]]
  expect_equal(hyperOf(x1$terms[[1]])$theta1$prior, "wishart2d")
  expect_null(hyperOf(x1$terms[[2]]))                          # the copy of the slope
  expect_match(hyperOf(inl$random$x2[[1]]$terms[[1]])$prec$prior, "^expression: sigma = exp\\(-theta/2\\); s = 1;")
  # logtnormal precision = 1 / 0.5^2 = 4
  expect_equal(hyperOf(inl$random$x2[[2]]$terms[[1]]), list(prec = list(prior = "logtnormal", param = c(0, 4))))
  gam <- list(prec = list(prior = "loggamma", param = c(1, 0.01)))
  expect_equal(hyperOf(inl$random$x3[[1]]$terms[[1]]), gam)    # `||`: the same prior on each effect
  expect_equal(hyperOf(inl$random$x3[[2]]$terms[[1]]), gam)
  expect_equal(inl$random$x3[[2]]$prior, gammaPrec(1, 0.01))
})

test_that("buildInlaFormula(): scaled() coefficient priors, given and default", {
  inl <- fxModel("priors", fit = FALSE)$inla
  # normalCoef(0, 2): precision 1 / 2^2 = 0.25; default normalCoef(0, 10): 0.01
  expect_equal(hyperOf(inl$scaled$yg_x1$term), list(beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.25))))
  expect_equal(hyperOf(inl$scaled$yg_x2$term), list(beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.01))))
  expect_equal(hyperOf(inl$scaled$yg_x3$term), hyperOf(inl$scaled$yg_x2$term))
  expect_equal(inl$scaled$yg_x1$prior, normalCoef(0, 2))
  expect_null(inl$scaled$yg_x2$prior)
})

test_that("buildInlaFormula(): wishartPrec(scale = \"data\") is resolved from the data and stored", {
  inl <- fxModel("priors", fit = FALSE)$inla
  S <- vcBlock(refFit("lmerX1"))
  pr <- inl$random$x1[[1]]$prior
  expect_equal(pr$args$source, "data")
  expect_equal(pr$args$scale, solve(S) / 4, tolerance = 1e-6)                  # V = S^-1 / df, df = k + 2
  expect_match(format(pr), "(from the data)", fixed = TRUE)
  h <- hyperOf(inl$random$x1[[1]]$terms[[1]])$theta1
  expect_equal(h$param, c(4, 4 * S[1, 1], 4 * S[2, 2], 4 * S[1, 2]), tolerance = 1e-6)   # R = df S
})

test_that("buildInlaFormula(): INLA reads the embedded hyper lists of every term", {
  skip_if_not_installed("INLA")
  inl <- fxModel("priors", fit = FALSE)$inla
  spec <- INLA:::inla.interpret.formula(inl$formula, data.same.len = FALSE, data = NULL)$random.spec
  names(spec) <- vapply(spec, `[[`, "", "label")
  plain <- function(h) list(prior = as.vector(h$prior), param = as.vector(h$param))
  expect_equal(spec$x1_1_id$model, "iid2d")
  expect_equal(plain(spec$x1_1_id$hyper$theta1),
               plain(hyperOf(inl$random$x1[[1]]$terms[[1]])$theta1))
  expect_equal(spec$x1_t_id$model, "copy")
  expect_equal(plain(spec$yg_1_site$hyper$theta), list(prior = "pc.prec", param = c(1, 0.01)))
  expect_match(spec$x2_1_id$hyper$theta$prior, "^expression:")
  expect_equal(plain(spec$x3_t_id$hyper$theta), list(prior = "loggamma", param = c(1, 0.01)))
  expect_equal(plain(spec$yg_FROM_x1_id_t$hyper$theta), list(prior = "normal", param = c(0, 0.25)))
  expect_equal(plain(spec$yg_FROM_x2_id_t$hyper$theta), list(prior = "normal", param = c(0, 0.01)))
  # the zero-equation field: fixed diffuse precision (zeroFieldHyper)
  expect_true(isTRUE(spec$x1_TO_yg_id_t$hyper$theta$fixed))
  expect_equal(as.vector(spec$x1_TO_yg_id_t$hyper$theta$initial), -6)
})

test_that("buildInlaFormula(): block priors are checked after factor expansion", {
  d <- fxSmallData()                                           # g has 3 levels: `1 + g` has k = 3 effects
  f <- function(prior) collectFormulas(list(y ~ t + scaled(x)),
                                       list(x ~ t + re(1 + g | id, prior = prior)))
  expect_error(buildInlaFormula(f(wishartPrec(scale = diag(2))), prepareData(f(wishartPrec(scale = diag(2))), d)),
               "scale of block `1 \\+ g \\| id` must be a 3 x 3 matrix")
  expect_error(buildInlaFormula(f(pcSd(1)), prepareData(f(pcSd(1)), d)),
               "Correlated block `1 \\+ g \\| id` \\(k = 3\\) needs wishartPrec\\(\\)")
  expect_silent(inl <- buildInlaFormula(f(wishartPrec(scale = diag(3))), prepareData(f(wishartPrec(scale = diag(3))), d)))
  expect_equal(hyperOf(inl$random$x[[1]]$terms[[1]])$theta1$prior, "wishart3d")
})

test_that("buildInlaFormula(): a data-scaled prior needs a family lme4 can fit", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ t + scaled(pt)), list(pt ~ t + re(1 + t | id, prior = wishartPrec(scale = "data"))),
                       family = list(pt = gammaFamily()))
  expect_error(buildInlaFormula(f, prepareData(f, d)),
               "needs a Gaussian, lognormal, binomial or Poisson response .* \\(formula of 'pt'\\)")
  f2 <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + re(1 + t || id, prior = wishartPrec(scale = "data"))))
  expect_error(buildInlaFormula(f2, prepareData(f2, d)), "needs a correlated block")
})

test_that("buildInlaFormula(): a singular lme4 fit warns and centres the prior on the regularised covariance", {
  skip_if_not(lme4::isSingular(refFit("lmerSmallX")), "the lme4 fit of x on fxSmallData() is not singular")
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + re(1 + t | id, prior = wishartPrec(scale = "data"))))
  expect_warning(inl <- buildInlaFormula(f, prepareData(f, d)), "singular; its covariance matrix was regularised")
  V <- inl$random$x[[1]]$prior$args$scale
  expect_true(all(is.finite(V)))
  expect_true(all(eigen(V, only.values = TRUE)$values > 0))
  S <- solve(4 * V)                                            # the regularised covariance (df = 4)
  expect_equal(S, refBlock(refFit("lmerSmallX"))$cov, tolerance = 1e-4)
  expect_equal(abs(stats::cov2cor(S)[1, 2]), 0.95)             # |correlation| clamped to 0.95
})

# ---------------------------------------------------------------------------
# the fitted "priors" fixture

test_that("jsm(): the \"priors\" fixture builds and fits without warnings", {
  # no hyper-override warning (no user `hyper` next to a prior) and no regularisation
  # warning (the lme4 fit of x1 on fxData() is not singular)
  expect_identical(fxWarnings("priors", fit = FALSE), character())
  expect_identical(fxWarnings("priors"), character())
})

test_that("hyperparameters(): family rows of the \"priors\" fixture are labelled by INLA's short names", {
  h <- hyperparameters(fxModel("priors"))
  fam <- h[h$type == "family", c("response", "parameter")]
  rownames(fam) <- NULL
  expect_equal(fam, data.frame(response = c("yg", "x1", "x2", "x3", "x3"),
                               parameter = c("residual precision", "residual precision", "size",
                                             "precision", "dof"),
                               stringsAsFactors = FALSE))
  expect_true(all(is.finite(h$mean)))
})

test_that("inlaFamily(): without INLA the name and link are accepted with a warning", {
  local_mocked_bindings(hasInla = function() FALSE)
  expect_warning(f <- inlaFamily("whatever", link = "anylink"), "INLA is not installed")
  expect_equal(f$link, "anylink")
  expect_null(familyLinks("gaussian"))
  expect_silent(gaussianFamily(link = "anylink"))   # links are validated only with INLA
})

