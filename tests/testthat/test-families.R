test_that("constructors give the INLA likelihood and control.family", {
  f <- betaFamily()
  expect_s3_class(f, "mvremFamily")
  expect_equal(f$family, "beta")
  expect_equal(controlFamily(f), list())                 # INLA defaults: nothing set

  f <- betaFamily(phi = gammaPrec(1, 0.1), link = "probit")
  expect_equal(controlFamily(f), list(
    control.link = list(model = "probit"),
    hyper = list(phi = list(prior = "loggamma", param = c(1, 0.1)))))
  expect_equal(controlFamily(gaussianFamily(prec = pcSd(2)))$hyper$prec,
               list(prior = "pc.prec", param = c(2, 0.01)))
  expect_equal(controlFamily(negBinomialFamily(size = gammaPrec(2, 1)))$hyper$size,
               list(prior = "loggamma", param = c(2, 1)))
  # raw INLA entries are passed through; control is merged
  f <- betaBinomialFamily(rho = list(prior = "gaussian", param = c(0, 1)),
                          control = list(variant = 1))
  expect_equal(controlFamily(f), list(hyper = list(rho = list(prior = "gaussian", param = c(0, 1))),
                                      variant = 1))
  expect_equal(zeroInflatedPoissonFamily(0)$family, "zeroinflatedpoisson0")
  expect_equal(zeroInflatedNegBinomialFamily()$family, "zeroinflatednbinomial1")
  expect_equal(studentTFamily()$family, "t")
  expect_warning(controlFamily(gammaFamily(prec = gammaPrec(),
                                           control = list(hyper = list(prec = list(param = c(1, 1)))))),
                 "overrides the `prior` of family 'gamma' \\(prec\\)")
  expect_output(print(betaFamily(phi = gammaPrec(1, 0.1))), "<mvrem family> beta, link = logit")
})

test_that("invalid family arguments", {
  expect_error(betaBinomialFamily(rho = gammaPrec()), "does not accept this prior constructor")
  expect_error(negBinomialFamily(size = pcSd(1)), "use gammaPrec\\(\\)")
  expect_error(gaussianFamily(prec = list(a = 1)), "must be NULL, a prior constructor")
  expect_error(gaussianFamily(prec = 3), "must be NULL, a prior constructor")
  expect_error(zeroInflatedPoissonFamily(type = 2), "0 or 1")
  expect_error(betaFamily(control = list(1)), "named list")
  expect_error(inlaFamily("logistic", hyper = list(1)), "named list")
  skip_if_not_installed("INLA")
  expect_error(betaFamily(link = "identity"), "not available for INLA likelihood 'beta'")
  expect_error(inlaFamily("betta"), "Unknown INLA likelihood")
  expect_error(inlaFamily("coxph"), "not supported yet")
  expect_equal(controlFamily(inlaFamily("logistic", link = "identity",
                                        hyper = list(prec = gammaPrec(1, 1)))),
               list(control.link = list(model = "identity"),
                    hyper = list(prec = list(prior = "loggamma", param = c(1, 1)))))
})

test_that("trials() and exposure() on the left-hand side", {
  p <- parseJsmFormula(k | trials(n) ~ t + (1 | id), part = "out")
  expect_equal(p$response, "k")
  expect_equal(p$lhs, "k")
  expect_equal(p$aux, list(trials = quote(n)))
  expect_equal(p$variables$aux, "n")
  expect_equal(deparse1(p$lme4Formula), "k ~ t + (1 | id)")
  expect_equal(deparse1(p$fixedFormula), "k ~ t")

  p <- parseJsmFormula(c | exposure(2 * pt) + trials(m) ~ t, part = "long")
  expect_equal(p$aux, list(exposure = quote(2 * pt), trials = quote(m)))
  expect_setequal(p$variables$aux, c("pt", "m"))
  expect_equal(parseJsmFormula(y ~ t, "out")$variables$aux, character())

  expect_error(parseJsmFormula(y | weights(w) ~ t, "out"), "Unsupported term .* weights\\(w\\)")
  expect_error(parseJsmFormula(y | trials(n) + trials(m) ~ t, "out"), "more than once")
  expect_error(parseJsmFormula(y | trials(x = n) ~ t, "out"), "one unnamed argument")
  expect_error(parseJsmFormula(y | trials(n, m) ~ t, "out"), "one unnamed argument")

  expect_error(collectFormulas(y | trials(n) ~ t, x ~ t, familyOut = c(y = "poisson")),
               "trials\\(\\) .* not used by family 'poisson' \\(it takes exposure\\(\\)\\)")
  expect_error(collectFormulas(y | exposure(n) ~ t, x ~ t), "not used by family 'gaussian'")
  f <- collectFormulas(y ~ t, x | exposure(e) ~ t, familyLong = c(x = "nbinomial"))
  expect_equal(f$variables$role[f$variables$variable == "e"], "aux")
})

d <- data.frame(id = rep(1:4, each = 3), t = rep(0:2, 4), x = rnorm(12),
                k = c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5), n = 5,
                b = rep(0:1, 6), pt = seq(1, 2.1, by = 0.1), cnt = 0:11)

test_that("prepareData checks responses against the family", {
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x), b ~ t),
                       list(x ~ t + (1 | id), cnt | exposure(pt) ~ t),
                       familyOut = c(k = "binomial", b = "binomial"),
                       familyLong = c(cnt = "poisson"))
  dat <- prepareData(d, d, f)
  expect_equal(names(dat$out$k), c(".key", "k", "t", "n", "id"))
  # missing trials / exposure -> row dropped
  dNA <- d
  dNA$n[2] <- NA
  dNA$pt[3] <- NA
  datNA <- prepareData(dNA, dNA, f)
  expect_equal(nrow(datNA$out$k), 11)
  expect_equal(nrow(datNA$long$cnt), 11)

  bad <- function(col, value, family = NULL) {
    dd <- d
    dd[[col]][1] <- value
    prepareData(dd, dd, f)
  }
  expect_error(bad("k", 6), "Response 'k' \\(family 'binomial'\\) must be integers between 0 and trials\\(\\)")
  expect_error(bad("k", 1.5), "between 0 and trials")
  expect_error(bad("b", 2), "must be 0 or 1")
  expect_error(bad("n", 0), "trials\\(\\) of 'k' must be positive integers")
  expect_error(bad("pt", -1), "exposure\\(\\) of 'cnt' must be positive")
  expect_error(bad("cnt", -1), "non-negative integers \\(1 value\\(s\\) not\\)")

  g <- function(fam, y) {
    dd <- data.frame(y = y, z = 0, t = seq_along(y))
    prepareData(dd, dd, collectFormulas(y ~ t, z ~ t, familyOut = list(y = fam)))
  }
  expect_error(g(betaFamily(), c(0.5, 1)), "must be in \\(0, 1\\)")
  expect_error(g(gammaFamily(), c(1, 0)), "must be positive")
  expect_error(g(lognormalFamily(), c(1, -1)), "must be positive")
  expect_error(g(gaussianFamily(), c("a", "b")), "must be numeric")
  expect_silent(g(inlaFamily("logistic"), c(-1, 1)))
})

test_that("INLA formula: families, control.family, inputs and stacks", {
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x)),
                       list(x ~ t + (1 | id), cnt | exposure(pt) ~ t),
                       familyOut = list(k = binomialFamily(link = "probit")),
                       familyLong = list(x = gaussianFamily(prec = pcSd(1)),
                                         cnt = zeroInflatedPoissonFamily()))
  inl <- buildInlaFormula(f, prepareData(d, d, f))
  expect_equal(inl$families, c(k_k = "binomial", x_x = "gaussian",
                               cnt_cnt = "zeroinflatedpoisson1", x_x_0eq = "gaussian",
                               cnt_cnt_0eq = "gaussian"))
  expect_equal(names(inl$controlFamily), names(inl$families))
  expect_equal(inl$controlFamily$k_k, list(control.link = list(model = "probit")))
  expect_equal(inl$controlFamily$x_x$hyper$prec$prior, "pc.prec")
  expect_equal(inl$controlFamily$cnt_cnt, list())
  expect_equal(inl$controlFamily$x_x_0eq,
               list(hyper = list(prec = list(initial = 10, fixed = TRUE))))
  expect_equal(inl$inputs, list(k = c(Ntrials = "k_Ntrials"), x = stats::setNames(character(), character()),
                                cnt = c(E = "cnt_E")))
  expect_equal(inl$stacks$k$main$k_Ntrials, d$n)
  expect_equal(inl$stacks$cnt$main$cnt_E, d$pt)
  expect_equal(names(inl$stacks$k$main)[1:2], c("k_k", "k_Ntrials"))

  # longitudinal response with exposure: zero-equation rows have NA
  f2 <- collectFormulas(list(y ~ t + scaled(cnt)), list(cnt | exposure(pt) ~ t + (1 | id)),
                        familyLong = c(cnt = "poisson"))
  d2 <- transform(d, y = x)
  inl2 <- buildInlaFormula(f2, prepareData(d2, d2, f2))
  expect_equal(inl2$stacks$cnt$main$cnt_E, d$pt)
  expect_true(all(is.na(inl2$stacks$cnt$y$cnt_E)))
  expect_equal(names(inl2$stacks$cnt$y), names(inl2$stacks$cnt$main))

  # naming of the input columns
  inl3 <- buildInlaFormula(f2, prepareData(d2, d2, f2), naming = namingScheme(exposure = "E_{resp}"))
  expect_equal(inl3$inputs$cnt, c(E = "E_cnt"))
  expect_error(buildInlaFormula(f2, prepareData(d2, d2, f2), naming = namingScheme(exposure = "{resp}_1")),
               "duplicated INLA names: cnt_1")

  txt <- formatInlaFormula(inl, nrows = 1)
  expect_true(any(grepl('# family = "binomial", control.family = list(control.link = list(model = "probit"))',
                        txt, fixed = TRUE)))
  expect_true(any(grepl('# zero equations: family = "gaussian", control.family = list(hyper',
                        txt, fixed = TRUE)))
  expect_true(any(grepl('# family = "zeroinflatedpoisson1"$', txt)))
  expect_equal(deparse1(eval(parse(text = txt))), deparse1(inl$formula))
})

test_that("summary reports the INLA family", {
  f <- collectFormulas(list(k | trials(n) ~ t), list(x ~ t), familyOut = c(k = "binomial"))
  s <- summarizeData(f, d, d)
  expect_equal(s$equations$family, c("gaussian", "binomial"))
  expect_true("n" %in% s$missing$variable)
})

test_that("LaTeX of the families", {
  tex <- function(fam, lhs = y ~ t) {
    p <- parseJsmFormula(lhs, "out")            # standalone: no association terms
    p$family <- fam
    body <- sub("^\\\\begin\\{aligned\\}\n(.*)\n\\\\end\\{aligned\\}$", "\\1", unclass(latexEquation(p)))
    strsplit(body, " \\\\\n", fixed = TRUE)[[1]]
  }
  expect_equal(tex(binomialFamily(), y | trials(n) ~ t)[1:2], c(
    "\\operatorname{logit}\\left(p^{(y)}_{i}\\right) &= \\beta^{(y)}_{0} + \\beta^{(y)}_{1}\\, t_{i}",
    "y_{i} &\\sim \\mathrm{Binomial}\\left(n_{i}, p^{(y)}_{i}\\right)"))
  l <- tex(binomialFamily(link = "probit"))            # no trials(): Bernoulli
  expect_true(startsWith(l[1], "\\Phi^{-1}\\left(p^{(y)}_{i}\\right) &="))
  expect_equal(l[2], "y_{i} &\\sim \\mathrm{Bernoulli}\\left(p^{(y)}_{i}\\right)")
  expect_equal(tex(poissonFamily(), y | exposure(pt) ~ t)[2],
               "y_{i} &\\sim \\mathrm{Poisson}\\left(\\mathit{pt}_{i}\\, \\lambda^{(y)}_{i}\\right)")
  expect_equal(tex(negBinomialFamily(size = gammaPrec(1, 2)))[2:3], c(
    "y_{i} &\\sim \\mathrm{NegBin}\\left(\\mu^{(y)}_{i}, \\kappa^{(y)}\\right)",
    "\\kappa^{(y)} &\\sim \\mathrm{Gamma}\\left(1, 2\\right)"))
  expect_equal(tex(zeroInflatedPoissonFamily(1))[2], paste0(
    "y_{i} &\\sim \\pi^{(y)}\\,\\delta_{0} + \\left(1 - \\pi^{(y)}\\right) ",
    "\\mathrm{Poisson}\\left(\\lambda^{(y)}_{i}\\right)"))
  expect_match(tex(zeroInflatedNegBinomialFamily(0))[2], "\\mathrm{NegBin}^{+}\\left(", fixed = TRUE)
  expect_equal(tex(gammaFamily())[2],
               "y_{i} &\\sim \\mathrm{Gamma}\\left(\\phi^{(y)}, \\phi^{(y)} / \\mu^{(y)}_{i}\\right)")
  expect_equal(tex(betaBinomialFamily(), y | trials(n) ~ t)[2],
               "y_{i} &\\sim \\mathrm{BetaBinomial}\\left(n_{i}, p^{(y)}_{i}, \\rho^{(y)}\\right)")
  l <- tex(studentTFamily())
  expect_true(endsWith(l[1], " + \\varepsilon^{(y)}_{i}"))
  expect_equal(l[2], paste0("\\varepsilon^{(y)}_{i} &\\sim t_{\\nu^{(y)}}\\left(0, ",
                            "\\left(\\tau^{(y)}_{\\varepsilon}\\right)^{-1}\\right)"))
  expect_match(tex(lognormalFamily())[1], "^\\\\log\\\\left\\(y_\\{i\\}\\\\right\\) &= .*varepsilon")
  l <- tex(gaussianFamily(link = "log", prec = halfNormalSd(2)))
  expect_equal(l[1:2], c(
    "\\log\\left(\\mu^{(y)}_{i}\\right) &= \\beta^{(y)}_{0} + \\beta^{(y)}_{1}\\, t_{i}",
    "y_{i} &\\sim \\mathcal{N}\\left(\\mu^{(y)}_{i}, \\left(\\tau^{(y)}_{\\varepsilon}\\right)^{-1}\\right)"))
  expect_match(l[3], "^\\\\left\\(\\\\tau\\^\\{\\(y\\)\\}_\\{\\\\varepsilon\\}\\\\right\\)\\^\\{-1/2\\} &\\\\sim")
  skip_if_not_installed("INLA")
  expect_equal(tex(inlaFamily("logistic", hyper = list(prec = gammaPrec(1, 1))))[1:3], c(
    "g_{y}\\left(\\mathrm{E}\\left[y_{i}\\right]\\right) &= \\beta^{(y)}_{0} + \\beta^{(y)}_{1}\\, t_{i}",
    "y_{i} &\\sim \\mathrm{logistic}",
    "\\theta^{(y)}_{\\mathrm{prec}} &\\sim \\mathrm{Gamma}\\left(1, 1\\right)"))
})
