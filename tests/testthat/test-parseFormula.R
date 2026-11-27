# Tests of R/parseFormula.R: parseJsmFormula(), collectFormulas() and their
# helpers (parseLhs, flattenTerms, parseRe, parseScaled, extraParams, buildFormula,
# collectVariables, resolveFamilies, checkFamilyList, checkAux, variableTable).
# Blocks are grouped by function, in the order of the source file's pipeline.

# ---- helpers ---------------------------------------------------------------

# Targets of the scaled() terms of a parsed outcome formula, in the order written.
scaledTargets <- function(p) vapply(p$scaled, `[[`, "", "target")

# collectFormulas() of the one-link model y ~ scaled(x1), x1 ~ t, with family arguments.
oneLink <- function(...) collectFormulas(list(y ~ scaled(x1)), list(x1 ~ t), ...)

# The message of the error raised by `expr` (fails if there is none).
errorMessage <- function(expr) conditionMessage(expect_error(expr))

# ---- parseJsmFormula(): input and left-hand side --------------------------

test_that("parseJsmFormula(): rejects non-formulae and one-sided formulae", {
  expect_error(parseJsmFormula("y ~ t"), "two-sided formula")
  expect_error(parseJsmFormula(list(y ~ t)), "two-sided formula")
  expect_error(parseJsmFormula(~ t), "two-sided formula")
})

test_that("parseJsmFormula(): `part` defaults to \"out\" and is validated", {
  expect_equal(parseJsmFormula(y ~ t)$part, "out")
  expect_equal(parseJsmFormula(y ~ t, "long")$part, "long")
  expect_error(parseJsmFormula(y ~ t, part = "middle"), "should be one of")
})

test_that("parseJsmFormula(): a plain left-hand side gives the name, no inputs, the formula as given", {
  p <- parseJsmFormula(y ~ t, "long")
  expect_equal(p$response, "y")
  expect_equal(p$lhs, "y")
  expect_true(is.name(p$responseExpr))
  expect_equal(p$aux, list())
  expect_identical(p$variables$aux, character())
  expect_identical(p$formula, y ~ t)
})

test_that("parseJsmFormula(): a transformed left-hand side keeps the expression and names its variable", {
  p <- parseJsmFormula(I(y / 2) ~ 1, "long")
  expect_equal(p$response, "y")
  expect_equal(p$lhs, "I(y/2)")
  expect_true(is.call(p$responseExpr))
  expect_equal(p$variables$response, "y")
  expect_equal(deparse1(p$lme4Formula), "I(y/2) ~ 1")
  expect_equal(deparse1(p$fixedFormula), "I(y/2) ~ 1")

  p2 <- parseJsmFormula(log(y + 1) ~ t, "long")
  expect_equal(p2$response, "y")
  expect_equal(p2$lhs, "log(y + 1)")
})

test_that("parseJsmFormula(): a left-hand side with zero or two variables is an error", {
  expect_error(parseJsmFormula(1 ~ t), "exactly one variable \\(found '1'\\)")
  expect_error(parseJsmFormula(log(2) ~ t), "exactly one variable")
  expect_error(parseJsmFormula(y + z ~ t), "exactly one variable \\(found 'y \\+ z'\\)")
  expect_error(parseJsmFormula(cbind(y, n - y) ~ t), "exactly one variable")
  expect_error(parseJsmFormula(I(y / z) ~ t), "exactly one variable")
  # two bars: the first `|` part is taken as the response expression
  expect_error(parseJsmFormula(y | trials(n) | exposure(E) ~ t), "exactly one variable")
})

test_that("parseJsmFormula(): trials() / exposure() / event() are family inputs, not formula terms", {
  p <- parseJsmFormula(k | trials(n) ~ t + (1 | id), "out")
  expect_equal(p$response, "k")
  expect_equal(p$lhs, "k")
  expect_equal(p$aux, list(trials = quote(n)))
  expect_equal(p$variables$aux, "n")
  expect_equal(deparse1(p$lme4Formula), "k ~ t + (1 | id)")
  expect_equal(deparse1(p$fixedFormula), "k ~ t")

  p2 <- parseJsmFormula(cnt | exposure(pt * 2) ~ 1, "out")
  expect_equal(p2$aux, list(exposure = quote(pt * 2)))
  expect_equal(p2$variables$aux, "pt")

  p3 <- parseJsmFormula(time | event(status) ~ 1 + sex, "out")
  expect_equal(p3$response, "time")
  expect_equal(p3$aux, list(event = quote(status)))
  expect_equal(p3$variables$aux, "status")
})

test_that("parseJsmFormula(): several inputs and input expressions of several variables", {
  p <- parseJsmFormula(log(k) | trials(n) + exposure(E) ~ 1, "out")
  expect_named(p$aux, c("trials", "exposure"))
  expect_equal(p$variables$aux, c("n", "E"))
  expect_equal(p$lhs, "log(k)")
  # the order written is kept, also in a long formula
  p2 <- parseJsmFormula(c | exposure(2 * pt) + trials(m) ~ t, "long")
  expect_equal(p2$aux, list(exposure = quote(2 * pt), trials = quote(m)))
  # any expression inside, also of several variables
  p3 <- parseJsmFormula(k | trials(n1 + n2) ~ 1, "out")
  expect_equal(p3$variables$aux, c("n1", "n2"))
})

test_that("parseJsmFormula(): an input must have exactly one unnamed argument", {
  expect_error(parseJsmFormula(y | trials(n = m) ~ t), "exactly one unnamed argument")
  expect_error(parseJsmFormula(y | trials(n, m) ~ t), "exactly one unnamed argument")
  expect_error(parseJsmFormula(y | trials() ~ t), "exactly one unnamed argument")
  expect_no_error(parseJsmFormula(y | trials(n) ~ t))
})

test_that("parseJsmFormula(): a repeated input is an error", {
  expect_error(parseJsmFormula(y | trials(n) + trials(m) ~ t), "trials\\(\\) is given more than once")
  expect_error(parseJsmFormula(y | exposure(a) + exposure(b) ~ t), "exposure\\(\\) is given more than once")
})

test_that("parseJsmFormula(): other terms after `|` on the left-hand side are errors", {
  expect_error(parseJsmFormula(y | offset(w) ~ t), "Unsupported term.*offset\\(w\\)")   # weights() is an input since E13
  expect_error(parseJsmFormula(y | n ~ t), "Unsupported term.*: n")
  expect_error(parseJsmFormula(y | -trials(n) ~ t), "Unsupported term")
  expect_error(parseJsmFormula(y | trials(n) - exposure(E) ~ t), "Unsupported term.*exposure\\(E\\)")
})

# ---- parseJsmFormula(): fixed part ---------------------------------------

test_that("parseJsmFormula(): an intercept-only formula has no covariates, groups or blocks", {
  p <- parseJsmFormula(y ~ 1, "long")
  expect_equal(deparse1(p$fixedFormula), "y ~ 1")
  expect_length(p$variables$covariate, 0)
  expect_length(p$variables$group, 0)
  expect_length(p$re, 0)
})

test_that("parseJsmFormula(): `0 +`, `-1` and `- 1` are kept in both formulae", {
  p0 <- parseJsmFormula(y ~ 0 + t + re(1 | id), "long")
  expect_equal(deparse1(p0$lme4Formula), "y ~ 0 + t + (1 | id)")
  expect_equal(deparse1(p0$fixedFormula), "y ~ 0 + t")

  pm <- parseJsmFormula(y ~ -1 + re(1 | id), "long")
  expect_equal(deparse1(pm$lme4Formula), "y ~ -1 + (1 | id)")
  expect_equal(deparse1(pm$fixedFormula), "y ~ -1")

  pt <- parseJsmFormula(y ~ t - 1 + re(1 | id), "long")
  expect_equal(deparse1(pt$lme4Formula), "y ~ t - 1 + (1 | id)")
  expect_equal(deparse1(pt$fixedFormula), "y ~ t - 1")
})

test_that("parseJsmFormula(): without fixed terms the fixed formula falls back to the intercept", {
  # only random effects
  pr <- parseJsmFormula(y ~ re(1 | id), "long")
  expect_equal(deparse1(pr$lme4Formula), "y ~ (1 | id)")
  expect_equal(deparse1(pr$fixedFormula), "y ~ 1")
  # only scaled(): both formulae are intercept-only, `0 +` is kept
  ps <- parseJsmFormula(y ~ scaled(x1), "out")
  expect_equal(deparse1(ps$lme4Formula), "y ~ 1")
  expect_equal(deparse1(ps$fixedFormula), "y ~ 1")
  ps0 <- parseJsmFormula(y ~ 0 + scaled(x1), "out")
  expect_equal(deparse1(ps0$fixedFormula), "y ~ 0")
})

test_that("parseJsmFormula(): factor, interaction and function covariates give their variables", {
  p <- parseJsmFormula(y ~ factor(g) * t + a:b + ns(age, df = 3) + log(dose + 1), "long")
  expect_equal(p$variables$covariate, c("g", "t", "a", "b", "age", "dose"))
  expect_equal(deparse1(p$fixedFormula), "y ~ factor(g) * t + a:b + ns(age, df = 3) + log(dose + 1)")
  # the response may appear on the right (e.g. a lag): it is a covariate there
  p2 <- parseJsmFormula(y ~ lagY, "long")
  expect_equal(p2$variables$covariate, "lagY")
})

test_that("parseJsmFormula(): t^2 and I(t^2) without a separate t term yield covariate t", {
  p <- parseJsmFormula(x1 ~ 1 + I(t^2) + re(1 + I(t^2) | id), "long")
  expect_equal(p$variables$covariate, "t")
  expect_equal(p$variables$group, "id")
  expect_equal(deparse1(p$lme4Formula), "x1 ~ 1 + I(t^2) + (1 + I(t^2) | id)")
  # plain t^2 (formula operator, not I()) also refers to t
  p2 <- parseJsmFormula(x1 ~ t^2 + re(1 | id), "long")
  expect_equal(p2$variables$covariate, "t")
})

test_that("parseJsmFormula(): the formula environment is kept", {
  mk <- function() { local <- 1; y ~ t + re(1 | id) }
  f <- mk()
  p <- parseJsmFormula(f, "long")
  expect_identical(environment(p$lme4Formula), environment(f))
  expect_identical(environment(p$fixedFormula), environment(f))
})

# ---- parseJsmFormula(): random effects -----------------------------------

test_that("parseJsmFormula(): plain bar terms and re() give the same blocks", {
  a <- parseJsmFormula(y ~ t + (1 + t | id), "long")
  b <- parseJsmFormula(y ~ t + re(1 + t | id), "long")
  expect_equal(a$re, b$re)
  expect_equal(deparse1(a$lme4Formula), deparse1(b$lme4Formula))
  blk <- a$re[[1]]
  expect_equal(deparse1(blk$term), "1 + t | id")
  expect_equal(blk$group, "id")
  expect_true(blk$correlated)
  expect_equal(blk$params, list())
  expect_null(blk$prior)
})

test_that("parseJsmFormula(): `||` blocks are uncorrelated and `site:id` groups give both variables", {
  p <- parseJsmFormula(y ~ t + re(0 + t || site:id) + (1 | site), "long")
  expect_length(p$re, 2)
  expect_false(p$re[[1]]$correlated)
  expect_equal(p$re[[1]]$group, "site:id")
  expect_equal(p$re[[2]]$group, "site")
  expect_true(p$re[[2]]$correlated)
  expect_equal(p$variables$group, c("site", "id"))
  expect_equal(p$variables$covariate, "t")
  expect_equal(deparse1(p$lme4Formula), "y ~ t + (0 + t || site:id) + (1 | site)")
})

test_that("parseJsmFormula(): random-effect covariates are covariates even without a fixed term", {
  p <- parseJsmFormula(y ~ 1 + (0 + age | id), "long")
  expect_equal(p$variables$covariate, "age")
  p2 <- parseJsmFormula(y ~ 1 + re(1 + log(t) | id), "long")
  expect_equal(p2$variables$covariate, "t")
})

test_that("parseJsmFormula(): the lme4 formula is accepted by lme4::lmer()", {
  skip_if_not_installed("lme4")
  p <- parseJsmFormula(Reaction ~ Days + re(Days | Subject, foo = 2), part = "long")
  fit <- lme4::lmer(p$lme4Formula, data = lme4::sleepstudy)
  expect_s4_class(fit, "merMod")
})

test_that("parseJsmFormula(): re() keeps `prior` apart from the extra f() arguments", {
  p <- parseJsmFormula(y ~ re(1 + t | id, prior = wishartPrec(), constr = TRUE,
                              hyper = list(theta1 = list(initial = 1))), "long")
  blk <- p$re[[1]]
  expect_s3_class(blk$prior, "mvremPrior")
  expect_equal(blk$prior, wishartPrec())
  expect_named(blk$params, c("constr", "hyper"))        # prior is not an f() argument
  expect_equal(blk$params$hyper, list(theta1 = list(initial = 1)))
  expect_equal(deparse1(p$lme4Formula), "y ~ (1 + t | id)")
})

test_that("parseJsmFormula(): re() with only a prior has an empty (unnamed) params list", {
  p <- parseJsmFormula(y ~ re(1 | id, prior = pcSd(1)), "long")
  expect_identical(p$re[[1]]$params, list())
  expect_equal(p$re[[1]]$prior, pcSd(1))
})

test_that("parseJsmFormula(): re() parameter names are not data variables", {
  p <- parseJsmFormula(y ~ t + re(1 | id, a = 1, hyper = list(prec = list(initial = 2))), "long")
  expect_equal(p$variables$covariate, "t")
  expect_false(any(c("a", "hyper", "prec") %in% unlist(p$variables)))
})

test_that("parseJsmFormula(): the first argument of re() must be an unnamed bar term", {
  expect_error(parseJsmFormula(y ~ re()), "first argument of re\\(\\) must be a random-effects term")
  expect_error(parseJsmFormula(y ~ re(t, a = 1)), "first argument of re\\(\\) must be a random-effects term")
  expect_error(parseJsmFormula(y ~ re(prior = pcSd(1), 1 | id)), "first argument of re\\(\\) must be a random-effects term")
  expect_error(parseJsmFormula(y ~ re(b = 1 | id)), "first argument of re\\(\\) must not be named")
})

test_that("parseJsmFormula(): extra arguments of re() must be named", {
  expect_error(parseJsmFormula(y ~ re(1 | id, pcSd(1))), "extra arguments of re\\(\\) must be named")
  expect_error(parseJsmFormula(y ~ re(1 | id, a = 1, 2)), "must be named \\(formula for 'y'\\)")
})

test_that("parseJsmFormula(): `prior` of re() must be a block prior constructor", {
  expect_error(parseJsmFormula(y ~ re(1 | id, prior = normalCoef())), "`prior` in re\\(\\) must be created by a prior constructor")
  expect_error(parseJsmFormula(y ~ re(1 | id, prior = list(prec = list(prior = "pc.prec")))),
               "`prior` in re\\(\\)")
  expect_error(parseJsmFormula(y ~ re(1 | id, prior = "pc.prec")), "`prior` in re\\(\\)")
})

test_that("parseJsmFormula(): nested re() / scaled() are errors", {
  expect_error(parseJsmFormula(y ~ t:re(1 | id)), "must be top-level terms")
  expect_error(parseJsmFormula(y ~ log(re(1 | id))), "must be top-level terms")
  expect_error(parseJsmFormula(y ~ I(2 * scaled(x1))), "must be top-level terms")
  expect_error(parseJsmFormula(y ~ exp(scaled(x1)), "out"), "must be top-level terms \\(formula for 'y'\\)")
})

test_that("parseJsmFormula(): extra arguments are evaluated in the formula's environment", {
  mk <- function() {
    myPrior <- pcSd(2, 0.05)
    myCoef <- normalCoef(1, 2)
    myHyper <- list(prec = list(initial = 3))
    y ~ t + re(1 | id, prior = myPrior, hyper = myHyper) + scaled(x1, prior = myCoef, k = myHyper)
  }
  p <- parseJsmFormula(mk(), "out")
  expect_equal(p$re[[1]]$prior, pcSd(2, 0.05))
  expect_equal(p$re[[1]]$params, list(hyper = list(prec = list(initial = 3))))
  expect_equal(p$scaled[[1]]$prior, normalCoef(1, 2))
  expect_equal(p$scaled[[1]]$params, list(k = list(prec = list(initial = 3))))
  # also through collectFormulas()
  f <- collectFormulas(list(mk()), list(x1 ~ t + (1 | id)))
  expect_equal(f$out$y$re[[1]]$prior, pcSd(2, 0.05))
})

test_that("parseJsmFormula(): an undefined variable in an extra argument is an error at parse time", {
  expect_error(parseJsmFormula(y ~ re(1 | id, prior = notDefinedAnywhere123)), "notDefinedAnywhere123")
})

# ---- parseJsmFormula(): scaled() -----------------------------------------

test_that("parseJsmFormula(): scaled() terms are collected in order with prior and params", {
  p <- parseJsmFormula(y ~ scaled(x2) + t + scaled(x1, prior = normalCoef(0, 1), hyper = list(beta = list(initial = 0.5))), "out")
  expect_equal(scaledTargets(p), c("x2", "x1"))
  expect_null(p$scaled[[1]]$prior)
  expect_identical(p$scaled[[1]]$params, list())
  expect_equal(p$scaled[[2]]$prior, normalCoef(0, 1))
  expect_named(p$scaled[[2]]$params, "hyper")
})

test_that("parseJsmFormula(): scaled() is dropped from both formulae; targets and params are not variables", {
  p <- parseJsmFormula(y ~ t + scaled(x1, k = 2), "out")
  expect_equal(deparse1(p$lme4Formula), "y ~ t")
  expect_equal(deparse1(p$fixedFormula), "y ~ t")
  expect_equal(p$variables$covariate, "t")
  expect_false(any(c("x1", "k") %in% unlist(p$variables)))
})

test_that("parseJsmFormula(): a long formula has no `scaled` element, an outcome without scaled() an empty list", {
  expect_null(parseJsmFormula(x ~ t, "long")$scaled)
  expect_identical(parseJsmFormula(y ~ t, "out")$scaled, list())
})

test_that("parseJsmFormula(): the first argument of scaled() must be an unnamed response name", {
  expect_error(parseJsmFormula(y ~ scaled()), "first argument of scaled\\(\\) must be the name")
  expect_error(parseJsmFormula(y ~ scaled("x1")), "first argument of scaled\\(\\) must be the name")
  expect_error(parseJsmFormula(y ~ scaled(log(x1))), "first argument of scaled\\(\\) must be the name")
  expect_error(parseJsmFormula(y ~ scaled(~x1)), "first argument of scaled\\(\\) must be the name")
  expect_error(parseJsmFormula(y ~ scaled(x = x1)), "first argument of scaled\\(\\) must not be named")
})

test_that("parseJsmFormula(): extra arguments of scaled() must be named, `prior` a normalCoef()", {
  expect_error(parseJsmFormula(y ~ scaled(x1, 2)), "extra arguments of scaled\\(\\) must be named")
  expect_error(parseJsmFormula(y ~ scaled(x1, prior = pcSd(1))), "`prior` in scaled\\(\\) must be normalCoef")
  expect_error(parseJsmFormula(y ~ scaled(x1, prior = c(0, 1))), "`prior` in scaled\\(\\) must be normalCoef")
  expect_error(parseJsmFormula(y ~ scaled(x1, prior = list(beta = list()))), "must be normalCoef")
})

test_that("parseJsmFormula(): scaled() in a long formula is an error", {
  expect_error(parseJsmFormula(x ~ t + scaled(x1), "long"), "only allowed in outcome and latent formulae \\(found in formula for 'x'\\)")
})

test_that("parseJsmFormula(): a covariate named like re / scaled is an ordinary covariate", {
  p <- parseJsmFormula(y ~ t + re, "long")
  expect_equal(p$variables$covariate, c("t", "re"))
  expect_length(p$re, 0)
  expect_equal(deparse1(p$fixedFormula), "y ~ t + re")
  p2 <- parseJsmFormula(y ~ log(scaled) + scaled(x1), "out")
  expect_equal(p2$variables$covariate, "scaled")
  expect_equal(p2$scaled[[1]]$target, "x1")
  # also as a random-effect covariate / inside I()
  p3 <- parseJsmFormula(y ~ I(re^2) + (re | id), "long")
  expect_equal(p3$variables$covariate, "re")
  expect_length(p3$re, 1)
  # nested CALLS are still rejected
  expect_error(parseJsmFormula(y ~ log(re(1 | id)), "long"), "must be top-level terms")
  expect_error(parseJsmFormula(y ~ (t + re(1 | id)), "long"), "must be top-level terms")
  expect_error(parseJsmFormula(y ~ t + I(scaled(x1)), "out"), "must be top-level terms")
})

test_that("parseJsmFormula(): a doubly parenthesised bar term is a random effect", {
  p <- parseJsmFormula(y ~ t + ((1 | id)), "long")
  expect_length(p$re, 1)
  expect_equal(p$variables$group, "id")
  expect_equal(p$variables$covariate, "t")
  expect_equal(p$re[[1]]$term, quote(1 | id))
  expect_equal(deparse1(p$lme4Formula), "y ~ t + (1 | id)")
  expect_equal(deparse1(p$fixedFormula), "y ~ t")
  # any number of brackets, also around `||` and inside re()
  p3 <- parseJsmFormula(y ~ (((1 + t || id))), "long")
  expect_length(p3$re, 1)
  expect_false(p3$re[[1]]$correlated)
  expect_equal(deparse1(p3$lme4Formula), "y ~ (1 + t || id)")
  p4 <- parseJsmFormula(y ~ re((1 | id), prior = pcSd(1)), "long")
  expect_equal(p4$re[[1]]$term, quote(1 | id))
  expect_equal(p4$re[[1]]$prior, pcSd(1))
  # the same block is found in a whole model
  f <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + ((1 | id))))
  expect_equal(f$variables$role[f$variables$variable == "id"], "group")
})

# ---- internal helpers ----------------------------------------------------

test_that("parseLhs(): splits the response expression from the inputs", {
  expect_equal(parseLhs(quote(y)), list(response = quote(y), aux = list()))
  r <- parseLhs(quote(log(k) | trials(n)))
  expect_equal(r$response, quote(log(k)))
  expect_equal(r$aux, list(trials = quote(n)))
  expect_equal(parseLhs(quote(t | event(s)))$aux, list(event = quote(s)))
})

test_that("flattenTerms(): splits additive terms and tracks signs", {
  tt <- flattenTerms(quote(a + b - c - (d + e) + -f))
  expect_equal(vapply(tt, function(x) deparse1(x$expr), ""), c("a", "b", "c", "(d + e)", "f"))
  expect_equal(vapply(tt, `[[`, 1L, "sign"), c(1L, 1L, -1L, -1L, -1L))
  expect_equal(vapply(flattenTerms(quote(+a)), `[[`, 1L, "sign"), 1L)
  expect_equal(vapply(flattenTerms(quote(-a)), `[[`, 1L, "sign"), -1L)
  expect_length(flattenTerms(quote(a:b)), 1)
})

test_that("buildFormula(): rebuilds signs and falls back to an intercept", {
  e <- globalenv()
  expect_equal(deparse1(buildFormula(quote(y), list(), e)), "y ~ 1")
  f <- buildFormula(quote(log(y)), flattenTerms(quote(-1 + t - u)), e)
  expect_equal(deparse1(f), "log(y) ~ -1 + t - u")
  expect_identical(environment(f), e)
})

test_that("collectVariables(): roles from the formula parts", {
  re <- list(makeReBlock(quote(1 + age | site:id), list()))
  v <- collectVariables(quote(log(y)), list(quote(t), quote(I(t^2))), re, list(trials = quote(n)))
  expect_equal(v, list(response = "y", covariate = c("t", "age"), group = c("site", "id"), aux = "n"))
  v0 <- collectVariables(quote(y), list(), list(), list())
  expect_identical(v0$aux, character())
  expect_length(v0$covariate, 0)
})

# ---- collectFormulas(): structure ------------------------------------------

test_that("collectFormulas(): single formulae are accepted; elements are named by the responses", {
  f <- collectFormulas(y ~ t + scaled(x1), x1 ~ t + (1 | id))
  expect_named(f, c("out", "long", "variables", "order"))
  expect_identical(f$order, c("y", "x1"))
  expect_named(f$out, "y")
  expect_named(f$long, "x1")
})

test_that("collectFormulas(): names of the input lists are ignored", {
  f <- collectFormulas(list(a = log(y) ~ t + scaled(x1)), list(b = x1 ~ t + (1 | id)))
  expect_named(f$out, "y")
  expect_named(f$long, "x1")
  expect_equal(f$out$y$lhs, "log(y)")
})

test_that("collectFormulas(): parses the reference model specification", {
  myHyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  res <- collectFormulas(
    formOut = list(
      z ~ 1 + t + (1 + t | id) + scaled(x1),
      y ~ 1 + t + re(1 + t | id, hyper = myHyper) +
        scaled(x1, hyper = list(beta = list(initial = 1))) + scaled(x2) + scaled(x3)
    ),
    formLong = list(
      x1 ~ 1 + t + re(0 + t | id, a = 1) + re(1 | id),
      x2 ~ 1 + t + re(1 + t | id),
      x3 ~ 1 + t + re(1 + t || id)
    )
  )
  expect_named(res$out, c("z", "y"))
  expect_named(res$long, c("x1", "x2", "x3"))
  # lme4 simplification
  expect_equal(deparse1(res$out$z$lme4Formula), "z ~ 1 + t + (1 + t | id)")
  expect_equal(deparse1(res$out$y$lme4Formula), "y ~ 1 + t + (1 + t | id)")
  expect_equal(deparse1(res$long$x1$lme4Formula), "x1 ~ 1 + t + (0 + t | id) + (1 | id)")
  expect_equal(deparse1(res$long$x3$lme4Formula), "x3 ~ 1 + t + (1 + t || id)")
  # re blocks
  expect_length(res$out$z$re, 1)
  expect_equal(res$out$z$re[[1]]$params, list())
  expect_equal(res$out$y$re[[1]]$params, list(hyper = myHyper))
  expect_length(res$long$x1$re, 2)
  expect_equal(res$long$x1$re[[1]]$params, list(a = 1))
  expect_true(res$long$x2$re[[1]]$correlated)
  expect_false(res$long$x3$re[[1]]$correlated)
  # the links of an outcome are its scaled() terms, in the order written
  expect_equal(scaledTargets(res$out$y), c("x1", "x2", "x3"))
  expect_equal(res$out$y$scaled[[1]]$params, list(hyper = list(beta = list(initial = 1))))
  expect_equal(res$out$y$scaled[[3]]$params, list())
  expect_equal(scaledTargets(res$out$z), "x1")
  expect_null(res$long$x1$scaled)
  # the original formula is kept as written
  expect_match(deparse1(res$out$y$formula), "re(1 + t | id, hyper = myHyper)", fixed = TRUE)
  expect_match(deparse1(res$out$y$formula), "scaled(x1", fixed = TRUE)
})

test_that("collectFormulas(): a response defined more than once (within or across parts) is an error", {
  expect_error(collectFormulas(list(y ~ scaled(x)), list(x ~ t, log(x) ~ t)), "more than once: x")
  expect_error(collectFormulas(list(y ~ scaled(x)), list(x ~ t, y ~ t)), "more than once: y")
  expect_error(collectFormulas(list(k | trials(n) ~ scaled(x), I(k / 2) ~ t), list(x ~ t)),
               "more than once: k")
  expect_error(collectFormulas(list(y ~ scaled(x), y ~ t, y ~ 1), list(x ~ t)), "more than once: y$")
})

# ---- collectFormulas(): associations ---------------------------------------

test_that("collectFormulas(): a scaled() target must be a longitudinal response", {
  expect_error(collectFormulas(y ~ scaled(x1) + scaled(q) + scaled(r), x1 ~ t),
               "for 'y' refers to unknown longitudinal response\\(s\\): q, r")
  # an outcome response is not a valid target
  expect_error(collectFormulas(list(y ~ scaled(x1), z ~ scaled(y)), x1 ~ t),
               "unknown longitudinal response\\(s\\): y")
})

test_that("collectFormulas(): at most one scaled() per target in a formula; the same target from two outcomes is fine", {
  expect_error(collectFormulas(y ~ scaled(x1) + t + scaled(x1, prior = normalCoef()), x1 ~ t),
               "given more than once for: x1")
  expect_error(collectFormulas(y ~ scaled(x1) + scaled(x1), x1 ~ t), "more than once for: x1")
  expect_no_error(collectFormulas(list(y ~ scaled(x1), z ~ scaled(x1)), x1 ~ t))
})

test_that("collectFormulas(): at least one scaled() term is required in the model", {
  expect_error(collectFormulas(y ~ t, x1 ~ t), "No scaled\\(\\) term")
  expect_error(collectFormulas(list(y ~ t, z ~ 1), list(x1 ~ t)), "No scaled\\(\\) term")
})

test_that("collectFormulas(): every longitudinal response must be linked; all unlinked ones are listed", {
  expect_error(collectFormulas(y ~ scaled(x1), list(x1 ~ t, x2 ~ t, x3 ~ t)),
               "not named in any scaled\\(\\) term: x2, x3")
})

test_that("collectFormulas(): an outcome without scaled() is allowed when another has one", {
  f <- collectFormulas(list(y ~ t + (1 | id), z ~ t + scaled(x1)), x1 ~ age + (1 | id))
  expect_identical(f$out$y$scaled, list())
  expect_identical(f$out$y$linkVariables, character())
  expect_equal(f$out$z$linkVariables, c("id", "age"))
})

test_that("collectFormulas(): outcomes are linked only to the targets of their scaled() terms, in order", {
  f <- collectFormulas(list(y ~ 1 + scaled(x2) + scaled(x1),
                            z ~ t + scaled(x2, hyper = list(beta = list(initial = 1))),
                            w ~ t + scaled(x3), v ~ t),
                       list(x1 ~ t + (1 | id), x2 ~ age, x3 ~ g))
  expect_equal(scaledTargets(f$out$y), c("x2", "x1"))
  expect_equal(scaledTargets(f$out$z), "x2")
  expect_equal(scaledTargets(f$out$w), "x3")
  expect_length(f$out$v$scaled, 0)
  expect_equal(f$out$z$scaled[[1]]$params, list(hyper = list(beta = list(initial = 1))))
  expect_equal(f$out$y$scaled[[1]]$params, list())
  # link variables come from the linked long formulae only
  expect_setequal(f$out$y$linkVariables, c("t", "id", "age"))
  expect_equal(f$out$z$linkVariables, "age")
  expect_equal(f$out$w$linkVariables, "g")
  expect_equal(f$out$v$linkVariables, character())
})

test_that("collectFormulas(): linkVariables = groups then covariates per target, no duplicates", {
  f <- collectFormulas(list(y ~ scaled(x2) + scaled(x1), z ~ age + (1 | id) + scaled(x2)),
                       list(x1 ~ t + (1 | id), x2 ~ age + (1 | site:id)))
  # targets in the order written (x2 then x1)
  expect_equal(f$out$y$linkVariables, c("site", "id", "age", "t"))
  # variables already in the outcome formula are dropped
  expect_equal(f$out$z$linkVariables, "site")
})

test_that("collectFormulas(): outcome inputs count as present; long inputs and responses are not link variables", {
  f <- collectFormulas(list(k | trials(t) ~ 1 + scaled(x)), list(x | exposure(pt) ~ t + (1 | id)),
                       family = list(k = "binomial", x = "poisson"))
  expect_equal(f$out$k$linkVariables, "id")
})

test_that("collectFormulas(): the outcome data must contain the link variables", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ 1 + scaled(x)), list(x ~ t + (1 | id)))
  expect_error(prepareData(f, dataOut = d[names(d) != "t"], dataLong = d), "not found in `dataOut`: t")
  expect_no_error(prepareData(f, d))
})

test_that("collectFormulas(): a covariate used only inside I() must be in the data", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ 1 + scaled(x)), list(x ~ I(t^2) + re(1 | id)))
  expect_error(prepareData(f, dataOut = d, dataLong = d[names(d) != "t"]), "not found in `dataLong`: t")
})

test_that("collectFormulas(): (1 || id) with a single intercept prepares data", {
  f <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + (1 || id)))
  expect_false(f$long$x$re[[1]]$correlated)
  expect_equal(f$long$x$re[[1]]$term, quote(1 || id))
  # lme4's double-bar expansion cannot parse an intercept-only `||`: it is passed as `|`
  expect_equal(deparse1(f$long$x$lme4Formula), "x ~ t + (1 | id)")
  d <- prepareData(f, fxSmallData())
  expect_identical(d$long$x, prepareData(collectFormulas(list(y ~ t + scaled(x)),
                                                         list(x ~ t + (1 | id))), fxSmallData())$long$x)
  # also via re(), with a composite group, and a block with slopes keeps `||`
  f2 <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + re(1 || id, prior = pcSd(1)) + (1 || site:id)))
  expect_equal(deparse1(f2$long$x$lme4Formula), "x ~ t + (1 | id) + (1 | site:id)")
  expect_no_error(prepareData(f2, fxSmallData()))
  f3 <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + (1 + t || id)))
  expect_equal(deparse1(f3$long$x$lme4Formula), "x ~ t + (1 + t || id)")
  expect_equal(deparse1(parseJsmFormula(y ~ (0 + t || id))$lme4Formula), "y ~ (0 + t || id)")
})

# ---- collectFormulas(): variables table ------------------------------------

test_that("collectFormulas(): variables table contents", {
  f <- collectFormulas(list(k | trials(n) ~ grp + scaled(x1)),
                       list(x1 ~ t + (1 + t | site:id)),
                       family = list(k = binomialFamily()))
  v <- f$variables
  expect_named(v, c("variable", "role", "part", "formula"))
  expect_equal(v$variable, c("k", "grp", "n", "x1", "t", "site", "id"))
  expect_equal(v$role, c("response", "covariate", "aux", "response", "covariate", "group", "group"))
  nOut <- 3; nLong <- 4                     # k, grp, n | x1, t, site, id
  expect_equal(v$part, rep(c("out", "long"), c(nOut, nLong)))
  expect_equal(v$formula, rep(c("k", "x1"), c(nOut, nLong)))
  expect_equal(rownames(v), as.character(seq_len(nrow(v))))
  # link variables are not in the table (they belong to the long formula)
  expect_equal(sum(v$variable == "t"), 1)
})

test_that("collectFormulas(): a long exposure() variable has role aux", {
  f <- collectFormulas(y ~ t + scaled(x), x | exposure(e) ~ t, familyLong = c(x = "nbinomial"))
  expect_equal(f$variables$role[f$variables$variable == "e"], "aux")
})

test_that("variableTable(): empty input and a variable repeated across formulae", {
  e <- variableTable(list())
  expect_equal(nrow(e), 0)
  expect_named(e, c("variable", "role", "part", "formula"))
  f <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t + (1 + t | id)))
  v <- f$variables
  # t appears once per formula
  expect_equal(v[v$variable == "t", "formula"], c("y", "x"))
  expect_true(all(vapply(v, is.character, NA)))
})

# ---- collectFormulas(): families -------------------------------------------

test_that("collectFormulas(): `family` serves both parts, unlisted responses are gaussian", {
  f <- collectFormulas(list(y ~ scaled(x1), z ~ 1), list(x1 ~ t),
                       family = list(y = betaFamily(), x1 = poissonFamily()))
  expect_equal(f$out$y$family, betaFamily())
  expect_equal(f$long$x1$family, poissonFamily())
  expect_equal(f$out$z$family, gaussianFamily())
  expect_s3_class(f$out$z$family, "mvremFamily")
})

test_that("collectFormulas(): familyOut / familyLong override `family` entry by entry", {
  f <- collectFormulas(list(y ~ scaled(x1), z ~ 1), list(x1 ~ t),
                       family = list(y = betaFamily(), x1 = "poisson", z = "gamma"),
                       familyOut = list(y = betaFamily(link = "probit")),
                       familyLong = c(x1 = "nbinomial"))
  expect_equal(f$out$y$family$link, "probit")
  expect_equal(f$out$z$family, gammaFamily())           # from `family`, not overridden
  expect_equal(f$long$x1$family, negBinomialFamily())
})

test_that("collectFormulas(): INLA names and family objects may be mixed; a name gives the constructor's defaults", {
  f <- collectFormulas(list(y ~ scaled(x1), z ~ 1), list(x1 ~ t),
                       family = c(y = "beta", x1 = "poisson"),
                       familyOut = list(z = betaFamily(link = "probit")))
  expect_equal(f$out$y$family, betaFamily())
  expect_equal(f$long$x1$family, poissonFamily())
  expect_equal(f$out$z$family$link, "probit")
  # element of a list given as a named length-1 character vector
  f2 <- oneLink(family = list(y = c(a = "beta")))
  expect_equal(f2$out$y$family, betaFamily())
})

test_that("collectFormulas(): empty / NULL family arguments give gaussian everywhere", {
  f <- oneLink(family = list(), familyOut = NULL, familyLong = list())
  expect_equal(f$out$y$family, gaussianFamily())
  expect_equal(f$long$x1$family, gaussianFamily())
})

test_that("collectFormulas(): family lists must hold family objects or single names", {
  expect_error(oneLink(family = c(y = NA_character_)), "`family` must be a named list of family objects")
  expect_error(oneLink(family = list(y = c("beta", "gamma"))), "`family` must be a named list of family objects")
  expect_error(oneLink(family = list(y = NA_character_)), "`family` must be a named list")
  expect_error(oneLink(familyOut = list(y = 1)), "`familyOut` must be a named list of family objects")
  expect_error(oneLink(family = gaussianFamily()), "`family` must be a named list of families")
})

test_that("collectFormulas(): family list names must be unique, non-empty responses of the part", {
  expect_error(oneLink(family = list(betaFamily())), "`family` must have unique, non-empty names")
  expect_error(oneLink(family = list(y = betaFamily(), betaFamily())), "non-empty names")
  expect_error(oneLink(familyOut = c(y = "beta", y = "gamma")), "unique, non-empty names")
  expect_error(oneLink(family = list(q = "beta")), "`family` names do not match any response: q")
  expect_error(oneLink(familyLong = list(y = "beta")), "`familyLong` names do not match.*: y \\(y belongs to the other part\\)")
  expect_error(oneLink(familyOut = list(x1 = "beta", q = "beta")),
               "`familyOut` names do not match any response of this part: x1, q \\(x1 belongs to the other part\\)")
  expect_error(oneLink(familyOut = list(q = "beta")), "of this part: q$")
})

test_that("collectFormulas(): INLA likelihood names are validated when INLA is installed", {
  skip_if_not_installed("INLA")
  expect_error(oneLink(familyOut = c(y = "betta")), "Unknown INLA likelihood")
  expect_error(oneLink(familyOut = c(y = "coxph")), "use weibullSurvFamily")
  expect_equal(oneLink(familyOut = c(y = "logistic"))$out$y$family$kind, "inla")
})
# Not tested: the no-INLA path (a warning "not validated" instead of the checks above).
# checkInlaFamily() calls requireNamespace("INLA") directly and the package has no
# `requireNamespace` binding of its own, so local_mocked_bindings() cannot reach it.

test_that("resolveFamilies(): defaults, shared entries and overrides", {
  r <- resolveFamilies(NULL, c("a", "b"), "familyOut", "x")
  expect_named(r, c("a", "b"))
  expect_true(all(vapply(r, inherits, NA, "mvremFamily")))
  r2 <- resolveFamilies(list(b = "beta"), c("a", "b"), "familyOut", "x",
                        shared = list(a = "gamma", b = "poisson", x = "poisson"))
  expect_equal(r2$a, gammaFamily())
  expect_equal(r2$b, betaFamily())
  expect_null(r2$x)                     # entries of the other part are not used here
  expect_identical(resolveFamilies(NULL, character(), "familyOut", "x"), stats::setNames(list(), character()))
})

test_that("checkFamilyList(): accepted forms and rejected ones", {
  expect_null(checkFamilyList(NULL, "family"))
  expect_equal(checkFamilyList(c(y = "beta"), "family"), list(y = "beta"))
  expect_equal(checkFamilyList(list(), "family"), list())
  expect_error(checkFamilyList(c("beta"), "familyOut"), "`familyOut` must have unique")
  expect_error(checkFamilyList(list(y = 1L), "familyOut"), "`familyOut` must be a named list")
})

# ---- checkAux() ------------------------------------------------------------

test_that("checkAux(): a left-hand side input must be used by the response's family", {
  expect_error(collectFormulas(y | exposure(E) ~ scaled(x), x ~ t),
               "exposure\\(\\) on the left-hand side of the formula for 'y' is not used by family 'gaussian'")
  expect_error(collectFormulas(y | exposure(E) ~ scaled(x), x ~ t, family = list(y = "binomial")),
               "not used by family 'binomial' \\(it takes trials\\(\\)\\)")
  expect_error(collectFormulas(y | trials(n) ~ scaled(x), x ~ t, family = list(y = poissonFamily())),
               "trials\\(\\) on the left-hand side.*'poisson' \\(it takes exposure\\(\\)\\)")
  expect_error(collectFormulas(y ~ scaled(x), x | trials(n) ~ t),
               "formula for 'x' is not used by family 'gaussian'")
  expect_error(collectFormulas(list(y | event(status) ~ 1 + scaled(x1)), list(x1 ~ t + (1 | id)),
                               family = list(y = betaFamily())),
               "event\\(\\) on the left-hand side.*not used by family 'beta'")
  expect_no_error(collectFormulas(y | trials(n) ~ scaled(x), x | exposure(E) ~ t,
                                  family = list(y = binomialFamily(), x = poissonFamily())))
})

test_that("checkAux(): a survival family requires event() on the left-hand side", {
  expect_error(collectFormulas(list(time ~ 1 + sex + scaled(x1)), list(x1 ~ t + (1 | id)),
                               family = list(time = exponentialSurvFamily())),
               "needs the event indicator.*time \\| event\\(status\\)")
  expect_error(collectFormulas(list(time ~ 1 + scaled(x1)), list(x1 ~ t + (1 | id)),
                               family = c(time = "weibullsurv")),
               "needs the event indicator")
})

test_that("checkAux(): direct calls on hand-made entries", {
  expect_invisible(checkAux(list(aux = list(), family = gaussianFamily(), response = "y", lhs = "y")))
  expect_error(checkAux(list(aux = list(), family = weibullSurvFamily(), response = "T", lhs = "T")),
               "needs the event indicator.*`T \\| event\\(status\\) ~ \\.\\.\\.`")
})

test_that("checkAux(): the list of accepted inputs is written once", {
  skip_if_not_installed("INLA")
  # inlaFamily() takes both trials() and exposure(); event() is not among them
  p <- list(aux = list(event = quote(s)), family = inlaFamily("logistic"), response = "y", lhs = "y")
  msg <- errorMessage(checkAux(p))
  expect_equal(lengths(regmatches(msg, gregexpr("it takes", msg, fixed = TRUE))), 1L)
  expect_match(msg, "(it takes trials(), exposure())", fixed = TRUE)
  # a family with one input: unchanged wording; none: no "(it takes" at all
  expect_match(errorMessage(checkAux(list(aux = list(exposure = quote(e)), family = binomialFamily(),
                                          response = "k", lhs = "k"))),
               "exposure() on the left-hand side of the formula for 'k' is not used by family 'binomial' (it takes trials()).",
               fixed = TRUE)
  expect_false(grepl("it takes", errorMessage(checkAux(list(aux = list(trials = quote(n)), family = gaussianFamily(),
                                                             response = "y", lhs = "y")))))
})

# ---- collectFormulas(): survival responses ---------------------------------

test_that("collectFormulas(): a survival outcome by family object or INLA name", {
  f <- collectFormulas(list(time | event(status) ~ 1 + sex + scaled(x1)), list(x1 ~ t + (1 | id)),
                       family = list(time = weibullSurvFamily()))
  expect_equal(f$out$time$family$family, "weibullsurv")
  expect_equal(f$variables$role[f$variables$variable == "status"], "aux")
  # the INLA name is the constructor with defaults
  f2 <- collectFormulas(list(time | event(status) ~ 1 + scaled(x1)), list(x1 ~ t + (1 | id)),
                        family = c(time = "weibullsurv"))
  expect_equal(f2$out$time$family, weibullSurvFamily())
  expect_equal(f2$out$time$family$aux, "event")
})

test_that("collectFormulas(): survival families are rejected for longitudinal responses", {
  expect_error(collectFormulas(list(x1 ~ 1 + scaled(T)), list(T | event(s) ~ 1 + t + (1 | id)),
                               family = list(T = weibullSurvFamily())),
               "longitudinal response\\(s\\) 'T'.*only for outcomes")
})

test_that("collectFormulas(): without INLA a family name is accepted with a warning", {
  local_mocked_bindings(hasInla = function() FALSE)
  expect_warning(f <- collectFormulas(list(y ~ t + scaled(x)), list(x ~ t), family = c(y = "whatever")),
                 "INLA is not installed.*'whatever'.*not validated")
  expect_equal(f$out$y$family$family, "whatever")
  expect_null(inlaFamilies())
})

test_that("parseJsmFormula(): offset() is rejected, as model.matrix() would drop it silently (B21)", {
  expect_error(parseJsmFormula(x1 ~ t + offset(log(e)) + (1 | id), "long"),
               "offset\\(\\) is not supported \\(formula for 'x1'\\).*`x1 \\| exposure\\(e\\) ~")
  expect_error(collectFormulas(list(y ~ t + offset(e) + scaled(x1), x1 ~ t + (1 | id))), "offset\\(\\)")
  # a variable named offset is an ordinary covariate (as for `re` / `scaled`, B2)
  p <- parseJsmFormula(x1 ~ t + offset + (1 | id), "long")
  expect_true("offset" %in% p$variables$covariate)
})
