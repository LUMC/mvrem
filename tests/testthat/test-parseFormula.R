test_that("collectFormulas parses the reference model specification", {
  myHyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  res <- collectFormulas(
    formOut = list(
      z ~ 1 + t + (1 + t | id),
      y ~ 1 + t + re(1 + t | id, hyper = myHyper) +
        scaled(x1, prior = "normal") + scaled(x2)
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
  expect_equal(res$long$x3$re[[1]]$group, "id")

  # scaled terms (outcome only)
  expect_equal(vapply(res$out$y$scaled, `[[`, "", "target"), c("x1", "x2"))
  expect_equal(res$out$y$scaled[[1]]$params, list(prior = "normal"))
  expect_equal(res$out$z$scaled, list())
  expect_null(res$long$x1$scaled)

  # original formula kept
  expect_match(deparse1(res$out$y$formula), "scaled(x1", fixed = TRUE)
})

test_that("lme4 formula is usable by lme4", {
  skip_if_not_installed("lme4")
  p <- parseJsmFormula(Reaction ~ Days + re(Days | Subject, foo = 2), part = "long")
  fit <- lme4::lmer(p$lme4Formula, data = lme4::sleepstudy)
  expect_s4_class(fit, "merMod")
})

test_that("minus terms are preserved", {
  p <- parseJsmFormula(y ~ t - 1 + re(1 | id), part = "long")
  expect_equal(deparse1(p$lme4Formula), "y ~ t - 1 + (1 | id)")
})

test_that("invalid specifications error", {
  expect_error(parseJsmFormula(x ~ t + scaled(z), part = "long"), "only allowed in outcome")
  expect_error(parseJsmFormula(y ~ re(t, a = 1)), "first argument of re")
  expect_error(parseJsmFormula(y ~ re(1 | id, 3)), "must be named")
  expect_error(parseJsmFormula(y ~ scaled(~x1)), "name of a longitudinal")
  expect_error(parseJsmFormula(y ~ log(re(1 | id))), "top-level")
  expect_error(parseJsmFormula(~ t), "two-sided")
  expect_error(parseJsmFormula(cbind(a, b) ~ t), "exactly one variable")
  expect_error(parseJsmFormula(log(2) ~ t), "exactly one variable")
  expect_error(collectFormulas(list(y ~ scaled(q)), list(x ~ t)), "unknown longitudinal")
  expect_error(collectFormulas(list(y ~ t), list(y ~ t)), "more than once")
  expect_error(collectFormulas(list(y ~ t, log(y) ~ t), list(x ~ t)), "more than once")
})

test_that("variables are collected by role", {
  p <- parseJsmFormula(log(y) ~ 1 + t + I(t^2) + poly(age, 2) + sex:t +
                         re(1 + log(t) | id, a = 1) + (0 + t | site:id) +
                         scaled(x1, k = 2),
                       part = "out")
  expect_equal(p$variables$response, "y")
  expect_equal(p$response, "y")
  expect_equal(p$lhs, "log(y)")
  expect_setequal(p$variables$covariate, c("t", "age", "sex"))
  expect_setequal(p$variables$group, c("id", "site"))
  expect_false("x1" %in% unlist(p$variables))   # scaled target is not a data variable
  expect_false("a" %in% unlist(p$variables))    # re() parameter names are not variables

  res <- collectFormulas(list(y ~ t + scaled(x1)), list(x1 ~ t + I(t^2) + re(1 | id)))
  v <- res$variables
  expect_equal(v[v$formula == "x1" & v$role == "group", "variable"], "id")
  expect_equal(nrow(v), 5)
})

test_that("t^2 without a separate t term yields covariate t", {
  p <- parseJsmFormula(x1 ~ 1 + I(t^2) + re(1 + I(t^2) | id), part = "long")
  expect_equal(p$variables$covariate, "t")
  expect_equal(p$variables$group, "id")
  expect_equal(deparse1(p$lme4Formula), "x1 ~ 1 + I(t^2) + (1 + I(t^2) | id)")

  # plain t^2 (formula operator, not I()) also refers to t
  p2 <- parseJsmFormula(x1 ~ t^2 + re(1 | id), part = "long")
  expect_equal(p2$variables$covariate, "t")

  # data must therefore contain t
  forms <- collectFormulas(list(y ~ 1 + scaled(x1)), list(x1 ~ I(t^2) + re(1 | id)))
  d <- data.frame(id = 1:3, t = 1:3, x1 = 0, y = 0)
  expect_error(prepareData(d, d[, -2], forms), "not found in `dataLong`: t")
})

test_that("likelihood families are validated and default to gaussian", {
  f <- collectFormulas(list(y ~ t + scaled(x1), z ~ t), list(x1 ~ t, n ~ t),
                       familyOut = c(y = "beta"), familyLong = c(n = "poisson"))
  expect_equal(vapply(f$out, `[[`, "", "family"), c(y = "beta", z = "gaussian"))
  expect_equal(vapply(f$long, `[[`, "", "family"), c(x1 = "gaussian", n = "poisson"))
  expect_equal(collectFormulas(y ~ t, x ~ t)$out$y$family, "gaussian")

  skip_if_not_installed("INLA")
  expect_true(all(c("gaussian", "beta", "poisson", "binomial") %in% inlaFamilies()))
  expect_error(collectFormulas(y ~ t, x ~ t, familyOut = c(y = "betta")), "Unknown INLA likelihood")
  expect_error(collectFormulas(y ~ t, x ~ t, familyOut = c(x = "beta")), "belongs to the other part")
  expect_error(collectFormulas(y ~ t, x ~ t, familyLong = c(q = "beta")), "do not match any response")
  expect_error(collectFormulas(y ~ t, x ~ t, familyOut = "beta"), "non-empty names")
  expect_error(collectFormulas(y ~ t, x ~ t, familyOut = c(y = "beta", y = "gamma")), "unique")
  expect_error(collectFormulas(y ~ t, x ~ t, familyOut = list(y = "beta")), "character vector")
})
