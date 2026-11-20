# E15: data-scaled scalar priors of single-effect and `||` blocks: pcSd(u = "data"),
# halfNormalSd(scale = "data"), halfCauchySd(scale = "data"), gammaPrec(rate = "data"),
# resolved per effect from an lme4 fit of the block's formula.

dataRules <- function() dataPriorRules

# x3 of fxData(): `||` block, gaussian
x3Block <- function(prior) dataBlock(stats::as.formula(bquote(x3 ~ t + re(1 + t || id, prior = .(prior)))),
                                     gaussianFamily())
lmerX3 <- function() {
  fit <- suppressMessages(lme4::lmer(x3 ~ t + (1 + t || id), fxData()$long))
  vc <- as.data.frame(lme4::VarCorr(fit))
  list(s = c(vc$sdcor[vc$var1 %in% "(Intercept)"], vc$sdcor[vc$var1 %in% "t"]), sigma = stats::sigma(fit))
}

# ---------------------------------------------------------------------------
# constructors, format, checks

test_that("scalar prior constructors accept \"data\" with the documented defaults", {
  expect_equal(pcSd(u = "data")$args, list(u = "data", alpha = dataRules()$pcAlpha))
  expect_equal(pcSd(u = "data", alpha = 0.1)$args$alpha, 0.1)
  expect_equal(pcSd(u = 1)$args$alpha, 0.01)                       # numeric u: unchanged default
  expect_equal(halfNormalSd("data")$args$scale, "data")
  expect_equal(halfCauchySd(scale = "data")$args$scale, "data")
  expect_equal(gammaPrec(rate = "data")$args, list(shape = 1, rate = "data"))
  expect_equal(gammaPrec(shape = 3, rate = "data")$args$shape, 3)
  expect_equal(gammaPrec()$args, list(shape = 1, rate = 5e-5))      # numeric rate: unchanged default
  for (bad in list("Data", "foo", c("data", "data"), NA_character_))
    expect_error(pcSd(u = bad), "`u` must be a single positive number or \"data\"")
  expect_error(halfNormalSd("x"), "`scale` must be a single positive number or \"data\"")
  expect_error(gammaPrec(shape = "data"), "`shape` must be a single positive number")
  expect_true(isScalarDataPrior(pcSd(u = "data")))
  expect_false(isScalarDataPrior(pcSd(u = 1)))
  expect_false(isScalarDataPrior(wishartPrec(scale = "data")))
  expect_false(isScalarDataPrior(NULL))
})

test_that("format() / print() of data-scaled priors say \"data\" and the rule", {
  expect_match(format(pcSd(u = "data")), "P(sigma > u) = 0.05 with u = data (6 x the lme4 SD of each effect)",
               fixed = TRUE)
  expect_match(format(halfNormalSd("data")), "sd = data (2 x the lme4 SD of each effect)", fixed = TRUE)
  expect_match(format(halfCauchySd("data")), "scale = data (2 x the lme4 SD of each effect)", fixed = TRUE)
  expect_match(format(gammaPrec(rate = "data")), "shape = 1, rate = data (shape x the lme4 variance", fixed = TRUE)
  expect_output(print(pcSd(u = "data")), "u = data")
  r <- dataScalarPrior(pcSd(u = "data"), 0.2)
  expect_match(format(r), "P(sigma > 1.2) = 0.05", fixed = TRUE)
  expect_match(format(r), "(from the data: lme4 SD 0.2)", fixed = TRUE)
})

test_that("checkBlockPrior() accepts data-scaled scalar priors for single and `||` blocks only", {
  for (pr in list(pcSd(u = "data"), halfNormalSd("data"), halfCauchySd("data"), gammaPrec(rate = "data"))) {
    expect_invisible(checkBlockPrior(pr, 1L, TRUE, quote(1 | id)))
    expect_invisible(checkBlockPrior(pr, 2L, FALSE, quote(1 + t || id)))
    expect_error(checkBlockPrior(pr, 2L, TRUE, quote(1 + t | id)), "needs wishartPrec")
  }
})

test_that("data-scaled priors are rejected for family hyperparameters and unresolved in priorHyper()", {
  expect_error(gaussianFamily(prec = pcSd(u = "data")), "data-scaled priors .* random-effect blocks only")
  expect_error(betaFamily(phi = gammaPrec(rate = "data")), "random-effect blocks only")
  expect_error(priorHyper(pcSd(u = "data"), 1L), "must be resolved")
  expect_error(priorHyper(gammaPrec(rate = "data"), 1L), "must be resolved")
})

# ---------------------------------------------------------------------------
# resolution

test_that("resolveEffectPriors(): `||` block, one prior per effect from the lme4 SD of that effect", {
  ref <- lmerX3()
  expect_true(all(ref$s > 0.1 * ref$sigma))                  # no floor in this fixture
  eff <- c("(Intercept)", "t")
  R <- dataRules()
  cases <- list(
    list(pcSd(u = "data"), function(s) list(u = R$pcU * s, alpha = R$pcAlpha)),
    list(halfNormalSd("data"), function(s) list(scale = R$halfNormalScale * s)),
    list(halfCauchySd("data"), function(s) list(scale = R$halfCauchyScale * s)),
    list(gammaPrec(rate = "data"), function(s) list(shape = 1, rate = 1 * s^2)))
  for (cs in cases) {
    b <- x3Block(cs[[1]])
    pr <- resolveEffectPriors(cs[[1]], b$p, b$df, eff, b$p$re[[1]])
    expect_length(pr, 2L)
    for (e in 1:2) {
      expect_s3_class(pr[[e]], "mvremPrior")
      expect_equal(pr[[e]]$dist, cs[[1]]$dist)
      expect_equal(pr[[e]]$args[names(cs[[2]](1))], cs[[2]](ref$s[e]), tolerance = 1e-4)
      expect_equal(pr[[e]]$args$source, "data")
      expect_equal(pr[[e]]$args$sdData, ref$s[e], tolerance = 1e-4)
      expect_false(isScalarDataPrior(pr[[e]]))
    }
  }
  # E[tau] of the gamma prior is 1 / s^2
  g <- resolveEffectPriors(gammaPrec(rate = "data"), b$p, b$df, eff, b$p$re[[1]])
  expect_equal(g[[2]]$args$shape / g[[2]]$args$rate, 1 / ref$s[2]^2, tolerance = 1e-4)
})

test_that("resolveEffectPriors(): non-data priors are repeated unchanged, without an lme4 fit", {
  b <- x3Block(pcSd(1))
  local_mocked_bindings(lme4BlockFit = function(...) stop("no lme4 fit expected"))
  expect_identical(resolveEffectPriors(pcSd(1), b$p, b$df, c("(Intercept)", "t"), b$p$re[[1]]),
                   list(pcSd(1), pcSd(1)))
  expect_identical(resolveEffectPriors(NULL, b$p, b$df, "t", b$p$re[[1]]), list(NULL))
})

test_that("resolveEffectPriors(): single-effect blocks, also two on one group", {
  f <- x1 ~ t + re(1 | id, prior = halfCauchySd("data")) + re(0 + t | id, prior = pcSd(u = "data"))
  b <- dataBlock(f, gaussianFamily())
  vc <- as.data.frame(lme4::VarCorr(suppressMessages(lme4::lmer(x1 ~ t + (1 | id) + (0 + t | id), fxData()$long))))
  s0 <- vc$sdcor[vc$var1 %in% "(Intercept)"]; s1 <- vc$sdcor[vc$var1 %in% "t"]
  p0 <- resolveEffectPriors(b$p$re[[1]]$prior, b$p, b$df, "(Intercept)", b$p$re[[1]])
  p1 <- resolveEffectPriors(b$p$re[[2]]$prior, b$p, b$df, "t", b$p$re[[2]])
  expect_equal(p0[[1]]$args$scale, dataRules()$halfCauchyScale * s0, tolerance = 1e-4)
  expect_equal(p1[[1]]$args$u, dataRules()$pcU * s1, tolerance = 1e-4)
})

test_that("resolveEffectPriors(): glmer for Poisson and binomial responses", {
  d <- fxData()$long
  b <- dataBlock(cnt | exposure(pt) ~ t + re(1 | id, prior = pcSd(u = "data")), poissonFamily())
  fit <- suppressWarnings(suppressMessages(lme4::glmer(cnt ~ t + (1 | id) + offset(log(pt)), d, family = poisson())))
  s <- max(sqrt(lme4::VarCorr(fit)$id[1]), dataRules()$floor)              # glmer: floor f x 1 on the link scale
  r <- suppressWarnings(resolveEffectPriors(pcSd(u = "data"), b$p, b$df, "(Intercept)", b$p$re[[1]]))
  expect_equal(r[[1]]$args$sdData, s, tolerance = 1e-3)
  b <- dataBlock(k | trials(n) ~ t + re(1 | id, prior = halfNormalSd("data")), binomialFamily())
  fit <- suppressWarnings(suppressMessages(lme4::glmer(cbind(k, n - k) ~ t + (1 | id), d, family = binomial())))
  s <- max(sqrt(lme4::VarCorr(fit)$id[1]), dataRules()$floor)
  r <- suppressWarnings(resolveEffectPriors(halfNormalSd("data"), b$p, b$df, "(Intercept)", b$p$re[[1]]))
  expect_equal(r[[1]]$args$scale, dataRules()$halfNormalScale * s, tolerance = 1e-3)
})

test_that("resolveEffectPriors(): a zero lme4 SD falls back to a fraction of the residual SD with a warning", {
  set.seed(11)
  d <- data.frame(id = rep(sprintf("S%02d", 1:30), each = 4), t = rep(0:3, 30))
  d$x <- 1 + 0.5 * d$t + rnorm(nrow(d), sd = 0.5)            # no subject effect at all
  d$y <- d$x + rnorm(nrow(d))
  fit <- suppressMessages(lme4::lmer(x ~ t + (1 + t || id), d))
  vc <- as.data.frame(lme4::VarCorr(fit))
  s <- c(vc$sdcor[vc$var1 %in% "(Intercept)"], vc$sdcor[vc$var1 %in% "t"])
  # floor = f x residual SD / rms of the effect's covariate (intercept: 1)
  floor <- dataRules()$floor * stats::sigma(fit) / c(1, sqrt(mean(d$t^2)))
  expect_true(all(s < floor))                                # precondition: singular fit
  b <- dataBlock(x ~ t + re(1 + t || id, prior = pcSd(u = "data")), gaussianFamily(), d = d)
  r <- withWarnings(resolveEffectPriors(pcSd(u = "data"), b$p, b$df, c("(Intercept)", "t"), b$p$re[[1]]))
  pct <- paste0(dataRules()$floor * 100, "% of the residual SD")
  expect_true(any(grepl(paste0("gives a small \\(or zero: singular fit\\) SD for `\\(Intercept\\)`, `t`.*floor .*", pct),
                        r$warnings)))
  expect_equal(vapply(r$value, function(x) x$args$sdData, 0), floor, tolerance = 1e-4)
  expect_equal(r$value[[2]]$args$u, dataRules()$pcU * floor[2], tolerance = 1e-4)
  # the same warning reaches buildInlaFormula()
  f <- collectFormulas(list(y ~ scaled(x)), list(x ~ t + re(1 + t || id, prior = pcSd(u = "data"))))
  expect_warning(buildInlaFormula(f, prepareData(f, d)), "small \\(or zero: singular fit\\) SD")
})

test_that("data-scaled scalar priors: other families are an error naming the block", {
  b <- dataBlock(x1 ~ t + re(1 | id, prior = pcSd(u = "data")), gaussianFamily())
  p <- b$p
  p$family <- betaFamily()
  expect_error(resolveEffectPriors(pcSd(u = "data"), p, b$df, "(Intercept)", p$re[[1]]),
               "pcSd\\(u = \"data\"\\) needs a Gaussian, lognormal, binomial or Poisson response .*numeric parameter .*block `1 \\| id` \\(formula of 'x1'\\)")
  p$family <- negBinomialFamily()
  expect_error(resolveEffectPriors(gammaPrec(rate = "data"), p, b$df, "(Intercept)", p$re[[1]]),
               "gammaPrec\\(rate = \"data\"\\) needs a Gaussian")
})

# ---------------------------------------------------------------------------
# INLA formula, LaTeX, notes

test_that("buildInlaFormula(): each effect of a `||` block gets its own resolved hyper in f()", {
  f <- collectFormulas(list(y ~ t + scaled(x3)), list(x3 ~ t + re(1 + t || id, prior = pcSd(u = "data"))))
  inl <- buildInlaFormula(f, prepareData(f, fxData()$long))
  ref <- lmerX3()
  rnd <- inl$random$x3
  expect_length(rnd, 2L)
  for (e in 1:2) {
    expect_equal(rnd[[e]]$prior$args$source, "data")
    expect_equal(rnd[[e]]$prior$args$u, dataRules()$pcU * ref$s[e], tolerance = 1e-4)
    h <- hyperOf(rnd[[e]]$terms[[1]])
    expect_equal(h$prec$prior, "pc.prec")
    expect_equal(h$prec$param, c(dataRules()$pcU * ref$s[e], dataRules()$pcAlpha), tolerance = 1e-4)
  }
})

test_that("buildInlaFormula(): gamma / half-normal / half-Cauchy hypers from the data", {
  ref <- lmerX3()
  R <- dataRules()
  build <- function(pr) {
    f <- collectFormulas(list(y ~ t + scaled(x3)), list(x3 ~ t + re(1 + t || id, prior = pr)))
    buildInlaFormula(f, prepareData(f, fxData()$long))$random$x3
  }
  h <- hyperOf(build(gammaPrec(rate = "data"))[[2]]$terms[[1]])
  expect_equal(h$prec$prior, "loggamma")
  expect_equal(h$prec$param, c(1, 1 * ref$s[2]^2), tolerance = 1e-4)
  h <- hyperOf(build(halfNormalSd("data"))[[1]]$terms[[1]])
  expect_equal(h$prec$prior, "logtnormal")
  expect_equal(h$prec$param, c(0, 1 / (R$halfNormalScale * ref$s[1])^2), tolerance = 1e-4)
  h <- hyperOf(build(halfCauchySd("data"))[[1]]$terms[[1]])
  expect_equal(h$prec$prior, halfCauchyExpression(R$halfCauchyScale * ref$s[1]), tolerance = 1e-4)
})

test_that("latexEquation(): \\hat{s}_{data} without data, the resolved numbers per effect with data", {
  f <- collectFormulas(list(y ~ t + scaled(x3)), list(x3 ~ t + re(1 + t || id, prior = pcSd(u = "data"))))
  without <- format(latexEquation(f))
  expect_match(without, r"(\mathrm{Exp}\left(\frac{3}{6\,\hat{s}_{k,\mathrm{data}}}\right))", fixed = TRUE)
  expect_match(without, r"(> 6\,\hat{s}_{k,\mathrm{data}}\right) = 0.05, \quad k = 0, 1)", fixed = TRUE)
  dat <- prepareData(f, fxData()$long)
  with <- latexEquation(f, dat)
  expect_false(grepl("hat{s}", with, fixed = TRUE))
  ref <- lmerX3()
  for (e in 1:2) {
    u <- numLatex(dataRules()$pcU * ref$s[e])
    expect_match(with, paste0(r"(\left(\tau^{(x_{3})}_{)", e - 1L, r"(,\mathrm{id}}\right)^{-1/2} > )", u,
                              r"(\right) = 0.05)"), fixed = TRUE)
  }
  expectWellFormed(with)
  expectWellFormed(latexEquation(f))
  # single effect: no k index on s
  f1 <- collectFormulas(list(y ~ t + scaled(x1)), list(x1 ~ t + re(1 | id, prior = halfNormalSd("data"))))
  l1 <- format(latexEquation(f1))
  expect_match(l1, r"(\mathcal{N}^{+}\left(0, \left(\left(2\,\hat{s}_{\mathrm{data}}\right)^{-2}\right)^{-1}\right))",
               fixed = TRUE)
  f2 <- collectFormulas(list(y ~ t + scaled(x1)), list(x1 ~ t + re(1 | id, prior = gammaPrec(rate = "data"))))
  expect_match(format(latexEquation(f2)), r"(\mathrm{Gamma}\left(1, \hat{s}_{\mathrm{data}}^{2}\right))", fixed = TRUE)
  f3 <- collectFormulas(list(y ~ t + scaled(x1)), list(x1 ~ t + re(1 | id, prior = halfCauchySd("data"))))
  expect_match(format(latexEquation(f3)), r"(\mathrm{C}^{+}\left(0, 2\,\hat{s}_{\mathrm{data}}\right))", fixed = TRUE)
})

test_that("priorNote(): the rule before resolution, the lme4 SD after", {
  expect_match(priorNote(pcSd(u = "data")), "with u 6 times the lme4 SD of the effect")
  expect_match(priorNote(gammaPrec(rate = "data")), "mean 1/s\\^2")
  r <- dataScalarPrior(halfCauchySd("data"), 0.25)
  expect_match(priorNote(r), "half-Cauchy prior on the SD with scale 0.5 - from the data: the lme4 SD of the effect is 0.25")
  expect_false(grepl("from the data", priorNote(halfCauchySd(1))))
})

# ---------------------------------------------------------------------------
# INLA fit

test_that("jsm(): a fit with data-scaled priors; tooltips show the resolved priors", {
  skip_if_not_installed("INLA")
  d <- fxData()$long
  m <- suppressWarnings(jsm(list(y ~ t + scaled(x3)),
                            list(x3 ~ t + re(1 + t || id, prior = pcSd(u = "data"))),
                            data = d, family = list(y = betaFamily())))
  skip_if(is.null(inlaFit(m)), "INLA fit failed")
  ref <- lmerX3()
  vc <- varianceComponents(m, "x3")$table
  sds <- vc$mean[vc$parameter == "sd" & vc$effect %in% c("(Intercept)", "t")]
  expect_length(sds, 2L)
  expect_true(all(sds > 0.5 * ref$s & sds < 2 * ref$s))          # near the lme4 SDs
  # the hyper reached INLA: pc.prec with the resolved parameters
  hy <- m$fit$all.hyper$random
  pcs <- Filter(function(h) isTRUE(as.character(h$hyper$theta$prior) == "pc.prec"), hy)
  expect_length(pcs, 2L)
  expect_equal(sort(vapply(pcs, function(h) as.numeric(h$hyper$theta$param[1]), 0)),
               sort(dataRules()$pcU * ref$s), tolerance = 1e-4)
  tt <- format(latexEquation(m, annotate = "tooltip"))
  expect_match(tt, "from the data: the lme4 SD of the effect is", fixed = TRUE)
  est <- estimatesTable(m)
  expect_equal(sum(est$type == "sd" & est$response == "x3"), 2L)     # estimatesTable() unaffected
  expect_silent(withNullDevice(plotModel(m)))                          # plotModel() unaffected
})
