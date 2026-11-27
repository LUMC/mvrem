# predict(method = "sample") (E12): predictions from joint posterior draws of the fit,
# without a new INLA run. Fixtures: helper-fixtures.R, helper-predictSample.R.
#
# INLA runs of this file: "multi" and "priors" fitted with config = TRUE, and the refit
# predictions they are compared with ("multi" type = "linear"; "priors" type = "mean";
# "multi" type = "mean" is the shared fxPredict()).

# ---------------------------------------------------------------------------
# validation against the refit method

test_that("sample and refit agree within Monte-Carlo error ('multi': binomial + trials, poisson + exposure, both types)", {
  for (type in c("mean", "linear")) {
    pr <- fxPredictRefit("multi", type)
    ps <- fxPredictSample("multi", type)
    expect_equal(ps[c("response", "id", "t")], pr[c("response", "id", "t")])
    d <- fxPredictDiff(pr, ps)
    # all four responses, new subjects (N1-N3) and an existing one at new times (S01)
    expect_setequal(d$response, c("k", "cnt", "x1", "x2"))
    expect_setequal(d$id, c("N1", "N2", "N3", "S01"))
    expect_lt(max(abs(d$z)), 4)                       # means: within 4 Monte-Carlo SEs
    expect_true(all(d$sdRatio > 0.9 & d$sdRatio < 1.1))
    # 95% bounds: within 0.3 SD. Their Monte-Carlo SE is about 0.04 SD (4000 draws), but the
    # refit is a separate INLA run whose tails vary between machines: the same seeded draws were
    # at most 0.11 SD from the refit in the dev container and 0.21 SD in GitLab CI (2026-11-21)
    b <- pmax(abs(d$dLo), abs(d$dHi))
    expect_lt(max(b), 0.3, label = paste0("max bound difference (", type, ", ",
                                          paste(d[which.max(b), c("response", "id", "t")], collapse = " "), ")"))
    # the time of both methods (refit: cached fxPredict() has none; it is a full INLA run)
    if (!is.null(attr(pr, "elapsed"))) expect_lt(attr(ps, "elapsed"), attr(pr, "elapsed"))
  }
})

test_that("sample and refit agree on the 'priors' model (outcome RE on site, `||` block, two blocks on id)", {
  pr <- fxPredictRefit("priors", "mean")
  ps <- fxPredictSample("priors", "mean")
  d <- fxPredictDiff(pr, ps)
  expect_setequal(d$response, c("yg", "x1", "x2", "x3"))
  # the refit is a separate INLA run, whose posterior differs from the fit's by more than the
  # Monte-Carlo error on this model (max |z| 2.4-5.7 between runs, also before any change of
  # mvrem: the threshold 4 of 'multi' failed in full test runs): the means are compared in
  # refit SDs (z / sqrt(n); 0.09 at |z| = 5.7) instead
  expect_lt(max(abs(d$z)) / sqrt(fxNsamples), 0.15)
  # the new subjects of x2 (negative binomial; random intercept with a half-Cauchy prior and
  # an independent slope with a half-normal one): the posterior of the slope precision is very
  # skewed and INLA integrates over only 2-5 points of the hyperparameters here, so the fit and
  # the refit disagree on it (mean 96 in the fit, 155 in the refit) and the spread of a new
  # subject varies between INLA runs by more than any fixed tolerance (SD ratio 0.84-0.91 in one
  # run, beyond 0.75 in another): the refit is no reference for their spread, only the means
  # (checked above) are compared; elsewhere as for 'multi'
  hard <- d$response == "x2" & d$id != "S01"
  expect_true(any(hard))
  # (wider than for 'multi', a coarse check: full test runs gave SD ratios beyond 0.9 / 1.1 and
  # bounds up to 0.37 SD apart, from the run-to-run variation of the refit above; 'multi' keeps
  # the tight checks)
  expect_true(all(d$sdRatio[!hard] > 0.8 & d$sdRatio[!hard] < 1.25))
  expect_lt(max(abs(c(d$dLo[!hard], d$dHi[!hard]))), 0.5)
})

test_that("the posterior draws reproduce INLA's marginals of the fixed effects and random effects", {
  m <- fxModelConfig("multi")
  withr::with_seed(3, dr <- posteriorDraws(m$fit, fxNsamples))
  sf <- m$fit$summary.fixed
  mine <- vapply(rownames(sf), function(n) dr$latent[dr$rows(n), ], numeric(fxNsamples))
  expect_lt(max(abs(colMeans(mine) - sf$mean) / (sf$sd / sqrt(fxNsamples))), 4)
  expect_equal(apply(mine, 2, sd), sf$sd, tolerance = 0.05, ignore_attr = TRUE)
  sr <- m$fit$summary.random$x1_1_id                # the correlated block: 2 x 40 effects
  b <- dr$latent[dr$rows("x1_1_id"), ]
  expect_equal(nrow(b), nrow(sr))
  expect_lt(max(abs(rowMeans(b) - sr$mean) / (sr$sd / sqrt(fxNsamples))), 4.5)
  expect_equal(apply(b, 1, sd), sr$sd, tolerance = 0.05)
  # hyperparameters on the user scale, named as summary.hyperpar
  expect_equal(rownames(dr$hyper), rownames(m$fit$summary.hyperpar))
  expect_equal(rowMeans(dr$hyper)[["Beta for k_FROM_x1_id_t"]],
               m$fit$summary.hyperpar["Beta for k_FROM_x1_id_t", "mean"], tolerance = 0.05)
  expect_true(all(dr$hyper["Rho1:2 for x1_1_id", ] > -1 & dr$hyper["Rho1:2 for x1_1_id", ] < 1))
})

# ---------------------------------------------------------------------------
# the result

test_that("predict(method = 'sample'): same class, columns and attributes as the refit method", {
  pr <- fxPredict()
  ps <- fxPredictSample("multi", "mean")
  m <- fxModelConfig("multi")
  expect_identical(class(ps), class(pr))
  expect_identical(names(ps), names(pr))
  expect_identical(vapply(ps, typeof, ""), vapply(pr, typeof, ""))
  expect_equal(attr(ps, "type"), "mean")
  expect_equal(attr(ps, "method"), "sample")
  expect_equal(attr(pr, "method"), "refit")
  expect_identical(attr(ps, "model"), m)                   # the original model, not refitted
  expectOrderedStats(ps)
  expect_true(all(ps$mode >= ps$q0.025 & ps$mode <= ps$q0.975))   # density-estimate mode
  k <- ps[ps$response == "k", ]
  expect_true(all(k$q0.025 > 0 & k$q0.975 < 1))                   # proportion scale
  expect_true(all(ps$q0.025[ps$response %in% c("cnt", "x2")] > 0)) # rate scale
  expect_equal(attr(fxPredictSample("multi", "linear"), "type"), "linear")
})

test_that("a new subject is one draw per subject: shared by its rows and by the linked outcomes", {
  m <- fxModelConfig("multi")
  # the same new subject twice at the same time gives identical summaries; a new subject's
  # longitudinal prediction is wider than the existing subject's
  nd <- data.frame(id = c("N9", "N9", "S01"), t = c(1, 1, 1), n = 5L, pt = 1, grp = "a")
  p <- predict(m, nd, method = "sample", nsamples = 500, seed = 2)
  for (r in unique(p$response)) {
    q <- p[p$response == r, ]
    expect_equal(q[1, statCols], q[2, statCols], ignore_attr = TRUE)
  }
  x1 <- p[p$response == "x1", ]
  expect_gt(x1$sd[1], 2 * x1$sd[3])
  # an existing subject at a fitted time: close to INLA's fitted value of that row
  fv <- fittedValues(m, "x1")
  f1 <- fv[fv$id == "S01" & fv$t == 1, ]
  expect_equal(x1$mean[3], f1$mean, tolerance = 0.02)
  expect_equal(x1$sd[3], f1$sd, tolerance = 0.1)
})

test_that("predict(method = 'sample'): seed reproduces the draws and restores the RNG state", {
  m <- fxModelConfig("multi")
  nd <- data.frame(id = c("N1", "S01"), t = 1)
  set.seed(11)
  before <- .Random.seed
  a <- predict(m, nd, response = "x1", method = "sample", nsamples = 200, seed = 5)
  expect_identical(.Random.seed, before)
  b <- predict(m, nd, response = "x1", method = "sample", nsamples = 200, seed = 5)
  expect_identical(a, b)
  c <- predict(m, nd, response = "x1", method = "sample", nsamples = 200, seed = 6)
  expect_false(identical(a$mean, c$mean))
  expect_identical(.Random.seed, before)
  predict(m, nd, response = "x1", method = "sample", nsamples = 50)
  expect_false(identical(.Random.seed, before))     # seed = NULL: draws from the current state
})

test_that("plot() of sampled predictions works as for the refit method", {
  ps <- fxPredictSample("multi", "mean")
  withNullDevice(expect_invisible(plot(ps, interactive = FALSE)))
  skip_if_not_installed("htmlwidgets")
  w <- plot(ps, interactive = TRUE)
  expect_s3_class(w, "htmlwidget")
})

# ---------------------------------------------------------------------------
# method = "auto" (Q24): sample with config, refit without

test_that("predict(method = 'auto') samples when the fit has INLA's configurations", {
  m <- fxModelConfig("multi")
  expect_false(is.null(m$fit$misc$configs))
  cap <- mockInlaRun()
  cap$model <- NULL
  p <- predict(m, fxNewdata, nsamples = 200, seed = 1)       # no `method`: auto
  expect_null(cap$model)                                       # no INLA run
  expect_equal(attr(p, "method"), "sample")
  expect_equal(attr(p, "methodChosen"), "auto")
  expect_equal(attr(p, "nsamples"), 200L)
  expect_identical(attr(p, "model"), m)
  # the same draws as an explicit method = "sample"
  pe <- predict(m, fxNewdata, method = "sample", nsamples = 200, seed = 1)
  expect_equal(as.data.frame(p), as.data.frame(pe), ignore_attr = TRUE)
  expect_equal(attr(pe, "methodChosen"), "sample")
  # printed: the method on one line after the table
  out <- capture.output(print(p))
  expect_match(out[length(out)], "^Predictions by sampling \\(200 posterior draws\\); chosen automatically: the fit has INLA's configurations\\.$")
  expect_match(utils::tail(capture.output(print(pe)), 1), "^Predictions by sampling \\(200 posterior draws\\)\\.$")
  # options of an INLA run select the refit; an explicit refit refits
  pr <- predict(m, fxNewdata, response = "x1", threads = 1)
  expect_equal(attr(pr, "method"), "refit")
  expect_equal(cap$options$threads, 1)
  cap$model <- NULL
  pr2 <- predict(m, fxNewdata, response = "x1", method = "refit")
  expect_false(is.null(cap$model))
  expect_equal(attr(pr2, "method"), "refit")
})

test_that("predict(method = 'auto') refits when the fit has no configurations (INLA mocked)", {
  m <- fxModel("multi")
  expect_null(m$fit$misc$configs)
  cap <- mockInlaRun()
  p <- predict(m, fxNewdata, response = "x1")
  expect_false(is.null(cap$model))                            # the (mocked) INLA run happened
  expect_equal(nrow(modelData(cap$model)$long$x1), fxCounts()$rows + nrow(fxNewdata))
  expect_equal(attr(p, "method"), "refit")
  expect_equal(attr(p, "methodChosen"), "auto")
  expect_null(attr(p, "nsamples"))
  expect_match(utils::tail(capture.output(print(p)), 1),
               "^Predictions by refitting with INLA; chosen automatically: the fit has no INLA configurations\\.$")
  # the cached refit prediction of the fixture was chosen the same way
  expect_equal(attr(fxPredict(), "method"), "refit")
})

test_that("predict(method = 'auto') falls back to the refit when sampling is not supported", {
  m <- fxModelConfig("multi")
  cap <- mockInlaRun()
  testthat::local_mocked_bindings(predictSampled = function(...) sampleUnsupported("`method = \"sample\"` does not support X; use `method = \"refit\"`."))
  expect_message(p <- predict(m, fxNewdata, response = "x1"), "sampling is not available for this model \\(`method = \"sample\"` does not support X\\); refitting")
  expect_equal(attr(p, "method"), "refit")
  expect_false(is.null(cap$model))
  # an explicit method = "sample" keeps the error
  expect_error(predict(m, fxNewdata, response = "x1", method = "sample"), "does not support X", class = "mvremSampleUnsupported")
})

# ---------------------------------------------------------------------------
# errors

test_that("predict(method = 'sample') without config: a clear error with the hints", {
  m <- fxModel("multi")
  expect_null(m$fit$misc$configs)
  expect_error(predict(m, fxNewdata, method = "sample"), "control.compute = list\\(config = TRUE\\)")
  expect_error(predict(m, fxNewdata, method = "sample"), "method = \"refit\"")
})

test_that("predict(method = 'sample'): argument checks; unseen levels and survival responses as for refit", {
  m <- fxModelConfig("multi")
  nd <- data.frame(id = "S01", t = 1)
  expect_error(predict(m, nd, response = "x1", method = "sample", nsamples = 1), "`nsamples` must be")
  expect_error(predict(m, nd, response = "x1", method = "sample", nsamples = 10.5), "`nsamples` must be")
  expect_error(predict(m, nd, response = "x1", method = "sample", seed = "a"), "`seed` must be")
  expect_error(predict(m, nd, response = "x1", method = "sample", threads = 2),
               "`threads` only apply to `method = \"refit\"`")
  expect_error(predict(m, nd, method = "nope"), "should be one of")
  nd2 <- data.frame(id = "S01", t = 1, n = 5, pt = 1, grp = "zzz")
  expect_error(predict(m, nd2, response = "cnt", method = "sample"), "Unknown level\\(s\\) of `grp`")
  ms <- fxModel("surv")
  expect_error(predict(ms, data.frame(id = "S01", t = 0, status = 1), response = "time", method = "sample"),
               "not available for survival responses")
})

test_that("inverseLink(): INLA's inverse links; links with parameters are rejected", {
  skip_if_not_installed("INLA")
  expect_equal(inverseLink(binomialFamily())(matrix(0, 1, 2)), matrix(0.5, 1, 2))
  expect_equal(inverseLink(poissonFamily())(matrix(log(2))), matrix(2))
  expect_equal(inverseLink(binomialFamily(link = "probit"))(matrix(1)), matrix(pnorm(1)), tolerance = 1e-6)
  expect_equal(inverseLink(gaussianFamily())(matrix(3)), matrix(3))
  expect_equal(inverseLink(inlaFamily("poisson"))(matrix(0)), matrix(1))         # INLA's default link
  expect_error(inverseLink(binomialFamily(link = "robit")), "does not support the link 'robit'")
})
