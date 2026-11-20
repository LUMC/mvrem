# Tests of Q3: every result follows the order in which the formulae were passed to jsm()
# (one list: as listed, outcomes and longitudinal responses interleaved; pair style:
# outcomes, then longitudinal). A `response =` argument selects, it does not reorder.
# One extra INLA run: the "multi" fixture passed as an interleaved list (cached below).

# the "multi" fixture as one list: x1, k, x2, cnt
fxInterleavedForm <- function() {
  s <- fxSpec("multi")
  list(s$formLong[[1]], s$formOut[[1]], s$formLong[[2]], s$formOut[[2]])
}
interleavedOrder <- c("x1", "k", "x2", "cnt")

fxInterleaved <- function(fit = TRUE) {
  key <- paste0("interleaved_", if (fit) "fit" else "nofit")
  if (fit) {
    testthat::skip_if_not_installed("INLA")
    testthat::skip_if_not_installed("Matrix")
  }
  if (is.null(.fx[[key]])) {
    s <- fxSpec("multi")
    m <- suppressWarnings(suppressMessages(jsm(
      form = fxInterleavedForm(), data = s$data, family = s$family, fit = fit, criteria = fit,
      inlaOptions = if (fit) list(control.compute = list(config = TRUE)) else list())))
    .fx[[key]] <- if (fit && is.null(inlaFit(m))) "failed" else m
  }
  if (identical(.fx[[key]], "failed")) testthat::skip("INLA fit of the interleaved model failed")
  .fx[[key]]
}

# the responses of a table / panel list in order of first appearance
firstSeen <- function(v) unique(v[!is.na(v) & nzchar(v)])

# ---- without a fit -------------------------------------------------------------------

test_that("the order passed is stored: one list as listed, pair outcomes first; old objects fall back", {
  m <- fxInterleaved(fit = FALSE)
  expect_identical(m$forms$order, interleavedOrder)
  expect_identical(modelResponses(m), interleavedOrder)
  expect_identical(modelResponses(modelFormulas(m)), interleavedOrder)
  expect_identical(modelResponses(fxModel("multi", fit = FALSE)), c("k", "cnt", "x1", "x2"))
  old <- m
  old$forms$order <- NULL                       # a model object made before Q3
  expect_identical(modelResponses(old), c("k", "cnt", "x1", "x2"))
  # the parts keep their own order; the INLA likelihoods / stacks stay outcomes first
  expect_identical(names(m$forms$out), c("k", "cnt"))
  expect_identical(names(m$forms$long), c("x1", "x2"))
  expect_identical(names(m$inla$responses), c("k", "cnt", "x1", "x2"))
})

test_that("dataSummary(), print(), plot panels of an unfitted model follow the order passed", {
  m <- fxInterleaved(fit = FALSE)
  expect_identical(dataSummary(m)$equations$response, interleavedOrder)
  expect_identical(firstSeen(dataSummary(m)$missing$response), interleavedOrder)
  out <- capture.output(print(m))
  expect_match(grep("^Observations:", out, value = TRUE), "x1: [0-9]+, k: [0-9]+, x2: [0-9]+, cnt: [0-9]+")
  d <- plotLinesData(m, NULL, NULL, NULL, 5, NULL, TRUE)
  expect_identical(vapply(d$panels, `[[`, "", "response"), interleavedOrder)
  d2 <- plotLinesData(m, c("cnt", "x2", "k"), NULL, NULL, 5, NULL, TRUE)
  expect_identical(vapply(d2$panels, `[[`, "", "response"), c("k", "x2", "cnt"))
  # plotModel(): markers left, outcomes right, each column in the order passed (unchanged)
  withNullDevice(expect_silent(plotModel(m)))
})

# ---- fitted --------------------------------------------------------------------------

test_that("the interleaved model is the 'multi' fixture: same estimates, only the order differs", {
  m <- fxInterleaved()
  ref <- fxModel("multi")
  fe <- fixedEffects(m)
  fr <- fixedEffects(ref)
  k <- match(paste(fe$response, fe$term), paste(fr$response, fr$term))
  expect_false(anyNA(k))
  expect_equal(fe$mean, fr$mean[k], tolerance = 0.02)
})

test_that("every result table follows the order passed; response = selects without reordering", {
  m <- fxInterleaved()
  sel <- c("cnt", "x2", "k")                           # given in another order
  selOrder <- c("k", "x2", "cnt")
  expect_identical(firstSeen(fixedEffects(m)$response), interleavedOrder)
  expect_identical(names(coef(m))[1], paste0("x1: ", fixedEffects(m)$term[1]))
  h <- hyperparameters(m)
  # grouped by response: x1 (residual precision, RE block), k (gammas), x2 (RE), cnt (gamma)
  expect_identical(firstSeen(h$response), interleavedOrder)
  expect_false(is.unsorted(match(h$response, interleavedOrder), na.rm = TRUE))
  expect_identical(h$parameter[h$response == "x1"][1], "residual precision")   # INLA's order within
  a <- associations(m)
  expect_identical(paste(a$outcome, a$longitudinal), c("k x1", "k x2", "cnt x2"))
  s <- summary(m)
  expect_identical(firstSeen(s$fixed$response), interleavedOrder)
  re <- randomEffects(m)
  expect_identical(firstSeen(re$response), c("x1", "x2"))
  expect_identical(firstSeen(randomEffects(m, c("x2", "x1"))$response), c("x1", "x2"))
  vc <- varianceComponents(m)
  expect_identical(firstSeen(vc$table$response), c("x1", "x2"))
  expect_length(vc$marginals, nrow(vc$table))
  expect_equal(vc$table$response[vc$table$block == "residual"], "x1")
  fv <- fittedValues(m)
  expect_identical(firstSeen(fv$response), interleavedOrder)
  expect_identical(firstSeen(fittedValues(m, sel)$response), selOrder)
  expect_identical(firstSeen(residuals(m)$response), interleavedOrder)
  expect_identical(firstSeen(residuals(m, sel)$response), selOrder)
  mf <- modelFit(m)
  expect_identical(mf$response, c(interleavedOrder, "overall"))
  et <- estimatesTable(m)
  expect_identical(firstSeen(et$response), interleavedOrder)
  expect_identical(firstSeen(estimatesTable(m, sel)$response), selOrder)
})

test_that("every plot follows the order passed", {
  skipNoWidgets()
  m <- fxInterleaved()
  w <- plot(m, interactive = TRUE)
  expect_identical(substr(panelTitles(w), 1, 2), c("x1", "k ", "x2", "cn"))
  wr <- plotResiduals(m, interactive = TRUE)
  expect_identical(vapply(plotResidualsData(m, NULL, "time", NULL, NULL, 5, NULL)$panels, `[[`, "", "response"),
                   interleavedOrder)
  expect_identical(vapply(plotResidualsData(m, c("cnt", "x1"), "time", NULL, NULL, 5, NULL)$panels, `[[`, "",
                          "response"), c("x1", "cnt"))
  expect_length(panelTitles(wr), 4L)
  expect_identical(panelTitles(plotFixedEffects(m, interactive = TRUE)),
                   paste0(interleavedOrder, ": fixed effects"))
  expect_identical(panelTitles(plotFixedEffects(m, response = c("cnt", "x1"), interactive = TRUE)),
                   paste0(c("x1", "cnt"), ": fixed effects"))
  expect_identical(panelTitles(plotFixedEffects(m, type = "density", interactive = TRUE)),
                   paste0(interleavedOrder, ": fixed effects"))
  expect_identical(substr(panelTitles(plotRandomEffects(m, interactive = TRUE)), 1, 2), c("x1", "x1", "x2"))
  expect_identical(substr(panelTitles(plotRandomEffects(m, c("x2", "x1"), interactive = TRUE)), 1, 2),
                   c("x1", "x1", "x2"))
  expect_identical(panelTitles(plotVarianceComponents(m, interactive = TRUE)),
                   c("x1 | id [(Intercept), t]", "correlations", "x2 | id [(Intercept)]", "residual"))
  # associations: forest labels and density curves in the same (associations()) order
  lab <- c("x1 -> k", "x2 -> k", "x2 -> cnt")
  af <- plotAssociations(m, interactive = TRUE)
  expect_identical(af$x$panels[[1]]$ylabels$text, lab)
  ad <- plotAssociations(m, type = "density", interactive = TRUE)
  cv <- ad$x$panels[[1]]$curves
  expect_identical(vapply(cv, `[[`, "", "label"), lab)
  a <- associations(m)
  means <- vapply(cv, function(c) sum(c$x * c$y) / sum(c$y), 0)    # density means ~ posterior means
  expect_equal(means, a$mean, tolerance = 0.1)
})

test_that("predict() and its plot follow the order passed; auto samples with config", {
  m <- fxInterleaved()
  p <- predict(m, fxNewdata, nsamples = 200, seed = 1)
  expect_identical(attr(p, "method"), "sample")
  expect_identical(firstSeen(p$response), interleavedOrder)
  nd <- list(cnt = fxNewdata, x1 = fxNewdata)              # list order differs from the model's
  pl <- predict(m, nd, nsamples = 200, seed = 1)
  expect_identical(firstSeen(pl$response), c("x1", "cnt"))
  ps <- predict(m, fxNewdata, response = c("cnt", "x2"), nsamples = 200, seed = 1)
  expect_identical(firstSeen(ps$response), c("x2", "cnt"))
  skipNoWidgets()
  w <- plot(p, interactive = TRUE)
  expect_identical(substr(panelTitles(w), 1, 2), c("x1", "k ", "x2", "cn"))
})
