# E13: observation weights `y | weights(w) ~ ...` (INLA `scale`): parsing, family checks,
# data checks, stacks, inla() arguments, naming, display, LaTeX, predict() and a
# simulation fit.

# fxSmallData() (a copy) with a weight column `w` (default 1 / 4 alternating).
weightData <- function(w = rep(c(1, 4), 6)) {
  d <- fxSmallData()
  d$w <- w
  d
}

# x | weights(w) ~ t + (1 | id), linked to y ~ t + scaled(x), unfitted.
weightModel <- function(lhs = "x | weights(w)", d = weightData(), family = NULL, out = y ~ t + scaled(x), ...) {
  jsm(list(out), list(stats::as.formula(paste(lhs, "~ t + (1 | id)"))), data = d,
      family = family, fit = FALSE, ...)
}

errMsg <- function(expr) conditionMessage(expect_error(expr))

# ---- parsing ----------------------------------------------------------------------

test_that("parseJsmFormula(): weights() is a family input, its variables have role aux", {
  p <- parseJsmFormula(x | weights(w) ~ t + (1 | id), "long")
  expect_identical(p$response, "x")
  expect_identical(p$lhs, "x")
  expect_identical(p$aux, list(weights = quote(w)))
  expect_identical(p$variables$aux, "w")
  expect_false("w" %in% p$variables$covariate)
  # any expression of the data, combined with another input
  p2 <- parseJsmFormula(log(x) | weights(1 / v^2) + trials(n) ~ t, "long")
  expect_identical(names(p2$aux), c("weights", "trials"))
  expect_identical(p2$aux$weights, quote(1 / v^2))
  expect_setequal(p2$variables$aux, c("v", "n"))
})

test_that("parseJsmFormula(): malformed weights() terms are errors", {
  expect_error(parseJsmFormula(x | weights(w) + weights(v) ~ t, "long"), "weights\\(\\) is given more than once")
  expect_error(parseJsmFormula(x | weights(a = w) ~ t, "long"), "weights\\(\\) takes exactly one unnamed argument")
  expect_error(parseJsmFormula(x | weights(w, v) ~ t, "long"), "exactly one unnamed argument")
  expect_error(parseJsmFormula(x | weight(w) ~ t, "long"), "weights\\(\\)")   # listed in the hint
})

# ---- families ---------------------------------------------------------------------

test_that("families: weights() is accepted exactly by the likelihoods that document `scale`", {
  takes <- function(f) isTRUE(f$weights)
  expect_true(takes(gaussianFamily()))
  expect_true(takes(studentTFamily()))
  expect_true(takes(gammaFamily()))
  expect_true(takes(betaFamily()))
  expect_false(takes(lognormalFamily()))
  expect_false(takes(poissonFamily()))
  expect_false(takes(negBinomialFamily()))
  expect_false(takes(binomialFamily()))
  expect_false(takes(betaBinomialFamily()))
  expect_false(takes(zeroInflatedPoissonFamily()))
  expect_false(takes(zeroInflatedNegBinomialFamily()))
  expect_false(takes(weibullSurvFamily()))
  expect_false(takes(exponentialSurvFamily()))
})

test_that("collectFormulas(): weights() with a family that does not use it is an error", {
  f <- list(y ~ t + scaled(x))
  msg <- errMsg(collectFormulas(f, list(x | weights(w) ~ t), family = list(x = poissonFamily())))
  expect_match(msg, "weights() on the left-hand side of the formula for 'x' is not used by family 'poisson' (it takes exposure()).",
               fixed = TRUE)
  msg <- errMsg(collectFormulas(f, list(x | weights(w) ~ t), family = list(x = lognormalFamily())))
  expect_match(msg, "not used by family 'lognormal'.", fixed = TRUE)
  expect_false(grepl("it takes", msg))
  expect_error(collectFormulas(f, list(x | weights(w) ~ t), family = list(x = "binomial")),
               "not used by family 'binomial' \\(it takes trials\\(\\)\\)")
  # an input the family lacks, next to weights(): weights() itself is fine
  expect_error(collectFormulas(f, list(x | weights(w) + exposure(e) ~ t)),
               "^exposure\\(\\) on the left-hand side.*'gaussian'\\.$")
})

test_that("collectFormulas(): weights() for outcome and long responses of every accepting family", {
  f <- collectFormulas(list(y | weights(v) ~ t + scaled(x)), list(x | weights(w) ~ t),
                       family = list(y = betaFamily(), x = studentTFamily()))
  expect_identical(f$out$y$aux$weights, quote(v))
  expect_identical(f$long$x$aux$weights, quote(w))
  expect_true(all(c("v", "w") %in% f$variables$variable[f$variables$role == "aux"]))
  expect_no_error(collectFormulas(list(y ~ t + scaled(x)), list(x | weights(w) ~ t),
                                  family = list(x = gammaFamily())))
})

test_that("inlaFamily(): weights() combined with trials() / exposure()", {
  skip_if_not_installed("INLA")
  f <- collectFormulas(list(k | trials(n) + weights(w) ~ t + scaled(x)), list(x ~ t),
                       family = list(k = inlaFamily("binomial")))
  expect_identical(names(f$out$k$aux), c("trials", "weights"))
  d <- weightData()
  inl <- buildInlaFormula(f, prepareData(f, d))
  expect_identical(inl$inputs$k, c(Ntrials = "k_Ntrials", scale = "k_scale"))
  st <- inl$stacks$k$main
  expect_identical(names(st)[1:3], c("k_k", "k_Ntrials", "k_scale"))
  a <- inlaArgumentList(inl)
  expect_equal(a$Ntrials[seq_len(nrow(d))], d$n)
  expect_equal(a$scale[seq_len(nrow(d))], d$w)
  expect_true(all(a$scale[-seq_len(nrow(d))] == 1))
})

# ---- data checks ------------------------------------------------------------------

test_that("prepareData(): weights must be positive and finite; missing ones need na.omit", {
  for (bad in list(c(0, rep(1, 11)), c(-1, rep(1, 11)), c(Inf, rep(1, 11))))
    expect_error(weightModel(d = weightData(bad)), "weights\\(\\) of 'x' must be positive and finite", info = bad[1])
  expect_error(weightModel(d = weightData(as.character(rep(1, 12)))), "weights\\(\\) of 'x' must be numeric")
  d <- weightData(c(NA, rep(1, 11)))
  expect_error(weightModel(d = d), "na.omit")
  m <- weightModel(d = d, na.action = na.omit)
  expect_equal(nrow(modelData(m)$long$x), 11L)
  expect_true("w" %in% names(modelData(m)$long$x))
  # an expression: evaluated on the data
  expect_no_error(weightModel("x | weights(1 / w)"))
})

test_that("summarizeData(): the weight variable is listed with its missing values", {
  d <- weightData(c(NA, rep(1, 11)))
  forms <- collectFormulas(list(y ~ t + scaled(x)), list(x | weights(w) ~ t + (1 | id)))
  s <- summarizeData(forms, d)
  miss <- s$missing[s$missing$response == "x", ]
  expect_equal(miss$missing[miss$variable == "w"], 1L)
  eq <- s$equations[s$equations$response == "x", ]
  expect_equal(c(eq$rows, eq$complete, eq$dropped), c(12L, 11L, 1L))
  expect_false("w" %in% s$missing$variable[s$missing$response == "y"])
})

# ---- stacks, inla() arguments, naming ---------------------------------------------

test_that("buildInlaFormula(): the scale column follows the responses, NA in zero-equation rows", {
  m <- weightModel()
  inl <- m$inla
  expect_identical(inl$inputs$x, c(scale = "x_scale"))
  expect_length(inl$inputs$y, 0L)
  main <- inl$stacks$x$main
  expect_identical(names(main)[1:3], c("x_x", "x_x_0eq", "x_scale"))
  expect_equal(main$x_scale, weightData()$w)
  zero <- inl$stacks$x$y
  expect_true(all(is.na(zero$x_scale)))
  expect_false("x_scale" %in% names(inl$stacks$y$main))
  # the formula does not use the column
  expect_false("x_scale" %in% all.vars(inl$formula))
})

test_that("inlaArgumentList(): `scale` is one value per row, the weights at their rows, 1 elsewhere", {
  m <- weightModel()
  st <- stackRows(m$inla)
  a <- inlaArgumentList(m$inla)
  expect_length(a$scale, nrow(st))
  rows <- st$.block == "x: observations"
  expect_equal(a$scale[rows], weightData()$w)
  expect_true(all(a$scale[!rows] == 1))
  expect_false("x_scale" %in% names(a$data))
  expect_null(a$Ntrials)
  # no weights(): no `scale` argument
  expect_null(inlaArgumentList(weightModel("x")$inla)$scale)
  # weights on both parts
  m2 <- jsm(list(y | weights(w) ~ t + scaled(x)), list(x | weights(1 / w) ~ t + (1 | id)), data = weightData(),
            fit = FALSE)
  st2 <- stackRows(m2$inla)
  a2 <- inlaArgumentList(m2$inla)
  expect_equal(a2$scale[st2$.block == "y: observations"], weightData()$w)
  expect_equal(a2$scale[st2$.block == "x: observations"], 1 / weightData()$w)
  expect_true(all(a2$scale[st2$.block == "x: zero equations for y"] == 1))
  expect_identical(inlaArguments(m2)$scale, a2$scale)
})

test_that("jsm(inlaOptions = ): `scale` is built by mvrem and cannot be set", {
  expect_error(weightModel(inlaOptions = list(scale = 2)), "must not set `scale`")   # checked up front
  m <- weightModel()
  m$options$inlaOptions <- list(scale = 2)
  expect_error(inlaArguments(m), "must not set `scale`")
})

test_that("namingScheme(): the weights template names the scale column", {
  expect_identical(namingScheme()$weights, "{resp}_scale")
  m <- weightModel(naming = namingScheme(weights = "W_{resp}"))
  expect_identical(m$inla$inputs$x, c(scale = "W_x"))
  expect_true("W_x" %in% names(m$inla$stacks$x$main))
  expect_equal(inlaArgumentList(m$inla)$scale[stackRows(m$inla)$.block == "x: observations"], weightData()$w)
  expect_error(namingScheme(weights = "{long}_s"), "Unknown placeholder.*`weights`.*\\{resp\\}")
  # a collision with another column is caught
  expect_error(weightModel(naming = namingScheme(weights = "{resp}_1")), "duplicated INLA names: x_1")
})

# ---- display ------------------------------------------------------------------------

test_that("formatInlaFormula(): an inputs comment line and the scale column in the data", {
  txt <- formatInlaFormula(weightModel(), nrows = 2)
  expect_true(any(grepl("   # inputs: scale = x_scale (rows of other likelihoods: 1)", txt, fixed = TRUE)))
  expect_true(any(grepl("x_scale", txt[grepl("x_x", txt)], fixed = TRUE)))
  expect_false(any(grepl("# inputs:", formatInlaFormula(weightModel("x")), fixed = TRUE)))
  expect_false(any(grepl("# inputs:", formatInlaFormula(weightModel(), nrows = 0), fixed = TRUE)))
  # still parses back to the formula
  expect_equal(deparse1(str2lang(paste(txt, collapse = "\n"))), deparse1(weightModel()$inla$formula))
})

test_that("latexEquation(): the residual precision of a weighted gaussian is w_i tau", {
  m <- weightModel()
  tex <- format(latexEquation(m))
  expect_match(tex, "\\varepsilon^{(x)}_{i} &\\sim \\mathcal{N}\\left(0, \\left(w_{i}\\, \\tau^{(x)}_{\\varepsilon}\\right)^{-1}\\right)",
               fixed = TRUE)
  # the outcome without weights keeps tau alone
  expect_match(tex, "\\varepsilon^{(y)}_{i} &\\sim \\mathcal{N}\\left(0, \\left(\\tau^{(y)}_{\\varepsilon}\\right)^{-1}\\right)",
               fixed = TRUE)
  expect_false(grepl("w_{i}", format(latexEquation(weightModel("x"))), fixed = TRUE))
  # an expression and the other families
  expect_match(format(latexEquation(weightModel("x | weights(1 / w)"))), "\\left(\\frac{1}{w_{i}}\\, \\tau^{(x)}", fixed = TRUE)
  ft <- format(latexEquation(weightModel(family = list(x = studentTFamily()))))
  expect_match(ft, "t_{\\nu^{(x)}}\\left(0, \\left(w_{i}\\, \\tau^{(x)}_{\\varepsilon}\\right)^{-1}\\right)", fixed = TRUE)
  d <- weightData()
  d$x <- exp(d$x) / (1 + exp(d$x))
  fb <- format(latexEquation(weightModel(d = d, family = list(x = betaFamily()))))
  expect_match(fb, "\\mathrm{Beta}\\left(\\mu^{(x)}_{i}\\, w_{i}\\, \\phi^{(x)}, \\left(1 - \\mu^{(x)}_{i}\\right) w_{i}\\, \\phi^{(x)}\\right)",
               fixed = TRUE)
  d$x <- exp(fxSmallData()$x)
  fg <- format(latexEquation(weightModel(d = d, family = list(x = gammaFamily()))))
  expect_match(fg, "\\mathrm{Gamma}\\left(w_{i}\\, \\phi^{(x)}, w_{i}\\, \\phi^{(x)} / \\mu^{(x)}_{i}\\right)", fixed = TRUE)
})

test_that("plotModel() / plot(): a weighted model draws", {
  m <- weightModel()
  expect_no_error(withNullDevice(plotModel(m)))
  expect_no_error(withNullDevice(plot(m, interactive = FALSE)))
})

# ---- fit: simulation ----------------------------------------------------------------

# A gaussian marker measured by two devices whose residual SDs differ by `ratio`
# (device B: ratio * sdA); w = the inverse variance ratio (A: 1, B: 1 / ratio^2), so the
# residual precision at w = 1 is 1 / sdA^2. One outcome per subject at t = 0.
simWeights <- function(seed = 13, n = 100, sdA = 0.2, ratio = 3) {
  set.seed(seed)
  ids <- sprintf("S%03d", seq_len(n))
  long <- data.frame(id = rep(ids, each = 5), t = rep(0:4 / 2, n), stringsAsFactors = FALSE)
  i <- match(long$id, ids)
  b0 <- rnorm(n, sd = 0.5)
  long$device <- factor(sample(c("A", "B"), nrow(long), replace = TRUE))
  long$x <- 1 + 0.5 * long$t + b0[i] + rnorm(nrow(long), sd = ifelse(long$device == "A", sdA, sdA * ratio))
  long$w <- ifelse(long$device == "A", 1, 1 / ratio^2)
  out <- data.frame(id = ids, t = 0, stringsAsFactors = FALSE)
  out$y <- 0.5 + 0.8 * (1 + b0) + rnorm(n, sd = 0.3)
  list(long = long, out = out)
}

weightFits <- function() {
  if (!is.null(.fx$weightFits)) return(.fx$weightFits)
  skip_if_not_installed("INLA")
  s <- simWeights()
  fit <- function(lhs) suppressMessages(jsm(list(y ~ 1 + scaled(x)), list(stats::as.formula(paste(lhs, "~ t + (1 | id)"))),
                                            dataOut = s$out, dataLong = s$long))
  .fx$weightFits <- list(w = fit("x | weights(w)"), u = fit("x"), sim = s)
  .fx$weightFits
}

test_that("fit: weights() recover the residual precision and shrink the fixed-effect SDs", {
  f <- weightFits()
  skip_if(is.null(inlaFit(f$w)) || is.null(inlaFit(f$u)), "INLA fit failed")
  hp <- function(m) { h <- hyperparameters(m); h[h$response == "x" & h$parameter == "residual precision", ] }
  hw <- hp(f$w); hu <- hp(f$u)
  # true precision at w = 1: 1 / 0.2^2 = 25 (weighted: 26.2 [22.7, 30.0]); unweighted: a
  # compromise of 25 and 25 / 9 (5.2 [4.5, 6.0])
  expect_true(hw$q0.025 < 25 && 25 < hw$q0.975)
  expect_true(hu$q0.975 < 10)
  fe <- function(m) { x <- fixedEffects(m); x[x$response == "x", ] }
  fw <- fe(f$w); fu <- fe(f$u)
  # SDs: (Intercept) 0.054 vs 0.059, t 0.018 vs 0.028
  expect_true(all(fw$sd < fu$sd))
  expect_lt(fw$sd[fw$term == "t"], 0.8 * fu$sd[fu$term == "t"])
  expect_true(fw$q0.025[fw$term == "t"] < 0.5 && 0.5 < fw$q0.975[fw$term == "t"])
  # the association is unaffected in substance
  expect_true(abs(associations(f$w)$mean - 0.8) < 0.15)
  # tables built from the fit work as usual
  et <- estimatesTable(f$w)
  expect_true(any(et$response == "x" & et$type == "residual"))
  expect_equal(et$mean[et$response == "x" & et$type == "residual"], 1 / sqrt(25), tolerance = 0.1)
})

test_that("fit: fitted values and residuals keep the response scale; the weight is a data column", {
  f <- weightFits()
  skip_if(is.null(inlaFit(f$w)), "INLA fit failed")
  fv <- fittedValues(f$w, "x")
  expect_equal(fv$observed, f$sim$long$x)
  expect_equal(fv$w, f$sim$long$w)
  r <- residuals(f$w, "x")
  expect_equal(r$residual, r$observed - r$mean)
  # tooltip annotation of the weight variable (an aux variable: key cov|w)
  tex <- format(latexEquation(f$w, annotate = "tooltip"))
  expect_match(tex, "variable w", fixed = TRUE)
  expect_no_error(withNullDevice(plotModel(f$w)))
})

test_that("predict(): the weight variable is optional in newdata (set to 1)", {
  f <- weightFits()
  skip_if(is.null(inlaFit(f$w)), "INLA fit failed")
  nd <- data.frame(id = c("S001", "NEW"), t = c(2.5, 0), stringsAsFactors = FALSE)
  p <- suppressMessages(predict(f$w, nd, response = "x"))
  expect_equal(nrow(p), 2L)
  expect_false("w" %in% names(p))
  expect_true(all(is.finite(p$mean)))
  refit <- attr(p, "model")
  newRows <- utils::tail(refit$data$long$x, 2)
  expect_equal(newRows$w, c(1, 1))
  # the inla() scale of the new rows: 1
  a <- inlaArguments(refit)
  st <- stackRows(refit$inla)
  obs <- which(st$.block == "x: observations")
  expect_equal(a$scale[utils::tail(obs, 2)], c(1, 1))
})

# ---- plots: point area ------------------------------------------------------------

test_that("weightSizes(): point area proportional to the weight, relative to the median, clamped", {
  expect_equal(weightSizes(c(1, 4, 1)), c(1, 2, 1))
  expect_equal(weightSizes(c(1, 1, 1 / 9)), c(1, 1, 1 / 3))
  expect_equal(weightSizes(c(1, 1, 1e-6, 1e6)), c(1, 1, 0.3, 3))
  expect_null(weightSizes(NULL))
})

test_that("plot() / plotResiduals() / predictions: weighted points carry their size and a note", {
  m <- weightModel()
  d <- linesData(m)
  px <- d$panels[[which(vapply(d$panels, `[[`, "", "response") == "x")]]
  py <- d$panels[[which(vapply(d$panels, `[[`, "", "response") == "y")]]
  expect_equal(px$note, "point area proportional to w")
  expect_null(py$note)
  wd <- weightData()
  for (s in px$subjects) {
    k <- which(wd$id == s$id)
    expect_equal(s$size, weightSizes(wd$w)[k][order(wd$t[k])])
  }
  expect_true(all(vapply(py$subjects, function(s) is.null(s$size), NA)))
  expect_no_error(withNullDevice(plot(m, interactive = FALSE)))
  skipNoWidgets()
  w <- plot(m, interactive = TRUE)
  wx <- w$x$panels[[which(vapply(d$panels, `[[`, "", "response") == "x")]]
  expect_equal(wx$note, "point area proportional to w")
  expect_true(all(vapply(Filter(function(l) l$kind == "obs", wx$lines), function(l) length(l$size) == length(l$x), NA)))
  if (hasNode()) {
    res <- runWidgetJs(list(w = w))$w
    expect_true(res$ok, label = paste(res$error, res$stack))
  }
  # fitted: residuals and the observed points behind predictions
  f <- weightFits()
  skip_if(is.null(inlaFit(f$w)), "INLA fit failed")
  rd <- plotResidualsData(f$w, "x", "time", NULL, NULL, 30, NULL)
  expect_equal(rd$panels[[1]]$note, "point area proportional to w")
  expect_true(all(vapply(rd$panels[[1]]$subjects, function(s) length(s$size) == length(s$x), NA)))
  expect_no_error(withNullDevice(plotResiduals(f$w, interactive = FALSE)))
  wr <- plotResiduals(f$w, "x", interactive = TRUE)
  expect_true(any(vapply(wr$x$panels[[1]]$lines, function(l) !is.null(l$size), NA)))
})
