# Extended tests of fitted models: jsm() arguments, the accessors of an
# mvremModel, refit() / update(), the result tables (fixed effects, hyperparameters,
# associations, random effects, variance components, fitted values, residuals,
# criteria, estimates table) and predict(). Fixtures: helper-fixtures.R; helpers:
# helper-fit.R.
#
# INLA runs of this file beyond the fixtures (4 fits + the cached predict()):
# refit() with a failing option (fails in INLA's argument checks), refit() with dic
# only, update() with criteria given through inlaOptions. Everything else reads the
# fixtures or mocks runInla() (mockInlaRun()).

# ---------------------------------------------------------------------------
# jsm(): arguments and the unfitted model

test_that("jsm(fit = FALSE): all components, empty options, no fit", {
  m <- fxModel("multi", fit = FALSE)
  expect_s3_class(m, "jsm")
  expect_s3_class(m, "mvremModel")
  expect_named(m, c("call", "forms", "data", "dataSummary", "inla", "time", "options", "fit", "inlaError"))
  expect_null(m$fit)
  expect_null(m$inlaError)
  expect_null(inlaFit(m))
  # the helper builds unfitted fixtures with criteria = FALSE, so inlaOptions is empty
  expect_equal(m$options, list(verbose = FALSE, threads = NULL, inlaOptions = list()))
})

test_that("jsm(): the construction accessors return the components", {
  m <- fxModel("multi", fit = FALSE)
  expect_identical(modelFormulas(m), m$forms)
  expect_named(modelFormulas(m), c("out", "long", "variables", "order"))
  expect_identical(modelData(m), m$data)
  expect_named(modelData(m)$out, c("k", "cnt"))
  expect_named(modelData(m)$long, c("x1", "x2"))
  expect_equal(vapply(modelData(m)$out, nrow, 0L), c(k = fxCounts()$rows, cnt = fxCounts()$rows))
  expect_identical(dataSummary(m), m$dataSummary)
  expect_s3_class(dataSummary(m), "mvremDataSummary")
  expect_identical(inlaSpec(m), m$inla)
  expect_identical(inlaFormula(m), m$inla$formula)
  expect_s3_class(inlaFormula(m), "formula")
})

test_that("jsm(fit = FALSE) equals the step-by-step pipeline", {
  m <- fxModel("beta", fit = FALSE)
  f <- fxForms("beta")
  inl <- buildInlaFormula(f, prepareData(f, fxSpec("beta")$data))
  expect_equal(deparse(inlaFormula(m)), deparse(inl$formula))
  expect_equal(unclass(formatInlaFormula(m)), unclass(formatInlaFormula(inl)))
  expect_equal(inlaData(m), inlaData(inl))
})

test_that("print() of an unfitted model lists equations, links, observations and the state", {
  m <- fxModel("multi", fit = FALSE)
  out <- capture.output(res <- print(m))
  expect_identical(res, m)
  expect_match(out, "k \\(binomial\\), cnt \\(poisson\\)", all = FALSE)
  expect_match(out, "x1 \\(gaussian\\), x2 \\(poisson\\)", all = FALSE)
  expect_match(out, "k <- x1, x2; cnt <- x2", all = FALSE)
  expect_match(out, paste0("x1: ", fxCounts()$rows, ", x2: ", fxCounts()$rows), all = FALSE)
  expect_match(out, "Not fitted \\(use refit", all = FALSE)
  expect_false(any(grepl("INLA failed", out)))
})

test_that("every result accessor of an unfitted model asks for a fit", {
  m <- fxModel("multi", fit = FALSE)
  fns <- list(fixedEffects = fixedEffects, hyperparameters = hyperparameters,
              associations = associations, randomEffects = randomEffects,
              varianceComponents = varianceComponents, fittedValues = fittedValues,
              residuals = residuals, modelFit = modelFit, estimatesTable = estimatesTable,
              summary = summary, coef = coef)
  for (f in names(fns)) expect_error(fns[[f]](m), "not fitted.*refit\\(\\)", info = f)
  expect_error(predict(m, data.frame(id = "S01", t = 1)), "not fitted")
  # the construction accessors work
  expect_s3_class(inlaFormula(m), "formula")
  expect_s3_class(inlaData(m), "data.frame")
  expect_type(inlaArguments(m), "list")
})

test_that("jsm(): dataOut / dataLong and familyOut follow prepareData() / collectFormulas()", {
  d <- fxSmallData()
  m <- jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), dataOut = d, dataLong = d,
           familyOut = c(y = "t"), fit = FALSE)
  expect_equal(modelFormulas(m)$out$y$family$family, "t")
  expect_equal(modelFormulas(m)$long$x$family$family, "gaussian")
  expect_error(jsm(list(y ~ t + scaled(x)), list(x ~ t), dataOut = d, fit = FALSE),
               "No data for the longitudinal")
  expect_error(jsm(list(y ~ t), list(x ~ t), data = d, fit = FALSE), "No scaled\\(\\) term")
})

test_that("jsm() options: threads, verbose and criteria reach the INLA arguments", {
  s <- fxSpec("beta")
  m <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE,
           threads = 2, verbose = TRUE, criteria = TRUE,
           inlaOptions = list(control.compute = list(dic = FALSE, config = TRUE)))
  a <- inlaArguments(m)
  expect_equal(a$num.threads, 2)
  expect_true(a$verbose)
  # criteria = TRUE is merged UNDER the user's control.compute
  expect_equal(a$control.compute, list(dic = FALSE, waic = TRUE, cpo = TRUE, config = TRUE))
  expect_equal(m$options$inlaOptions$control.compute$waic, TRUE)
  m2 <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE, criteria = TRUE)
  expect_equal(inlaArguments(m2)$control.compute, list(dic = TRUE, waic = TRUE, cpo = TRUE))
  m3 <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE,
            inlaOptions = list(control.compute = list(dic = TRUE)))
  expect_equal(inlaArguments(m3)$control.compute, list(dic = TRUE))
})

test_that("jsm(inlaOptions = ): a named list; the arguments mvrem builds are refused", {
  s <- fxSpec("beta")
  build <- function(opt) jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE,
                             inlaOptions = opt)
  expect_error(inlaArguments(build(list(1))), "must be a named list")
  for (nm in c("formula", "family", "data", "control.family", "Ntrials", "E"))
    expect_error(inlaArguments(build(stats::setNames(list(1), nm))), paste0("must not set `", nm, "`"), info = nm)
  expect_error(inlaArguments(build(list(Ntrials = 1, E = 2))), "must not set `Ntrials`, `E`")
  skip_if_not_installed("INLA")
  # with fit = TRUE the error comes before INLA runs (an error, not the "INLA failed" warning)
  expect_error(jsm(s$formOut, s$formLong, data = s$data, inlaOptions = list(formula = y ~ 1)),
               "must not set `formula`")
})

test_that("jsm(na.action = ): a function or its name; dropped rows are reported", {
  s <- fxSpec("beta")
  d <- s$data
  d$t[c(3, 10)] <- NA
  expect_error(jsm(s$formOut, s$formLong, data = d, fit = FALSE), "Missing values in `t`")
  m1 <- jsm(s$formOut, s$formLong, data = d, family = s$family, na.action = stats::na.omit, fit = FALSE)
  m2 <- jsm(s$formOut, s$formLong, data = d, family = s$family, na.action = "na.omit", fit = FALSE)
  expect_equal(nrow(modelData(m1)$out$y), nrow(d) - 2)
  expect_equal(nrow(modelData(m1)$long$x1), nrow(d) - 2)
  expect_equal(modelData(m1), modelData(m2))
  expect_match(capture.output(print(m1)), "rows dropped for missing values", all = FALSE)
  expect_match(capture.output(print(m1)), "y: 2, x1: 2", all = FALSE)
  expect_error(jsm(s$formOut, s$formLong, data = d, na.action = stats::na.fail, fit = FALSE))
})

test_that("jsm(): rows with a missing response are kept and get the link of their likelihood", {
  s <- fxSpec("beta")
  d <- s$data
  d$y[c(2, 7, 30)] <- NA
  d$x1[c(3, 11)] <- NA
  m <- jsm(s$formOut, s$formLong, data = d, family = s$family, fit = FALSE)
  expect_equal(nrow(modelData(m)$out$y), nrow(d))
  eq <- dataSummary(m)$equations
  expect_equal(eq$complete[eq$response == "y"], nrow(d) - 3)
  expect_equal(eq$complete[eq$response == "x1"], nrow(d) - 2)
  expect_match(capture.output(print(m)), "rows with missing values, kept", all = FALSE)
  expect_match(capture.output(print(m)), "y: 3, x1: 2", all = FALSE)      # in the order passed
  # control.predictor$link gives every row the likelihood of its block, so INLA's fitted
  # value of a missing y is on the mean scale (checked on a fit in the predict() tests)
  a <- inlaArguments(m)
  st <- inlaData(m)
  expect_equal(unique(a$control.predictor$link[st$.block == "y: observations"]), 1L)
  expect_equal(unique(a$control.predictor$link[st$.block == "x1: observations"]), 2L)
  expect_true(all(is.na(a$data$MX[which(st$.block == "y: observations")[c(2, 7, 30)], "y_y"])))
})

test_that("jsm(naming = ) and options(mvrem.naming) rename the INLA objects", {
  s <- fxSpec("beta")
  nm <- namingScheme(responseObject = "YY", reIndex = "re_{resp}_{term}_{group}",
                     scaledIndex = "sc_{out}_{long}")
  m <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, naming = nm, fit = FALSE)
  a <- inlaArguments(m)
  expect_true("YY" %in% names(a$data))
  expect_false("MX" %in% names(a$data))
  ftxt <- paste(deparse(inlaFormula(m)), collapse = " ")
  expect_match(ftxt, "re_x1_1_id", fixed = TRUE)
  expect_match(ftxt, "sc_y_x1", fixed = TRUE)
  withr::local_options(mvrem.naming = nm)
  m2 <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE)
  expect_equal(deparse(inlaFormula(m2)), deparse(inlaFormula(m)))
})

test_that("jsm(naming = ): a custom scheme builds the same model structure", {
  ref <- fxModel("beta", fit = FALSE)
  s <- fxSpec("beta")
  nm <- namingScheme(responseObject = "YY", reIndex = "re_{resp}_{term}_{group}",
                     reLevels = "lv_{resp}_{terms}_{group}", zeroIndex = "z_{long}_{out}",
                     scaledIndex = "sc_{out}_{long}", fixed = "b_{resp}_{term}")
  m <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, naming = nm, fit = FALSE)
  expect_equal(unname(inlaSpec(m)$families), unname(inlaSpec(ref)$families))
  expect_equal(names(inlaSpec(m)$families), c("y_y", "x1_x1", "x1_x1_0eq"))   # response templates kept
  expect_equal(dim(inlaData(m)), dim(inlaData(ref)))
  expect_equal(as.vector(table(inlaData(m)$.block)), as.vector(table(inlaData(ref)$.block)))
  expect_true(all(c("b_y_t", "re_x1_1_id", "z_x1_y", "sc_y_x1") %in% names(inlaData(m))))
  expect_false(any(c("y_t", "x1_1_id") %in% names(inlaData(m))))
  expect_equal(inlaSpec(m)$fixed$y$columns, inlaSpec(ref)$fixed$y$columns, ignore_attr = "names")
  expect_equal(unname(inlaSpec(m)$levels$lv_x1_1_t_id), unname(inlaSpec(ref)$levels$x1_1_t_id))
})

test_that("effectLabel(): the intercept gets naming$interceptLabel, other effects their column", {
  expect_equal(effectLabel(c("(Intercept)", "t", "grpb"), namingScheme()), c("(Intercept)", "t", "grpb"))
  expect_equal(effectLabel(c("(Intercept)", "t"), namingScheme(interceptLabel = "b0")), c("b0", "t"))
  expect_length(effectLabel(character(), namingScheme()), 0L)
  expect_equal(namingScheme()$interceptLabel, "(Intercept)")
})

test_that("jsm(naming = namingScheme(interceptLabel = )) is stored and renames no INLA object", {
  ref <- fxModel("beta", fit = FALSE)
  s <- fxSpec("beta")
  m <- jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE,
           naming = namingScheme(interceptLabel = "b0"))
  expect_equal(inlaSpec(m)$naming$interceptLabel, "b0")
  expect_equal(inlaSpec(ref)$naming$interceptLabel, "(Intercept)")
  expect_equal(deparse(inlaFormula(m)), deparse(inlaFormula(ref)))
  expect_equal(names(inlaData(m)), names(inlaData(ref)))
  expect_equal(unname(inlaSpec(m)$fixed$y$columns)[1], "(Intercept)")   # the design column keeps R's name
})

test_that("jsm() without INLA gives a clear error unless fit = FALSE", {
  local_mocked_bindings(hasInla = function() FALSE)
  s <- fxSpec("beta")
  expect_s3_class(jsm(s$formOut, s$formLong, data = s$data, family = s$family, fit = FALSE), "jsm")
  expect_error(jsm(s$formOut, s$formLong, data = s$data, family = s$family), "INLA package is required")
  expect_error(refit(fxModel("beta", fit = FALSE)), "INLA package is required")
})

# ---------------------------------------------------------------------------
# accessors of the INLA side: inlaData(), inlaArguments()

test_that("inlaData() stacks all blocks in likelihood order; n limits each part", {
  m <- fxModel("multi", fit = FALSE)
  st <- inlaData(m)
  expect_equal(levels(st$.block),
               c("k: observations", "cnt: observations", "x1: observations",
                 "x1: zero equations for k", "x2: observations", "x2: zero equations for k",
                 "x2: zero equations for cnt"))
  expect_equal(nrow(st), sum(vapply(unlist(inlaSpec(m)$stacks, recursive = FALSE), nrow, 0L)))
  expect_false(is.unsorted(as.integer(st$.block)))
  expect_true(all(c("k_k", "cnt_cnt", "x1_x1", "x1_x1_0eq", "x2_x2_0eq", ".block") %in% names(st)))
  expect_equal(as.vector(table(inlaData(m, n = 2)$.block)), rep(2L, 7))
  # the default method takes a buildInlaFormula() result
  expect_equal(inlaData(inlaSpec(m), n = 3), inlaData(m, n = 3))
  expect_error(inlaData(list(a = 1)), "must be a model or the result of buildInlaFormula")
  expect_error(inlaData(1:3), "must be a model or the result of buildInlaFormula")
})

test_that("inlaData(n = 0) gives an empty stack with the columns and block levels of the full one", {
  m <- fxModel("multi", fit = FALSE)
  st0 <- inlaData(m, n = 0)
  st <- inlaData(m)
  expect_s3_class(st0, "data.frame")
  expect_equal(nrow(st0), 0L)
  expect_equal(names(st0), names(st))
  expect_s3_class(st0$.block, "factor")
  expect_equal(levels(st0$.block), levels(st$.block))
  expect_equal(nrow(inlaData(inlaSpec(m), n = 0)), 0L)
})

test_that("inlaArguments(): one likelihood per response and zero equation, aligned rows", {
  m <- fxModel("multi", fit = FALSE)
  a <- inlaArguments(m)
  expect_named(a, c("formula", "family", "data", "control.family", "control.predictor", "Ntrials", "E",
                    "verbose"), ignore.order = TRUE)
  expect_equal(a$family, unname(inlaSpec(m)$families))
  expect_equal(a$family, c("binomial", "poisson", "gaussian", "poisson", "gaussian", "gaussian"))
  expect_length(a$control.family, length(a$family))
  expect_equal(colnames(a$data$MX), names(inlaSpec(m)$families))
  n <- nrow(inlaData(m))
  expect_equal(nrow(a$data$MX), n)
  expect_length(a$control.predictor$link, n)
  expect_equal(a$control.predictor$compute, TRUE)
  expect_true(all(a$control.predictor$link %in% seq_along(a$family)))
  expect_false(a$verbose)
  expect_null(a$num.threads)
  expect_null(a$control.compute)
  # each row has exactly one response column filled (zero-equation rows: the 0 of `_0eq`)
  expect_true(all(rowSums(!is.na(a$data$MX)) == 1L))
  # the input columns are not passed as data, nor .block
  expect_false(any(c("k_Ntrials", "cnt_E", ".block") %in% names(a$data)))
})

test_that("inlaArguments(): Ntrials / E hold the input of the row's likelihood, 1 elsewhere", {
  m <- fxModel("multi", fit = FALSE)
  a <- inlaArguments(m)
  d <- fxData()$long
  st <- inlaData(m)
  kRows <- st$.block == "k: observations"
  cntRows <- st$.block == "cnt: observations"
  expect_length(a$Ntrials, nrow(st))
  expect_length(a$E, nrow(st))
  expect_false(anyNA(a$Ntrials))
  expect_false(anyNA(a$E))
  expect_equal(a$Ntrials[kRows], d$n)
  expect_equal(a$E[cntRows], d$pt)
  expect_true(all(a$Ntrials[!kRows] == 1))
  expect_true(all(a$E[!cntRows] == 1))
})

test_that("inlaArguments() of a survival model: the response object is a list with an inla.surv", {
  skip_if_not_installed("INLA")
  a <- inlaArguments(fxModel("surv", fit = FALSE))
  expect_type(a$data$MX, "list")
  expect_s3_class(a$data$MX$time_time, "inla.surv")
  expect_equal(a$family[1], "weibullsurv")
  expect_null(a$Ntrials)
  expect_null(a$E)
})

test_that("inlaArguments() are the arguments INLA received for the fit", {
  m <- fxModel("multi")
  a <- inlaArguments(m)
  received <- inlaFit(m)$.args                      # inla() stores its call's arguments
  expect_equal(received$family, a$family)
  expect_equal(received$data$MX, a$data$MX)
  expect_equal(received$Ntrials, a$Ntrials)
  expect_equal(received$E, a$E)
  expect_equal(received$control.predictor$link, a$control.predictor$link)
  expect_equal(deparse(received$formula), deparse(a$formula))
  # the fixture is fitted with criteria = TRUE (helper): merged into control.compute
  expect_equal(a$control.compute, list(dic = TRUE, waic = TRUE, cpo = TRUE))
  expect_equal(received$control.compute[c("dic", "waic", "cpo")], a$control.compute)
  expect_equal(m$options$inlaOptions, list(control.compute = a$control.compute))
})

test_that("print() of a fitted model says so", {
  m <- fxModel("beta")
  expect_s3_class(inlaFit(m), "inla")
  out <- capture.output(print(m))
  expect_match(out, "Fitted with INLA", all = FALSE)
  expect_false(any(grepl("Not fitted", out)))
})

# ---------------------------------------------------------------------------
# refit(), update(), INLA failure

test_that("refit(): a failure keeps the options and the error; the model stays usable", {
  skip_if_not_installed("INLA")
  m <- fxModel("beta")
  bad <- list(control.inla = list(strategy = "no-such-strategy"))   # rejected by INLA's argument checks
  r <- NULL
  suppressMessages(expect_warning(r <- refit(m, inlaOptions = bad), "^INLA failed: "))
  expect_null(inlaFit(r))
  expect_s3_class(r$inlaError, "error")
  expect_equal(r$options$inlaOptions, bad)
  expect_error(fixedEffects(r), "not fitted \\(INLA failed: ")
  expect_error(summary(r), "not fitted")
  expect_match(capture.output(print(r)), "Not fitted: INLA failed", all = FALSE)
  expect_s3_class(inlaFormula(r), "formula")
  expect_equal(modelData(r), modelData(m))
})

test_that("refit() with other options fits again and replaces the options; dic only", {
  skip_if_not_installed("INLA")
  m <- fxModel("beta")
  r2 <- suppressMessages(refit(m, inlaOptions = list(control.compute = list(dic = TRUE)), threads = 1))
  skip_if(is.null(inlaFit(r2)), "INLA fit failed")
  expect_null(r2$inlaError)
  expect_equal(r2$options, list(verbose = FALSE, threads = 1,
                                inlaOptions = list(control.compute = list(dic = TRUE))))
  expect_false(is.null(inlaFit(r2)$dic))
  expect_null(inlaFit(r2)$waic)
  # criteria incomplete: modelFit() explains how to get them
  expect_error(modelFit(r2), "without model criteria.*criteria = TRUE")
  expect_equal(fixedEffects(r2)$mean, fixedEffects(m)$mean, tolerance = 1e-3)
  expect_equal(modelFit(m)$dic[3], inlaFit(r2)$dic$dic, tolerance = 1e-3)
})

test_that("refit(): unknown arguments are an error", {
  m <- fxModel("beta", fit = FALSE)
  expect_error(refit(m, foo = 1, bar = 2), "Unknown argument\\(s\\) to refit\\(\\): foo, bar")
})

test_that("update() repeats the call with changed formulas, family or options", {
  m <- fxModel("beta", fit = FALSE)
  m2 <- update(m, formOut = list(y ~ 1 + scaled(x1)), formLong = list(x1 ~ t + (1 | id)))
  expect_s3_class(m2, "jsm")
  expect_null(inlaFit(m2))
  expect_equal(unname(inlaSpec(m2)$fixed$y$columns), "(Intercept)")
  expect_length(inlaSpec(m2)$random$x1, 1L)
  expect_equal(inlaSpec(m2)$random$x1[[1]]$effects, "(Intercept)")
  m3 <- update(m, family = list(y = "beta", x1 = studentTFamily()))
  expect_equal(modelFormulas(m3)$long$x1$family$family, "t")
  expect_equal(modelFormulas(m)$long$x1$family$family, "gaussian")     # the original is unchanged
  m4 <- update(m, threads = 3)
  expect_equal(m4$options$threads, 3)
})

test_that("update(evaluate = FALSE) gives the changed call; a model without call is an error", {
  m <- fxModel("beta", fit = FALSE)
  cl <- update(m, threads = 3, evaluate = FALSE)
  expect_true(is.call(cl))
  expect_equal(cl$threads, 3)
  expect_false(cl$fit)                                      # the rest of the call is kept
  expect_true(update(m, fit = TRUE, evaluate = FALSE)$fit)
  m0 <- m
  m0$call <- NULL
  expect_error(update(m0, fit = FALSE), "no stored call")
})

# ---------------------------------------------------------------------------
# fixed effects, coef(), summary()

test_that("fixedEffects(): one row per design column, factor levels expanded, INLA's rows", {
  m <- fxModel("multi")
  fe <- fixedEffects(m)
  expect_named(fe, c("response", "term", statCols, "mode"))
  expect_false("kld" %in% names(fe))
  expect_equal(paste(fe$response, fe$term),
               c("k (Intercept)", "k t", "cnt (Intercept)", "cnt grpb", "cnt grpc",
                 "x1 (Intercept)", "x1 t", "x2 (Intercept)", "x2 t"))
  expectOrderedStats(fe)
  inlaNames <- unlist(lapply(inlaSpec(m)$fixed, function(f) names(f$columns)), use.names = FALSE)
  expect_equal(inlaNames[4], "cnt_grpb")
  expect_equal(fe$mean, unname(inlaFit(m)$summary.fixed[inlaNames, "mean"]))
  expect_equal(fe$q0.975, unname(inlaFit(m)$summary.fixed[inlaNames, "0.975quant"]))
})

test_that("coef(): the fixed-effect means named response: term", {
  m <- fxModel("multi")
  fe <- fixedEffects(m)
  expect_equal(unname(coef(m)), fe$mean)
  expect_equal(names(coef(m)), paste0(fe$response, ": ", fe$term))
})

test_that("summary() collects the tables; print(digits =) rounds them", {
  m <- fxModel("multi")
  s <- summary(m)
  expect_s3_class(s, "summary.mvremModel")
  expect_named(s, c("call", "model", "fixed", "hyper", "associations"))
  out <- capture.output(res <- print(s, digits = 2))
  expect_identical(res, s)
  expect_match(out, "^Fixed effects", all = FALSE)
  expect_match(out, "^Associations", all = FALSE)
  expect_match(out, "^Hyperparameters", all = FALSE)
  g <- format(signif(associations(m)$mean[1], 2), scientific = FALSE, drop0trailing = TRUE)
  expect_match(out, paste0("^ +k +x1 +", gsub(".", "\\.", g, fixed = TRUE), " "), all = FALSE)
  hyp <- out[seq(grep("^Hyperparameters", out), length(out))]
  expect_false(any(grepl("gamma", hyp)))                   # associations are listed separately
  expect_false(any(grepl("mode", out)))
})

# ---------------------------------------------------------------------------
# hyperparameters(), associations()

test_that("hyperparameters() of 'beta': grouped by response (order passed), INLA's order within, all recognised", {
  m <- fxModel("beta")
  h <- hyperparameters(m)
  expect_named(h, c("response", "parameter", "type", statCols, "mode"))
  expect_equal(nrow(h), nrow(inlaFit(m)$summary.hyperpar))
  expect_equal(h$parameter,
               c("phi", "gamma (x1 -> y)", "residual precision", "precision of (Intercept) | id",
                 "precision of t | id", "correlation of (Intercept), t | id"))
  expect_equal(h$response, c("y", "y", "x1", "x1", "x1", "x1"))
  expect_equal(h$type, c("family", "association", "family", "random", "random", "random"))
  expect_true(all(h$mean[grepl("precision", h$parameter)] > 0))
  expect_true(h$q0.025[6] >= -1 && h$q0.975[6] <= 1)
  # the internal hyperTable() keeps INLA's order: likelihood rows, then the latent terms
  ht <- hyperTable(m)
  expect_equal(ht$.inla, rownames(inlaFit(m)$summary.hyperpar))
  expect_equal(ht$type, c("family", "family", "random", "random", "random", "association"))
  expectOrderedStats(h)
})

test_that("hyperparameters() of 'multi': labels per response; binomial / poisson have no family row", {
  m <- fxModel("multi")
  h <- hyperparameters(m)
  expect_equal(nrow(h), nrow(inlaFit(m)$summary.hyperpar))
  expect_false(anyNA(h$type))
  expect_setequal(h$parameter,
                  c("residual precision", "precision of (Intercept) | id", "precision of t | id",
                    "correlation of (Intercept), t | id", "precision of (Intercept) | id",
                    "gamma (x1 -> k)", "gamma (x2 -> k)", "gamma (x2 -> cnt)"))
  expect_equal(h$response, rep(c("k", "cnt", "x1", "x2"), c(2, 1, 4, 1)))   # grouped, order passed
  expect_equal(h$type[4], "family")                       # x1 gaussian: first of x1
  expect_equal(h$response[4], "x1")
  expect_equal(sum(h$type == "family"), 1L)
  expect_equal(h$response[h$parameter == "precision of (Intercept) | id"], c("x1", "x2"))
  expect_equal(h$response[h$type == "association"], c("k", "k", "cnt"))
})

test_that("hyperparameters() of 'priors': every family and block hyperparameter is labelled", {
  m <- fxModel("priors")
  h <- hyperparameters(m)
  expect_equal(nrow(h), nrow(inlaFit(m)$summary.hyperpar))
  expect_false(anyNA(h$type))
  expect_equal(paste(h$response, h$parameter)[h$type == "family"],
               c("yg residual precision", "x1 residual precision", "x2 size", "x3 precision", "x3 dof"))
  expect_equal(paste(h$response, h$parameter)[h$type == "random"],
               c("yg precision of (Intercept) | site",
                 "x1 precision of (Intercept) | id", "x1 precision of t | id", "x1 correlation of (Intercept), t | id",
                 "x2 precision of (Intercept) | id", "x2 precision of t | id",
                 "x3 precision of (Intercept) | id", "x3 precision of t | id"))
  expect_equal(h$parameter[h$type == "association"],
               c("gamma (x1 -> yg)", "gamma (x2 -> yg)", "gamma (x3 -> yg)"))
  expect_true(all(h$mean[grepl("precision|size|dof", h$parameter)] > 0))
  expect_equal(fxWarnings("priors"), character())        # the fit raised no warning (e.g. of the data prior)
})

test_that("hyperparameters() of a survival outcome: alpha, no cure-model betas", {
  hs <- hyperparameters(fxModel("surv"))
  expect_equal(hs$parameter[hs$type == "family"], c("alpha", "residual precision"))
  expect_equal(hs$response[hs$type == "family"], c("time", "x1"))
  expect_false(any(grepl("^beta[0-9]", hs$parameter)))
  expect_false(anyNA(hs$type))
  expect_gt(hs$mean[hs$parameter == "alpha"], 0)
})

test_that("hyperTable() carries link / group / effect columns that the accessors drop", {
  m <- fxModel("multi")
  ht <- hyperTable(m)
  expect_true(all(hyperInternalColumns %in% names(ht)))
  expect_equal(hyperInternalColumns, c("link", "group", "effect", "effect2", ".inla"))
  expect_equal(ht$link[ht$type == "association"], c("x1", "x2", "x2"))
  expect_true(all(is.na(ht$link[ht$type != "association"])))
  expect_equal(ht$group[ht$type == "random"], rep("id", 4))
  expect_equal(ht$effect[ht$type == "random"], c("(Intercept)", "t", "(Intercept)", "(Intercept)"))
  expect_equal(ht$effect2[ht$type == "random"], c(NA, NA, "t", NA))
  expect_true(all(is.na(ht[ht$type == "family", c("group", "effect", "effect2")])))
  expect_false(any(hyperInternalColumns %in% names(hyperparameters(m))))
  expect_false(any(hyperInternalColumns %in% names(associations(m))))
})

test_that("naming$interceptLabel renames the intercept in hyperparameters() and estimatesTable()", {
  m <- fxModel("beta")
  mb <- m                                   # a copy: the label is read from the model's naming
  mb$inla$naming$interceptLabel <- "b0"
  h <- hyperparameters(mb)
  expect_equal(h$parameter[h$type == "random"],
               c("precision of b0 | id", "precision of t | id", "correlation of b0, t | id"))
  expect_equal(h[statCols], hyperparameters(m)[statCols])
  et <- estimatesTable(mb)
  expect_equal(et$parameter[et$type %in% c("fixed", "sd", "correlation")],
               c("b0", "t", "b0", "t", "SD b0 | id", "SD t | id", "corr b0, t | id"))
  expect_false(any(grepl("(Intercept)", et$parameter, fixed = TRUE)))
  expect_equal(fixedEffects(mb)$term[1], "(Intercept)")             # the term keeps R's name
  expect_equal(hyperparameters(m)$parameter[4], "precision of (Intercept) | id")
})

test_that("associations() equal the association rows of hyperparameters()", {
  for (nm in c("beta", "multi", "surv", "priors")) {
    m <- fxModel(nm)
    a <- associations(m)
    h <- hyperparameters(m)
    ha <- h[h$type == "association", ]
    expect_named(a, c("outcome", "longitudinal", statCols, "mode"))
    expect_equal(a$outcome, ha$response, info = nm)
    expect_equal(paste0("gamma (", a$longitudinal, " -> ", a$outcome, ")"), ha$parameter, info = nm)
    expect_equal(a[statCols], `rownames<-`(ha[statCols], NULL), info = nm)
    expectOrderedStats(a)
  }
})

test_that("associations(): one row per scaled() term, in formula order", {
  a <- associations(fxModel("multi"))
  expect_equal(paste(a$outcome, a$longitudinal), c("k x1", "k x2", "cnt x2"))
  ap <- associations(fxModel("priors"))
  expect_equal(paste(ap$outcome, ap$longitudinal), c("yg x1", "yg x2", "yg x3"))
  expect_equal(nrow(associations(fxModel("surv"))), 1L)
})

test_that("associations() recover the simulated signs (fxData(): k, y and time driven by x1 positively, yg by x3 negatively)", {
  expect_gt(associations(fxModel("multi"))$mean[1], 0)
  ab <- associations(fxModel("beta"))                      # simulated gamma 0.8 (plogis(-1 + 0.8 eta1))
  expect_true(ab$q0.025 < 0.8 && 0.8 < ab$q0.975)
  as <- associations(fxModel("surv"))                      # hazard exp(-2.5 + 0.8 (1 + b0))
  expect_true(as$q0.025 < 0.8 && 0.8 < as$q0.975)
  ap <- associations(fxModel("priors"))                    # yg = 0.5 eta1 + 0.3 eta2 - 0.4 x3 + ...
  expect_gt(ap$mean[ap$longitudinal == "x1"], 0)
  expect_lt(ap$mean[ap$longitudinal == "x3"], 0)
})

# ---------------------------------------------------------------------------
# random effects

test_that("randomEffects(): one row per subject and effect; correlated block read from one field", {
  m <- fxModel("multi")
  n <- fxCounts()$subjects
  ids <- sort(unique(fxData()$long$id))
  re <- randomEffects(m)
  expect_named(re, c("response", "group", "term", "level", statCols, "mode"))
  expect_equal(as.vector(table(re$response, re$term)[c("x1", "x2"), ]), c(n, n, n, 0L))
  expectOrderedStats(re)
  r1 <- randomEffects(m, "x1")
  expect_equal(r1$level[r1$term == "t"], ids)
  expect_equal(r1$level[r1$term == "(Intercept)"], ids)
  r2 <- randomEffects(m, "x2")
  expect_equal(nrow(r2), n)
  expect_equal(r2$level, ids)
  expect_true(all(r2$group == "id" & r2$term == "(Intercept)"))
  # the correlated block of x1 is one INLA field of 2 n rows, x2's single effect its own field
  sr <- inlaFit(m)$summary.random
  bx1 <- inlaSpec(m)$random$x1[[1]]
  expect_length(bx1$effects, 2L)   # a correlated block: all effects in one field
  expect_equal(nrow(sr[[bx1$index[1]]]), 2 * n)
  expect_length(inlaSpec(m)$random$x2[[1]]$effects, 1L)
})

test_that("randomEffects(response = ): selects, in the order of the model; empty for an outcome without REs, errors", {
  m <- fxModel("multi")
  re <- randomEffects(m)
  expect_equal(randomEffects(m, c("x2", "x1")), re)          # x1 before x2, as passed
  expect_equal(unique(randomEffects(m, c("x2", "k", "x1"))$response), c("x1", "x2"))
  expect_equal(nrow(randomEffects(m, "k")), 0L)
  expect_named(randomEffects(m, "k"), c("response", "group", "term", "level"))
  expect_error(randomEffects(m, c("x1", "zz")), "Unknown response\\(s\\): zz")
})

test_that("randomEffects() of 'priors': outcome RE on site, two blocks on one group, a `||` block", {
  m <- fxModel("priors")
  n <- fxCounts()$subjects
  re <- randomEffects(m)
  tab <- table(re$response, re$term, re$group)
  expect_equal(unname(tab["yg", "(Intercept)", "site"]), fxCounts()$sites)
  expect_equal(unname(tab["yg", , "id"]), c(0L, 0L))
  expect_equal(unname(tab[c("x1", "x2", "x3"), "(Intercept)", "id"]), rep(n, 3))
  expect_equal(unname(tab[c("x1", "x2", "x3"), "t", "id"]), rep(n, 3))
  expect_equal(sort(unique(re$level[re$response == "yg"])), levels(fxData()$long$site))
  # only x1's block is read from a correlated field (one field of 2 n rows); x2's two re()
  # terms and x3's `||` block (lme4 semantics: one independent block per effect) are
  # single-effect blocks with a field of n rows each
  sr <- inlaFit(m)$summary.random
  rnd <- inlaSpec(m)$random
  expect_length(rnd$x1[[1]]$effects, 2L)
  expect_equal(nrow(sr[[rnd$x1[[1]]$index[1]]]), 2 * n)
  for (r in c("x2", "x3")) {
    expect_length(rnd[[r]], 2L)
    expect_equal(vapply(rnd[[r]], function(b) b$effects, ""), c("(Intercept)", "t"), info = r)
    for (b in rnd[[r]]) expect_equal(nrow(sr[[b$index[1]]]), n)
  }
})

test_that("randomEffects() labels the levels of a model refitted with new subjects", {
  p <- fxPredict()
  re0 <- randomEffects(fxModel("multi"), "x1")
  re2 <- randomEffects(attr(p, "model"), "x1")
  s01 <- function(re) re$mean[re$level == "S01" & re$term == "(Intercept)"]
  expect_equal(s01(re2), s01(re0), tolerance = 0.05)
  expect_equal(nrow(re2), nrow(re0) + 2 * 3)                 # 3 new subjects, 2 effects
  expect_setequal(re2$level[re2$term == "t"], c(unique(fxData()$long$id), "N1", "N2", "N3"))
})

# B17: the subject ids of the fitted model may be text, a factor (any level order) or
# integers, and the rows may come in any order; every result row must carry the right
# subject. Each variant of the "beta" data is fitted once (~5 s) and compared with the
# fixture and with INLA's own field.
idVariants <- function() {
  d <- fxData()$long
  idNum <- as.integer(sub("S", "", d$id))
  set.seed(31)
  shuffle <- function(x) x[sample(nrow(x)), ]
  list(
    integer = transform(d, id = 41L - idNum),                                   # 40, 39, ... in order of appearance
    factor = transform(d, id = factor(d$id, levels = rev(sort(unique(d$id))))),  # levels in reverse order
    text = shuffle(transform(d, id = ifelse(idNum <= 20, paste0("Z", d$id), paste0("A", d$id))))  # sort != appearance, rows shuffled
  )
}
idMap <- function(variant, d0 = fxData()$long) {   # original id -> id of the variant, as character
  v <- unique(data.frame(from = d0$id, to = as.character(variant$id[match(d0$id, sub("^[ZA]", "", as.character(variant$id)))]),
                         stringsAsFactors = FALSE))
  stats::setNames(v$to, v$from)
}

test_that("randomEffects(), fittedValues(): subject ids as integers, factor or text, rows in any order", {
  s <- fxSpec("beta")
  m0 <- fxModel("beta")
  re0 <- randomEffects(m0, "x1")
  fv0 <- fittedValues(m0, "x1")
  key0 <- paste(fv0$id, fv0$t)
  for (nm in names(idVariants())) {
    d <- idVariants()[[nm]]
    m <- suppressMessages(jsm(s$formOut, s$formLong, data = d, family = s$family))
    skip_if(is.null(inlaFit(m)), paste("INLA fit of the", nm, "variant failed"))
    map <- if (nm == "integer") stats::setNames(as.character(41L - as.integer(sub("S", "", unique(fxData()$long$id)))), unique(fxData()$long$id)) else
      if (nm == "factor") stats::setNames(unique(fxData()$long$id), unique(fxData()$long$id)) else idMap(d)
    re <- randomEffects(m, "x1")
    # one row per subject and effect, levels are the variant's values as character
    expect_setequal(re$level, as.character(unique(d$id)))
    expect_equal(nrow(re), nrow(re0))
    # the same effects as the fixture, attached to the same subjects (tolerance: INLA's
    # optimisation with the rows in another order)
    for (term in c("(Intercept)", "t")) {
      a <- re0[re0$term == term, ]; b <- re[re$term == term, ]
      expect_equal(b$mean[match(map[a$level], b$level)], a$mean, tolerance = 0.01, info = paste(nm, term))   # INLA with the rows in another order
    }
    # ... and they are INLA's own values for that label
    sr <- inlaFit(m)$summary.random$x1_1_id
    lab <- paste0("x1|1|", re$level[re$term == "(Intercept)"])
    expect_equal(re$mean[re$term == "(Intercept)"], sr$mean[match(lab, as.character(sr$ID))], info = nm)
    # the intercept effect follows the subject's own data (association, not just labels)
    dev <- tapply(d$x1 - (1 + 0.5 * d$t), as.character(d$id), mean)   # deviation from the simulated line
    b0 <- re$mean[re$term == "(Intercept)"]
    expect_gt(stats::cor(b0, dev[re$level[re$term == "(Intercept)"]]), 0.9)
    # fitted values keep the rows of the data as given
    fv <- fittedValues(m, "x1")
    expect_equal(fv$observed, d$x1)
    expect_equal(as.character(fv$id), as.character(d$id))
    key <- paste(names(map)[match(as.character(fv$id), map)], fv$t)
    expect_equal(fv$mean, fv0$mean[match(key, key0)], tolerance = 0.01, info = nm)
    # residual plot data and random-effect plots carry the variant's ids
    pd <- plotResidualsData(m, "x1", "time", NULL, NULL, 300, NULL)
    expect_setequal(unlist(lapply(pd$panels, subjectIds)), as.character(unique(d$id)))
  }
})

# ---------------------------------------------------------------------------
# variance components

test_that("varianceComponents(): one row per SD / correlation + residual, sd and variance scales", {
  m <- fxModel("multi")
  vs <- varianceComponents(m)
  vv <- varianceComponents(m, scale = "variance")
  expect_s3_class(vs, "mvremVarCov")
  expect_named(vs, c("table", "cov", "cor", "marginals", "scale"))
  expect_equal(vs$scale, "sd")
  expect_equal(vv$scale, "variance")
  expect_equal(vs$table$parameter, c("sd", "sd", "correlation", "sd", "sd"))
  expect_equal(vv$table$parameter, c("variance", "variance", "correlation", "variance", "variance"))
  expect_equal(vs$table$block, c(rep("x1 | id [(Intercept), t]", 3), "x2 | id [(Intercept)]", "residual"))
  expect_equal(vs$table$response, c("x1", "x1", "x1", "x2", "x1"))  # residual: gaussian x1 only
  expect_equal(vs$table$effect, c("(Intercept)", "t", "(Intercept)", "(Intercept)", "residual"))
  expect_equal(vs$table$effect2, c(NA, NA, "t", NA, NA))
  expectOrderedStats(vs$table)
  expectOrderedStats(vv$table)
  expect_length(vs$marginals, nrow(vs$table))
  expect_true(all(vapply(vs$marginals, function(d) all(d[, 2] >= 0), NA)))
  # the median is invariant under the monotone transforms (quantiles come from inla.tmarginal
  # of the precision marginal: numerical, so not exact)
  isSd <- vs$table$parameter == "sd"
  expect_equal(vv$table$q0.5[isSd], vs$table$q0.5[isSd]^2, tolerance = 0.05)
  expect_equal(vv$table[!isSd, statCols], vs$table[!isSd, statCols])
})

test_that("varianceComponents() agrees with the precisions and correlations of hyperparameters()", {
  m <- fxModel("multi")
  vs <- varianceComponents(m)
  h <- hyperparameters(m)
  prec <- h$q0.5[h$response == "x2" & h$parameter == "precision of (Intercept) | id"]
  expect_equal(vs$table$q0.5[4], 1 / sqrt(prec), tolerance = 0.05)
  res <- h$q0.5[h$response == "x1" & h$parameter == "residual precision"]
  expect_equal(vs$table$q0.5[5], 1 / sqrt(res), tolerance = 0.05)
  expect_equal(vs$table$mean[3], h$mean[h$parameter == "correlation of (Intercept), t | id"])
})

test_that("varianceComponents(): cov / cor matrices are plug-ins of the mean SDs and correlations", {
  m <- fxModel("multi")
  vs <- varianceComponents(m)
  expect_named(vs$cov, c("x1 | id [(Intercept), t]", "x2 | id [(Intercept)]"))
  expect_named(vs$cor, names(vs$cov))
  expect_equal(vs$cov[[2]], matrix(vs$table$mean[4]^2, 1, 1, dimnames = list("(Intercept)", "(Intercept)")))
  expect_equal(dimnames(vs$cov[[1]]), list(c("(Intercept)", "t"), c("(Intercept)", "t")))
  expect_true(isSymmetric(vs$cov[[1]]))
  expect_true(isSymmetric(vs$cor[[1]]))
  expect_equal(unname(diag(vs$cor[[1]])), c(1, 1))
  expect_equal(vs$cor[[1]][1, 2], vs$table$mean[3])
  expect_equal(unname(sqrt(diag(vs$cov[[1]]))), vs$table$mean[1:2])
  expect_equal(vs$cov[[1]][1, 2], vs$table$mean[1] * vs$table$mean[2] * vs$table$mean[3])
})

test_that("varianceComponents(response = ): filters, empty result, errors, print", {
  m <- fxModel("multi")
  v2 <- varianceComponents(m, "x2")
  expect_equal(v2$table$response, "x2")
  expect_named(v2$cov, "x2 | id [(Intercept)]")
  expect_length(v2$marginals, 1L)
  vk <- varianceComponents(m, "k")
  expect_equal(nrow(vk$table), 0L)
  expect_length(vk$cov, 0L)
  expect_length(vk$marginals, 0L)
  expect_output(print(vk), "Variance components")
  expect_output(print(varianceComponents(m)), "x1 \\| id")
  expect_error(varianceComponents(m, "zz"), "Unknown response")
  expect_error(varianceComponents(m, scale = "precision"), "should be one of")
})

test_that("varianceComponents() of 'priors': one block per re() term, residual rows of gaussian / t", {
  vc <- varianceComponents(fxModel("priors"))
  tb <- vc$table
  expect_equal(tb$block, c("yg | site [(Intercept)]", rep("x1 | id [(Intercept), t]", 3),
                           "x2 | id [(Intercept)]", "x2 | id [t]", "x3 | id [(Intercept)]", "x3 | id [t]",
                           rep("residual", 3)))
  expect_equal(tb$response[tb$block == "residual"], c("yg", "x1", "x3"))     # not x2 (negbin)
  expect_equal(sum(tb$parameter == "correlation"), 1L)
  expect_equal(tb$response[tb$parameter == "correlation"], "x1")
  expect_named(vc$cov, c("yg | site [(Intercept)]", "x1 | id [(Intercept), t]", "x2 | id [(Intercept)]",
                         "x2 | id [t]", "x3 | id [(Intercept)]", "x3 | id [t]"))
  expect_equal(unname(vapply(vc$cov, nrow, 0L)), c(1L, 2L, 1L, 1L, 1L, 1L))
  expect_equal(unname(vapply(vc$cor, nrow, 0L)), c(1L, 2L, 1L, 1L, 1L, 1L))
  expect_equal(dimnames(vc$cov[["x2 | id [t]"]]), list("t", "t"))
  expectOrderedStats(tb[tb$parameter != "correlation", ])
})

test_that("varianceComponents() recovers the simulated SDs of 'beta'", {
  # fxData(): x1 has random intercepts with SD 0.5 and slopes with SD 0.2, residual SD 0.2
  vc <- varianceComponents(fxModel("beta"))
  sds <- vc$table[vc$table$parameter == "sd" & vc$table$block != "residual", ]
  expect_equal(sds$effect, c("(Intercept)", "t"))
  expect_equal(sds$mean, c(0.5, 0.2), tolerance = 0.35)
  expect_equal(vc$table$mean[vc$table$block == "residual"], 0.2, tolerance = 0.2)
})

# ---------------------------------------------------------------------------
# fitted values, residuals, model criteria

test_that("fittedValues(): one row per data row, observed on the proportion / rate scale", {
  m <- fxModel("multi")
  d <- fxData()$long
  nr <- fxCounts()$rows
  fv <- fittedValues(m)
  expect_equal(as.vector(table(factor(fv$response, fxCounts()$responses))), rep(nr, 4))
  expect_equal(unique(fv$response), fxCounts()$responses)
  expect_false(".key" %in% names(fv))
  expectOrderedStats(fv)
  k <- fv[fv$response == "k", ]
  expect_equal(k$observed, d$k / d$n)
  expect_true(all(k$mean > 0 & k$mean < 1))
  cn <- fv[fv$response == "cnt", ]
  expect_equal(cn$observed, d$cnt / d$pt)
  expect_true(all(is.na(cn$k)))                             # columns of other equations: NA
  expect_equal(fv$observed[fv$response == "x2"], d$x2)
  expect_equal(fv$observed[fv$response == "x1"], d$x1)
})

test_that("fittedValues(type = 'linear'): no observed column; medians map through the link", {
  m <- fxModel("multi")
  nr <- fxCounts()$rows
  k <- fittedValues(m, "k")
  lin <- fittedValues(m, "k", type = "linear")
  expect_false("observed" %in% names(lin))
  expect_equal(nrow(lin), nr)
  expect_equal(lin$q0.5, stats::qlogis(k$q0.5), tolerance = 0.01)
  cn <- fittedValues(m, "cnt")
  lc <- fittedValues(m, "cnt", type = "linear")
  expect_equal(lc$q0.5, log(cn$q0.5), tolerance = 0.01)
  expect_equal(fittedValues(m, "x1", type = "linear")$mean, fittedValues(m, "x1")$mean, tolerance = 1e-6)
})

test_that("fittedValues(response = ): selects, in the order of the model; errors, empty selection", {
  m <- fxModel("multi")
  nr <- fxCounts()$rows
  expect_equal(fittedValues(m, c("x2", "k"))$response, rep(c("k", "x2"), each = nr))
  expect_error(fittedValues(m, "zz"), "Unknown response")
  expect_error(fittedValues(m, type = "response"), "should be one of")
  expect_equal(NROW(fittedValues(m, character())), 0L)
})

test_that("fittedValues() of a survival outcome: lambda per subject, observed NA", {
  m <- fxModel("surv")
  fv <- fittedValues(m, "time")
  expect_equal(nrow(fv), nrow(fxData()$surv))
  expect_equal(fv$id, fxData()$surv$id)
  expect_true(all(is.na(fv$observed)))
  expect_true(all(fv$mean > 0))
  lin <- fittedValues(m, "time", type = "linear")
  expect_equal(lin$q0.5, log(fv$q0.5), tolerance = 0.01)
})

test_that("residuals(): fittedValues() plus a residual column on the response scale", {
  m <- fxModel("multi")
  r <- residuals(m, c("k", "cnt"))
  expect_equal(unique(r$response), c("k", "cnt"))
  expect_equal(r[setdiff(names(r), "residual")], fittedValues(m, c("k", "cnt")))
  expect_equal(names(r)[length(r)], "residual")
  expect_equal(sign(r$residual), sign(r$observed - r$mean))
  expect_false(anyNA(r$residual))
  rx <- residuals(m, "x1")
  expect_lt(abs(mean(rx$residual)), 0.05)                  # gaussian with intercept: centred
  expect_equal(nrow(residuals(m)), 4 * fxCounts()$rows)
  expect_error(residuals(m, "zz"), "Unknown response")
})

test_that("residuals(): survival responses are skipped, with a warning when named", {
  ms <- fxModel("surv")
  expect_no_warning(rs <- residuals(ms))                   # all responses: silently without time
  expect_equal(unique(rs$response), "x1")
  expect_error(residuals(ms, "time"), "No residuals: the selected responses are survival")
  expect_warning(r2 <- residuals(ms, c("time", "x1")), "No residuals for survival response\\(s\\) time")
  expect_equal(unique(r2$response), "x1")
  expect_equal(residualResponses(ms, NULL), "x1")
  expect_warning(expect_equal(residualResponses(ms, c("time", "x1")), "x1"), "time; skipped")
})

test_that("modelFit(): one row per response + overall, INLA's totals in the overall row", {
  m <- fxModel("multi")
  mf <- modelFit(m)
  expect_named(mf, c("response", "dic", "waic", "logScore", "cpoFailures", "mlik"))
  expect_equal(mf$response, c(fxCounts()$responses, "overall"))
  expect_true(all(is.na(mf$mlik[1:4])))
  expect_true(is.finite(mf$mlik[5]))
  expect_equal(mf$mlik[5], unname(inlaFit(m)$mlik[1]))
  expect_equal(mf$dic[5], inlaFit(m)$dic$dic)
  expect_equal(mf$waic[5], inlaFit(m)$waic$waic)
  expect_true(all(is.finite(mf$dic) & is.finite(mf$waic) & is.finite(mf$logScore)))
  expect_true(all(mf$cpoFailures >= 0))
  expect_true(all(is.finite(modelFit(fxModel("surv"))$logScore)))
  expect_equal(modelFit(fxModel("beta"))$response, c("y", "x1", "overall"))
})

test_that("modelFit() works with the criteria given through inlaOptions instead of criteria = TRUE", {
  skip_if_not_installed("INLA")
  m <- fxModel("beta")
  mi <- suppressMessages(update(m, criteria = FALSE,
                                inlaOptions = list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))))
  skip_if(is.null(inlaFit(mi)), "INLA fit failed")
  expect_false(mi$call$criteria)
  expect_equal(inlaArguments(mi)$control.compute, list(dic = TRUE, waic = TRUE, cpo = TRUE))
  mf <- modelFit(mi)
  expect_named(mf, c("response", "dic", "waic", "logScore", "cpoFailures", "mlik"))
  expect_equal(mf, modelFit(m), tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# estimatesTable()

test_that("estimatesTable() is consistent with the tables it joins", {
  m <- fxModel("multi")
  et <- estimatesTable(m)
  expect_s3_class(et, "mvremEstimates")
  expect_named(et, c("response", "type", "parameter", statCols))
  expectOrderedStats(et[et$type != "correlation", ])
  fe <- fixedEffects(m)
  a <- associations(m)
  vc <- varianceComponents(m)$table
  expect_equal(sum(et$type == "fixed"), nrow(fe))
  expect_equal(sum(et$type == "association"), nrow(a))
  expect_equal(sum(et$type %in% c("sd", "correlation", "residual")), nrow(vc))
  expect_equal(sum(et$type == "family"), 0L)               # binomial / poisson: no parameters
  expect_equal(sum(et$type == "other"), 0L)
  expect_equal(nrow(et), nrow(fe) + nrow(a) + nrow(vc))
  expect_equal(sort(et$mean[et$type == "fixed"]), sort(fe$mean))
  expect_equal(et$mean[et$type == "association"], a$mean)
  expect_equal(et$mean[et$type == "residual"], vc$mean[vc$block == "residual"])
  expect_equal(et$mean[et$type == "sd"], vc$mean[vc$parameter == "sd" & vc$block != "residual"])
  expect_true(all(nzchar(et$parameter)))
})

test_that("estimatesTable(): rows grouped by response (in the order passed), types in a fixed order", {
  m <- fxModel("multi")
  et <- estimatesTable(m)
  expect_equal(unique(et$response), fxCounts()$responses)
  ord <- c("fixed", "association", "sd", "correlation", "residual", "family", "other")
  for (r in unique(et$response)) expect_false(is.unsorted(match(et$type[et$response == r], ord)), info = r)
  expect_equal(et$parameter[et$response == "cnt"], c("(Intercept)", "grpb", "grpc", "gamma (x2 -> cnt)"))
  expect_equal(et$parameter[et$response == "x2"], c("(Intercept)", "t", "SD (Intercept) | id"))
  expect_equal(et$parameter[et$response == "x1"],
               c("(Intercept)", "t", "SD (Intercept) | id", "SD t | id", "corr (Intercept), t | id", "residual SD"))
})

test_that("estimatesTable() of 'priors' and 'surv': family rows size, dof, phi, alpha; none unrecognised", {
  ep <- estimatesTable(fxModel("priors"))
  expect_equal(paste(ep$response, ep$parameter)[ep$type == "family"], c("x2 size", "x3 dof"))
  expect_equal(ep$response[ep$type == "residual"], c("yg", "x1", "x3"))
  expect_equal(ep$parameter[ep$response == "yg" & ep$type == "sd"], "SD (Intercept) | site")
  expect_equal(ep$parameter[ep$response == "x2" & ep$type == "sd"], c("SD (Intercept) | id", "SD t | id"))
  expect_equal(sum(ep$type == "correlation"), 1L)
  expect_equal(sum(ep$type == "other"), 0L)
  h <- hyperparameters(fxModel("priors"))
  expect_equal(ep$mean[ep$parameter == "size"], h$mean[h$parameter == "size"])
  eb <- estimatesTable(fxModel("beta"))
  expect_equal(paste(eb$response, eb$parameter)[eb$type == "family"], "y phi")
  es <- estimatesTable(fxModel("surv"))
  expect_equal(es$parameter[es$response == "time"], c("(Intercept)", "gamma (x1 -> time)", "alpha"))
  expect_equal(es$type[es$response == "time"], c("fixed", "association", "family"))
})

test_that("estimatesTable(response = ) selects responses", {
  m <- fxModel("multi")
  ex <- estimatesTable(m, "x2")
  expect_true(all(ex$response == "x2"))
  expect_equal(ex[statCols], estimatesTable(m)[estimatesTable(m)$response == "x2", statCols] |> `rownames<-`(NULL))
  expect_error(estimatesTable(m, "zz"), "Unknown response")
  expect_equal(nrow(estimatesTable(m, character())), 0L)
})

test_that("estimatesTable(response = ) follows the order passed whatever the order given (Q3)", {
  m <- fxModel("multi")                                   # passed: k, cnt (outcomes), x1, x2
  et <- estimatesTable(m)
  for (sel in list(c("x2", "k"), c("k", "x2"), c("x2", "x1", "cnt"), c("x1", "k", "x2", "cnt"))) {
    e <- estimatesTable(m, sel)
    expect_equal(unique(e$response), intersect(fxCounts()$responses, sel), info = toString(sel))
    expect_equal(e[statCols], et[et$response %in% sel, statCols] |> `rownames<-`(NULL), info = toString(sel))
  }
})

test_that("format(): the text style is what print() shows; digits; unknown style", {
  et <- estimatesTable(fxModel("beta"))
  expect_type(format(et), "character")
  expect_length(format(et), 1L)
  expect_equal(format(et), format(et, "text"))
  expect_output(res <- print(et), "corr \\(Intercept\\), t \\| id")
  expect_identical(res, et)
  expect_output(print(et), "mean \\(sd\\)")
  expect_error(format(et, "html"), "should be one of")
  # digits: the first row's mean with 2 significant digits, not 4
  txt <- strsplit(format(et, "text", digits = 2), "\n")[[1]]
  expect_match(txt[2], trimws(formatC(signif(et$mean[1], 2), digits = 2, format = "g")), fixed = TRUE)
  expect_false(grepl(trimws(formatC(signif(et$mean[1], 4), digits = 4, format = "g")), txt[2], fixed = TRUE))
  expect_match(format(et, "text", digits = 2), "gamma \\(x1 -> y\\) +[0-9.]+ \\([0-9.]+\\) +\\[[0-9.-]+, [0-9.]+\\]")
  expect_snapshot(print(et, digits = 2))
})

test_that("format(style = 'markdown'): a pipe table, `|` of labels escaped, response on its first row", {
  et <- estimatesTable(fxModel("beta"))
  md <- strsplit(format(et, "markdown"), "\n")[[1]]
  expect_equal(md[2], "|---|---|---|---|")
  expect_match(md[3], "^\\| y \\| \\(Intercept\\) \\| ")
  expect_match(md[4], "^\\|  \\| t \\| ")
  expect_true(any(grepl("SD (Intercept) \\| id", md, fixed = TRUE)))
  expect_true(all(lengths(regmatches(md, gregexpr("(?<!\\\\)\\|", md, perl = TRUE))) == 5L))
  expect_equal(sum(grepl("^\\| x1 \\|", md)), 1L)
  expect_snapshot(cat(format(et, "markdown", digits = 2)))
})

test_that("format(style = 'latex'): a tabular, -> as $\\to$, special characters escaped", {
  et <- estimatesTable(fxModel("beta"))
  tex <- format(et, "latex")
  expect_match(tex, "^\\\\begin\\{tabular\\}\\{llrr\\}")
  expect_match(tex, "\\\\end\\{tabular\\}$")
  expect_match(tex, "gamma (x1 $\\to$ y)", fixed = TRUE)
  expect_match(tex, "SD (Intercept) | id & ", fixed = TRUE)
  expect_match(tex, "95\\% interval", fixed = TRUE)
  expect_snapshot(cat(format(et, "latex", digits = 2)))
})

test_that("format(style = 'array'): the table as display math, cells in \\text{}", {
  et <- estimatesTable(fxModel("beta"))
  arr <- strsplit(format(et, "array", digits = 2), "\n")[[1]]
  expect_equal(arr[1], "\\begin{array}{llrr}")
  expect_equal(arr[length(arr)], "\\end{array}")
  expect_match(arr[5], "^\\\\text\\{y\\} & \\\\text\\{\\(Intercept\\)\\} & -?[0-9.]+\\\\ \\(")
  expect_match(arr[6], "^ & \\\\text\\{t\\} & ")                   # response only on its first row
  expect_true(any(grepl("\\text{gamma (x1 $\\to$ y)}", arr, fixed = TRUE)))
  rows <- arr[5:(4 + nrow(et))]
  expect_true(all(grepl("\\\\$", rows)))
  expect_snapshot(cat(format(et, "array", digits = 2)))
})

test_that("format() / print() of an empty estimates table", {
  et <- estimatesTable(fxModel("multi"), character())
  expect_equal(nrow(et), 0L)
  txt <- format(et, "text")
  expect_match(txt, "response type parameter mean \\(sd\\) 95% interval")
  expect_match(txt, "no estimates")
  expect_output(print(et), "no estimates")
  # header lines only, no empty data row
  expect_length(strsplit(format(et, "markdown"), "\n")[[1]], 2L)
  tex <- strsplit(format(et, "latex"), "\n")[[1]]
  expect_length(tex, 6L)                                     # begin, hline, header, hline, hline, end
  expect_false(any(grepl("^ & ", tex)))
  arr <- strsplit(format(et, "array"), "\n")[[1]]
  expect_length(arr, 6L)
  expect_false(any(grepl("^ & ", arr)))
})

# ---------------------------------------------------------------------------
# predict()

test_that("predict(): one data frame, responses picked by its columns, one block per response", {
  p <- fxPredict()                        # "multi" on fxNewdata: N1, N2, N3 new; S01 at t = 2, 3
  nd <- fxNewdata
  expect_s3_class(p, "mvremPredictions")
  expect_s3_class(p, "data.frame")
  expect_equal(attr(p, "type"), "mean")
  expect_named(p, c("response", names(nd), statCols, "mode"))
  expect_equal(p$response, rep(fxCounts()$responses, each = nrow(nd)))
  expect_equal(p$id, rep(nd$id, 4))
  expect_equal(p$t, rep(nd$t, 4))
  expectOrderedStats(p)
  expect_true(all(p$mean[p$response == "k"] > 0 & p$mean[p$response == "k"] < 1))   # proportion scale
  expect_true(all(p$mean[p$response %in% c("cnt", "x2")] > 0))                      # rate scale
})

test_that("predict(): new subjects follow the population line with wide intervals; an existing one its own", {
  p <- fxPredict()
  m <- fxModel("multi")
  x1 <- p[p$response == "x1", ]
  new <- x1$id != "S01"
  # a new subject's random effects sit at their prior mean: fitted mean = fixed effects only
  fe <- fixedEffects(m)
  b <- fe$mean[fe$response == "x1"]                        # (Intercept), t
  expect_equal(x1$mean[new], b[1] + b[2] * x1$t[new], tolerance = 0.01)
  expect_gt(mean(x1$sd[new]), 2 * mean(x1$sd[!new]))
  # an existing subject continues its own line: fixed + its random effects (from the original fit)
  re <- randomEffects(m, "x1")
  u <- re$mean[re$level == "S01"]                          # (Intercept), t
  expect_equal(x1$mean[!new], (b[1] + u[1]) + (b[2] + u[2]) * x1$t[!new], tolerance = 0.01)
})

test_that("predict(): the refitted model holds the new rows with a missing response; the original is unchanged", {
  p <- fxPredict()
  m <- fxModel("multi")
  nr <- fxCounts()$rows
  nn <- nrow(fxNewdata)
  m2 <- attr(p, "model")
  expect_s3_class(m2, "mvremModel")
  expect_equal(nrow(modelData(m2)$long$x1), nr + nn)
  expect_true(all(is.na(utils::tail(modelData(m2)$long$x1$x1, nn))))
  expect_true(all(is.na(utils::tail(modelData(m2)$out$k$k, nn))))
  expect_s3_class(modelData(m2)$out$cnt$grp, "factor")
  expect_equal(levels(modelData(m2)$out$cnt$grp), c("a", "b", "c"))
  expect_equal(levels(modelData(m2)$long$x1$id), c(sort(unique(fxData()$long$id)), "N1", "N2", "N3"))
  expect_equal(dataSummary(m2)$equations$complete, rep(nr, 4))
  expect_equal(dataSummary(m2)$equations$rows, rep(nr + nn, 4))
  # the predictions are the fitted values of the last rows, on the mean scale (missing response)
  fv <- fittedValues(m2, "x1")
  expect_equal(nrow(fv), nr + nn)
  expect_true(all(is.na(utils::tail(fv$observed, nn))))
  expect_equal(utils::tail(fv$mean, nn), p$mean[p$response == "x1"])
  expect_equal(utils::tail(fittedValues(m2, "k")[statCols], nn), p[p$response == "k", statCols], ignore_attr = TRUE)
  expect_equal(nrow(fittedValues(m, "x1")), nr)
  expect_equal(nrow(modelData(m)$long$x1), nr)
  # the linear predictor of an exposure response is log(rate): medians map
  lc <- fittedValues(m2, "cnt", type = "linear")
  expect_equal(log(p$q0.5[p$response == "cnt"]), utils::tail(lc$q0.5, nn), tolerance = 0.01)
})

test_that("predict(): list newdata; the response column is ignored; type = 'linear' (INLA mocked)", {
  m <- fxModel("multi")
  cap <- mockInlaRun()
  nd <- list(x1 = data.frame(id = "S02", t = c(0, 3), x1 = 999),
             cnt = data.frame(id = "S02", t = 1, pt = 2, grp = factor("b", levels = c("a", "b", "c"))))
  p <- predict(m, nd, type = "linear")
  expect_s3_class(p, "mvremPredictions")
  expect_equal(attr(p, "type"), "linear")
  expect_equal(p$response, c("cnt", "x1", "x1"))          # the order of the model, not of the list
  expect_setequal(names(p), c("response", "id", "t", "x1", "pt", "grp", statCols, "mode"))
  expect_true(all(is.na(p$pt[2:3])) && is.na(p$x1[1]))
  m2 <- cap$model
  expect_true(all(is.na(utils::tail(modelData(m2)$long$x1$x1, 2))))   # 999 not used
  expect_equal(nrow(modelData(m2)$long$x1), fxCounts()$rows + 2)
  expect_equal(nrow(modelData(m2)$out$cnt), fxCounts()$rows + 1)
  expect_equal(nrow(modelData(m2)$long$x2), fxCounts()$rows)            # not predicted: unchanged
  expect_equal(utils::tail(modelData(m2)$out$cnt$grp, 1), factor("b", levels = c("a", "b", "c")))
})

test_that("predict(verbose = , threads = , inlaOptions = ) reach the new INLA run (mocked)", {
  m <- fxModel("multi")
  cap <- mockInlaRun()
  p <- predict(m, data.frame(id = "S01", t = 1), response = "x1", verbose = TRUE, threads = 1,
               inlaOptions = list(control.compute = list(dic = FALSE, config = TRUE)))
  expect_equal(nrow(p), 1L)
  expect_true(cap$options$verbose)
  expect_equal(cap$options$threads, 1)
  expect_true(cap$args$verbose)
  expect_equal(cap$args$num.threads, 1)
  # merged over the model's options (the fixture has dic / waic / cpo from criteria = TRUE)
  expect_equal(cap$args$control.compute, list(dic = FALSE, waic = TRUE, cpo = TRUE, config = TRUE))
  expect_equal(attr(p, "model")$options, cap$options)
  expect_equal(m$options$verbose, FALSE)                    # the original model's options are kept
  # without them: the model's own options
  predict(m, data.frame(id = "S01", t = 1), response = "x1")
  expect_equal(cap$options, m$options)
})

test_that("predict() accepts a tibble as newdata (INLA mocked)", {
  skip_if_not_installed("tibble")
  m <- fxModel("multi")
  cap <- mockInlaRun()
  nd <- tibble::tibble(id = c("S01", "S02"), t = c(1, 2))
  p <- predict(m, nd, response = "x1")
  expect_s3_class(p, "mvremPredictions")
  expect_equal(class(p), c("mvremPredictions", "data.frame"))
  expect_equal(nrow(p), 2L)
  expect_equal(p$id, c("S01", "S02"))
  expect_equal(class(modelData(cap$model)$long$x1), "data.frame")
  pl <- predict(m, list(x1 = nd))
  expect_equal(pl$t, c(1, 2))
})

test_that("predict(): argument checks before refitting", {
  m <- fxModel("multi")
  expect_error(predict(m, data.frame(id = "S01", t = 1), response = "nope"), "Unknown response\\(s\\): nope")
  expect_error(predict(m, list(x1 = data.frame(id = "S01", t = 1), foo = data.frame(id = "S01", t = 1))),
               "Unknown response\\(s\\): foo")
  expect_error(predict(m, list(x1 = 1)), "must be a data frame or a list of data frames named by response")
  expect_error(predict(m, list(data.frame(id = "S01", t = 1))), "must be a data frame or a list")
  expect_error(predict(m, 1:3), "must be a data frame or a list")
  expect_error(predict(m, data.frame(id = "S01", t = 1)[0, ], response = "x1"), "for 'x1' has no rows")
  expect_error(predict(m, data.frame(id = "S01", t = 1), response = "k"), "for 'k' lacks the variable\\(s\\) `n`")
  expect_error(predict(m, data.frame(t = 1), response = "x1"), "lacks the variable\\(s\\) `id`")
  expect_error(predict(m, data.frame(id = "S01", t = NA_real_), response = "x1"), "Missing values in `t`")
  expect_error(predict(m, data.frame(id = "S01", t = 1), type = "response"), "should be one of")
  expect_error(predict(m, data.frame(foo = 1)), "has the variables of no equation")
  expect_error(predict(m, data.frame(foo = 1)), "needed: k: ")
  expect_error(predict(m, data.frame(id = "S01", t = 1), foo = 1), "Unknown argument\\(s\\) to predict\\(\\): foo")
  ms <- fxModel("surv")
  expect_error(predict(ms, list(time = data.frame(id = "S01", t = 0, status = 1))),
               "not available for survival responses \\(time\\)")
  expect_error(predict(ms, data.frame(id = "S01", t = 0, status = 1), response = "time"),
               "not available for survival responses")
})

test_that("predict() with a level of a factor covariate not seen in the fit is an error", {
  m <- fxModel("multi")
  nd <- data.frame(id = "S01", t = 1, n = 5, pt = 1, grp = "zzz")
  expect_error(predict(m, nd, response = "cnt"),
               "Unknown level\\(s\\) of `grp` in `newdata` for response cnt: zzz \\(the fit has: a, b, c\\)")
  expect_error(predict(m, nd), "Unknown level\\(s\\) of `grp`")               # picked by columns
  nd$grp <- factor("zzz")
  expect_error(predict(m, nd, response = "cnt"), "`grp`.*zzz")
})

test_that("predict() with an unseen value of a character covariate or of a link variable is an error", {
  d <- fxSmallData()                       # site: character A / B; g: factor u / v / w
  # the checks run before INLA: an unfitted model with a stub fit suffices
  mc <- stubFitted(jsm(list(y ~ site + scaled(x)), list(x ~ t + (1 | id)), data = d, fit = FALSE))
  expect_type(modelData(mc)$out$y$site, "character")
  expect_error(predict(mc, data.frame(id = "S1", site = "Z", t = 1), response = "y"),
               "Unknown level\\(s\\) of `site` in `newdata` for response y: Z \\(the fit has: A, B\\)")
  # g is a link variable of cnt (a covariate of x), not in cnt's own formula
  ml <- stubFitted(jsm(list(cnt | exposure(pt) ~ 1 + scaled(x)), list(x ~ g + (1 | id)), data = d,
                       family = list(cnt = poissonFamily()), fit = FALSE))
  expect_setequal(modelFormulas(ml)$out$cnt$linkVariables, c("id", "g"))
  expect_error(predict(ml, data.frame(id = "S1", pt = 1, g = "zz"), response = "cnt"),
               "Unknown level\\(s\\) of `g` in `newdata` for response cnt: zz \\(the fit has: u, v, w\\)")
  expect_error(predict(ml, data.frame(id = "S1", pt = 1, g = "zz")), "`g`.*zz")
})

test_that("predict() allows a new subject id (a grouping variable may take new values)", {
  m <- fxModel("multi")
  cap <- mockInlaRun()
  p <- predict(m, data.frame(id = c("ZZZ", "S01"), t = 1), response = "x1")
  expect_equal(p$id, c("ZZZ", "S01"))
  expect_true("ZZZ" %in% levels(modelData(cap$model)$long$x1$id))
  expect_equal(nlevels(modelData(cap$model)$long$x1$id), fxCounts()$subjects + 1L)
})

test_that("predict() with list newdata of different columns keeps the summaries last", {
  m <- fxModel("multi")
  nd <- list(x1 = data.frame(id = "S02", t = 0), cnt = data.frame(id = "S02", t = 1, pt = 2, grp = "b"))
  p <- suppressMessages(predict(m, nd))
  expect_equal(names(p), c("response", "id", "t", "pt", "grp", statCols, "mode"))
  expect_equal(p$response, c("cnt", "x1"))          # the order of the model, not of the list
  expect_equal(p$pt, c(2, NA))                      # NA where a response's newdata lacks the column
  expect_equal(as.character(p$grp), c("b", NA))
  expectOrderedStats(p)
})
