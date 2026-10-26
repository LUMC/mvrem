d <- data.frame(id = sprintf("S%03d", rep(1:3, each = 4)), t = rep(c(0, 0.5, 1, 1.5), 3),
                g = factor(rep(c("u", "v"), 6)),
                x1 = rnorm(12), x2 = rnorm(12), x3 = rnorm(12), z = rnorm(12), y = runif(12))
wish <- list(theta1 = list(prior = "wishart2d", param = c(4, 1, 0.0011, 0)))
forms <- collectFormulas(
  formOut = list(
    y ~ 1 + t + re(1 + t | id, hyper = wish) +
      scaled(x1, hyper = list(beta = list(initial = 2.2, param = c(0, 0.022)))) +
      scaled(x2)
  ),
  formLong = list(
    x1 ~ 1 + t + (1 + t | id),
    x2 ~ 1 + t + re(1 | id),
    x3 ~ 1 + I(t^2) + g + re(1 + t || id)
  )
)
inl <- buildInlaFormula(forms, prepareData(d, d, forms))
term <- function(x) as.list(x)

test_that("response columns follow the naming scheme", {
  expect_equal(inl$responses, c(y = "y_y", x1 = "x1_x1", x2 = "x2_x2", x3 = "x3_x3"))
  expect_equal(inl$zeroResponses, c(x1 = "x1_x1_0eq", x2 = "x2_x2_0eq", x3 = "x3_x3_0eq"))
  expect_equal(deparse1(inl$formula[[2]]), "MX")
})

test_that("fixed terms: <resp>_<term>, intercept as 1", {
  expect_equal(inl$fixed$y$columns, c(y_1 = "(Intercept)", y_t = "t"))
  expect_equal(inl$fixed$x3$columns, c(x3_1 = "(Intercept)", x3_I.t.2. = "I(t^2)", x3_gv = "gv"))
})

test_that("correlated random effects: iid2d + copy weighted by the fixed column", {
  b <- inl$random$y[[1]]
  expect_equal(b$index, c("y_1_id", "y_t_id"))
  expect_equal(b$weight, c(NA, "y_t"))
  expect_equal(b$levels, "y_1_t_id")
  t1 <- term(b$terms[[1]])
  expect_equal(deparse1(t1[[2]]), "y_1_id")
  expect_equal(t1$model, "iid2d")
  expect_equal(deparse1(t1$values), "LEVELS$y_1_t_id")
  expect_equal(deparse1(t1$n), "length(LEVELS$y_1_t_id)")
  expect_equal(t1$hyper, wish)
  expect_equal(deparse1(b$terms[[2]]), 'f(y_t_id, y_t, copy = "y_1_id")')
  expect_equal(inl$levels$y_1_t_id,
               c(paste0("y|1|S00", 1:3), paste0("y|t|S00", 1:3)))
})

test_that("single and uncorrelated random effects: one iid per effect", {
  expect_equal(deparse1(inl$random$x2[[1]]$terms[[1]]),
               'f(x2_1_id, model = "iid", values = LEVELS$x2_1_id, n = length(LEVELS$x2_1_id))')
  expect_length(inl$random$x3, 2)
  expect_equal(vapply(inl$random$x3, `[[`, "", "index"), c("x3_1_id", "x3_t_id"))
  expect_equal(deparse1(inl$random$x3[[2]]$terms[[1]]),
               'f(x3_t_id, x3_t, model = "iid", values = LEVELS$x3_t_id, n = length(LEVELS$x3_t_id))')
  expect_equal(inl$levels$x3_t_id, paste0("x3|t|S00", 1:3))
})

test_that("zero-equation terms: one per (long, out) pair, keyed by group and covariates", {
  expect_named(inl$zero, c("x1_y", "x2_y"))
  z <- inl$zero$x1_y
  expect_equal(z$vars, c("id", "t"))
  expect_equal(deparse1(z$term),
               'f(x1_TO_y_id_t, x1_WEIGHT, model = "iid", values = LEVELS$x1_TO_y_id_t, hyper = list())')
  expect_equal(head(inl$levels$x1_TO_y_id_t, 2), c("x1|TO|y|S001|0", "x1|TO|y|S001|0.5"))
  expect_length(inl$levels$x1_TO_y_id_t, 12)
  expect_equal(inl$zero$x2_y$index, "x2_TO_y_id_t")  # t: fixed covariate of x2
})

test_that("scaled terms: copy of the zero-equation field, user hyper merged", {
  expect_named(inl$scaled, c("y_x1", "y_x2"))
  s <- term(inl$scaled$y_x1$term)
  expect_equal(deparse1(s[[2]]), "y_FROM_x1_id_t")
  expect_equal(s$copy, "x1_TO_y_id_t")
  expect_equal(s$hyper$beta, list(fixed = FALSE, prior = "normal", param = c(0, 0.022),
                                  initial = 2.2))
})

test_that("formula environment holds LEVELS", {
  expect_identical(get("LEVELS", environment(inl$formula)), inl$levels)
})

test_that("naming scheme is configurable", {
  nm <- namingScheme(fixed = "b.{resp}.{term}", intercept = "Int", zeroWeight = "w.{long}",
                     responseObject = "Y", levelsObject = "LV", valueSep = "/",
                     zeroValue = "{long}>{out}/{values}")
  r <- buildInlaFormula(forms, prepareData(d, d, forms), naming = nm)
  expect_equal(names(r$fixed$y$columns), c("b.y.Int", "b.y.t"))
  expect_equal(r$random$y[[1]]$weight, c(NA, "b.y.t"))
  expect_equal(r$zero$x1_y$weight, "w.x1")
  expect_equal(r$levels$x1_TO_y_id_t[1], "x1>y/S001/0")
  expect_equal(deparse1(r$formula[[2]]), "Y")
  expect_match(deparse1(r$formula), "values = LV$y_Int_t_id", fixed = TRUE)

  oldOpt <- options(mvrem.naming = nm)
  on.exit(options(oldOpt))
  expect_equal(names(buildInlaFormula(forms, prepareData(d, d, forms))$fixed$y$columns),
               c("b.y.Int", "b.y.t"))

  expect_error(namingScheme(foo = "x"), "Unknown naming template")
  expect_error(buildInlaFormula(forms, prepareData(d, d, forms),
                                naming = namingScheme(reIndex = "{resp}_{group}")),
               "duplicated")
})

test_that("complete formula is accepted by INLA", {
  skip_if_not_installed("INLA")
  i <- INLA:::inla.interpret.formula(inl$formula, data.same.len = FALSE, data = NULL)
  models <- vapply(i$random.spec, `[[`, "", "model")
  expect_equal(as.vector(table(models)[c("iid", "iid2d", "copy")]), c(5L, 2L, 4L))
})

test_that("stacks hold the formula's columns for observation and zero-equation rows", {
  y <- inl$stacks$y$main
  expect_equal(names(y), c("y_y", "y_1", "y_t", "y_1_id", "y_t_id",
                           "y_FROM_x1_id_t", "y_FROM_x2_id_t"))
  expect_equal(y$y_y, d$y)
  expect_equal(y$y_t_id[1:2], c("y|t|S001", "y|t|S001"))
  expect_equal(y$y_FROM_x1_id_t[2], "x1|TO|y|S001|0.5")

  expect_named(inl$stacks$x1, c("main", "y"))
  m <- inl$stacks$x1$main
  z <- inl$stacks$x1$y
  expect_equal(names(m), names(z))
  expect_equal(names(m), c("x1_x1", "x1_x1_0eq", "x1_1", "x1_t", "x1_1_id", "x1_t_id",
                           "x1_WEIGHT", "x1_TO_y_id_t"))
  expect_true(all(is.na(m$x1_x1_0eq)) && all(is.na(m$x1_WEIGHT)) && all(is.na(m$x1_TO_y_id_t)))
  expect_true(all(is.na(z$x1_x1)))
  expect_true(all(z$x1_x1_0eq == 0) && all(z$x1_WEIGHT == -1))
  # zero-equation labels are exactly the LEVELS of the field / scaled copies
  expect_setequal(z$x1_TO_y_id_t, inl$levels$x1_TO_y_id_t)
  expect_equal(z$x1_TO_y_id_t, y$y_FROM_x1_id_t)
  # RE labels are all in LEVELS
  expect_true(all(c(m$x1_1_id, m$x1_t_id) %in% inl$levels$x1_1_t_id))

  # uncorrelated slope without fixed counterpart still gets its weight column
  expect_true("x3_t" %in% names(inl$stacks$x3$main))
  expect_false("x3_t" %in% names(inl$fixed$x3$columns))
  expect_equal(inl$stacks$x3$main$x3_I.t.2., d$t^2)
})

test_that("formatInlaFormula shows terms with data comments", {
  txt <- formatInlaFormula(inl, nrows = 2)
  expect_s3_class(txt, "mvremFormulaText")
  expect_equal(txt[1], "MX ~ -1 +")
  expect_true(any(grepl("^   # ----- y: y ~ 1 \\+ t", txt)))
  expect_true("y_1 +" %in% txt)
  expect_true('f(y_t_id, y_t, copy = "y_1_id") +' %in% txt)
  expect_false(grepl("\\+$", txt[length(txt)]))
  # nrows rows per part: header + 2 observation + 2 zero-equation rows for x1
  i <- grep("----- x1:", txt)
  expect_equal(sum(grepl("x1\\|1\\|S001", txt[i:(i + 8)])), 4)

  # comments do not change the formula
  f <- eval(parse(text = paste(txt, collapse = "\n")))
  expect_equal(deparse1(f), deparse1(inl$formula))

  txt0 <- formatInlaFormula(inl, nrows = 0)
  expect_false(any(grepl("S001", txt0)))
  expect_output(print(txt0), "MX ~ -1 +", fixed = TRUE)
})

test_that("families per response column; zero equations are gaussian", {
  f <- collectFormulas(list(y ~ t + (1 | id) + scaled(x1)), list(x1 ~ t + (1 | id)),
                       familyOut = c(y = "beta"))
  r <- buildInlaFormula(f, prepareData(d, d, f))
  expect_equal(r$families, c(y_y = "beta", x1_x1 = "gaussian", x1_x1_0eq = "gaussian"))
})

test_that("formatInlaFormula rounds data to `digits` significant digits", {
  y1 <- format(inl$stacks$y$main$y_y[1], digits = 3)
  txt <- formatInlaFormula(inl, nrows = 1)
  expect_true(any(grepl(y1, txt, fixed = TRUE)))
  expect_false(any(grepl(format(inl$stacks$y$main$y_y[1], digits = 7), txt, fixed = TRUE)))
  expect_true(any(grepl(format(inl$stacks$y$main$y_y[1], digits = 6),
                        formatInlaFormula(inl, nrows = 1, digits = 6), fixed = TRUE)))
})

test_that("zero equations: one row per unique combination of matching variables", {
  # outcome with two rows per (id, t): both must use the same u, defined once
  dOut <- rbind(d, d)
  f <- collectFormulas(list(y ~ t + (1 | id) + scaled(x1)), list(x1 ~ t + (1 | id)))
  r <- buildInlaFormula(f, prepareData(dOut, d, f))
  z <- r$stacks$x1$y
  expect_equal(nrow(r$stacks$y$main), 24L)
  expect_equal(nrow(z), 12L)
  expect_false(anyDuplicated(z$x1_TO_y_id_t) > 0)
  expect_setequal(z$x1_TO_y_id_t, r$levels$x1_TO_y_id_t)
  # each outcome row points to the zero equation with its own id and t
  y <- r$stacks$y$main
  i <- match(y$y_FROM_x1_id_t, z$x1_TO_y_id_t)
  expect_false(anyNA(i))
  expect_equal(z$x1_t[i], dOut$t)
  expect_equal(sub("^x1\\|TO\\|y\\|([^|]+)\\|.*$", "\\1", z$x1_TO_y_id_t[i]), dOut$id)
})

test_that("matching by covariate values works with missing covariates", {
  dm <- data.frame(id = sprintf("S%03d", rep(1:3, each = 4)), t = rep(c(0, 0.5, 1, 1.5), 3),
                   x1 = rnorm(12), y = runif(12))
  # longitudinal data: response missing (row 2), covariate missing (row 5),
  # subject S003 without any longitudinal observation
  dLong <- dm
  dLong$x1[2] <- NA
  dLong$t[5] <- NA
  dLong <- dLong[dLong$id != "S003", ]
  # outcome data: covariate t (used only for matching) missing in row 7
  dOut <- dm
  dOut$t[7] <- NA

  f <- collectFormulas(list(y ~ 1 + (1 | id) + scaled(x1)), list(x1 ~ 1 + t + (1 + t | id)))
  expect_equal(f$out$y$linkVariables, "t")
  dat <- prepareData(dOut, dLong, f)
  r <- buildInlaFormula(f, dat)

  y <- r$stacks$y$main
  m <- r$stacks$x1$main
  z <- r$stacks$x1$y

  # incomplete rows dropped: outcome row 7 (t missing), long rows 2 and 5 + S003
  expect_equal(dat$out$y$.key, paste0("K", setdiff(1:12, 7)))
  expect_equal(nrow(y), 11L)
  expect_equal(nrow(m), 6L)
  expect_false(any(grepl("S003", m$x1_1_id)))

  # one zero equation per (id, t) of the outcome, including combinations without
  # a longitudinal observation (S001 at t = 0.5, S002 at t = 0, all of S003)
  expect_equal(nrow(z), 11L)
  expect_setequal(z$x1_TO_y_id_t, r$levels$x1_TO_y_id_t)
  expect_true(all(c("x1|TO|y|S001|0.5", "x1|TO|y|S002|0", "x1|TO|y|S003|1.5") %in%
                    z$x1_TO_y_id_t))
  expect_false(any(grepl("\\|NA$|\\|NA\\|", z$x1_TO_y_id_t)))

  # every outcome row is matched to the zero equation with its own id and t ...
  i <- match(y$y_FROM_x1_id_t, z$x1_TO_y_id_t)
  expect_false(anyNA(i))
  expect_equal(z$x1_t[i], dOut$t[-7])
  expect_equal(sub("^x1\\|1\\|", "", z$x1_1_id[i]), dOut$id[-7])
  # ... whose covariates and random-effect labels are those of x1 at that (id, t)
  expect_equal(z$x1_t, as.numeric(sub("^.*\\|", "", z$x1_TO_y_id_t)))
  expect_true(all(c(z$x1_1_id, z$x1_t_id) %in% r$levels$x1_1_t_id))  # S003 has RE levels
  expect_true("x1|1|S003" %in% r$levels$x1_1_t_id)

  # na.fail stops instead of dropping rows
  expect_error(prepareData(dOut, dLong, f, na.action = "na.fail"), "missing values")
})
