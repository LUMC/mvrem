# E38: latent equations, `latent(eta) ~ ... + scaled(x1)` linked to outcomes by scaled(eta).
# Fixtures: helper-latent.R.

# ---------------------------------------------------------------------------
# parsing / collectFormulas()

test_that("collectFormulas(): latent(name) on the left makes a latent equation, in either style", {
  fc <- collectFormulas(ltForms())
  expect_named(fc$latent, "eta")
  expect_equal(names(fc$out), c("z", "y"))
  expect_equal(names(fc$long), c("x1", "x2"))
  expect_equal(fc$order, c("z", "y", "eta", "x1", "x2"))
  expect_equal(fc$latent$eta$part, "latent")
  expect_identical(fc$latent$eta$variables$response, character())
  expect_equal(scaledTargets(fc$latent$eta), c("x1", "x2"))
  # the default options
  expect_s3_class(fc$latent$eta$spec, "mvremLatentSpec")
  expect_true(fc$latent$eta$spec$center)
  # pair style: a latent formula may be in either list
  fp <- collectFormulas(formOut = list(z ~ 1 + scaled(eta), latent(eta) ~ 1 + scaled(x1)),
                        formLong = list(x1 ~ 1 + (1 | id)))
  expect_named(fp$latent, "eta")
  expect_equal(fp$order, c("z", "eta", "x1"))
  # models without latent equations are unchanged
  expect_null(collectFormulas(list(y ~ scaled(x1), x1 ~ t))$latent)
})

test_that("collectFormulas(): an outcome's link variables include those of the latent's targets", {
  fc <- collectFormulas(list(y ~ 1 + scaled(eta), latent(eta) ~ 1 + (1 | id) + scaled(x1),
                             x1 ~ 1 + t + (1 | site)))
  expect_setequal(fc$out$y$linkVariables, c("id", "t", "site"))
  expect_setequal(linkedVariables(fc$latent$eta, fc$long), c("id", "t", "site"))
})

test_that("collectFormulas(): latent options come from family = list(eta = latentSpec())", {
  fc <- collectFormulas(ltForms(), family = list(eta = latentSpec(centerTargets = TRUE, logPrec = -15),
                                                 y = "gaussian"))
  expect_true(fc$latent$eta$spec$centerTargets)
  expect_equal(fc$latent$eta$spec$logPrec, -15)
  expect_equal(fc$out$y$family$family, "gaussian")
})

test_that("collectFormulas(): errors of latent equations", {
  expect_error(collectFormulas(list(y ~ scaled(eta), latent(eta) ~ 1 + scaled(x), x ~ 1),
                               family = list(eta = "gaussian")),
               "'eta' is a latent equation.*latentSpec")
  expect_error(collectFormulas(list(y ~ scaled(eta), latent(eta) ~ 1 + scaled(x), x ~ 1),
                               family = list(x = latentSpec())),
               "latentSpec\\(\\) in `family` for x, which is not a latent equation")
  expect_error(collectFormulas(list(y ~ scaled(x), latent(eta) ~ 1 + scaled(x), x ~ 1)),
               "Latent equation\\(s\\) not named in any scaled\\(\\) term of an outcome: eta")
  expect_error(collectFormulas(list(y ~ scaled(e1), latent(e1) ~ 1 + scaled(e2),
                                    latent(e2) ~ 1 + scaled(x), x ~ 1)),
               "can only be linked to longitudinal responses")
  expect_error(collectFormulas(list(y ~ scaled(eta), latent(eta, b) ~ 1 + scaled(x), x ~ 1)),
               "latent\\(\\) takes exactly one unnamed argument")
  expect_error(collectFormulas(list(y ~ scaled(eta), latent(eta) | weights(w) ~ 1 + scaled(x), x ~ 1)),
               "has no observations.*weights\\(\\)")
  expect_error(collectFormulas(list(y ~ scaled(eta), latent(eta) ~ 1 + scaled(x), x ~ 1),
                               familyOut = list(eta = latentSpec())),
               "`familyOut` names the latent equation")
  expect_error(latentSpec(logPrec = "a"), "`logPrec` must be NULL or one finite number")
  expect_error(latentSpec(center = NA), "center")
})

test_that("parseJsmFormula(): scaled(x, coef = c) fixes the coefficient", {
  p <- parseJsmFormula(y ~ 1 + scaled(x, coef = 2), "out")
  expect_equal(p$scaled[[1L]]$coef, 2)
  expect_length(p$scaled[[1L]]$params, 0L)
  expect_error(parseJsmFormula(y ~ scaled(x, coef = "a"), "out"), "`coef` in scaled\\(x\\) must be one finite number")
  expect_error(parseJsmFormula(y ~ scaled(x, coef = 1, prior = normalCoef()), "out"), "both `coef`")
})

# ---------------------------------------------------------------------------
# buildInlaFormula()

test_that("buildInlaFormula(): a latent equation is one shared field over its evaluation points", {
  m <- ltModel(fit = FALSE)
  inl <- m$inla
  pts <- inl$latentPoints$eta
  expect_equal(nrow(pts), nrow(unique(ltData()[c("id", "t")])))
  expect_named(pts, c("id", "t"))
  lf <- inl$latent$eta
  expect_equal(lf$index, "eta_LATENT_id_t")
  expect_equal(lf$logPrec, -6)
  expect_true(lf$center)
  args <- as.list(lf$term)
  expect_true(isTRUE(args$constr))
  expect_equal(eval(args$hyper), list(prec = list(initial = -6, fixed = TRUE)))
  # the stack: zero rows only, labels of the field, weight -1
  st <- inl$stacks$eta[[latentStackPart]]
  expect_equal(nrow(st), nrow(pts))
  expect_true(all(st$eta_eta_0eq == 0))
  expect_true(all(st$eta_WEIGHT == -1))
  expect_equal(st$eta_LATENT_id_t, latentLabels(pts, "eta", c("id", "t"), inl$naming))
  # both outcomes copy the same field, z with its coefficient fixed at 1
  expect_equal(as.character(as.list(inl$scaled$z_eta$term)$copy), "eta_LATENT_id_t")
  expect_equal(as.character(as.list(inl$scaled$y_eta$term)$copy), "eta_LATENT_id_t")
  expect_equal(eval(as.list(inl$scaled$z_eta$term)$hyper), list(beta = list(initial = 1, fixed = TRUE)))
  # markers: a zero field towards the latent, at its points
  expect_equal(nrow(inl$stacks$x1$eta), nrow(pts))
  expect_true(isTRUE(inl$zero$x1_eta$toLatent))
  expect_equal(as.character(as.list(inl$scaled$eta_x1$term)$copy), "x1_TO_eta_id_t")
  # one likelihood for the latent's zero rows (fixed precision), after the markers'
  expect_equal(names(inl$families), c("z_z", "y_y", "x1_x1", "x2_x2", "x1_x1_0eq", "x2_x2_0eq", "eta_eta_0eq"))
  expect_equal(inl$controlFamily$eta_eta_0eq, zeroControlFamily)
})

test_that("buildInlaFormula(): a subject without marker data gets marker labels in LEVELS (latent without RE)", {
  d <- ltData()
  dl <- d[d$id != "P60", ]
  f <- list(z ~ 1 + scaled(eta, coef = 1), latent(eta) ~ 1 + t + scaled(x1) + scaled(x2),
            x1 ~ 1 + t + (1 + t | id), x2 ~ 1 + t + (1 | id))
  inl <- jsm(f, data = list(z = d, x1 = dl, x2 = dl), fit = FALSE)$inla
  expect_true("P60" %in% inl$groups$id)
  # the RE labels of the markers' zero rows at the latent's points (a correlated block has one
  # LEVELS entry for all its effects, x1_1_t_id)
  for (r in c("x1", "x2")) {
    st <- inl$stacks[[r]]$eta
    for (i in grep("^x[12]_(1|t)_id$", names(st), value = TRUE))
      expect_true(all(st[[i]] %in% unlist(inl$levels)), info = i)
  }
  expect_true("x1|1|P60" %in% inl$levels$x1_1_t_id)
  expect_true("x2|1|P60" %in% inl$levels$x2_1_id)
})

test_that("buildInlaFormula(): latentSpec(logPrec, centerTargets) and jsm(zeroLogPrec)", {
  m <- jsm(ltForms(), data = ltData(), fit = FALSE, zeroLogPrec = -8,
           family = list(eta = latentSpec(center = FALSE, centerTargets = TRUE, logPrec = -15)))
  inl <- m$inla
  expect_equal(inl$zeroLogPrec, -8)
  expect_null(as.list(inl$latent$eta$term)$constr)
  expect_equal(eval(as.list(inl$latent$eta$term)$hyper)$prec$initial, -15)
  z <- as.list(inl$zero$x1_eta$term)
  expect_true(isTRUE(z$constr))
  expect_equal(eval(z$hyper)$prec$initial, -15)
  # a marker field towards an outcome takes the model's value
  m2 <- jsm(list(y ~ scaled(eta) + scaled(x2), latent(eta) ~ 1 + scaled(x1), x1 ~ t + (1 | id),
                 x2 ~ t + (1 | id)), data = ltData(), fit = FALSE, zeroLogPrec = -8)
  expect_equal(eval(as.list(m2$inla$zero$x2_y$term)$hyper)$prec$initial, -8)
  expect_null(as.list(m2$inla$zero$x2_y$term)$constr)
  # update() keeps it, also when it came from options(mvrem.zeroLogPrec)
  expect_equal(update(m2, fit = FALSE)$inla$zeroLogPrec, -8)
  m3 <- withr::with_options(list(mvrem.zeroLogPrec = -9),
                            jsm(ltForms(), data = ltData(), fit = FALSE))
  expect_equal(m3$call$zeroLogPrec, -9)
  expect_error(jsm(ltForms(), data = ltData(), fit = FALSE, zeroLogPrec = NA), "`zeroLogPrec` must be one finite")
})

test_that("buildInlaFormula(): a latent random effect cannot take a data-scaled prior", {
  f <- ltForms()
  f[[3L]] <- latent(eta) ~ 1 + t + re(1 | id, prior = pcSd(u = "data")) + scaled(x1) + scaled(x2)
  expect_error(jsm(f, data = ltData(), fit = FALSE), "cannot take a prior scaled from the data")
})

test_that("formatInlaFormula(), latexEquation(), print(): latent equations", {
  m <- ltModel(fit = FALSE)
  txt <- formatInlaFormula(m, nrows = 1)
  expect_true(any(grepl("latent equations: family = \"gaussian\"", txt, fixed = TRUE)))
  expect_true(any(grepl("f(eta_LATENT_id_t, eta_WEIGHT, model = \"iid\"", txt, fixed = TRUE)))
  tex <- unclass(latexEquation(m))
  expect_match(tex, "u^{(\\mathit{eta})}_{\\mathrm{id},t} &= ", fixed = TRUE)
  expect_match(tex, "\\gamma_{y,\\mathit{eta}}\\,u^{(\\mathit{eta})}_{\\mathrm{id}[i],t_{i}}", fixed = TRUE)
  expect_match(tex, "u^{(x_{1} \\to \\mathit{eta})}_{\\mathrm{id},t}", fixed = TRUE)
  expect_match(tex, "\\textstyle\\sum u^{(\\mathit{eta})} = 0", fixed = TRUE)
  out <- capture.output(print(m))
  expect_true(any(grepl("^Latent: +eta \\(no observations\\)", out)))
  expect_true(any(grepl("eta <- x1, x2", out, fixed = TRUE)))
  expect_false(any(grepl("set only by the priors", out)))   # z fixes the scale
  free <- jsm(ltForms(coef = FALSE), data = ltData(), fit = FALSE)
  expect_true(any(grepl("scale and sign of eta set only by the priors", capture.output(print(free)))))
  # no row for the latent in the data summary
  expect_false("eta" %in% m$dataSummary$equations$response)
  expect_true(all(c("eta", "eta") %in% m$dataSummary$zero$long))
})

test_that("plotModel(): a middle column for latent equations", {
  m <- ltModel(fit = FALSE)
  nodes <- diagramNodes(m$forms)
  expect_equal(vapply(nodes, `[[`, "", "side"), c(x1 = "L", x2 = "L", eta = "M", z = "R", y = "R"))
  links <- diagramLinks(c(m$forms$out, latentParts(m$forms)))
  expect_equal(paste(links$from, links$to), c("eta z", "eta y", "x1 eta", "x2 eta"))
  expect_equal(links$coef, c(1, NA, NA, NA))
  expect_silent(withNullDevice(plotModel(m)))
})

# ---------------------------------------------------------------------------
# fitted model

test_that("jsm(): the latent model recovers the simulated coefficients (scale fixed by coef = 1)", {
  m <- ltModel()
  truth <- attr(ltData(), "truth")
  a <- associations(m)
  expect_equal(a$level, c("latent -> outcome", "marker -> latent", "marker -> latent"))
  within <- function(row, value) row$q0.025 - 0.05 <= value && value <= row$q0.975 + 0.05
  expect_true(within(a[a$outcome == "y", ], truth$beta))
  expect_true(within(a[a$longitudinal == "x1", ], truth$gamma1))
  expect_true(within(a[a$longitudinal == "x2", ], truth$gamma2))
  fe <- fixedEffects(m)
  expect_true(within(fe[fe$response == "eta" & fe$term == "t", ], truth$slope))
  vc <- varianceComponents(m)$table
  expect_true(within(vc[vc$response == "eta" & vc$block != "residual", ], truth$sdV))
})

test_that("fittedValues(): a latent equation only when named, at its evaluation points", {
  m <- ltModel()
  fv <- fittedValues(m)
  expect_false("eta" %in% fv$response)
  le <- fittedValues(m, "eta")
  expect_equal(unique(le$response), "eta")
  expect_equal(nrow(le), nrow(m$inla$latentPoints$eta))
  expect_true(all(c("id", "t", "mean", "sd", "q0.025", "q0.975") %in% names(le)))
  expect_false("observed" %in% names(le))
  # centred (sum-to-zero constraint)
  expect_lt(abs(mean(le$mean)), 1e-3)
  expect_error(residuals(m, "eta"), "latent equation")
  expect_false("eta" %in% residuals(m)$response)
})

test_that("modelFit(): per observed response, all observations, and INLA's totals", {
  m <- ltModel()
  mf <- modelFit(m)
  expect_equal(mf$response, c("z", "y", "x1", "x2", "observations", "overall"))
  expect_equal(mf$waic[5], sum(mf$waic[1:4]), tolerance = 1e-8)
  expect_equal(mf$waic[6], inlaFit(m)$waic$waic)
  expect_true(is.finite(mf$mlik[6]))
})

test_that("predict(): a latent model refits ('sample' is not supported yet)", {
  m <- ltModel()
  nd <- data.frame(id = "NEW", t = c(0, 0.5), x1 = NA_real_, x2 = NA_real_)
  expect_error(predict(m, nd, response = "y", method = "sample"), "latent equations")
  p <- suppressMessages(predict(m, nd, response = "y", threads = 1))
  expect_equal(attr(p, "method"), "refit")
  expect_equal(nrow(p), 2L)
  expect_true(all(is.finite(p$mean)))
})
