d <- data.frame(id = sprintf("S%03d", rep(1:3, each = 4)), t = rep(c(0, 0.5, 1, 1.5), 3),
                g = factor(rep(c("u", "v"), 6)),
                x1 = rnorm(12), x2 = rnorm(12), x3 = rnorm(12), y = rnorm(12))
build <- function(formOut, formLong) {
  f <- collectFormulas(formOut, formLong)
  buildInlaFormula(f, prepareData(d, d, f))
}
hyperOf <- function(term) as.list(term)$hyper
lines <- function(x) {
  body <- sub("\n\\\\end\\{aligned\\}$", "", sub("^\\\\begin\\{aligned\\}\n", "", unclass(x)))
  strsplit(body, " \\\\\n", fixed = TRUE)[[1]]
}

test_that("constructors create mvremPrior objects and validate arguments", {
  p <- pcSd(u = 2)
  expect_s3_class(p, "mvremPrior")
  expect_equal(p$args, list(u = 2, alpha = 0.01))
  expect_equal(gammaPrec()$args, list(shape = 1, rate = 5e-5))
  expect_null(wishartPrec()$args$df)
  expect_error(pcSd(u = 1, alpha = 1.5), "`alpha` must be in \\(0, 1\\)")
  expect_error(pcSd(), "argument \"u\" is missing")
  expect_error(halfNormalSd(-1), "`scale` must be a single positive number")
  expect_error(halfCauchySd(c(1, 2)), "single positive")
  expect_error(gammaPrec(rate = 0), "`rate`")
  expect_error(wishartPrec(scale = matrix(c(1, 2, 0, 1), 2)), "symmetric positive-definite")
  expect_error(wishartPrec(scale = matrix(c(1, 2, 2, 1), 2)), "symmetric positive-definite")
  expect_error(wishartPrec(df = -1), "`df`")
})

test_that("print shows the distribution and the INLA mapping", {
  expect_output(print(pcSd(u = 1)),
                "PC prior on SD: P\\(sigma > 1\\) = 0.01, i.e. sigma ~ Exp\\(4.605\\).*pc.prec")
  expect_output(print(wishartPrec(df = 4, scale = diag(c(1, 0.1)))),
                "scale = diag\\(1, 0.1\\)", fixed = FALSE)
  expect_output(print(wishartPrec()), "df = k \\+ 2, scale = 1 \\* I")
  expect_output(print(halfNormalSd(0.5)), "logtnormal, param = c\\(0, 4\\)")
})

test_that("re() stores the prior separately from the f() parameters", {
  p <- parseJsmFormula(x1 ~ t + re(1 | id, prior = pcSd(u = 1), constr = TRUE), "long")
  expect_equal(p$re[[1]]$prior, pcSd(u = 1))
  expect_equal(p$re[[1]]$params, list(constr = TRUE))
  p <- parseJsmFormula(x1 ~ t + re(1 | id, prior = pcSd(u = 1)) + (0 + t | id), "long")
  expect_equal(p$re[[1]]$params, list())
  expect_null(p$re[[2]]$prior)
  pInt <- halfNormalSd(2)
  p <- parseJsmFormula(x1 ~ t + re(1 | id, prior = pInt), "long")
  expect_equal(p$re[[1]]$prior, pInt)
})

test_that("invalid prior arguments in re() / scaled() are errors", {
  expect_error(parseJsmFormula(x1 ~ t + re(1 | id, prior = list(prior = "pc.prec")), "long"),
               "must be created by a prior constructor")
  expect_error(collectFormulas(y ~ scaled(x1, prior = pcSd(1)), x1 ~ t + (1 | id)),
               "not supported in scaled\\(\\)")
})

test_that("scalar priors give the INLA hyper of every iid term of the block", {
  inl <- build(y ~ t + re(1 | id, prior = halfCauchySd(25)),
               list(x1 ~ t + re(1 | id, prior = gammaPrec(1, 0.01)),
                    x2 ~ t + re(1 + t || id, prior = pcSd(u = 0.5, alpha = 0.05)),
                    x3 ~ t + re(0 + t | id, prior = halfNormalSd(0.5))))
  expect_equal(hyperOf(inl$random$x1[[1]]$terms[[1]]),
               list(prec = list(prior = "loggamma", param = c(1, 0.01))))
  pc <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.05)))
  expect_equal(hyperOf(inl$random$x2[[1]]$terms[[1]]), pc)
  expect_equal(hyperOf(inl$random$x2[[2]]$terms[[1]]), pc)
  expect_equal(hyperOf(inl$random$x3[[1]]$terms[[1]]),
               list(prec = list(prior = "logtnormal", param = c(0, 4))))
  hc <- hyperOf(inl$random$y[[1]]$terms[[1]])$prec
  expect_named(hc, "prior")
  expect_match(hc$prior, "^expression: sigma = exp\\(-theta/2\\); s = 25;")
  expect_equal(deparse1(inl$random$x1[[1]]$terms[[1]]),
               paste0('f(x1_1_id, model = "iid", values = LEVELS$x1_1_id, n = length(LEVELS$x1_1_id), ',
                      'hyper = list(prec = list(prior = "loggamma", param = c(1, 0.01))))'))
})

test_that("wishartPrec gives wishart<k>d with R = scale^-1; copies get no hyper", {
  inl <- build(y ~ t, list(
    x1 ~ t + re(1 + t | id, prior = wishartPrec(df = 4, scale = diag(c(1, 0.1)))),
    x2 ~ t + I(t^2) + re(1 + t + I(t^2) | id, prior = wishartPrec(df = 5, scale = 0.5)),
    x3 ~ t + re(1 + t | id, prior = wishartPrec())))
  b <- inl$random$x1[[1]]
  expect_equal(hyperOf(b$terms[[1]]), list(theta1 = list(prior = "wishart2d", param = c(4, 1, 10, 0))))
  expect_null(hyperOf(b$terms[[2]]))
  expect_equal(hyperOf(inl$random$x2[[1]]$terms[[1]]),
               list(theta1 = list(prior = "wishart3d", param = c(5, 2, 2, 2, 0, 0, 0))))
  expect_equal(hyperOf(inl$random$x3[[1]]$terms[[1]])$theta1$param, c(4, 1, 1, 0))  # df = k + 2
})

test_that("wishart parameters: diagonal of R, then R12, R13, ..., R23, ...", {
  V <- matrix(c(1, .3, .1, .3, .5, .2, .1, .2, .4), 3)
  R <- solve(V)
  h <- mvrem:::priorHyper(wishartPrec(df = 5, scale = V), 3L)
  expect_equal(h$theta1$param, c(5, R[1, 1], R[2, 2], R[3, 3], R[1, 2], R[1, 3], R[2, 3]))
  h6 <- mvrem:::priorHyper(wishartPrec(scale = diag(6)), 6L)
  expect_equal(h6$theta1$prior, "wishartkd")
  expect_length(h6$theta1$param, 1 + 6 * 7 / 2)
  expect_equal(h6$theta1$param[1], 8)
})

test_that("priors are checked against the block when building", {
  expect_error(build(y ~ t, x1 ~ t + re(1 + t | id, prior = pcSd(1))),
               "Correlated block `1 \\+ t \\| id` \\(k = 2\\) needs wishartPrec\\(\\)")
  expect_error(build(y ~ t, x1 ~ t + re(1 | id, prior = wishartPrec())),
               "wishartPrec\\(\\) needs a correlated block")
  expect_error(build(y ~ t, x1 ~ t + re(1 + t || id, prior = wishartPrec())),
               "needs a correlated block")
  expect_error(build(y ~ t, x1 ~ t + re(1 + t | id, prior = wishartPrec(df = 3))),
               "df must be > k \\+ 1 = 3")
  # factor: k = 2 effects from the data, so a 3 x 3 scale is wrong
  expect_error(build(y ~ t, x1 ~ t + re(1 + g | id, prior = wishartPrec(scale = diag(3)))),
               "must be a 2 x 2 matrix")
  expect_silent(build(y ~ t, x1 ~ t + re(1 + g | id, prior = wishartPrec(scale = diag(2)))))
})

test_that("other f() arguments are merged; hyper overriding the prior warns", {
  inl <- build(y ~ t, x1 ~ t + re(1 | id, prior = pcSd(u = 1), constr = TRUE,
                                  hyper = list(prec = list(initial = 4))))
  tm <- as.list(inl$random$x1[[1]]$terms[[1]])
  expect_true(tm$constr)
  expect_equal(tm$hyper, list(prec = list(prior = "pc.prec", param = c(1, 0.01), initial = 4)))
  expect_warning(
    inl <- build(y ~ t, x1 ~ t + re(1 | id, prior = pcSd(u = 1),
                                    hyper = list(prec = list(prior = "loggamma", param = c(1, 1))))),
    "`hyper` overrides the `prior` of block `1 \\| id` \\(prec\\)")
  expect_equal(hyperOf(inl$random$x1[[1]]$terms[[1]]),
               list(prec = list(prior = "loggamma", param = c(1, 1))))
})

test_that("blocks without prior keep INLA's default", {
  inl <- build(y ~ t, x1 ~ t + (1 + t | id))
  expect_null(hyperOf(inl$random$x1[[1]]$terms[[1]]))
})

test_that("formula with priors is accepted by INLA", {
  skip_if_not_installed("INLA")
  inl <- build(y ~ t + re(1 | id, prior = halfCauchySd(5)) + scaled(x1),
               x1 ~ t + re(1 + t | id, prior = wishartPrec(df = 4, scale = diag(c(1, 0.1)))))
  i <- INLA:::inla.interpret.formula(inl$formula, data.same.len = FALSE, data = NULL)
  expect_true("iid2d" %in% vapply(i$random.spec, `[[`, "", "model"))
})

test_that("LaTeX: a prior line after the block's distribution", {
  f <- collectFormulas(
    y ~ t + re(1 | id, prior = halfCauchySd(5)) + scaled(x1),
    list(x1 ~ t + re(1 + t | id, prior = wishartPrec(df = 4, scale = diag(c(1, 0.1)))),
         x2 ~ t + re(1 + t || id, prior = gammaPrec(1, 0.01)),
         x3 ~ t + re(1 | id, prior = pcSd(u = 1)) + re(0 + t | id, prior = halfNormalSd(0.5))))
  dat <- prepareData(d, d, f)
  l1 <- lines(latexEquation(f$long$x1, dat$long$x1))
  i <- grep("\\mathbf{Q}^{(x_{1})}_{\\mathrm{id}}\\right)^{-1}", l1, fixed = TRUE)
  expect_equal(l1[i + 1L], paste0("\\mathbf{Q}^{(x_{1})}_{\\mathrm{id}} &\\sim \\mathcal{W}_{2}\\left(4, ",
                                  "\\begin{pmatrix}1 & 0 \\\\ 0 & 0.1\\end{pmatrix}\\right)"))
  l2 <- lines(latexEquation(f$long$x2))
  expect_equal(l2[length(l2)],
               "\\tau^{(x_{2})}_{k,\\mathrm{id}} &\\sim \\mathrm{Gamma}\\left(1, 0.01\\right), \\quad k = 0, 1")
  tex3 <- latexEquation(f$long$x3)
  expect_match(tex3, paste0("\\left(\\tau^{(x_{3})}_{0,\\mathrm{id}}\\right)^{-1/2} &\\sim ",
                            "\\mathrm{Exp}\\left(4.61\\right), \\quad P\\left(\\left(\\tau^{(x_{3})}_{0,",
                            "\\mathrm{id}}\\right)^{-1/2} > 1\\right) = 0.01"), fixed = TRUE)
  expect_match(tex3, "\\left(\\tau^{(x_{3})}_{1,\\mathrm{id}}\\right)^{-1/2} &\\sim \\mathcal{N}^{+}\\left(0, \\left(4\\right)^{-1}\\right)",
               fixed = TRUE)
  whole <- latexEquation(f, dat)
  expect_match(whole, "\\left(\\tau^{(y)}_{0,\\mathrm{id}}\\right)^{-1/2} &\\sim \\mathrm{C}^{+}\\left(0, 5\\right)",
               fixed = TRUE)
  expect_false(grepl("sigma", whole))
  # scalar wishart scale and default df
  f2 <- collectFormulas(y ~ t, x1 ~ t + re(1 + t | id, prior = wishartPrec(scale = 0.5)))
  expect_match(latexEquation(f2$long$x1), "\\mathcal{W}_{2}\\left(4, 0.5\\,\\mathbf{I}_{2}\\right)",
               fixed = TRUE)
  # with data the prior is checked against the block, as in buildInlaFormula()
  f3 <- collectFormulas(y ~ t, x1 ~ t + re(1 + t | id, prior = pcSd(1)))
  expect_error(latexEquation(f3$long$x1, prepareData(d, d, f3)$long$x1), "needs wishartPrec")
})

test_that("numbers in LaTeX", {
  expect_equal(mvrem:::numLatex(5e-5), "5 \\cdot 10^{-5}")
  expect_equal(mvrem:::numLatex(1e-6), "10^{-6}")
  expect_equal(mvrem:::numLatex(4.60517), "4.61")
  expect_equal(mvrem:::numLatex(250), "250")
})
