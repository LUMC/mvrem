# Tests of R/latex.R (latexEquation() and its helpers) and of the LaTeX helpers of
# R/priors.R (priorLatex(), coefPriorLatex(), numLatex(), matrixLatex()). Helpers in
# helper-latex.R. The LaTeX of the families and of the priors is covered here (not in
# test-families.R / test-priors.R).

# ---------------------------------------------------------------------------
# Helper sanity

test_that("helper-latex.R: brace, wrapper and \\htmlData helpers work", {
  expect_true(bracesBalanced("a_{b} \\{0, 1\\} \\\\ c"))
  expect_false(bracesBalanced("a_{b"))
  expect_false(bracesBalanced("}{"))
  expect_true(leftRightBalanced("\\left(a\\right) \\rightarrow"))
  expect_equal(unwrapCommand("\\htmlData{a=1}{x_{\\htmlData{b=2}{i}}} + y", "\\htmlData"), "x_{i} + y")
  expect_equal(unwrapCommand("\\underset{\\substack{1 \\\\ 2}}{\\beta_{0}}", "\\underset"), "\\beta_{0}")
  w <- htmlDataWrappers("\\htmlData{a=1}{x_{\\htmlData{b=2}{i}}} + y")
  expect_equal(w$arg, c("a=1", "b=2"))
  expect_equal(w$sym, c("x_{\\htmlData{b=2}{i}}", "i"))
  expect_equal(nrow(htmlDataWrappers("x + y")), 0L)
})

# ---------------------------------------------------------------------------
# Families: link + distribution line (familyLatex())

test_that("familyLatex(): every family writes its link and its distribution line", {
  cases <- list(
    list(y ~ t, betaFamily(),
         c(r"(\operatorname{logit}\left(\mu^{(y)}_{i}\right) &= )",
           r"(y_{i} &\sim \mathrm{Beta}\left(\mu^{(y)}_{i}\, \phi^{(y)}, \left(1 - \mu^{(y)}_{i}\right) \phi^{(y)}\right))")),
    list(k | trials(n) ~ t, binomialFamily(),
         c(r"(\operatorname{logit}\left(p^{(k)}_{i}\right) &= )",
           r"(k_{i} &\sim \mathrm{Binomial}\left(n_{i}, p^{(k)}_{i}\right))")),
    list(k ~ t, binomialFamily(link = "probit"),
         c(r"(\Phi^{-1}\left(p^{(k)}_{i}\right) &= )",
           r"(k_{i} &\sim \mathrm{Bernoulli}\left(p^{(k)}_{i}\right))")),
    list(k ~ t, binomialFamily(link = "cloglog"),
         r"(\operatorname{cloglog}\left(p^{(k)}_{i}\right) &= )"),
    list(cnt | exposure(pt) ~ t, poissonFamily(),
         c(r"(\log\left(\lambda^{(\mathit{cnt})}_{i}\right) &= )",
           r"(\mathit{cnt}_{i} &\sim \mathrm{Poisson}\left(\mathit{pt}_{i}\, \lambda^{(\mathit{cnt})}_{i}\right))")),
    list(cnt ~ t, poissonFamily(),
         r"(\mathit{cnt}_{i} &\sim \mathrm{Poisson}\left(\lambda^{(\mathit{cnt})}_{i}\right))"),
    list(cnt | exposure(pt) ~ t, negBinomialFamily(),
         c(r"(\log\left(\mu^{(\mathit{cnt})}_{i}\right) &= )",
           r"(\mathrm{NegBin}\left(\mathit{pt}_{i}\, \mu^{(\mathit{cnt})}_{i}, \kappa^{(\mathit{cnt})}\right))")),
    list(y ~ t, gammaFamily(),
         c(r"(\log\left(\mu^{(y)}_{i}\right) &= )",
           r"(y_{i} &\sim \mathrm{Gamma}\left(\phi^{(y)}, \phi^{(y)} / \mu^{(y)}_{i}\right))")),
    list(y ~ t, lognormalFamily(),
         c(r"(\log\left(y_{i}\right) &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + \varepsilon^{(y)}_{i})",
           r"(\varepsilon^{(y)}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(y)}_{\varepsilon}\right)^{-1}\right))")),
    list(y ~ t, studentTFamily(),
         c(r"(y_{i} &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + \varepsilon^{(y)}_{i})",
           r"(\varepsilon^{(y)}_{i} &\sim t_{\nu^{(y)}}\left(0, \left(\tau^{(y)}_{\varepsilon}\right)^{-1}\right))")),
    list(y ~ t, gaussianFamily(),
         c(r"(y_{i} &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + \varepsilon^{(y)}_{i})",
           r"(\varepsilon^{(y)}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(y)}_{\varepsilon}\right)^{-1}\right))")),
    list(y ~ t, gaussianFamily(link = "log"),
         c(r"(\log\left(\mu^{(y)}_{i}\right) &= )",
           r"(y_{i} &\sim \mathcal{N}\left(\mu^{(y)}_{i}, \left(\tau^{(y)}_{\varepsilon}\right)^{-1}\right))")),
    list(k | trials(n) ~ t, betaBinomialFamily(),
         r"(k_{i} &\sim \mathrm{BetaBinomial}\left(n_{i}, p^{(k)}_{i}, \rho^{(k)}\right))"),
    list(k ~ t, betaBinomialFamily(),
         r"(k_{i} &\sim \mathrm{BetaBinomial}\left(1, p^{(k)}_{i}, \rho^{(k)}\right))"),
    list(cnt ~ t, zeroInflatedPoissonFamily(type = 0),
         r"(\pi^{(\mathit{cnt})}\,\delta_{0} + \left(1 - \pi^{(\mathit{cnt})}\right) \mathrm{Poisson}^{+}\left(\lambda^{(\mathit{cnt})}_{i}\right))"),
    list(cnt ~ t, zeroInflatedPoissonFamily(type = 1),
         r"(\pi^{(\mathit{cnt})}\,\delta_{0} + \left(1 - \pi^{(\mathit{cnt})}\right) \mathrm{Poisson}\left(\lambda^{(\mathit{cnt})}_{i}\right))"),
    list(cnt | exposure(pt) ~ t, zeroInflatedNegBinomialFamily(type = 0),
         r"(\mathrm{NegBin}^{+}\left(\mathit{pt}_{i}\, \mu^{(\mathit{cnt})}_{i}, \kappa^{(\mathit{cnt})}\right))"),
    list(cnt ~ t, zeroInflatedNegBinomialFamily(type = 1),
         r"(\left(1 - \pi^{(\mathit{cnt})}\right) \mathrm{NegBin}\left(\mu^{(\mathit{cnt})}_{i}, \kappa^{(\mathit{cnt})}\right))"),
    list(time | event(status) ~ 1, weibullSurvFamily(),
         c(r"(\log\left(\lambda^{(\mathit{time})}_{i}\right) &= \beta^{(\mathit{time})}_{0})",
           r"(\mathrm{Weibull}\left(\alpha^{(\mathit{time})}, \lambda^{(\mathit{time})}_{i}\right), \quad S(t) = \exp\left(-\lambda^{(\mathit{time})}_{i}\, t^{\alpha^{(\mathit{time})}}\right), \quad \mathit{status}_{i} \in \{0, 1\})")),
    list(time | event(status) ~ 1, exponentialSurvFamily(),
         r"(\mathit{time}_{i} &\sim \mathrm{Exp}\left(\lambda^{(\mathit{time})}_{i}\right), \quad S(t) = \exp\left(-\lambda^{(\mathit{time})}_{i}\, t\right), \quad \mathit{status}_{i} \in \{0, 1\})"),
    list(y ~ t, inlaFamily("sn"),
         c(r"(g_{y}\left(\mathrm{E}\left[y_{i}\right]\right) &= )", r"(y_{i} &\sim \mathrm{sn})"))
  )
  for (cs in cases) {
    tex <- famTex(cs[[1]], cs[[2]])
    for (s in cs[[3]]) expect_match(tex, s, fixed = TRUE, info = deparse(cs[[1]]))
    expectWellFormed(tex)
  }
})

test_that("familyLatex(): only families with a residual term write an epsilon", {
  for (fam in list(gaussianFamily(), lognormalFamily(), studentTFamily()))
    expect_match(famTex(y ~ t, fam), r"(\varepsilon^{(y)}_{i})", fixed = TRUE, info = fam$family)
  for (fam in list(betaFamily(), poissonFamily(), binomialFamily(), gammaFamily(), weibullSurvFamily(),
                   gaussianFamily(link = "log")))                   # log link: N(mu, tau^-1), no epsilon term
    expect_false(grepl(r"(\varepsilon^{(y)}_{i})", famTex(y ~ t, fam), fixed = TRUE), info = fam$family)
  # zero inflation of type 1 is not truncated
  expect_false(grepl("^{+}", famTex(cnt ~ t, zeroInflatedPoissonFamily(1)), fixed = TRUE))
})

test_that("familyLatex(): a family without a prior writes exactly one distribution line", {
  fams <- list(betaFamily(), binomialFamily(), poissonFamily(), negBinomialFamily(), gammaFamily(),
               zeroInflatedPoissonFamily(0), betaBinomialFamily(), inlaFamily("sn"))
  for (fam in fams) {
    l <- texLines(famTex(y ~ t, fam))
    expect_equal(sum(grepl("&\\sim", l, fixed = TRUE)), 1L, info = fam$family)
    expect_match(l[1], "&= ", fixed = TRUE)
  }
})

test_that("familyLatex(): a transformed response is written as its expression", {
  expect_match(famTex(log(y) ~ t, lognormalFamily()), r"(\log\left(\log\left(y_{i}\right)\right) &=)", fixed = TRUE)
  tex <- famTex(log(y) ~ t, betaFamily())
  expect_match(tex, r"(\operatorname{logit}\left(\mu^{(y)}_{i}\right) &=)", fixed = TRUE)
  expect_match(tex, r"(\log\left(y_{i}\right) &\sim \mathrm{Beta}\left()", fixed = TRUE)
})

test_that("familyLatex(): a family hyperparameter prior adds one line right after the distribution", {
  cases <- list(
    list(betaFamily(phi = gammaPrec(1, 0.1)), r"(\phi^{(y)} &\sim \mathrm{Gamma}\left(1, 0.1\right))"),
    list(gaussianFamily(prec = pcSd(1)),                   # Exp rate = -log(0.01) / 1 = 4.61
         r"(\left(\tau^{(y)}_{\varepsilon}\right)^{-1/2} &\sim \mathrm{Exp}\left(4.61\right), \quad P\left(\left(\tau^{(y)}_{\varepsilon}\right)^{-1/2} > 1\right) = 0.01)"),
    list(studentTFamily(prec = halfNormalSd(2)),           # precision 1 / 2^2
         r"(\left(\tau^{(y)}_{\varepsilon}\right)^{-1/2} &\sim \mathcal{N}^{+}\left(0, \left(0.25\right)^{-1}\right))"),
    list(gammaFamily(prec = halfCauchySd(3)),
         r"(\left(\phi^{(y)}\right)^{-1/2} &\sim \mathrm{C}^{+}\left(0, 3\right))"),
    list(negBinomialFamily(size = gammaPrec(2, 0.5)), r"(\kappa^{(y)} &\sim \mathrm{Gamma}\left(2, 0.5\right))"),
    list(inlaFamily("sn", hyper = list(skew = gammaPrec())),
         r"(\theta^{(y)}_{\mathrm{skew}} &\sim \mathrm{Gamma}\left(1, 5 \cdot 10^{-5}\right))")
  )
  for (cs in cases) {
    l <- texLines(famTex(y ~ t, cs[[1]]))
    iDist <- grep("^[^&]*&\\\\sim ", l)[1]
    iPrior <- grep(cs[[2]], l, fixed = TRUE)
    expect_equal(iPrior, iDist + 1L, info = cs[[1]]$family)
    expect_equal(sum(grepl("&\\sim", l, fixed = TRUE)), 2L, info = cs[[1]]$family)
  }
})

test_that("familyLatex(): a raw INLA hyper list gets no prior line", {
  l <- texLines(famTex(y ~ t, betaFamily(phi = list(prior = "loggamma", param = c(1, 0.1)))))
  expect_equal(sum(grepl("&\\sim", l, fixed = TRUE)), 1L)
  expect_false(any(grepl("Gamma", l, fixed = TRUE)))
})

test_that("latexEquation(): a parsed formula without family: outcome = generic link, long = gaussian", {
  out <- latexEquation(parseJsmFormula(y ~ t + scaled(a), part = "out"))
  expect_match(out, r"(g_{y}\left(\mathrm{E}\left[y_{i}\right]\right) &=)", fixed = TRUE)
  expect_false(grepl("\\sim", out, fixed = TRUE))
  long <- latexEquation(parseJsmFormula(y ~ t, part = "long"))
  expect_match(long, r"(\varepsilon^{(y)}_{i} &\sim \mathcal{N})", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Random-effect blocks and their priors (priorLatex())

test_that("latexEquation(): a `|` block with several effects is a multivariate normal with Q", {
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + (1 + t | id), part = "long")))
  mvn <- texLine(l, r"(\left(b^{(y)}_{0,l}, b^{(y)}_{1,l}\right)^\top &\sim)")
  expect_match(mvn, r"(\mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(y)}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}})",
               fixed = TRUE)
  expect_false(any(grepl("b^{(y)}_{0,l} &", l, fixed = TRUE)))
  # three effects: a 3-tuple
  tex <- latexEquation(parseJsmFormula(y ~ t + re(1 + t + I(t^2) | id, prior = wishartPrec()), part = "long"))
  expect_match(tex, r"(\left(b^{(y)}_{0,l}, b^{(y)}_{1,l}, b^{(y)}_{2,l}\right)^\top)", fixed = TRUE)
})

test_that("latexEquation(): single-effect and `||` blocks are univariate normals", {
  l <- texLines(latexEquation(parseJsmFormula(y ~ (1 | id), part = "long")))
  expect_match(texLine(l, "b^{(y)}_{0,l} &"), r"(b^{(y)}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(y)}_{0,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}})",
               fixed = TRUE)
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + (1 + t || id), part = "long")))
  expect_match(texLine(l, "b^{(y)}_{0,l} &"), r"(\left(\tau^{(y)}_{0,\mathrm{id}}\right)^{-1})", fixed = TRUE)
  expect_match(texLine(l, "b^{(y)}_{1,l} &"), r"(\left(\tau^{(y)}_{1,\mathrm{id}}\right)^{-1})", fixed = TRUE)
  expect_false(any(grepl("mathbf{Q}", l, fixed = TRUE)))
})

test_that("latexEquation(): without a block prior there is no prior line", {
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + (1 | id) + (0 + t | id) + (1 + t | g), part = "long")))
  dist <- grep("&\\sim", l, fixed = TRUE, value = TRUE)
  expect_length(dist, 4L)                                    # residual + two univariate + one MVN
  expect_false(any(grepl("Gamma|Exp|mathcal\\{W\\}|\\^\\{\\+\\}", dist)))
})

test_that("latexEquation(): a scalar block prior is one line right after the block's distribution", {
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 | id, prior = gammaPrec(2, 3)) +
                                                re(0 + t | id, prior = pcSd(0.5)), part = "long")))
  iB0 <- grep("b^{(y)}_{0,l} &", l, fixed = TRUE)
  iB1 <- grep("b^{(y)}_{1,l} &", l, fixed = TRUE)
  expect_match(l[iB0 + 1L], r"(\tau^{(y)}_{0,\mathrm{id}} &\sim \mathrm{Gamma}\left(2, 3\right))", fixed = TRUE)
  # pcSd(0.5): Exp rate = -log(0.01) / 0.5 = 9.21
  expect_match(l[iB1 + 1L], r"(\left(\tau^{(y)}_{1,\mathrm{id}}\right)^{-1/2} &\sim \mathrm{Exp}\left(9.21\right), \quad P\left(\left(\tau^{(y)}_{1,\mathrm{id}}\right)^{-1/2} > 0.5\right) = 0.01)",
               fixed = TRUE)
  expect_false(any(grepl("quad k =", l, fixed = TRUE)))      # single effects: no k suffix
})

test_that("latexEquation(): a `||` block shares one prior line over its effects (k suffix)", {
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 + t || id, prior = halfNormalSd(1)), part = "long")))
  expect_match(l[length(l)], r"(\left(\tau^{(y)}_{k,\mathrm{id}}\right)^{-1/2} &\sim \mathcal{N}^{+}\left(0, \left(1\right)^{-1}\right), \quad k = 0, 1)",
               fixed = TRUE)
  expect_equal(sum(grepl("mathcal{N}^{+}", l, fixed = TRUE)), 1L)
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 + t || id, prior = halfCauchySd(2)), part = "long")))
  expect_match(l[length(l)], r"(\left(\tau^{(y)}_{k,\mathrm{id}}\right)^{-1/2} &\sim \mathrm{C}^{+}\left(0, 2\right), \quad k = 0, 1)",
               fixed = TRUE)
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 + t || id, prior = gammaPrec(1, 0.01)), part = "long")))
  expect_match(l[length(l)], r"(\tau^{(y)}_{k,\mathrm{id}} &\sim \mathrm{Gamma}\left(1, 0.01\right), \quad k = 0, 1)", fixed = TRUE)
})

test_that("priorLatex(): Wishart line: default df = k + 2, matrix / scalar / data scales", {
  wishartLine <- function(prior) {
    l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 + t | id, prior = prior), part = "long")))
    texLine(l, r"(\mathbf{Q}^{(y)}_{\mathrm{id}} &\sim)")
  }
  expect_match(wishartLine(wishartPrec()), r"(\mathcal{W}_{2}\left(4, \mathbf{I}_{2}\right))", fixed = TRUE)
  expect_match(wishartLine(wishartPrec(df = 5, scale = diag(c(1, 0.1)))),
               r"(\mathcal{W}_{2}\left(5, \begin{pmatrix}1 & 0 \\ 0 & 0.1\end{pmatrix}\right))", fixed = TRUE)
  expect_match(wishartLine(wishartPrec(scale = 2)), r"(\left(4, 2\,\mathbf{I}_{2}\right))", fixed = TRUE)
  expect_match(wishartLine(wishartPrec(scale = "data")), r"(\left(4, \hat{\mathbf{V}}_{\mathrm{data}}\right))", fixed = TRUE)
  # k = 3: W_3 with df 5
  tex <- latexEquation(parseJsmFormula(y ~ t + re(1 + t + I(t^2) | id, prior = wishartPrec()), part = "long"))
  expect_match(tex, r"(\mathcal{W}_{3}\left(5, \mathbf{I}_{3}\right))", fixed = TRUE)
  # the Wishart line follows the MVN line
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + re(1 + t | id, prior = wishartPrec()), part = "long")))
  expect_equal(grep("mathcal{W}", l, fixed = TRUE), grep("^\\\\left\\(b", l) + 1L)
})

test_that("priorLatex(): unit cases: NULL, k suffix only for several effects, Wishart", {
  expect_identical(priorLatex(NULL, 1L, tau = "t", Q = "Q", idx = 0L), character())
  expect_equal(priorLatex(gammaPrec(1, 2), 1L, tau = "T", Q = NULL, idx = 3L), r"(T &\sim \mathrm{Gamma}\left(1, 2\right))")
  expect_match(priorLatex(gammaPrec(1, 2), 3L, tau = "T", Q = NULL, idx = c(0L, 2L, 5L)), r"(\quad k = 0, 2, 5$)")
  expect_equal(priorLatex(wishartPrec(df = 6, scale = 0.5), 3L, tau = NULL, Q = "Q", idx = 0:2),
               r"(Q &\sim \mathcal{W}_{3}\left(6, 0.5\,\mathbf{I}_{3}\right))")
})

test_that("matrixLatex(): identity, scaled identity and pmatrix", {
  expect_equal(matrixLatex(1, 2L), r"(\mathbf{I}_{2})")
  expect_equal(matrixLatex(0.5, 3L), r"(0.5\,\mathbf{I}_{3})")
  expect_equal(matrixLatex(matrix(c(1, 2, 2, 5e-5), 2), 2L),
               r"(\begin{pmatrix}1 & 2 \\ 2 & 5 \cdot 10^{-5}\end{pmatrix})")
})

test_that("latexEquation(): wishartPrec(scale = 'data') is symbolic without data, the lme4 matrix with data", {
  f <- fxForms("beta")
  f$long$x1$re[[1]]$prior <- wishartPrec(scale = "data")
  dat <- prepareData(f, fxData()$long)
  without <- latexEquation(f)
  with <- latexEquation(f, dat)
  expect_match(without, r"(\hat{\mathbf{V}}_{\mathrm{data}})", fixed = TRUE)
  expect_false(grepl("hat{", with, fixed = TRUE))
  # V = S^-1 / df with S the lme4 covariance of the block and df = k + 2 = 4
  S <- as.matrix(lme4::VarCorr(lme4::lmer(x1 ~ t + (1 + t | id), fxData()$long))$id)
  V <- solve(S) / 4
  sig3 <- function(v) numText(v, 3)
  pm <- paste0("\\begin{pmatrix}", sig3(V[1, 1]), " & ", sig3(V[1, 2]), " \\\\ ", sig3(V[2, 1]), " & ",
               sig3(V[2, 2]), "\\end{pmatrix}")
  expect_match(with, paste0(r"(\mathbf{Q}^{(x_{1})}_{\mathrm{id}} &\sim \mathcal{W}_{2}\left(4, )", pm, r"(\right))"), fixed = TRUE)
  # the single formula with its data frame gives the same matrix
  expect_match(latexEquation(f$long$x1, dat$long$x1), pm, fixed = TRUE)
  expectWellFormed(with)
})

test_that("latexEquation(): block priors are checked against the block only with data", {
  d <- fxSmallData()
  p <- parseJsmFormula(y ~ t + re(1 + t | id, prior = gammaPrec()), part = "long")
  expect_match(latexEquation(p), r"(\tau^{(y)}_{k,\mathrm{id}} &\sim \mathrm{Gamma})", fixed = TRUE)
  expect_error(latexEquation(p, d), "needs wishartPrec")
  p <- parseJsmFormula(y ~ t + re(1 | id, prior = wishartPrec()), part = "long")
  expect_error(latexEquation(p, d), "correlated block")
})

# ---------------------------------------------------------------------------
# Association priors (normalCoef() on scaled())

test_that("coefPriorLatex(): normalCoef() -> N(mean, (1 / sd^2)^{-1}); NULL -> nothing", {
  expect_equal(coefPriorLatex(normalCoef(-0.5, 0.1), "g"), r"(g &\sim \mathcal{N}\left(-0.5, \left(100\right)^{-1}\right))")
  expect_equal(coefPriorLatex(normalCoef(), "g"), r"(g &\sim \mathcal{N}\left(0, \left(0.01\right)^{-1}\right))")
  expect_identical(coefPriorLatex(NULL, "g"), character())
})

test_that("latexEquation(): the gamma prior line follows the outcome's distribution; none for the default", {
  p <- parseJsmFormula(y ~ t + scaled(x1, prior = normalCoef(1, 2)) + scaled(x2), part = "out")
  p$family <- betaFamily()
  l <- texLines(latexEquation(p))
  iBeta <- grep("Beta", l, fixed = TRUE)
  expect_equal(grep(r"(\gamma_{y,x_{1}} &\sim)", l, fixed = TRUE), iBeta + 1L)
  expect_match(l[iBeta + 1L], r"(\gamma_{y,x_{1}} &\sim \mathcal{N}\left(1, \left(0.25\right)^{-1}\right))", fixed = TRUE)
  expect_false(any(grepl(r"(\gamma_{y,x_{2}} &\sim)", l, fixed = TRUE)))
  # in a whole model: after the outcome's family line, before the u lines
  f <- collectFormulas(y ~ t + scaled(x, prior = normalCoef(0, 1)), x ~ t + (1 | id), family = list(y = betaFamily()))
  l <- texLines(latexEquation(f))
  iPrior <- grep(r"(\gamma_{y,x} &\sim)", l, fixed = TRUE)
  expect_equal(iPrior, grep("Beta", l, fixed = TRUE) + 1L)
  expect_true(all(grep("^u\\^", l) > iPrior))
})

# ---------------------------------------------------------------------------
# numLatex()

test_that("numLatex(): significant digits and powers of ten", {
  expect_equal(numLatex(5e-5), r"(5 \cdot 10^{-5})")
  expect_equal(numLatex(0), "0")
  expect_equal(numLatex(-2.5), "-2.5")
  expect_equal(numLatex(1), "1")
  expect_equal(numLatex(250), "250")
  expect_equal(numLatex(0.001234567), "0.00123")
  expect_equal(numLatex(0.001234567, 2), "0.0012")
  expect_equal(numLatex(123456), r"(1.23 \cdot 10^{5})")
  expect_equal(numLatex(123456, 2), r"(1.2 \cdot 10^{5})")
  expect_equal(numLatex(1e6), "10^{6}")
  expect_equal(numLatex(1e-10), "10^{-10}")
  expect_equal(numLatex(-3e-7), r"(-3 \cdot 10^{-7})")
  expect_equal(numLatex(4.60517), "4.61")
  expect_equal(numLatex(4.60517, 5), "4.6052")
  for (v in c(5e-5, -1e6, 123456, 1e-10, 0.25)) expect_true(bracesBalanced(numLatex(v)))
  expect_false(grepl("e", numLatex(-1e6), fixed = TRUE))
  expect_match(numLatex(-1e6), "^-")
})

# ---------------------------------------------------------------------------
# Names and expressions

test_that("latexName(): one letter, letter + digits, other names; escapeLatex()", {
  expect_equal(latexName("t"), "t")
  expect_equal(latexName("t", "i"), "t_{i}")
  expect_equal(latexName("x1"), "x_{1}")
  expect_equal(latexName("x12", "i"), "x_{12,i}")
  expect_equal(latexName("B12", "i"), "B_{12,i}")
  expect_equal(latexName("age"), r"(\mathit{age})")
  expect_equal(latexName("x1a", "i"), r"(\mathit{x1a}_{i})")
  expect_equal(latexName("x_1", "i"), r"(\mathit{x\_1}_{i})")
  expect_equal(latexName("my.var_2", "i"), r"(\mathit{my.var\_2}_{i})")
  expect_equal(latexName("xy", NULL), r"(\mathit{xy})")
  expect_equal(escapeLatex("a_b&c%d$e#f{g}"), r"(a\_b\&c\%d\$e\#f\{g\})")
})

test_that("latexEquation(): variable and group names with underscores, dots and digits", {
  tex <- latexEquation(parseJsmFormula(my.var_2 ~ x_a + a.b + z9 + B12 + (1 | my_id), part = "long"))
  expect_match(tex, r"(\mathit{my.var\_2}_{i} &= \beta^{(\mathit{my.var\_2})}_{0} + \beta^{(\mathit{my.var\_2})}_{1}\, \mathit{x\_a}_{i})", fixed = TRUE)
  expect_match(tex, r"(\beta^{(\mathit{my.var\_2})}_{2}\, \mathit{a.b}_{i} + \beta^{(\mathit{my.var\_2})}_{3}\, z_{9,i} + \beta^{(\mathit{my.var\_2})}_{4}\, B_{12,i})", fixed = TRUE)
  expect_match(tex, r"(b^{(\mathit{my.var\_2})}_{0,\mathrm{my\_id}[i]})", fixed = TRUE)
  expect_match(tex, r"(n_{\mathrm{my\_id}})", fixed = TRUE)
  expectWellFormed(tex)
  # an interaction grouping term keeps its colon
  tex <- latexEquation(parseJsmFormula(y ~ t + (1 | site:id), part = "long"))
  expect_match(tex, r"(b^{(y)}_{0,\mathrm{site:id}[i]})", fixed = TRUE)
  expect_match(tex, r"(n_{\mathrm{site:id}})", fixed = TRUE)
})

test_that("exprLatex(): operators, functions, constants and the observation index", {
  expect_equal(exprLatex(quote(t)), "t_{i}")
  expect_equal(exprLatex(quote(t), obs = NULL), "t")
  expect_equal(exprLatex(quote(I(t^2))), "{t_{i}}^{2}")
  expect_equal(exprLatex(quote((t))), r"(\left(t_{i}\right))")
  expect_equal(exprLatex(quote(-t)), "-t_{i}")
  expect_equal(exprLatex(quote(+t)), "t_{i}")
  expect_equal(exprLatex(quote(t - 1)), "t_{i} - 1")
  expect_equal(exprLatex(quote(t + x1)), "t_{i} + x_{1,i}")
  expect_equal(exprLatex(quote(a:b)), r"(a_{i} \cdot b_{i})")
  expect_equal(exprLatex(quote(a * b)), r"(a_{i} \cdot b_{i})")
  expect_equal(exprLatex(quote(sqrt(t) / exp(u))), r"(\frac{\sqrt{t_{i}}}{\exp\left(u_{i}\right)})")
  expect_equal(exprLatex(quote(abs(t))), r"(\left|t_{i}\right|)")
  for (fn in c("log", "exp", "sin", "cos", "tan"))
    expect_equal(exprLatex(str2lang(paste0(fn, "(t)"))), paste0("\\", fn, r"(\left(t_{i}\right))"))
  expect_equal(exprLatex(quote(log(t + 1)), obs = NULL), r"(\log\left(t + 1\right))")
  expect_equal(exprLatex(quote(my_fun(t, 2))), r"(\operatorname{my\_fun}\left(t_{i}, 2\right))")
  expect_equal(exprLatex(quote(ns(t, 3))), r"(\operatorname{ns}\left(t_{i}, 3\right))")
  expect_equal(exprLatex(2.5), "2.5")
  expect_equal(exprLatex(quote(sqrt(age))), r"(\sqrt{\mathit{age}_{i}})")
  for (e in list(quote(log(sqrt(t) / (1 + t^2))), quote(abs(a - b) * exp(-t))))
    expect_true(bracesBalanced(exprLatex(e)))
})

test_that("columnLatex(): factor levels need data; unparsable columns are \\mathrm", {
  df <- data.frame(g = factor(c("u", "v")), h = c("p_q", "50%"), t = 0:1)
  expect_equal(columnLatex("gv", df), r"(\mathbb{I}\left[g_{i} = \text{v}\right])")
  expect_equal(columnLatex("gv:t", df), r"(\mathbb{I}\left[g_{i} = \text{v}\right] \cdot t_{i})")
  expect_equal(columnLatex("hp_q", df), r"(\mathbb{I}\left[h_{i} = \text{p\_q}\right])")   # character column
  expect_equal(columnLatex("h50%", df), r"(\mathbb{I}\left[h_{i} = \text{50\%}\right])")
  expect_equal(columnLatex("gv", NULL), r"(\mathit{gv}_{i})")
  expect_equal(columnLatex("gv", df, obs = NULL), r"(\mathbb{I}\left[g = \text{v}\right])")
  expect_equal(columnLatex("poly(t, 2)1", df), r"(\mathrm{poly(t, 2)1})")
})

test_that("formulaColumns(): factor expansion with data, one column per term without", {
  df <- fxSmallData()
  expect_equal(formulaColumns(~ 0 + g, df), c("gu", "gv", "gw"))
  expect_equal(formulaColumns(~ 0 + g, NULL), "g")
  expect_equal(formulaColumns(~ 1, NULL), "(Intercept)")
})

# ---------------------------------------------------------------------------
# Effect numbering

# First line (the equation) of a longitudinal formula.
eqLine <- function(f, data = NULL) texLines(latexEquation(parseJsmFormula(f, part = "long"), data))[1]

test_that("latexEquation(): fixed-effect indices: intercept 0, `0 +` / `- 1` start at 1, empty", {
  expect_equal(eqLine(y ~ 1), r"(y_{i} &= \beta^{(y)}_{0} + \varepsilon^{(y)}_{i})")
  expect_equal(eqLine(y ~ 0 + t), r"(y_{i} &= \beta^{(y)}_{1}\, t_{i} + \varepsilon^{(y)}_{i})")
  expect_equal(eqLine(y ~ t - 1), r"(y_{i} &= \beta^{(y)}_{1}\, t_{i} + \varepsilon^{(y)}_{i})")
  expect_equal(eqLine(y ~ 0), r"(y_{i} &= \varepsilon^{(y)}_{i})")
  # without data a factor is one term
  expect_equal(eqLine(y ~ 0 + h), r"(y_{i} &= \beta^{(y)}_{1}\, h_{i} + \varepsilon^{(y)}_{i})")
})

test_that("latexEquation(): random-effect indices: shared with the fixed column, else the next free ones", {
  # the slope shares k = 1 with the fixed slope
  expect_equal(eqLine(y ~ t + (0 + t | id)),
               r"(y_{i} &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + b^{(y)}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(y)}_{i})")
  # RE blocks are written in the order given
  expect_match(eqLine(y ~ t + re(0 + t | id) + re(1 | id)),
               r"(b^{(y)}_{1,\mathrm{id}[i]}\, t_{i} + b^{(y)}_{0,\mathrm{id}[i]} + \varepsilon)", fixed = TRUE)
  # no fixed intercept: the RE intercept gets the next free index after the fixed ones
  expect_equal(eqLine(y ~ t - 1 + (1 | id)),
               r"(y_{i} &= \beta^{(y)}_{1}\, t_{i} + b^{(y)}_{2,\mathrm{id}[i]} + \varepsilon^{(y)}_{i})")
  expect_equal(eqLine(y ~ 0 + t + (1 | id) + (0 + s | id)),
               r"(y_{i} &= \beta^{(y)}_{1}\, t_{i} + b^{(y)}_{2,\mathrm{id}[i]} + b^{(y)}_{3,\mathrm{id}[i]}\, s_{i} + \varepsilon^{(y)}_{i})")
  expect_equal(eqLine(y ~ 1 + (0 + s | id) + (0 + w | g)),
               r"(y_{i} &= \beta^{(y)}_{0} + b^{(y)}_{1,\mathrm{id}[i]}\, s_{i} + b^{(y)}_{2,\mathrm{g}[i]}\, w_{i} + \varepsilon^{(y)}_{i})")
  # the free index is reused across blocks for the same column
  expect_equal(eqLine(y ~ 1 + (0 + s | id) + (0 + s | g)),
               r"(y_{i} &= \beta^{(y)}_{0} + b^{(y)}_{1,\mathrm{id}[i]}\, s_{i} + b^{(y)}_{1,\mathrm{g}[i]}\, s_{i} + \varepsilon^{(y)}_{i})")
})

test_that("latexEquation(): with data, factor columns expand and REs on them share the indices", {
  d <- fxSmallData()                                         # g has the levels u, v, w
  l <- eqLine(y ~ 0 + g + g:t + (0 + g | id), d)
  expect_match(l, r"(\beta^{(y)}_{1}\, \mathbb{I}\left[g_{i} = \text{u}\right] + \beta^{(y)}_{2}\, \mathbb{I}\left[g_{i} = \text{v}\right] + \beta^{(y)}_{3}\, \mathbb{I}\left[g_{i} = \text{w}\right])", fixed = TRUE)
  expect_match(l, r"(\beta^{(y)}_{4}\, \mathbb{I}\left[g_{i} = \text{u}\right] \cdot t_{i})", fixed = TRUE)
  expect_match(l, r"(\beta^{(y)}_{6}\, \mathbb{I}\left[g_{i} = \text{w}\right] \cdot t_{i})", fixed = TRUE)
  expect_match(l, r"(b^{(y)}_{1,\mathrm{id}[i]}\, \mathbb{I}\left[g_{i} = \text{u}\right] + b^{(y)}_{2,\mathrm{id}[i]}\, \mathbb{I}\left[g_{i} = \text{v}\right] + b^{(y)}_{3,\mathrm{id}[i]}\, \mathbb{I}\left[g_{i} = \text{w}\right])", fixed = TRUE)
  # an RE slope without a fixed counterpart: the next index after the expanded fixed columns
  l <- eqLine(y ~ 0 + g + (1 + t || id), d)
  expect_match(l, r"(b^{(y)}_{4,\mathrm{id}[i]} + b^{(y)}_{5,\mathrm{id}[i]}\, t_{i})", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Whole models

test_that("latexEquation(): line order: long equations, zero equations, outcomes, distributions, u lines", {
  f <- fxForms("multi")
  l <- texLines(latexEquation(f, prepareData(f, fxData()$long)))
  iLong <- c(grep("^x_\\{1,i\\} &= ", l), grep(r"(^\\log\\left\(\\lambda\^\{\(x_\{2\}\)\}_\{i\}\\right\) &= )", l))
  iZero <- grep("^0 &= ", l)
  iOut <- grep(r"(^(\\operatorname\{logit\}\\left\(p\^\{\(k\)\}|\\log\\left\(\\lambda\^\{\(\\mathit\{cnt\}\)\}))", l)
  iU <- grep("^u\\^", l)
  iDist <- setdiff(grep("&\\sim", l, fixed = TRUE), iU)
  expect_length(iLong, 2L)
  expect_length(iZero, 3L)                                   # one per scaled() term
  expect_length(iOut, 2L)
  expect_length(iU, 3L)
  expect_equal(sort(c(iLong, iZero, iOut, iDist, iU)), seq_along(l))   # every line classified
  expect_true(max(iLong) < min(iZero) && max(iZero) < min(iOut) && max(iOut) < min(iDist) && max(iDist) < min(iU))
  # long distributions before the outcome ones
  expect_match(l[iDist[1]], r"(^\\varepsilon\^\{\(x_\{1\}\)\})")
  expect_lt(grep(r"(^x_\{2,i\} &\\sim \\mathrm\{Poisson\})", l), grep(r"(^k_\{i\} &\\sim \\mathrm\{Binomial\})", l))
  expect_lt(grep("Binomial", l), grep(r"(^\\mathit\{cnt\}_\{i\} &\\sim \\mathrm\{Poisson\})", l))
  expect_match(l[iU[3]], r"(\tau_{x_{2} \to \mathit{cnt}})", fixed = TRUE)
  expect_false(grepl("sigma", paste(l, collapse = "")))      # precision parametrisation only
})

test_that("latexEquation(): a zero equation is the long predictor without index minus u over D", {
  f <- fxForms("multi")
  tex <- latexEquation(f, prepareData(f, fxData()$long))
  l <- texLines(tex)
  z <- texLine(l, r"(u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id},t}, \quad)")
  expect_equal(z, paste0(r"(0 &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t + b^{(x_{2})}_{0,\mathrm{id}})",
                         r"( - u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{cnt}})"))
  zero <- grep("^0 &=", l, value = TRUE)
  expect_false(any(grepl("[i]", zero, fixed = TRUE)))
  expect_false(any(grepl("varepsilon", zero, fixed = TRUE)))
  expect_false(grepl("+ -", tex, fixed = TRUE))
  expectWellFormed(tex)
})

test_that("latexEquation(): an outcome uses u at its own values; factor covariates expand with data", {
  f <- fxForms("multi")
  l <- texLines(latexEquation(f, prepareData(f, fxData()$long)))
  cnt <- texLine(l, r"(\log\left(\lambda^{(\mathit{cnt})}_{i}\right) &=)")
  expect_match(cnt, r"(\mathbb{I}\left[\mathit{grp}_{i} = \text{b}\right])", fixed = TRUE)
  expect_match(cnt, r"(\gamma_{\mathit{cnt},x_{2}}\,u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id}[i],t_{i}})", fixed = TRUE)
  expect_false(grepl("x_{1}", cnt, fixed = TRUE))            # cnt is not linked to x1
  k <- texLine(l, r"(\operatorname{logit}\left(p^{(k)}_{i}\right) &=)")
  expect_match(k, r"(\gamma_{k,x_{1}}\,u^{(x_{1} \to k)}_{\mathrm{id}[i],t_{i}} + \gamma_{k,x_{2}}\,u^{(x_{2} \to k)}_{\mathrm{id}[i],t_{i}})", fixed = TRUE)
})

test_that("latexEquation(): u is indexed by the long groups and covariates; factors expand with data", {
  # site (A / B, a character column) is the factor-like covariate: fxSmallData()'s g is collinear with t
  f <- collectFormulas(y ~ 1 + scaled(x), x ~ t + site + (1 | id))
  dat <- prepareData(f, fxSmallData())
  l <- texLines(latexEquation(f, dat))
  z <- texLine(l, "0 &=")
  expect_match(z, r"(\mathbb{I}\left[\mathit{site} = \text{B}\right])", fixed = TRUE)
  expect_match(z, r"(u^{(x \to y)}_{\mathrm{id},t,\mathit{site}}, \quad \left(\mathrm{id}, t, \mathit{site}\right) \in \mathcal{D}_{y})", fixed = TRUE)
  expect_match(texLine(l, "y_{i} &="), r"(u^{(x \to y)}_{\mathrm{id}[i],t_{i},\mathit{site}_{i}})", fixed = TRUE)
  # without data the factor is one term
  expect_match(texLine(texLines(latexEquation(f)), "0 &="), r"(\beta^{(x)}_{2}\, \mathit{site} + b)", fixed = TRUE)
  # a long formula without random effects: u indexed by its covariates only
  f0 <- collectFormulas(y ~ t + scaled(x), x ~ t)
  l0 <- texLines(latexEquation(f0))
  expect_match(texLine(l0, "0 &="), r"(- u^{(x \to y)}_{t}, \quad \left(t\right) \in \mathcal{D}_{y})", fixed = TRUE)
  expect_match(texLine(l0, "y_{i} &="), r"(\gamma_{y,x}\,u^{(x \to y)}_{t_{i}})", fixed = TRUE)
})

test_that("latexEquation(): a standalone parsed outcome shows only its written scaled() terms", {
  f <- fxForms("multi")
  tex <- latexEquation(f$out$cnt)
  expect_match(tex, r"(\gamma_{\mathit{cnt},x_{2}}\,u^{(x_{2} \to \mathit{cnt})}_{i})", fixed = TRUE)
  expect_false(grepl("x_{1}", tex, fixed = TRUE))
  expect_false(grepl("0 &=", tex, fixed = TRUE))
  tex <- latexEquation(f$out$k)
  expect_equal(lengths(regmatches(tex, gregexpr(r"(\gamma_)", tex, fixed = TRUE))), 2L)
  # parseJsmFormula() result (no family): generic link, one line
  l <- texLines(latexEquation(parseJsmFormula(y ~ t + scaled(a) + scaled(b2), part = "out")))
  expect_length(l, 1L)
  expect_match(l, r"(\gamma_{y,a}\,u^{(a \to y)}_{i} + \gamma_{y,b_{2}}\,u^{(b_{2} \to y)}_{i})", fixed = TRUE)
})

test_that("latexEquation(): a survival outcome in a whole model with data", {
  m <- fxModel("surv", fit = FALSE)
  l <- texLines(latexEquation(m))
  expect_match(l[1], "^x_\\{1,i\\} &= ")
  expect_match(texLine(l, "0 &="), r"(u^{(x_{1} \to \mathit{time})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{time}})",
               fixed = TRUE)
  expect_match(texLine(l, r"(\log\left(\lambda^{(\mathit{time})}_{i}\right) &=)"),
               r"(\beta^{(\mathit{time})}_{0} + \gamma_{\mathit{time},x_{1}}\,u^{(x_{1} \to \mathit{time})}_{\mathrm{id}[i],t_{i}})", fixed = TRUE)
  expect_match(texLine(l, "Weibull"), r"(\mathit{status}_{i} \in \{0, 1\})", fixed = TRUE)
  expect_match(l[length(l)], "^u\\^")
})

test_that("latexEquation(): two survival outcomes (Weibull, exponential) linked to one marker", {
  f <- collectFormulas(list(time | event(status) ~ 1 + scaled(x1), w | event(ev) ~ 1 + scaled(x1)),
                       list(x1 ~ t + (1 | id)),
                       family = list(time = weibullSurvFamily(), w = exponentialSurvFamily()))
  l <- texLines(latexEquation(f))
  expect_match(texLine(l, "Weibull"), r"(\mathit{time}_{i} &\sim \mathrm{Weibull}\left(\alpha^{(\mathit{time})}, \lambda^{(\mathit{time})}_{i}\right))", fixed = TRUE)
  expect_match(texLine(l, r"(\mathrm{Exp}\left()"), r"(w_{i} &\sim \mathrm{Exp}\left(\lambda^{(w)}_{i}\right), \quad S(t) =)", fixed = TRUE)
  expect_match(texLine(l, r"(\log\left(\lambda^{(w)}_{i}\right) &=)"), r"(\gamma_{w,x_{1}}\,u^{(x_{1} \to w)})", fixed = TRUE)
  expect_length(grep("^0 &=", l), 2L)
})

test_that("latexEquation(): snapshot of the whole equation of each fixture", {
  # cran = TRUE: also verified when NOT_CRAN is unset (plain test_file()); the output is
  # deterministic (seeded fixtures, 3 significant digits of the lme4-derived Wishart scale)
  for (nm in c("beta", "multi", "surv", "priors"))
    expect_snapshot(cat(format(latexEquation(fxModel(nm, fit = FALSE)))), cran = TRUE)
})

# ---------------------------------------------------------------------------
# format / print / errors

test_that("format() and print() of mvremLatex", {
  tex <- latexEquation(parseJsmFormula(y ~ t, part = "long"))
  expect_s3_class(tex, "mvremLatex")
  expect_match(tex, "^\\\\begin\\{aligned\\}\n")
  expect_match(tex, "\n\\\\end\\{aligned\\}$")
  expect_identical(format(tex), unclass(tex))
  expect_identical(format(tex, delim = "$$"), paste0("$$", unclass(tex), "$$"))
  expect_identical(format(tex, delim = "\\["), paste0("\\[", unclass(tex), "\\["))
  expect_output(print(tex, delim = "$$"), "^\\$\\$\\\\begin\\{aligned\\}")
  expect_output(print(tex), "end{aligned}", fixed = TRUE)
  capture.output(vis <- withVisible(print(tex)))
  expect_false(vis$visible)
  expect_identical(capture.output(r <- print(tex)), strsplit(unclass(tex), "\n")[[1]])
  expect_identical(r, tex)
})

test_that("latexEquation(): invalid input and mismatched data error", {
  f <- fxForms("beta")
  dat <- prepareData(f, fxData()$long)
  expect_error(latexEquation(list(a = 1)), "parsed formula")
  expect_error(latexEquation(f$long$x1, dat), "data frame")
  expect_error(latexEquation(f, fxData()$long), "prepareData")
  expect_identical(latexEquation(f), latexEquation(f, NULL))
})

test_that("latexEquation.mvremModel(): annotate = 'none' (the default) needs no fit and is latexEquation(forms, data)", {
  m <- fxModel("beta", fit = FALSE)
  expect_error(latexEquation(m, annotate = "tooltip"), "not fitted")
  expect_error(latexEquation(m, annotate = "inline"), "not fitted")
  expect_error(latexEquation(m, annotate = "bogus"), "should be one of")
  expect_identical(latexEquation(m, annotate = "none"), latexEquation(m$forms, m$data))
  expect_identical(latexEquation(m), latexEquation(m, annotate = "none"))
  expect_null(latexState$est)
})

# ---------------------------------------------------------------------------
# Internal annotation helpers (no fit needed)

test_that("annotKey(): one key format for producer and consumers, vectorised, NULL parts empty", {
  expect_equal(annotKey("prec", "x1", "id", "t"), "prec|x1|id|t")
  expect_equal(annotKey("b", "x1", "id", c("(Intercept)", "t")), c("b|x1|id|(Intercept)", "b|x1|id|t"))
  expect_equal(annotKey("family", NULL, "prec"), "family||prec")     # standalone formula: no response
  expect_equal(distKey("normal"), "dist|normal")
  expect_equal(blockKey("x1", "id", c("(Intercept)", "t")), "block|x1|id|(Intercept)|t")
})

test_that("annot(): plain symbol without state or for unknown keys; tooltip / inline wrappers", {
  expect_identical(annot("x", "fixed|y|t"), "x")
  keys <- list("k" = estEntry("fixed effect", mean = 1.23456, sd = 0.1, lo = 1, hi = 2),
               "n" = estEntry("note only", note = "just a note"),
               "r" = estEntry("range", mean = 1, lo = 0, hi = 2, band = "range"))
  withEstState(list(mode = "tooltip", digits = 3, keys = keys), {
    expect_equal(annot("x", "k"), r"(\htmlData{par=fixed effect,mean=1.23,sd=0.1,lo=1,hi=2}{x})")
    expect_equal(annot("x", "unknown"), "x")
    expect_equal(annot(c("a", "b"), c("k", "zz")), c(r"(\htmlData{par=fixed effect,mean=1.23,sd=0.1,lo=1,hi=2}{a})", "b"))
    expect_equal(annot("x", "n"), r"(\htmlData{par=note only,note=just a note}{x})")
  })
  withEstState(list(mode = "inline", digits = 2, keys = keys), {
    expect_equal(annot("x", "k"), r"(\underset{\substack{\scriptstyle 1.2 \\ \scriptstyle (0.1)}}{x})")
    expect_equal(annot("x", "n"), "x")                         # no mean
    expect_equal(annot("x", "r"), "x")                         # range, not a 95% interval
  })
  expect_null(latexState$est)
})

test_that("htmlDataString(): sanitised values, vectors joined by ';', NA, band, digits", {
  s <- htmlDataString(estEntry("a, b = {c}; d", mean = c(1.23456, NA), sd = 0.1, lo = -1, hi = 2,
                               band = "80%", note = "x=1, {y}"), 2)
  expect_equal(s, "par=a b c d,mean=1.2;NA,sd=0.1,lo=-1,hi=2,band=80,note=x 1 y")
  expectValidHtmlData(s)
  s <- htmlDataString(estEntry("p", mean = 1e-7, lo = 0, hi = 1, item = c("a,b", "c"), band = "95%"), 3)
  expect_equal(s, "par=p,item=a b;c,mean=1e-07,lo=0,hi=1")
  expectValidHtmlData(s)
  expect_equal(htmlDataString(estEntry("p", mean = 123456), 3), "par=p,mean=1.23e+05")
})

test_that("inlineEstimate(): mean and sd on two lines, one line per item for matrices", {
  expect_equal(inlineEstimate(estEntry("p", mean = 0.123456, sd = 0.01)), r"(\scriptstyle 0.123 \\ \scriptstyle (0.01))")
  expect_equal(inlineEstimate(estEntry("p", mean = 2, item = "SD")), r"(\scriptstyle \text{SD}\ 2)")
  expect_equal(inlineEstimate(estEntry("p", mean = c(1, 2), sd = c(0.1, 0.2), item = c("SD a_b", "corr a:b"))),
               r"(\scriptstyle \text{SD a\_b}\ 1\,(0.1) \\ \scriptstyle \text{corr a:b}\ 2\,(0.2))")
})

test_that("priorNote(): text for every prior kind, raw lists and the default", {
  expect_equal(priorNote(NULL), "prior: INLA default")
  expect_match(priorNote(gammaPrec(1, 5e-5)), "Gamma with shape 1 and rate 5e-05")
  expect_match(priorNote(pcSd(0.5)), "P(SD > 0.5) is 0.01 (SD is Exp(9.21))", fixed = TRUE)
  expect_match(priorNote(halfNormalSd(2)), "half-normal prior on the SD with scale 2")
  expect_match(priorNote(halfCauchySd(3)), "half-Cauchy prior on the SD with scale 3")
  expect_match(priorNote(wishartPrec()), "with k + 2 degrees of freedom", fixed = TRUE)
  expect_match(priorNote(wishartPrec(df = 5)), "with 5 degrees of freedom")
  expect_match(priorNote(wishartPrec(scale = "data")), "lme4 fit")
  expect_match(priorNote(normalCoef(1, 2)), "mean 1 and SD 2")
  expect_equal(priorNote(list(prior = "loggamma", param = c(1, 0.1), fixed = TRUE)),
               "prior: INLA loggamma with param 1 0.1 (fixed)")
  expect_equal(priorNote(list(initial = 1)), "prior: INLA hyper list")
  expect_null(priorNote("x"))
  for (pr in list(gammaPrec(), pcSd(1), halfNormalSd(1), halfCauchySd(1), wishartPrec(), normalCoef()))
    expectValidHtmlData(htmlDataString(estEntry("p", note = priorNote(pr)), 3))
})

test_that("distributionNotes(): a note for every distribution and prior name written in the equations", {
  notes <- distributionNotes()
  used <- paste0("dist|", c("normal", "t", "gamma", "beta", "poisson", "negbin", "binomial", "bernoulli",
                            "betabinomial", "zero", "truncated", "inla", "weibull", "exp", "event",
                            "priorGamma", "priorPc", "priorHalfNormal", "priorHalfCauchy", "priorWishart",
                            "priorNormal"))
  expect_setequal(names(notes), used)
  # the keys of distName() are tolower() of the written names
  expect_true(all(paste0("dist|", tolower(c("Beta", "Poisson", "NegBin", "Binomial", "Bernoulli",
                                             "BetaBinomial", "Gamma", "Weibull", "Exp"))) %in% names(notes)))
  for (k in names(notes)) {
    expect_false(is.null(notes[[k]]$par), info = k)
    expect_false(is.null(notes[[k]]$note), info = k)
    expectValidHtmlData(htmlDataString(notes[[k]], 3))
  }
})

test_that("familyLatex(): every family's distribution name is annotated when its key exists", {
  keys <- distributionNotes()
  fams <- list(betaFamily(), binomialFamily(), poissonFamily(), negBinomialFamily(), gammaFamily(),
               zeroInflatedPoissonFamily(0), zeroInflatedNegBinomialFamily(1), betaBinomialFamily(),
               weibullSurvFamily(), exponentialSurvFamily(), inlaFamily("sn"), studentTFamily(), gaussianFamily())
  withEstState(list(mode = "tooltip", digits = 3, keys = keys), {
    for (fam in fams) {
      f <- if (fam$kind %in% c("weibull", "exponential")) time | event(status) ~ 1 else y ~ t
      tex <- famTex(f, fam)
      expect_match(tex, r"(\htmlData{par=)", fixed = TRUE, info = fam$family)
      expectWellFormed(tex)
      expect_identical(unwrapCommand(tex, "\\htmlData"), unclass(withEstState(NULL, famTex(f, fam))))
    }
  })
})

# ---------------------------------------------------------------------------
# Fitted models: annotate = "tooltip" / "inline"

# formatC() as htmlDataString() writes numbers (independent of the code under test)
tipNum <- function(v, d = 3L) numText(v, d)

test_that("latexEquation(annotate = 'tooltip'): balanced, valid \\htmlData arguments, same equations underneath", {
  for (nm in c("beta", "multi", "surv", "priors")) {
    m <- fxModel(nm)
    tex <- latexEquation(m, annotate = "tooltip")
    plain <- latexEquation(m)
    expect_s3_class(tex, "mvremLatex")
    expectWellFormed(tex)
    w <- htmlDataWrappers(tex)
    expect_equal(nrow(w), lengths(regmatches(tex, gregexpr("\\htmlData{", tex, fixed = TRUE))))
    for (a in unique(w$arg)) expectValidHtmlData(a)
    # every fixed effect of the model is annotated once per occurrence of its symbol
    fe <- fixedEffects(m)
    for (i in seq_len(nrow(fe)))
      expect_true(any(startsWith(w$arg, paste0("par=fixed effect of ", effectLabel(fe$term[i], m$inla$naming),
                                                " in ", fe$response[i], ","))), info = paste(nm, fe$term[i]))
    expect_identical(unwrapCommand(tex, "\\htmlData"), unclass(plain))
    expect_null(latexState$est)                             # state restored
    expect_false(grepl("htmlData", plain, fixed = TRUE))
  }
})

test_that("latexEquation(annotate = 'tooltip'): posterior summaries of the parameters with `digits`", {
  m <- fxModel("beta")
  fe <- fixedEffects(m)
  for (d in c(3L, 2L, 5L)) {
    tex <- latexEquation(m, annotate = "tooltip", digits = d)
    for (i in seq_len(nrow(fe)))
      expect_match(tex, paste0("\\htmlData{par=fixed effect of ", fe$term[i], " in ", fe$response[i], ",mean=",
                               tipNum(fe$mean[i], d), ",sd=", tipNum(fe$sd[i], d), ",lo=", tipNum(fe$q0.025[i], d),
                               ",hi=", tipNum(fe$q0.975[i], d), "}{\\beta^{"), fixed = TRUE)
  }
  tex <- latexEquation(m, annotate = "tooltip")
  a <- associations(m)
  expect_match(tex, paste0("\\htmlData{par=association gamma (x1 -> y),mean=", tipNum(a$mean[1]), ",sd=", tipNum(a$sd[1]),
                           ",lo=", tipNum(a$q0.025[1]), ",hi=", tipNum(a$q0.975[1]), ",note="), fixed = TRUE)
  expect_match(tex, "normal prior on the coefficient with mean 0 and SD 10 (default)}{\\gamma_{y,x_{1}}}", fixed = TRUE)
  h <- hyperparameters(m)
  phi <- h[h$response == "y" & h$parameter == "phi", ]
  expect_match(tex, paste0("\\htmlData{par=phi of y,mean=", tipNum(phi$mean), ",sd=", tipNum(phi$sd), ",lo=",
                           tipNum(phi$q0.025), ",hi=", tipNum(phi$q0.975), ",note=prior: INLA default}{\\phi^{(y)}}"), fixed = TRUE)
  # the Q matrix: SDs and the correlation of varianceComponents()
  vc <- varianceComponents(m)$table
  blk <- vc[vc$block == "x1 | id [(Intercept), t]", ]          # SD (Intercept), SD t, correlation
  expect_match(tex, paste0("item=SD (Intercept);SD t;corr (Intercept):t,mean=", paste(tipNum(blk$mean), collapse = ";"),
                           ",sd=", paste(tipNum(blk$sd), collapse = ";")), fixed = TRUE)
  # precision tooltips also give the SD and the (default) prior
  eps <- vc[vc$response == "x1" & vc$parameter == "sd" & vc$effect == "residual", ]
  rp <- h[h$response == "x1" & h$parameter == "residual precision", ]
  expect_match(tex, paste0("\\htmlData{par=residual precision of x1,mean=", tipNum(rp$mean), ",sd=", tipNum(rp$sd),
                           ",lo=", tipNum(rp$q0.025), ",hi=", tipNum(rp$q0.975), ",note=SD ", tipNum(eps$mean), " (",
                           tipNum(eps$sd), "). prior: INLA default}{\\tau^{(x_{1})}_{\\varepsilon}}"), fixed = TRUE)
  # a `|` block with several effects writes Q, not the single precisions: the prec keys exist
  # in the lookup (SD + prior note) but wrap no symbol
  sdInt <- blk[blk$parameter == "sd" & blk$effect == "(Intercept)", ]
  k <- estimateLookup(m)
  expect_equal(k[["prec|x1|id|(Intercept)"]]$note, paste0("SD ", tipNum(sdInt$mean), " (", tipNum(sdInt$sd), "). prior: INLA default"))
  expect_false(any(startsWith(htmlDataWrappers(tex)$arg, "par=precision of random effect")))
})

test_that("latexEquation(annotate = 'tooltip'): notes for data, responses, fitted means and the zero equation", {
  m <- fxModel("beta")
  tex <- latexEquation(m, annotate = "tooltip")
  # t takes the values 0, 0.5, 1, 1.5 (mean 0.75); 40 subjects x 4 times = 160 outcome rows
  expect_match(tex, "\\htmlData{par=variable t,item=values,mean=0.75,lo=0,hi=1.5,band=range}{t_{i}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=levels of id in the data of x1,note=40 levels}{n_{\\mathrm{id}}}", fixed = TRUE)
  expect_match(tex, "note=160 unique combinations}{\\mathcal{D}_{y}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=response y (outcome),item=observed,mean=", fixed = TRUE)
  expect_match(tex, "band=range,note=160 observations}{y_{i}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=response x1 (longitudinal),", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=fitted mean of y (one per observation: fittedValues()),item=posterior means,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=zero-equation field u (x1 -> y),note=linear predictor of x1 at the 160 combinations of (id t) in the data of y. Its precision is fixed}{u^{(x_{1} \\to y)}_{\\mathrm{id},t}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=precision of u (x1 -> y),note=fixed at exp(-6) (not estimated)}{\\tau_{x_{1} \\to y}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=random effect (Intercept) | id of x1,item=SD,mean=", fixed = TRUE)
  expect_match(tex, "note=one value per level: randomEffects()}{b^{(x_{1})}_{0,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=residual of x1,item=SD,mean=", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=normal distribution,note=", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=beta distribution,note=", fixed = TRUE)
})

test_that("latexEquation(annotate = 'tooltip'): specifics of the 'multi' and 'surv' fixtures", {
  tex <- latexEquation(fxModel("multi"), annotate = "tooltip")
  expect_match(tex, "\\htmlData{par=variable grp,note=3 levels: a b c}{\\mathit{grp}_{i}}", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=variable n,item=values,", fixed = TRUE)          # trials
  expect_match(tex, "\\htmlData{par=variable pt,item=values,", fixed = TRUE)         # exposure
  expect_match(tex, "\\htmlData{par=binomial distribution,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=Poisson distribution,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=fixed effect of grpc in cnt,", fixed = TRUE)
  tex <- latexEquation(fxModel("surv"), annotate = "tooltip")
  expect_match(tex, "par=fitted hazard scale lambda of time", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=event indicator,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=Weibull distribution,", fixed = TRUE)
  expect_match(tex, "\\htmlData{par=alpha of time,", fixed = TRUE)
})

test_that("effectLabel(): the intercept is labelled by namingScheme(interceptLabel = ), default `(Intercept)`", {
  expect_equal(effectLabel(c("(Intercept)", "t", "gb"), namingScheme()), c("(Intercept)", "t", "gb"))
  expect_equal(effectLabel(c("(Intercept)", "t"), namingScheme(interceptLabel = "intercept")), c("intercept", "t"))
  expect_equal(namingScheme()$interceptLabel, "(Intercept)")
  # jsm(naming = ) stores the scheme in the model (no fit needed)
  s <- fxSpec("beta")
  m <- suppressMessages(do.call(jsm, c(s, list(fit = FALSE, naming = namingScheme(interceptLabel = "intercept")))))
  expect_equal(m$inla$naming$interceptLabel, "intercept")
  expect_equal(fxModel("beta", fit = FALSE)$inla$naming$interceptLabel, "(Intercept)")
})

test_that("latexEquation(annotate = 'tooltip'): the intercept notes follow namingScheme(interceptLabel = )", {
  m <- fxModel("beta")
  w <- htmlDataWrappers(latexEquation(m, annotate = "tooltip"))
  expect_true(any(startsWith(w$arg, "par=fixed effect of (Intercept) in y,")))
  expect_true(any(startsWith(w$arg, "par=random effect (Intercept) | id of x1,")))
  expect_true(any(grepl("item=SD (Intercept);SD t;corr (Intercept):t,", w$arg, fixed = TRUE)))
  expect_equal(estimateLookup(m)[["prec|x1|id|(Intercept)"]]$par, "precision of random effect (Intercept) | id in x1")
  # the same fit with another label (a copy of the fixture; the label lives in the stored scheme)
  m2 <- m
  m2$inla$naming$interceptLabel <- "intercept"
  w2 <- htmlDataWrappers(latexEquation(m2, annotate = "tooltip"))
  expect_true(any(startsWith(w2$arg, "par=fixed effect of intercept in y,")))
  expect_true(any(startsWith(w2$arg, "par=random effect intercept | id of x1,")))
  expect_true(any(grepl("item=SD intercept;SD t;corr intercept:t,", w2$arg, fixed = TRUE)))
  expect_equal(estimateLookup(m2)[["prec|x1|id|(Intercept)"]]$par, "precision of random effect intercept | id in x1")
  expect_match(hyperTable(m2)$parameter, "precision of intercept \\| id", all = FALSE)
  expect_false(any(grepl("(Intercept)", w2$arg, fixed = TRUE)))
  expect_setequal(gsub("(Intercept)", "intercept", w$arg, fixed = TRUE), w2$arg)
  expect_identical(unwrapCommand(latexEquation(m2, annotate = "tooltip"), "\\htmlData"), unclass(latexEquation(m)))
})

test_that("estimateLookup(): prec keys come from hyperTable()'s group / effect columns", {
  m <- fxModel("priors")
  h <- hyperTable(m)
  k <- estimateLookup(m)
  prec <- h[h$type == "random" & is.na(h$effect2), ]
  expect_gt(nrow(prec), 0L)
  expect_true(all(!is.na(prec$group) & !is.na(prec$effect)))
  precKeys <- paste0("prec|", prec$response, "|", prec$group, "|", prec$effect)
  expect_true(all(precKeys %in% names(k)))
  expect_setequal(grep("^prec\\|", names(k), value = TRUE), precKeys)       # correlations give no prec key
  # an outcome RE on `site`, two blocks on one group, a `||` block
  expect_true(all(c("prec|yg|site|(Intercept)", "prec|x2|id|(Intercept)", "prec|x2|id|t",
                    "prec|x3|id|(Intercept)", "prec|x3|id|t") %in% names(k)))
  expect_equal(k[["prec|yg|site|(Intercept)"]]$par, "precision of random effect (Intercept) | site in yg")
  expect_equal(k[["prec|x2|id|t"]]$par, "precision of random effect t | id in x2")
  # the tooltips wrap the tau of the block's distribution line
  w <- fxTips("priors")
  symOf <- function(par) unique(w$sym[startsWith(w$arg, paste0("par=", par, ","))])
  expect_equal(symOf("precision of random effect (Intercept) | site in yg"), r"(\tau^{(\mathit{yg})}_{0,\mathrm{site}})")
  expect_equal(symOf("precision of random effect (Intercept) | id in x2"), r"(\tau^{(x_{2})}_{0,\mathrm{id}})")
  expect_equal(symOf("precision of random effect t | id in x2"), r"(\tau^{(x_{2})}_{1,\mathrm{id}})")
  expect_equal(symOf("precision of random effect t | id in x3"), r"(\tau^{(x_{3})}_{1,\mathrm{id}})")
})

# --- one test per prior kind on the "priors" fixture: its line and its tooltip note ---

# The unique \htmlData argument whose par is `par`.
tipArg <- function(w, par) {
  a <- unique(w$arg[startsWith(w$arg, paste0("par=", par, ","))])
  testthat::expect_length(a, 1L)
  if (length(a) == 1L) a else NA_character_
}

test_that("latexEquation(): wishartPrec(scale = 'data') on x1: the stored resolved matrix and its note", {
  m <- fxModel("priors")
  V <- m$inla$random$x1[[1]]$prior$args$scale
  expect_true(is.matrix(V) && all(dim(V) == 2L))
  expect_equal(m$inla$random$x1[[1]]$prior$args$source, "data")
  # the stored matrix is S^-1 / df of the lme4 fit of the x1 formula (df = k + 2 = 4)
  S <- as.matrix(lme4::VarCorr(lme4::lmer(x1 ~ t + (1 + t | id), fxData()$long))$id)
  expect_equal(V, unname(solve(S) / 4), tolerance = 1e-6)
  l <- fxTex("priors")
  pm <- paste0("\\begin{pmatrix}", tipNum(V[1, 1]), " & ", tipNum(V[1, 2]), " \\\\ ", tipNum(V[2, 1]), " & ",
               tipNum(V[2, 2]), "\\end{pmatrix}")
  expect_match(texLine(l, r"(\mathbf{Q}^{(x_{1})}_{\mathrm{id}} &\sim)"), paste0(r"(\mathcal{W}_{2}\left(4, )", pm, r"(\right))"), fixed = TRUE)
  w <- fxTips("priors")
  expect_match(tipArg(w, "covariance of the random effects | id of x1"),
               "note=Wishart prior on the precision matrix with k + 2 degrees of freedom (E[Q] is df times the scale matrix - here the precision matrix of an lme4 fit of the formula divided by df)", fixed = TRUE)
  # a `|` block with several effects writes Q, not its single precisions: the note is in the lookup only
  expect_match(estimateLookup(m)[["prec|x1|id|t"]]$note, "lme4 fit of the formula divided by df", fixed = TRUE)
  expect_false(any(startsWith(w$arg, "par=precision of random effect t | id in x1,")))
  expect_true(any(startsWith(w$arg, "par=Wishart prior,note=")))
  expect_equal(unique(w$sym[startsWith(w$arg, "par=Wishart prior,")]), r"(\mathcal{W}_{2})")
})

test_that("latexEquation(): pcSd(1) on the yg residual and the site block: lines and notes", {
  m <- fxModel("priors")
  l <- fxTex("priors")
  # Exp rate = -log(0.01) / 1 = 4.61
  pc <- r"(^{-1/2} &\sim \mathrm{Exp}\left(4.61\right), \quad P\left(\left()"
  expect_match(texLine(l, r"(\left(\tau^{(\mathit{yg})}_{\varepsilon}\right)^{-1/2} &\sim)"),
               paste0(r"(\left(\tau^{(\mathit{yg})}_{\varepsilon}\right))", pc, r"(\tau^{(\mathit{yg})}_{\varepsilon}\right)^{-1/2} > 1\right) = 0.01)"), fixed = TRUE)
  expect_match(texLine(l, r"(\left(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1/2} &\sim)"),
               paste0(r"(\left(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right))", pc, r"(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1/2} > 1\right) = 0.01)"), fixed = TRUE)
  expect_equal(grep(r"(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1/2} &\sim)", l, fixed = TRUE),
               grep(r"(b^{(\mathit{yg})}_{0,l} &\sim)", l, fixed = TRUE) + 1L)
  w <- fxTips("priors")
  note <- "PC prior on the SD: P(SD > 1) is 0.01 (SD is Exp(4.61))"
  expect_match(tipArg(w, "residual precision of yg"), paste0(",note=SD [0-9.]+ \\([0-9.]+\\)\\. ", gsub("([().])", "\\\\\\1", note), "$"))
  expect_match(tipArg(w, "precision of random effect (Intercept) | site in yg"), note, fixed = TRUE)
  # a single-effect block writes its precision, not Q: the block entry is in the lookup only
  blk <- estimateLookup(m)[["block|yg|site|(Intercept)"]]
  expect_equal(blk$item, "SD (Intercept)")
  expect_equal(blk$note, note)
  expect_false(any(startsWith(w$arg, "par=covariance of the random effects | site of yg,")))
  expect_true(any(startsWith(w$arg, "par=penalised-complexity prior,note=")))
  expect_equal(unique(w$sym[startsWith(w$arg, "par=penalised-complexity prior,")]), r"(\mathrm{Exp})")
})

test_that("latexEquation(): halfCauchySd(1) on the x2 intercept block: line and note", {
  m <- fxModel("priors")
  l <- fxTex("priors")
  expect_match(texLine(l, r"(\mathrm{C}^{+})"), r"(\left(\tau^{(x_{2})}_{0,\mathrm{id}}\right)^{-1/2} &\sim \mathrm{C}^{+}\left(0, 1\right))", fixed = TRUE)
  expect_equal(grep(r"(\mathrm{C}^{+})", l, fixed = TRUE), grep(r"(b^{(x_{2})}_{0,l} &\sim)", l, fixed = TRUE) + 1L)
  w <- fxTips("priors")
  expect_match(tipArg(w, "precision of random effect (Intercept) | id in x2"), "half-Cauchy prior on the SD with scale 1$")
  expect_match(tipArg(w, "precision of random effect (Intercept) | id in x2"), ",note=SD [0-9.]+ \\([0-9.]+\\)\\. half-Cauchy")
  expect_equal(unique(w$sym[startsWith(w$arg, "par=half-Cauchy prior,")]), r"(\mathrm{C}^{+})")
})

test_that("latexEquation(): halfNormalSd(0.5) on the x2 slope block: line and note", {
  m <- fxModel("priors")
  l <- fxTex("priors")
  # precision 1 / 0.5^2 = 4
  expect_match(texLine(l, r"(\mathcal{N}^{+})"), r"(\left(\tau^{(x_{2})}_{1,\mathrm{id}}\right)^{-1/2} &\sim \mathcal{N}^{+}\left(0, \left(4\right)^{-1}\right))", fixed = TRUE)
  expect_equal(grep(r"(\mathcal{N}^{+})", l, fixed = TRUE), grep(r"(b^{(x_{2})}_{1,l} &\sim)", l, fixed = TRUE) + 1L)
  w <- fxTips("priors")
  expect_match(tipArg(w, "precision of random effect t | id in x2"), "half-normal prior on the SD with scale 0.5$")
  expect_equal(unique(w$sym[startsWith(w$arg, "par=half-normal prior,")]), r"(\mathcal{N}^{+})")
})

test_that("latexEquation(): gammaPrec() on the x3 `||` block and on the x2 size / x3 precision: lines and notes", {
  m <- fxModel("priors")
  l <- fxTex("priors")
  expect_match(texLine(l, r"(\tau^{(x_{3})}_{k,\mathrm{id}} &\sim)"), r"(\mathrm{Gamma}\left(1, 0.01\right), \quad k = 0, 1)", fixed = TRUE)
  expect_equal(grep(r"(\tau^{(x_{3})}_{k,\mathrm{id}} &\sim)", l, fixed = TRUE), grep(r"(b^{(x_{3})}_{1,l} &\sim)", l, fixed = TRUE) + 1L)
  expect_match(texLine(l, r"(\kappa^{(x_{2})} &\sim)"), r"(\mathrm{Gamma}\left(1, 0.01\right))", fixed = TRUE)
  expect_match(texLine(l, r"(\tau^{(x_{3})}_{\varepsilon} &\sim)"), r"(\mathrm{Gamma}\left(1, 5 \cdot 10^{-5}\right))", fixed = TRUE)
  expect_equal(grep(r"(\tau^{(x_{3})}_{\varepsilon} &\sim)", l, fixed = TRUE), grep(r"(t_{\nu^{(x_{3})}})", l, fixed = TRUE) + 1L)
  w <- fxTips("priors")
  expect_match(tipArg(w, "size of x2"), "note=prior on the precision: Gamma with shape 1 and rate 0.01 (INLA loggamma)", fixed = TRUE)
  expect_match(tipArg(w, "precision of x3"), "Gamma with shape 1 and rate 5e-05 (INLA loggamma)", fixed = TRUE)
  expect_match(tipArg(w, "precision of x3"), ",note=SD [0-9.]+ \\([0-9.]+\\)\\. prior on the precision")
  expect_match(tipArg(w, "dof of x3"), "note=prior: INLA default$")
  for (e in c("(Intercept)", "t"))
    expect_match(tipArg(w, paste0("precision of random effect ", e, " | id in x3")), "Gamma with shape 1 and rate 0.01 (INLA loggamma)", fixed = TRUE)
  expect_equal(unique(w$sym[startsWith(w$arg, "par=gamma prior,")]), r"(\mathrm{Gamma})")
})

test_that("latexEquation(): normalCoef(0, 2) on scaled(x1): the gamma line and its note; defaults marked", {
  m <- fxModel("priors")
  l <- fxTex("priors")
  # precision 1 / 2^2 = 0.25
  g <- texLine(l, r"(\gamma_{\mathit{yg},x_{1}} &\sim)")
  expect_equal(g, r"(\gamma_{\mathit{yg},x_{1}} &\sim \mathcal{N}\left(0, \left(0.25\right)^{-1}\right))")
  expect_false(any(grepl(r"(\gamma_{\mathit{yg},x_{2}} &\sim)", l, fixed = TRUE)))
  expect_false(any(grepl(r"(\gamma_{\mathit{yg},x_{3}} &\sim)", l, fixed = TRUE)))
  # after the outcome's distributions, before the u lines
  expect_true(grep(r"(\gamma_{\mathit{yg},x_{1}} &\sim)", l, fixed = TRUE) > grep(r"(b^{(\mathit{yg})}_{0,l} &\sim)", l, fixed = TRUE))
  expect_true(all(grep("^u\\^", l) > grep(r"(\gamma_{\mathit{yg},x_{1}} &\sim)", l, fixed = TRUE)))
  w <- fxTips("priors")
  expect_match(tipArg(w, "association gamma (x1 -> yg)"), "note=normal prior on the coefficient with mean 0 and SD 2$")
  expect_match(tipArg(w, "association gamma (x2 -> yg)"), "note=normal prior on the coefficient with mean 0 and SD 10 (default)", fixed = TRUE)
  expect_equal(unique(w$sym[startsWith(w$arg, "par=normal prior,")]), r"(\mathcal{N})")
})

test_that("latexEquation(annotate = 'inline'): \\underset for parameters only; same equations underneath; digits", {
  for (nm in c("beta", "multi", "surv", "priors")) {
    m <- fxModel(nm)
    tex <- latexEquation(m, annotate = "inline")
    expectWellFormed(tex)
    expect_match(tex, "\\underset{\\substack{\\scriptstyle ", fixed = TRUE)
    expect_false(grepl("htmlData", tex, fixed = TRUE))
    expect_identical(unwrapCommand(tex, "\\underset"), unclass(latexEquation(m)))
    # covariates, u, the normal distribution name and D are not annotated
    expect_false(grepl("}{t_{i}}", tex, fixed = TRUE))
    expect_false(grepl("}{u^{", tex, fixed = TRUE))
    expect_false(grepl("}{\\mathcal{N}}", tex, fixed = TRUE))
    expect_false(grepl("}{\\mathcal{D}_", tex, fixed = TRUE))
    expect_null(latexState$est)
  }
  m <- fxModel("beta")
  fe <- fixedEffects(m)
  r <- fe[fe$response == "y" & fe$term == "t", ]
  a <- associations(m)
  for (d in c(2L, 4L)) {
    tex <- latexEquation(m, annotate = "inline", digits = d)
    expect_match(tex, paste0("\\underset{\\substack{\\scriptstyle ", numLatex(r$mean, d), " \\\\ \\scriptstyle (",
                             numLatex(r$sd, d), ")}}{\\beta^{(y)}_{1}}"), fixed = TRUE)
    expect_match(tex, paste0("\\underset{\\substack{\\scriptstyle ", numLatex(a$mean[1], d), " \\\\ \\scriptstyle (",
                             numLatex(a$sd[1], d), ")}}{\\gamma_{y,x_{1}}}"), fixed = TRUE)
  }
  # a random effect shows its SD; the Q matrix one line per item
  tex <- latexEquation(m, annotate = "inline")
  expect_match(tex, "\\\\substack\\{\\\\scriptstyle \\\\text\\{SD\\}\\\\ [0-9.]+ \\\\\\\\ \\\\scriptstyle \\([0-9.]+\\)\\}\\}\\{b\\^")
  expect_match(tex, "\\\\text\\{SD \\(Intercept\\)\\}\\\\ [0-9.]+\\\\,\\([0-9.]+\\) \\\\\\\\ \\\\scriptstyle \\\\text\\{SD t\\}")
  expect_match(tex, "\\\\text\\{corr \\(Intercept\\):t\\}\\\\ -?[0-9.]+\\\\,\\([0-9.]+\\)\\}\\}\\{\\\\mathbf\\{Q\\}")
})

test_that("estimateLookup(): keys of every symbol of the 'beta' model", {
  m <- fxModel("beta")
  k <- estimateLookup(m)
  expected <- c("fixed|y|(Intercept)", "fixed|y|t", "fixed|x1|(Intercept)", "fixed|x1|t", "assoc|y|x1",
                "family|y|phi", "family|x1|prec", "prec|x1|id|(Intercept)", "prec|x1|id|t", "resp|y", "mean|y",
                "resp|x1", "mean|x1", "b|x1|id|(Intercept)", "b|x1|id|t", "nlev|x1|id", "block|x1|id|(Intercept)|t",
                "eps|x1", "cov|t", "u|x1|y", "utau|x1|y", "D|y|x1")
  expect_true(all(expected %in% names(k)))
  expect_true(all(names(distributionNotes()) %in% names(k)))
  b <- k[["block|x1|id|(Intercept)|t"]]
  expect_length(b$mean, 3L)                                  # two SDs + one correlation
  expect_equal(b$item, c("SD (Intercept)", "SD t", "corr (Intercept):t"))
  expect_equal(k[["eps|x1"]]$item, "SD")
  expect_equal(k[["resp|y"]]$band, "range")
  expect_equal(k[["fixed|y|t"]]$band, "95%")
  fe <- fixedEffects(m)
  expect_equal(k[["fixed|y|t"]]$mean, fe$mean[fe$response == "y" & fe$term == "t"])
})

# ---------------------------------------------------------------------------
# KaTeX validation (node + vendored KaTeX; one node call per test)

test_that("latexEquation(): equations of unfitted models and special names are valid KaTeX", {
  skipWithoutKatex()
  texs <- list()
  for (nm in c("beta", "multi", "surv", "priors")) texs[[nm]] <- latexEquation(fxModel(nm, fit = FALSE))
  texs$names <- latexEquation(parseJsmFormula(my.var_2 ~ x_a + a.b + z9 + B12 + (1 | site:my_id), part = "long"))
  fg <- collectFormulas(y ~ 1 + scaled(x), x ~ t + site + (1 | id))
  texs$factorZero <- latexEquation(fg, prepareData(fg, fxSmallData()))
  texs$factorZeroNoData <- latexEquation(fg)
  texs$wishart3 <- latexEquation(parseJsmFormula(y ~ t + re(1 + t + I(t^2) | id, prior = wishartPrec()), part = "long"))
  texs$uncorrelated <- latexEquation(parseJsmFormula(y ~ t + re(1 + t || id, prior = halfCauchySd(2)), part = "long"))
  texs$survivalPair <- latexEquation(collectFormulas(
    list(time | event(status) ~ 1 + scaled(x1), w | event(ev) ~ 1 + scaled(x1)), list(x1 ~ t + (1 | id)),
    family = list(time = weibullSurvFamily(), w = exponentialSurvFamily())))
  texs$families <- famTex(cnt | exposure(pt) ~ t + I(t^2), zeroInflatedNegBinomialFamily(0, size = gammaPrec(2, 0.5)))
  out <- katexCheck(texs)
  expect_equal(unname(out), rep("OK", length(texs)))
})

for (nm in c("beta", "multi", "surv", "priors")) {
  test_that(paste0("latexEquation(annotate = ): annotated equations of the fitted '", nm, "' fixture are valid KaTeX"), {
    skipWithoutKatex()
    m <- fxModel(nm)
    texs <- lapply(stats::setNames(nm = c("none", "tooltip", "inline")), function(a) latexEquation(m, annotate = a))
    out <- katexCheck(texs)
    expect_equal(unname(out), rep("OK", length(texs)))
  })
}

# ---------------------------------------------------------------------------
# Escaping of special characters (B16)

test_that("escapeLatex(): \\text{} escapes all LaTeX special characters (factor levels, inline items)", {
  # The escaped strings: text mode (\text{}) and math mode (\mathrm{} etc.) differ for ^ ~ \
  # because KaTeX accepts \textasciicircum / \textasciitilde / \textbackslash only in text mode.
  expect_equal(escapeLatex("b^2~c\\d_e", "text"),
               r"(b\textasciicircum{}2\textasciitilde{}c\textbackslash{}d\_e)")
  expect_equal(escapeLatex("b^2~c\\d_e"), r"(b\hat{}2\tilde{}c\backslash{}d\_e)")
  expect_equal(escapeLatex("a_b&c%d$e#f{g}", "text"), r"(a\_b\&c\%d\$e\#f\{g\})")
  expect_equal(escapeLatex(c("a", NA, "")), c("a", NA, ""))
  expect_equal(escapeLatex(character()), character())
  # Not idempotent by design: call sites pass raw strings only (a second call escapes the escapes)
  expect_equal(escapeLatex(escapeLatex("a_b")), r"(a\backslash{}\_b)")

  d <- data.frame(y = seq_len(6) / 10, g = factor(rep(c("a", "b^2", "c\\d~e"), 2)))
  tex <- latexEquation(parseJsmFormula(y ~ g, part = "long"), d)
  expect_false(grepl("\\text{b^2}", tex, fixed = TRUE))
  expect_false(grepl("\\text{c\\d~e}", tex, fixed = TRUE))
  expect_match(tex, r"(\text{b\textasciicircum{}2})", fixed = TRUE)
  expect_match(tex, r"(\text{c\textbackslash{}d\textasciitilde{}e})", fixed = TRUE)
  expectWellFormed(tex)
  item <- inlineEstimate(estEntry("p", mean = c(1, 2), sd = c(0.1, 0.2), item = c("SD I(t^2)", "SD t")))
  expect_match(item, r"(\text{SD I(t\textasciicircum{}2)})", fixed = TRUE)
  one <- inlineEstimate(estEntry("p", mean = 1, sd = 0.1, item = "SD a\\b"))
  expect_match(one, r"(\text{SD a\textbackslash{}b})", fixed = TRUE)

  skipWithoutKatex()
  out <- katexCheck(list(levels = tex,
                         items = paste0("\\underset{\\substack{", item, "}}{b}"),
                         item1 = paste0("\\underset{\\substack{", one, "}}{b}"),
                         names = paste0("\\mathrm{", escapeLatex("b^2~c\\d_e"), "} + \\mathit{",
                                        escapeLatex("b^2~c\\d_e"), "}_{i} + \\operatorname{",
                                        escapeLatex("b^2~c\\d_e"), "}(x)")))
  expect_equal(unname(out), rep("OK", 4L))
})
