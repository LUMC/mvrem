# Helpers of test-families-priors.R: response checks on a standalone parsed formula,
# the half-Cauchy INLA expression as an R function, and the data-scaled Wishart prior
# (blocks, memoised lme4 reference fits and their regularised covariance matrices).

# checkResponse() on a standalone parsed outcome formula `lhs` with family `fam`
checkResp <- function(fam, lhs, df) {
  p <- parseJsmFormula(lhs, "out")
  p$family <- fam
  checkResponse(p, df)
}

# the halfCauchySd() INLA expression as an R function of theta = log precision
halfCauchyLogDens <- function(s) {
  body <- sub("^expression:", "", halfCauchyExpression(s))
  eval(parse(text = paste0("function(theta) {", body, "}")))
}

# the `hyper` argument of an f() term (embedded as the call that builds it)
hyperOf <- function(term) eval(as.list(term)$hyper)

# evaluate `expr`, returning list(value, warnings = the messages of the warnings raised)
withWarnings <- function(expr) {
  ws <- character()
  val <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}

# ---------------------------------------------------------------------------
# data-scaled Wishart prior

# a VarCorr block as a plain numeric matrix (no dimnames / attributes)
plainMatrix <- function(m) matrix(as.numeric(m), nrow(m))

dataPrior <- function() wishartPrec(scale = "data")
eff2 <- c("(Intercept)", "t")       # the effects of a `1 + t | id` block

# parsed long formula + prepared data of a data-scaled block on the fixture data
# (fxData()$long: 40 subjects, enough for lme4; fxSmallData() has only 4)
dataBlock <- function(formLong, family, d = fxData()$long) {
  resp <- all.vars(formLong[[2]])[1]
  fam <- stats::setNames(list(family), resp)
  formOut <- stats::as.formula(paste0("y ~ t + scaled(", resp, ")"))
  f <- collectFormulas(list(formOut), list(formLong), family = fam)
  list(p = f$long[[resp]], df = prepareData(f, d)$long[[resp]])
}

# lme4 reference fits, each fitted at most once per test run (cached in `.fxFam`); the
# singular-fit messages and warnings of lme4 are not under test here
.fxFam <- new.env(parent = emptyenv())
refFit <- function(name) {
  if (is.null(.fxFam[[name]])) {
    d <- fxData()$long
    quietly <- function(x) suppressWarnings(suppressMessages(x))
    .fxFam[[name]] <- switch(name,
      lmerX1 = quietly(lme4::lmer(x1 ~ t + (1 + t | id), d)),
      lmerLogPt = quietly(lme4::lmer(log(pt) ~ t + (1 + t | id), d)),
      glmerBinomProbit = quietly(lme4::glmer(cbind(k, n - k) ~ t + (1 + t | id), d,
                                             family = stats::binomial(link = "probit"))),
      glmerPois = quietly(lme4::glmer(cnt ~ t + (1 + t | id) + offset(log(pt)), d, family = stats::poisson())),
      lmerX1X3 = quietly(lme4::lmer(x1 ~ t + (1 + t | id) + (0 + x3 | id), d)),
      lmerSmallX = quietly(lme4::lmer(x ~ t + (1 + t | id), fxSmallData())),
      stop("unknown reference fit ", name))
  }
  .fxFam[[name]]
}

# the VarCorr block of a reference fit whose effects are `effects` (in that order)
vcBlock <- function(fit, effects = eff2) {
  for (m in lme4::VarCorr(fit))
    if (nrow(m) == length(effects) && setequal(rownames(m), effects))
      return(plainMatrix(as.matrix(m)[effects, effects, drop = FALSE]))
  stop("no random-effect block with effects ", paste(effects, collapse = ", "))
}

# The covariance blockCovariance() should give for a reference fit: its VarCorr block,
# regularised as blockCovariance() does (a singular block changes and warns).
refBlock <- function(fit, effects = eff2) {
  r <- withWarnings(regularizeCovariance(vcBlock(fit, effects), "reference"))
  list(cov = r$value, regularised = length(r$warnings) > 0L)
}

# blockCovariance() of the data-scaled block `b` (a dataBlock()) compared with the
# reference `ref` (a refBlock()): same matrix, the regularisation warning exactly when the
# reference block was regularised, no other warnings than lme4's own (passed on)
expectBlockCov <- function(b, ref, effects = eff2) {
  r <- withWarnings(blockCovariance(b$p, b$df, effects, b$p$re[[1]]))
  expect_equal(r$value, ref$cov, tolerance = 1e-4)
  expect_equal(any(grepl("singular; its covariance matrix was regularised", r$warnings)), ref$regularised)
  expect_true(all(grepl("regularised|^lme4 fit for", r$warnings)))
  r$value
}
