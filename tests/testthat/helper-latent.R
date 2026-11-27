# Fixtures of test-latent.R (E38, latent equations). Read-only for the tests.
#
#   ltData()          simulated data (60 subjects x 5 visits): markers x1 (correlated
#                     intercept + slope), x2 (intercept); latent
#                     eta = 0.5 t + v_id + 0.8 m1 - 0.5 m2 (m_k = marker k's predictor,
#                     v_id ~ N(0, 0.3^2)); outcomes z = 1 + eta + e (sd 0.3),
#                     y = -0.5 + 0.7 eta + e (sd 0.3); attribute `truth`
#   ltForms(coef)     the formulae: z links with scaled(eta, coef = 1) (scale fixed) when
#                     `coef`, else freely
#   ltModel(fit)      jsm() of ltForms(TRUE) on ltData(), criteria = TRUE (fitted once per run)

.lt <- new.env(parent = emptyenv())

ltData <- function() {
  if (!is.null(.lt$data)) return(.lt$data)
  set.seed(38)
  n <- 60
  ids <- sprintf("P%02d", seq_len(n))
  d <- data.frame(id = rep(ids, each = 5), t = rep(seq(0, 1, by = 0.25), n), stringsAsFactors = FALSE)
  i <- match(d$id, ids)
  b0 <- rnorm(n, sd = 0.6); b1 <- rnorm(n, sd = 0.4); c0 <- rnorm(n, sd = 0.5); v <- rnorm(n, sd = 0.3)
  m1 <- 0.2 + 0.5 * d$t + b0[i] + b1[i] * d$t
  m2 <- -0.1 - 0.3 * d$t + c0[i]
  eta <- 0.5 * d$t + v[i] + 0.8 * m1 - 0.5 * m2
  d$x1 <- m1 + rnorm(nrow(d), sd = 0.2)
  d$x2 <- m2 + rnorm(nrow(d), sd = 0.2)
  d$z <- 1 + eta + rnorm(nrow(d), sd = 0.3)
  d$y <- -0.5 + 0.7 * eta + rnorm(nrow(d), sd = 0.3)
  attr(d, "truth") <- list(gamma1 = 0.8, gamma2 = -0.5, beta = 0.7, slope = 0.5, sdV = 0.3)
  .lt$data <- d
  d
}

ltForms <- function(coef = TRUE) {
  list(
    if (coef) z ~ 1 + scaled(eta, coef = 1) else z ~ 1 + scaled(eta),
    y ~ 1 + scaled(eta),
    latent(eta) ~ 1 + t + (1 | id) + scaled(x1) + scaled(x2),
    x1 ~ 1 + t + (1 + t | id),
    x2 ~ 1 + t + (1 | id))
}

ltModel <- function(fit = TRUE) {
  if (!fit) return(jsm(ltForms(), data = ltData(), fit = FALSE))
  skip_if_not_installed("INLA")
  if (is.null(.lt$model))
    .lt$model <- suppressMessages(jsm(ltForms(), data = ltData(), criteria = TRUE, threads = 1))
  skip_if(is.null(inlaFit(.lt$model)), "INLA fit of the latent fixture failed")
  .lt$model
}
