# Helpers of test-fit.R (fitted models, result tables, predict()).

statCols <- c("mean", "sd", "q0.025", "q0.5", "q0.975")

# Posterior summaries present, SDs non-negative, quantiles ordered.
expectOrderedStats <- function(df, tol = 1e-8) {
  expect_true(all(statCols %in% names(df)))
  expect_true(all(df$sd >= 0))
  expect_true(all(df$q0.025 <= df$q0.5 + tol))
  expect_true(all(df$q0.5 <= df$q0.975 + tol))
}

# Sizes of the fixture data (fxData()$long), from which expected counts are derived.
fxCounts <- function() {
  d <- fxData()$long
  list(rows = nrow(d), subjects = length(unique(d$id)), sites = nlevels(d$site),
       responses = c("k", "cnt", "x1", "x2"))              # the "multi" fixture, in the order passed
}

# Mock INLA out for the rest of the calling test: runInla() records the model it would have
# run (its options and the inla() arguments) and returns the model's current fit, so that
# refit() / predict() go through their whole path without an INLA run. The returned
# environment holds `options`, `args` and `model` of the last (mocked) run.
mockInlaRun <- function(env = parent.frame()) {
  cap <- new.env()
  testthat::local_mocked_bindings(runInla = function(x) {
    cap$options <- x$options
    cap$args <- inlaArgumentList(x$inla, x$options)
    cap$model <- x
    x$fit
  }, .package = "mvrem", .env = env)
  cap
}

# A copy of an unfitted model marked as fitted with an empty fit object: only for argument
# checks that run before anything of the fit is used.
stubFitted <- function(m) {
  m$fit <- structure(list(), class = "inla")
  m
}
