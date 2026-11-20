# Leftovers of the old jsm() / results suites: plot helpers that test-plots.R (the
# plots area) does not cover yet. Everything else moved to test-fit.R (fixtures:
# helper-fixtures.R).

test_that("plotTimeVariable() / plotGroupVariable() pick the defaults of plot()", {
  m <- fxModel("multi", fit = FALSE)
  forms <- c(m$forms$out, m$forms$long)
  data <- c(m$data$out, m$data$long)
  expect_equal(plotTimeVariable(forms, data, NULL), "t")
  expect_equal(plotGroupVariable(forms, data, NULL), "id")
})

test_that("plotModel(detail = ) rejects an unknown detail", {
  expect_error(plotModel(fxModel("multi", fit = FALSE), detail = "other"), "should be one of")
})

test_that("ellipsePoints(): the 95% ellipse of a bivariate normal with the given SDs", {
  e <- ellipsePoints(2, 1, 0)
  expect_equal(max(e$x), 2 * sqrt(qchisq(0.95, 2)), tolerance = 1e-3)
  expect_equal(max(e$y), sqrt(qchisq(0.95, 2)), tolerance = 1e-3)
})
