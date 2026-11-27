# E19: the time variable of a model (jsm(time = ), timeVariable(), the plots' default).

timeData <- function() {
  d <- fxData()$long
  d$age <- 50 + as.numeric(factor(d$id))            # a second numeric covariate
  d
}
# the model of these tests (age in both formulae, t only in the marker's), unfitted, with
# further jsm() arguments (time = , naming = , ...)
timeModel <- function(d = timeData(), ...)
  jsm(list(y ~ age + scaled(x1), x1 ~ t + age + (1 | id)), data = d, family = list(y = betaFamily()),
      fit = FALSE, ...)

test_that("jsm(time = ): stored, printed and used by the plots", {
  d <- timeData()
  m <- timeModel(d, time = "t")
  tv <- timeVariable(m)
  expect_equal(as.character(tv), "t")
  expect_equal(attr(tv, "source"), "given")
  expect_output(print(m), "Time variable: t\n")
  expect_identical(linesData(m, n = 5)$time, "t")
  # the plot's own `time` still overrides the model's
  expect_identical(linesData(m, time = "age", n = 5)$time, "age")
})

test_that("jsm() without time: the numeric covariate of most formulae, with a message when ambiguous", {
  d <- timeData()
  # age is in both formulae, t only in one: age wins, and the message names t
  expect_message(
    m <- timeModel(d),
    "Time variable: 'age'.*other candidates: t.*jsm\\(time = \\)")
  expect_equal(as.character(timeVariable(m)), "age")
  expect_equal(attr(timeVariable(m), "source"), "detected")
  expect_output(print(m), "Time variable: age \\(detected; set with jsm\\(time = \\)\\)")
  # update() declares it
  m2 <- update(m, time = "t")
  expect_equal(as.character(timeVariable(m2)), "t")
  # a single candidate: no message
  expect_silent(m3 <- jsm(list(y ~ t + scaled(x1), x1 ~ t + (1 | id)), data = d,
                          family = list(y = betaFamily()), fit = FALSE))
  expect_equal(as.character(timeVariable(m3)), "t")
  # no numeric covariate: NA, nothing printed, the plots still ask for `time`
  m4 <- jsm(list(y ~ 1 + scaled(x1), x1 ~ 1 + (1 | id)), data = d,
            family = list(y = betaFamily()), fit = FALSE)
  expect_true(is.na(timeVariable(m4)))
  expect_equal(attr(timeVariable(m4), "source"), "none")
  expect_false(any(grepl("Time variable", capture.output(print(m4)))))
  expect_error(linesData(m4, n = 5), "No numeric covariate")
})

test_that("jsm(time = ): clear errors", {
  d <- timeData()
  f <- list(y ~ age + scaled(x1), x1 ~ t + age + (1 | id))
  fam <- list(y = betaFamily())
  expect_error(jsm(f, data = d, family = fam, time = "id", fit = FALSE),
               "`time` = 'id' must be a numeric variable")
  expect_error(jsm(f, data = d, family = fam, time = "y", fit = FALSE),
               "`time` = 'y' is not a numeric covariate of any formula \\(candidates: age, t\\)")
  expect_error(jsm(f, data = d, family = fam, time = "nope", fit = FALSE),
               "'nope' is not a variable of the model's data")
  expect_error(jsm(f, data = d, family = fam, time = c("t", "age"), fit = FALSE),
               "must be the name of one variable")
})

test_that("timeVariable(): models built before E19 detect it", {
  m <- fxModel("beta", fit = FALSE)
  old <- m
  old$time <- NULL
  expect_equal(timeVariable(old), timeVariable(m))
  expect_equal(as.character(timeVariable(m)), "t")
})

test_that("timeVariable(): a survival model takes the link variable of the outcome", {
  m <- fxModel("surv", fit = FALSE)
  expect_equal(as.character(timeVariable(m)), "t")
})

test_that("timeCandidates(): a numeric grouping variable is not a time candidate (B27)", {
  d <- fxData()$long
  d$id <- as.integer(factor(d$id))
  expect_message(m <- jsm(list(y ~ t + scaled(x1)), list(x1 ~ t + (1 | id)), data = d,
                          family = list(y = betaFamily()), fit = FALSE), NA)
  expect_equal(m$time$candidates, "t")
  expect_error(jsm(list(y ~ t + scaled(x1)), list(x1 ~ t + (1 | id)), data = d, time = "id", fit = FALSE),
               "must be a numeric variable")
})
