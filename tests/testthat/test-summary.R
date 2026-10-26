dm <- data.frame(id = sprintf("S%03d", rep(1:3, each = 4)), t = rep(c(0, 0.5, 1, 1.5), 3),
                 site = rep(c("a", "b"), 6), x1 = rnorm(12), x2 = rnorm(12), y = runif(12))
dLong <- dm
dLong$x1[2] <- NA
dLong$t[5] <- NA
dLong <- dLong[dLong$id != "S003", ]
dOut <- dm
dOut$t[7] <- NA
f <- collectFormulas(list(y ~ 1 + (1 | id) + (1 | site) + scaled(x1)),
                     list(x1 ~ 1 + t + (1 + t | id), x2 ~ t),
                     familyOut = c(y = "beta"))
s <- summarizeData(f, dOut, list(x1 = dLong, x2 = dm))

test_that("rows supplied, complete and dropped per equation", {
  expect_s3_class(s, "mvremDataSummary")
  e <- s$equations
  expect_equal(e$part, c("long", "long", "out"))
  expect_equal(e$response, c("x1", "x2", "y"))
  expect_equal(e$family, c("gaussian", "gaussian", "beta"))
  expect_equal(e$rows, c(8L, 12L, 12L))
  expect_equal(e$complete, c(6L, 12L, 11L))
  expect_equal(e$dropped, c(2L, 0L, 1L))
  expect_equal(e$groups, c("id: 2", "", "id: 3, site: 2"))
})

test_that("complete rows agree with prepareData(na.action = na.omit)", {
  dat <- prepareData(dOut, list(x1 = dLong, x2 = dm), f)
  expect_equal(s$equations$complete,
               c(nrow(dat$long$x1), nrow(dat$long$x2), nrow(dat$out$y)))
})

test_that("missing values per variable", {
  m <- s$missing
  get <- function(r, v) m$missing[m$response == r & m$variable == v]
  expect_equal(get("x1", "x1"), 1L)
  expect_equal(get("x1", "t"), 1L)
  expect_equal(get("x1", "id"), 0L)
  expect_equal(get("y", "t"), 1L)          # link variable of scaled(x1) is included
  expect_equal(get("y", "y"), 0L)
  expect_equal(sum(m$missing[m$response == "x2"]), 0L)
})

test_that("zero equations per scaled() pair", {
  z <- s$zero
  expect_equal(nrow(z), 1L)
  expect_equal(unlist(z[1, c("long", "out", "vars")], use.names = FALSE), c("x1", "y", "id, t"))
  expect_equal(z$rows, 11L)
  # duplicated outcome rows do not add zero equations
  expect_equal(summarizeData(f, rbind(dOut, dOut), list(x1 = dLong, x2 = dm))$zero$rows, 11L)
})

test_that("print shows the tables; models without scaled()", {
  expect_output(print(s), "Equations")
  expect_output(print(s), "Zero equations")
  expect_output(print(s), "Missing values per variable")
  s0 <- summarizeData(collectFormulas(y ~ t, x ~ t), data.frame(y = 1:3, t = 1:3, x = 1:3),
                      data.frame(y = 1:3, t = 1:3, x = 1:3))
  expect_equal(nrow(s0$zero), 0L)
  expect_output(print(s0), "none")
  expect_error(summarizeData(f, dOut, list(x1 = dLong)), "No data frame")
})
