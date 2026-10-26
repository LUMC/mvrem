forms <- collectFormulas(
  formOut = list(z ~ t + (1 | id), log(y) ~ t + re(1 | id) + scaled(x1)),
  formLong = list(x1 ~ t + re(1 + t | id), x2 ~ I(t^2) + re(1 | id))
)
d <- data.frame(id = rep(1:3, each = 2), t = rep(0:1, 3),
                z = 1:6, y = runif(6), x1 = rnorm(6), x2 = rnorm(6),
                unused = letters[1:6])

test_that("single data frame is shared by all formulae", {
  res <- prepareData(d, d, forms)
  expect_named(res$out, c("z", "y"))
  expect_named(res$long, c("x1", "x2"))
  for (df in c(res$out, res$long)) expect_equal(df$.key, paste0("K", 1:6))
})

test_that("only the formula's variables are selected", {
  res <- prepareData(d, d, forms)
  expect_named(res$out$y, c(".key", "y", "t", "id"))
  expect_named(res$long$x2, c(".key", "x2", "t", "id"))
  expect_equal(res$out$y$y, d$y)                  # untransformed response
  expect_equal(res$long$x1$t, d$t)
})

test_that("grouping factors are built from the bar notation", {
  f <- collectFormulas(list(y ~ t + (1 | site:id) + scaled(x1)), list(x1 ~ t + re(1 + t || id)))
  ds <- cbind(d, site = rep(c("a", "b"), 3))
  res <- prepareData(ds, ds, f)
  expect_named(res$out$y, c(".key", "y", "t", "site", "id", "site:id"))
  expect_s3_class(res$out$y[["site:id"]], "factor")
  expect_equal(as.character(res$out$y[["site:id"]]), paste(ds$site, ds$id, sep = ":"))
  expect_s3_class(res$long$x1$id, "factor")      # grouping variable replaced by factor
  expect_equal(levels(res$long$x1$id), c("1", "2", "3"))
})

test_that("formulae without random effects are supported", {
  f <- collectFormulas(list(y ~ t + scaled(x1)), list(x1 ~ t + (1 | id)))
  res <- prepareData(d, d, f)
  expect_named(res$out$y, c(".key", "y", "t", "id"))   # id: link variable of scaled(x1)
  expect_equal(nrow(res$out$y), 6L)
})

test_that("na.action controls rows with missing values", {
  dm <- d
  dm$x1[2] <- NA
  dm$t[5] <- NA
  dm$unused[1] <- NA                               # not a model variable: no effect

  res <- prepareData(d, dm, forms)                 # default na.omit
  expect_equal(res$long$x1$.key, paste0("K", c(1, 3, 4, 6)))
  expect_equal(res$long$x2$.key, paste0("K", c(1:4, 6)))
  expect_equal(res$out$y$.key, paste0("K", 1:6))

  expect_equal(prepareData(d, dm, forms, na.action = "na.omit"), res)
  expect_error(prepareData(d, dm, forms, na.action = stats::na.fail), "missing values")
  expect_error(prepareData(d, dm, forms, na.action = "na.fail"), "missing values")
  expect_no_error(prepareData(d, d, forms, na.action = stats::na.fail))
})

test_that("named list of data frames is matched by response", {
  dl <- list(x2 = d[1:4, ], x1 = d[3:6, ])
  res <- prepareData(d, dl, forms)
  expect_named(res$long, c("x1", "x2"))
  expect_equal(res$long$x1$x1, d$x1[3:6])
  expect_equal(res$long$x1$.key, paste0("K", 3:6))
})

test_that("list names use the response variable, not the transformed lhs", {
  res <- prepareData(list(z = d, y = d[1:3, ]), d, forms)
  expect_equal(nrow(res$out$y), 3)
  expect_error(prepareData(list(z = d, `log(y)` = d), d, forms), "No data frame .* y")
})

test_that("keys use row names when provided", {
  dn <- d
  rownames(dn) <- paste0("obs", 6:1)
  res <- prepareData(dn, d, forms)
  expect_equal(res$out$y$.key, paste0("Kobs", 6:1))
  expect_equal(res$long$x1$.key, paste0("K", 1:6))

  # tibble-like / non-unique row names are rejected
  bad <- d
  attr(bad, "row.names") <- c("a", "a", "b", "c", "d", "e")
  expect_error(prepareData(bad, d, forms), "must be unique")
})

test_that("row names are kept in .key for every data frame and subset", {
  dn <- d
  rownames(dn) <- c("a1", "b2", "c3", "d4", "e5", "f6")
  res <- prepareData(dn, list(x1 = dn[c(2, 5), ], x2 = dn), forms)

  for (df in c(res$out, res$long)) expect_type(df$.key, "character")
  expect_equal(res$out$z$.key, paste0("K", rownames(dn)))
  expect_equal(res$long$x1$.key, c("Kb2", "Ke5"))
  # same row -> same key across formulae
  i <- match(res$long$x1$.key, res$out$y$.key)
  expect_equal(res$out$y$t[i], res$long$x1$t)
  expect_equal(res$out$y$y[i], dn$y[c(2, 5)])
})

test_that("invalid data input errors", {
  expect_error(prepareData(d, list(x1 = d), forms), "No data frame .* x2")
  expect_error(prepareData(d, list(x1 = d, x2 = d, x9 = d), forms), "do not correspond")
  expect_error(prepareData(d, list(d, d), forms), "unique names")
  expect_error(prepareData(d, list(x1 = d, x2 = 1), forms), "not data frames")
  expect_error(prepareData(d, 1, forms), "must be a data frame")
  expect_error(prepareData(d[, -2], d, forms), "not found in `dataOut`: t")
  expect_error(prepareData(cbind(d, .key = 1), d, forms), "reserved")
})
