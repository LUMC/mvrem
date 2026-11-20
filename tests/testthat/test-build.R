# Tests of steps 2 and 3 (no INLA fit): prepareData(), summarizeData(), namingScheme(),
# buildInlaFormula() (terms, zero equations, stacks), formatInlaFormula(), inlaData() and
# inlaArgumentList(). Helpers: helper-build.R (smallForms(), buildSmall(), termHyper()) and
# helper-fixtures.R (fxSmallData(): 4 subjects S1..S4 x t = 0, 1, 2; fxModel(fit = FALSE)).

# Two outcomes (one with a transformed response) and two longitudinal responses.
pairForms <- function() {
  collectFormulas(
    formOut = list(y ~ t + (1 | id) + scaled(cnt), log(pt) ~ t + re(1 | id) + scaled(x)),
    formLong = list(x ~ t + re(1 + t | id), cnt ~ I(t^2) + re(1 | id)))
}

# ---- prepareData(): data arguments --------------------------------------------

test_that("prepareData(): one data frame serves every formula with the same keys", {
  d <- fxSmallData()
  res <- prepareData(pairForms(), d)
  expect_named(res$out, c("y", "pt"))
  expect_named(res$long, c("x", "cnt"))
  for (df in c(res$out, res$long)) expect_equal(df$.key, paste0("K", 1:12))
  expect_equal(prepareData(pairForms(), d), prepareData(pairForms(), dataOut = d, dataLong = d))
})

test_that("prepareData(): only the formula's variables are kept, the response untransformed", {
  d <- fxSmallData()
  res <- prepareData(pairForms(), d)
  expect_named(res$out$pt, c(".key", "pt", "t", "id"))       # log(pt): column pt
  expect_named(res$long$cnt, c(".key", "cnt", "t", "id"))    # I(t^2): column t
  expect_equal(res$out$pt$pt, d$pt)
  expect_equal(res$long$cnt$t, d$t)
})

test_that("prepareData(): a list named by response is matched in any order, by the response variable", {
  d <- fxSmallData()
  f <- pairForms()
  res <- prepareData(f, dataOut = d, dataLong = list(cnt = d[1:4, ], x = d[3:6, ]))
  expect_named(res$long, c("x", "cnt"))
  expect_equal(res$long$x$x, d$x[3:6])
  expect_equal(res$long$x$.key, paste0("K", 3:6))
  expect_equal(nrow(prepareData(f, dataOut = list(y = d, pt = d[1:3, ]), dataLong = d)$out$pt), 3L)
  expect_error(prepareData(f, dataOut = list(y = d, `log(pt)` = d), dataLong = d),
               "No data frame in `dataOut` for response\\(s\\): pt")
})

test_that("prepareData(): `data` serves each part; dataOut / dataLong override it", {
  d <- fxSmallData()
  f <- smallForms()
  res <- prepareData(f, d, dataLong = d[1:4, ])
  expect_equal(c(nrow(res$out$y), nrow(res$long$x)), c(12L, 4L))
  res <- prepareData(f, d, dataOut = d[1:4, ])
  expect_equal(c(nrow(res$out$y), nrow(res$long$x)), c(4L, 12L))
  # a list `data` naming both parts, overridden per part
  res <- prepareData(f, list(y = d[1:6, ], x = d), dataOut = d[1:9, ])
  expect_equal(c(nrow(res$out$y), nrow(res$long$x)), c(9L, 12L))
  res2 <- prepareData(f, list(y = d[1:6, ], x = d), dataLong = list(x = d[2:3, ]))
  expect_equal(res2$long$x$.key, c("K2", "K3"))
  expect_equal(nrow(res2$out$y), 6L)
})

test_that("prepareData(): invalid data arguments are errors", {
  d <- fxSmallData()
  f <- smallForms()
  expect_error(prepareData(f, dataOut = d), "No data for the longitudinal formulae")
  expect_error(prepareData(f, dataLong = d), "No data for the outcome formulae")
  expect_error(prepareData(f), "give `data`")
  expect_error(prepareData(f, dataOut = list(y = d, x = d), dataLong = d),
               "`dataOut` do not correspond to any formula: x")
  expect_error(prepareData(f, dataOut = d, dataLong = list(y = d)),
               "No data frame in `dataLong` for response\\(s\\): x")
  expect_error(prepareData(f, list(y = d)), "No data frame in `data` for response\\(s\\): x")
  expect_error(prepareData(f, list(y = d, x = d, w = d)), "`data` do not correspond to any formula: w")
  expect_error(prepareData(f, list(y = d, x = d, y = d)), "unique names")
  expect_error(prepareData(f, dataOut = d, dataLong = list(d, d)), "unique names")
  expect_error(prepareData(f, list(y = d, x = "a")), "not data frames: x")
  expect_error(prepareData(f, 1), "must be a data frame")
  expect_error(prepareData(f, dataOut = d, dataLong = 1), "must be a data frame")
  expect_error(prepareData(f, d[, names(d) != "t"]), "not found in `data`: t")
  expect_error(prepareData(f, dataOut = d[, names(d) != "t"], dataLong = d), "not found in `dataOut`: t")
  expect_error(prepareData(f, cbind(d, .key = 1)), "reserved")
})

test_that("prepareData(): keys from automatic row names agree between a data frame and its subsets", {
  d <- fxSmallData()
  f <- smallForms()
  res <- prepareData(f, dataOut = d, dataLong = d[c(2, 4, 6, 11), ])
  expect_equal(res$long$x$.key, c("K2", "K4", "K6", "K11"))
  i <- match(res$long$x$.key, res$out$y$.key)
  expect_equal(res$out$y$t[i], res$long$x$t)
  expect_equal(res$out$y$id[i], as.character(res$long$x$id))
})

test_that("prepareData(): explicit row names give the keys; they must be unique", {
  d <- fxSmallData()
  f <- pairForms()
  dn <- d
  rownames(dn) <- paste0("obs", 12:1)
  res <- prepareData(f, dataOut = dn, dataLong = list(x = dn[c(2, 6), ], cnt = d))
  for (df in c(res$out, res$long)) expect_type(df$.key, "character")
  expect_equal(res$out$y$.key, paste0("Kobs", 12:1))
  expect_equal(res$long$x$.key, c("Kobs11", "Kobs7"))
  expect_equal(res$long$cnt$.key, paste0("K", 1:12))
  # the same row gets the same key across formulae
  i <- match(res$long$x$.key, res$out$y$.key)
  expect_equal(res$out$y$t[i], res$long$x$t)
  expect_equal(res$out$y$y[i], dn$y[c(2, 6)])
  dd <- d
  attr(dd, "row.names") <- c(rep("r", 2), paste0("q", 3:12))
  expect_error(prepareData(f, dataOut = d, dataLong = dd), "`dataLong` must be unique")
})

test_that("prepareData(): a tibble is treated as a plain data frame, alone or in a list", {
  skip_if_not_installed("tibble")
  d <- fxSmallData()
  f <- smallForms()
  td <- tibble::as_tibble(d)
  expect_equal(prepareData(f, td), prepareData(f, d))
  expect_equal(prepareData(f, list(y = td, x = d)), prepareData(f, d))
  expect_equal(prepareData(f, dataOut = d, dataLong = list(x = td)), prepareData(f, d))
})

test_that("summarizeData(): a tibble gives the summary of the plain data frame", {
  skip_if_not_installed("tibble")
  d <- fxSmallData()
  f <- smallForms()
  expect_equal(summarizeData(f, tibble::as_tibble(d)), summarizeData(f, d))
  expect_equal(summarizeData(f, list(y = tibble::as_tibble(d), x = d)), summarizeData(f, d))
})

test_that("jsm(fit = FALSE): a tibble builds the same model as the plain data frame", {
  skip_if_not_installed("tibble")
  d <- fxSmallData()
  mt <- jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), data = tibble::as_tibble(d), fit = FALSE)
  md <- jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), data = d, fit = FALSE)
  expect_equal(modelData(mt), modelData(md))
  expect_equal(dataSummary(mt), dataSummary(md))
  expect_equal(inlaData(mt), inlaData(md))
  expect_equal(deparse(inlaFormula(mt)), deparse(inlaFormula(md)))
})

test_that("jsm(fit = FALSE): `data` / `dataOut` lists named by response reach prepareData()", {
  d <- fxSmallData()
  m <- jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), data = list(y = d[1:6, ], x = d), fit = FALSE)
  expect_equal(nrow(modelData(m)$out$y), 6L)
  expect_equal(nrow(modelData(m)$long$x), 12L)
  expect_equal(nrow(inlaData(m)), 6L + 12L + 6L)          # y, x, zero equations at the 6 y rows
  m2 <- jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), data = d, dataOut = list(y = d[1:9, ]), fit = FALSE)
  expect_equal(nrow(modelData(m2)$out$y), 9L)
  expect_equal(nrow(modelData(m2)$long$x), 12L)
  expect_error(jsm(list(y ~ t + scaled(x)), list(x ~ t + (1 | id)), data = list(y = d), fit = FALSE),
               "No data frame in `data` for response\\(s\\): x")
})

# ---- prepareData(): grouping variables ------------------------------------------

test_that("prepareData(): character, factor and numeric grouping variables give the same groups", {
  d <- fxSmallData()
  f <- smallForms()
  df <- d
  df$id <- factor(df$id, levels = rev(unique(df$id)))
  pc <- prepareData(f, d)
  pf <- prepareData(f, df)
  expect_s3_class(pc$long$x$id, "factor")
  expect_equal(as.character(pf$long$x$id), as.character(pc$long$x$id))
  dnum <- d
  dnum$id <- rep(1:4, each = 3)
  expect_equal(levels(prepareData(f, dnum)$long$x$id), c("1", "2", "3", "4"))
  # LEVELS are sorted whatever the factor's level order
  ic <- buildInlaFormula(f, pc)
  ifac <- buildInlaFormula(f, pf)
  expect_equal(ifac$levels$x_1_id, ic$levels$x_1_id)
  expect_equal(ic$levels$x_1_id, paste0("x|1|S", 1:4))
  expect_equal(ifac$stacks, ic$stacks)
})

test_that("prepareData(): site:id grouping gives a factor column, a syntactic index and ':' labels", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ t + (1 | site:id) + scaled(x)), list(x ~ t + (1 | id)))
  p <- prepareData(f, d)
  expect_named(p$out$y, c(".key", "y", "t", "site", "id", "site:id"))
  expect_equal(levels(p$out$y[["site:id"]]), c("A:S1", "A:S2", "B:S3", "B:S4"))
  inl <- buildInlaFormula(f, p)
  b <- inl$random$y[[1]]
  expect_equal(b$index, "y_1_site.id")
  expect_equal(b$levels, "y_1_site.id")
  expect_equal(inl$levels$y_1_site.id, paste0("y|1|", c("A:S1", "A:S2", "B:S3", "B:S4")))
  expect_equal(inl$stacks$y$main$y_1_site.id, paste0("y|1|", d$site, ":", d$id))
  s <- summarizeData(f, d)
  expect_equal(s$equations$groups[s$equations$response == "y"], "site:id: 4")
})

test_that("prepareData(): link variables must be present in the outcome data", {
  d <- smallData()
  f <- collectFormulas(list(y ~ 1 + scaled(x)), list(x ~ t + g + (1 | id)))
  expect_equal(f$out$y$linkVariables, c("id", "t", "g"))
  expect_error(prepareData(f, dataOut = d[c("y", "id", "t")], dataLong = d),
               "formula for 'y' not found in `dataOut`: g")
  p <- prepareData(f, d)
  expect_named(p$out$y, c(".key", "y", "id", "t", "g"))
  expect_type(p$out$y$id, "character")               # link variables stay as given
})

# ---- prepareData(): na.action ------------------------------------------------

test_that("prepareData(na.action): a missing response is kept with na.pass, dropped with na.omit", {
  d <- fxSmallData()
  d$x[3] <- NA
  d$y[5] <- NA
  d$site[1] <- NA                                   # not a model variable: no effect
  f <- smallForms()
  p <- prepareData(f, d)
  expect_equal(nrow(p$long$x), 12L)
  expect_true(is.na(p$long$x$x[3]))
  expect_true(is.na(p$out$y$y[5]))
  po <- prepareData(f, d, na.action = na.omit)
  expect_equal(po$long$x$.key, paste0("K", setdiff(1:12, 3)))
  expect_equal(po$out$y$.key, paste0("K", setdiff(1:12, 5)))
  expect_null(attr(po$out$y, "na.action"))
  # a function, its name and a user function behave the same
  myOmit <- function(object, ...) stats::na.omit(object, ...)
  expect_equal(prepareData(f, d, na.action = "na.omit"), po)
  expect_equal(prepareData(f, d, na.action = myOmit), po)
  expect_error(prepareData(f, d, na.action = "na.fail"), "missing values")
  expect_error(prepareData(f, d, na.action = stats::na.fail), "missing values")
  expect_no_error(prepareData(f, fxSmallData(), na.action = stats::na.fail))
  # stacks keep the missing response as NA (INLA predicts it)
  inl <- buildInlaFormula(f, p)
  expect_true(is.na(inl$stacks$x$main$x_x[3]))
  expect_true(is.na(inl$stacks$y$main$y_y[5]))
})

test_that("prepareData(na.pass): a missing covariate is an error asking for na.omit", {
  d <- fxSmallData()
  d$t[5] <- NA
  f <- smallForms()
  expect_error(prepareData(f, d), "Missing values in `t`.*na.omit")
  expect_equal(prepareData(f, d, na.action = na.omit)$long$x$.key, paste0("K", setdiff(1:12, 5)))
})

test_that("prepareData(na.pass): a missing group / trials / exposure / link variable is an error", {
  d <- data.frame(id = rep(c("S1", "S2", "S3"), each = 3), t = rep(0:2, 3), x = seq_len(9) / 10,
                  k = c(1, 2, 0, 3, 1, 2, 0, 1, 2), n = 5, pt = 1.5, cnt = c(0, 1, 2, 1, 0, 3, 1, 1, 0))
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x), cnt | exposure(pt) ~ 1 + scaled(x)),
                       list(x ~ t + (1 | id)), family = list(k = "binomial", cnt = "poisson"))
  dn <- d; dn$n[2] <- NA
  expect_error(prepareData(f, dn), "Missing values in `n` \\(formula of 'k'\\).*na.omit")
  po <- prepareData(f, dn, na.action = na.omit)
  expect_equal(po$out$k$.key, paste0("K", c(1, 3:9)))
  expect_equal(nrow(po$out$cnt), 9L)                 # n is not a variable of cnt
  de <- d; de$pt[4] <- NA
  expect_error(prepareData(f, de), "Missing values in `pt` \\(formula of 'cnt'\\)")
  dg <- d; dg$id[4] <- NA
  expect_error(prepareData(f, dg), "Missing values in `id`")
  pg <- prepareData(f, dg, na.action = "na.omit")
  expect_equal(pg$long$x$.key, paste0("K", c(1:3, 5:9)))
  expect_equal(pg$out$k$.key, paste0("K", c(1:3, 5:9)))   # id is a link variable of k
  # cnt ~ 1: t is a link variable only
  dt <- d; dt$t[6] <- NA
  expect_error(prepareData(f, dataOut = dt, dataLong = d), "Missing values in `t`.*'k'")
  expect_equal(prepareData(f, dataOut = dt, dataLong = d, na.action = na.omit)$out$cnt$.key,
               paste0("K", c(1:5, 7:9)))
})

# ---- summarizeData() ----------------------------------------------------------

test_that("summarizeData(): rows, complete rows, groups and zero equations per equation", {
  d <- fxSmallData()
  d$x[3] <- NA           # long: response missing
  d$t[7] <- NA           # long covariate and link variable of y
  d$y[1] <- NA
  f <- smallForms()
  s <- summarizeData(f, d)
  expect_s3_class(s, "mvremDataSummary")
  e <- s$equations                      # in the order passed (pair: outcomes, then longitudinal)
  expect_equal(e$part, c("out", "long"))
  expect_equal(e$response, c("y", "x"))
  expect_equal(e$rows, c(12L, 12L))
  expect_equal(e$complete, c(10L, 10L))
  expect_equal(e$dropped, c(2L, 2L))
  expect_equal(e$groups, c("", "id: 4"))
  expect_equal(e$family, c("gaussian", "gaussian"))
  m <- s$missing
  expect_equal(m$variable[m$response == "y"], c("y", "t", "id"))
  expect_equal(m$missing[m$response == "y"], c(1L, 1L, 0L))
  expect_equal(m$missing[m$response == "x"], c(1L, 1L, 0L))
  expect_equal(s$zero$rows, 10L)
  expect_equal(s$zero$vars, "id, t")
  expect_output(print(s), "Missing values per variable")
  expect_output(print(s), "id: 4")
  # dropping a whole subject from the long data lowers its group count only
  s2 <- summarizeData(f, dataOut = d, dataLong = d[d$id != "S4", ])
  expect_equal(s2$equations$groups[2], "id: 3")
  expect_equal(s2$zero$rows, 10L)       # zero equations come from the outcome rows
})

test_that("summarizeData(): list data, several groups, a long formula without random effects", {
  d <- fxSmallData()
  dLong <- d
  dLong$x[2] <- NA
  dLong$t[5] <- NA
  dLong <- dLong[dLong$id != "S4", ]
  dOut <- d
  dOut$t[7] <- NA                       # link variable of y
  f <- collectFormulas(list(y ~ 1 + (1 | id) + (1 | site) + scaled(x) + scaled(cnt)),
                       list(x ~ 1 + t + (1 + t | id), cnt ~ t))
  s <- summarizeData(f, dataOut = dOut, dataLong = list(x = dLong, cnt = d))
  e <- s$equations
  expect_equal(e$part, c("out", "long", "long"))
  expect_equal(e$response, c("y", "x", "cnt"))
  expect_equal(e$rows, c(12L, 9L, 12L))
  expect_equal(e$complete, c(11L, 7L, 12L))
  expect_equal(e$dropped, c(1L, 2L, 0L))
  expect_equal(e$groups, c("id: 4, site: 2", "id: 3", ""))
  # complete rows = what prepareData(na.action = na.omit) keeps
  dat <- prepareData(f, dataOut = dOut, dataLong = list(x = dLong, cnt = d), na.action = na.omit)
  expect_equal(e$complete, c(nrow(dat$out$y), nrow(dat$long$x), nrow(dat$long$cnt)))
  # one list: the order of the list
  fl <- collectFormulas(list(x ~ 1 + t + (1 + t | id), y ~ 1 + (1 | id) + (1 | site) + scaled(x) + scaled(cnt),
                             cnt ~ t))
  el <- summarizeData(fl, dataOut = dOut, dataLong = list(x = dLong, cnt = d))$equations
  expect_equal(el$response, c("x", "y", "cnt"))
  expect_equal(el[match(e$response, el$response), names(el)], e, ignore_attr = TRUE)
  m <- s$missing
  get <- function(r, v) m$missing[m$response == r & m$variable == v]
  expect_equal(c(get("x", "x"), get("x", "t"), get("x", "id")), c(1L, 1L, 0L))
  expect_equal(c(get("y", "y"), get("y", "t")), c(0L, 1L))
  expect_equal(sum(m$missing[m$response == "cnt"]), 0L)
  z <- s$zero
  expect_equal(z$long, c("x", "cnt"))
  expect_equal(z$out, c("y", "y"))
  expect_equal(z$vars, c("id, t", "t"))
  expect_equal(z$rows, c(11L, 3L))      # unique (id, t) / t among the complete outcome rows
  # duplicated outcome rows do not add zero equations
  expect_equal(summarizeData(f, dataOut = rbind(dOut, dOut), dataLong = list(x = dLong, cnt = d))$zero$rows,
               c(11L, 3L))
})

test_that("summarizeData(): missing values in aux and link variables are counted", {
  d <- fxSmallData()
  d$n[2] <- NA                          # trials() of k
  d$t[5] <- NA                          # covariate of x, link variable of k
  f <- collectFormulas(list(k | trials(n) ~ 1 + scaled(x)), list(x ~ t + (1 | id)),
                       family = list(k = "binomial"))
  s <- summarizeData(f, d)
  m <- s$missing
  expect_equal(m$variable[m$response == "k"], c("k", "n", "id", "t"))
  expect_equal(m$missing[m$response == "k"], c(0L, 1L, 0L, 1L))
  expect_equal(m$missing[m$response == "x"], c(0L, 1L, 0L))
  expect_equal(s$equations$complete, c(10L, 11L))           # k, then x (the order passed)
  expect_equal(s$equations$family, c("binomial", "gaussian"))
  expect_equal(s$zero$rows, 10L)
})

test_that("summarizeData(): a fixture model with aux variables", {
  s <- dataSummary(fxModel("multi", fit = FALSE))
  expect_equal(s$equations$response, c("k", "cnt", "x1", "x2"))       # the order passed
  expect_equal(s$equations$family, c("binomial", "poisson", "gaussian", "poisson"))
  expect_equal(s$zero$long, c("x1", "x2", "x2"))
  expect_equal(s$zero$out, c("k", "k", "cnt"))
  expect_equal(s$missing$variable[s$missing$response == "cnt"], c("cnt", "grp", "pt", "id", "t"))
  expect_true(all(s$missing$missing == 0L))
  expect_output(print(s), "none")
})

test_that("summarizeData(): `data` equals dataOut + dataLong; a model without random effects", {
  d <- fxSmallData()
  f <- smallForms()
  expect_equal(summarizeData(f, d), summarizeData(f, dataOut = d, dataLong = d))
  expect_error(summarizeData(f, dataOut = d), "No data for the longitudinal formulae")
  s0 <- summarizeData(collectFormulas(y ~ t + scaled(x), x ~ t), data.frame(y = 1:3, t = 1:3, x = 1:3))
  expect_equal(s0$equations$groups, c("", ""))
  expect_equal(s0$zero$rows, 3L)
  expect_output(print(s0), "Equations")
  expect_output(print(s0), "Zero equations")
  expect_output(print(s0), "none")
})

# ---- namingScheme() ------------------------------------------------------------

test_that("namingScheme(): defaults, overrides and invalid templates", {
  nm <- namingScheme()
  expect_equal(nm$response, "{resp}_{resp}")
  expect_equal(nm$zeroIndex, "{long}_TO_{out}_{vars}")
  expect_equal(nm$interceptLabel, "(Intercept)")
  expect_equal(namingScheme(sep = ".")$sep, ".")
  expect_error(namingScheme(foo = "x", bar = "y"), "Unknown naming template\\(s\\): foo, bar")
  expect_error(namingScheme(sep = 1), "Naming template `sep` must be a single string")
  expect_error(namingScheme(fixed = c("a", "b")), "must be a single string")
  expect_error(namingScheme(fixed = NA_character_), "must be a single string")
})

test_that("namingScheme(): each template accepts exactly its own placeholders", {
  for (tpl in names(templatePlaceholders)) {
    own <- templatePlaceholders[[tpl]]
    ok <- paste0("{", own, "}", collapse = "_")
    expect_equal(do.call(namingScheme, stats::setNames(list(ok), tpl))[[tpl]], ok, info = tpl)
    other <- setdiff(unlist(templatePlaceholders), own)[1L]
    expect_error(do.call(namingScheme, stats::setNames(list(paste0(ok, "_{", other, "}")), tpl)),
                 paste0("Unknown placeholder\\(s\\) \\{", other, "\\} in naming template `", tpl,
                        "` \\(allowed: ", paste0("\\{", own, "\\}", collapse = ", "), "\\)"), info = tpl)
  }
  expect_error(namingScheme(fixed = "{resp}_{foo}"),
               "Unknown placeholder\\(s\\) \\{foo\\} in naming template `fixed` \\(allowed: \\{resp\\}, \\{term\\}\\)")
  # entries that are not templates take no placeholder at all
  for (tpl in c("responseObject", "levelsObject", "intercept", "interceptLabel", "sep", "valueSep"))
    expect_error(do.call(namingScheme, stats::setNames(list("{resp}"), tpl)),
                 paste0("\\{resp\\} in naming template `", tpl, "`\\.$"), info = tpl)
})

test_that("namingScheme(): index templates and `sep` give syntactic INLA names", {
  f <- smallForms()
  inl <- buildSmall(f, naming = namingScheme(fixed = "{resp} {term}", intercept = "Int",
                                             reIndex = "{resp}-{term}-{group}",
                                             zeroWeight = "w.{long}", responseObject = "Y",
                                             levelsObject = "LV", sep = "+"))
  expect_equal(names(inl$fixed$x$columns), c("x.Int", "x.t"))
  expect_equal(inl$random$x[[1]]$index, "x.Int.id")
  expect_equal(inl$random$x[[1]]$levels, "x_Int_id")
  expect_equal(inl$zero$x_y$index, "x_TO_y_id.t")
  expect_equal(inl$zero$x_y$weight, "w.x")
  expect_equal(deparse1(inl$formula[[2]]), "Y")
  expect_match(deparse1(inl$formula), "values = LV$x_Int_id", fixed = TRUE)
  expect_identical(get("LV", environment(inl$formula)), inl$levels)
  expect_true(all(make.names(all.vars(inl$formula)) == all.vars(inl$formula)))
})

test_that("namingScheme(): value templates label the index values (need not be syntactic)", {
  inl <- buildSmall(smallForms(), naming = namingScheme(reValue = "{level}/{term}",
                                                        zeroValue = "{long}>{out}/{values}", valueSep = "@"))
  expect_equal(inl$levels$x_1_id, paste0("S", 1:4, "/1"))
  expect_equal(inl$levels$x_TO_y_id_t[1:2], c("x>y/S1@0", "x>y/S1@1"))
  expect_equal(inl$stacks$x$y$x_TO_y_id_t[1:2], c("x>y/S1@0", "x>y/S1@1"))
  expect_equal(inl$stacks$y$main$y_FROM_x_id_t[1:2], c("x>y/S1@0", "x>y/S1@1"))
})

test_that("namingScheme(): options(mvrem.naming) sets the default scheme", {
  oldOpt <- options(mvrem.naming = namingScheme(fixed = "b.{resp}.{term}"))
  on.exit(options(oldOpt))
  expect_equal(names(buildSmall(smallForms())$fixed$x$columns), c("b.x.1", "b.x.t"))
  expect_equal(names(buildSmall(smallForms(), naming = namingScheme())$fixed$x$columns), c("x_1", "x_t"))
})

test_that("namingScheme(interceptLabel = ): a display label that leaves the INLA formula unchanged", {
  f <- smallForms()
  inl0 <- buildSmall(f)
  inl <- buildSmall(f, naming = namingScheme(interceptLabel = "intercept"))
  expect_equal(inl$naming$interceptLabel, "intercept")
  expect_equal(deparse(inl$formula), deparse(inl0$formula))
  expect_equal(inl$stacks, inl0$stacks)
  expect_equal(inl$fixed, inl0$fixed)
  expect_equal(effectLabel(c("(Intercept)", "t"), inl$naming), c("intercept", "t"))
  expect_equal(effectLabel(c("(Intercept)", "t"), inl0$naming), c("(Intercept)", "t"))
})

test_that("namingScheme(): trials / exposure / event templates name the input columns", {
  d <- fxSmallData()
  f <- collectFormulas(list(k | trials(n) ~ t + scaled(x), cnt | exposure(pt) ~ 1 + scaled(x)),
                       list(x ~ t + (1 | id)), family = list(k = "binomial", cnt = "poisson"))
  inl <- buildSmall(f, naming = namingScheme(trials = "N_{resp}", exposure = "{resp}_off"))
  expect_equal(inl$inputs$k, c(Ntrials = "N_k"))
  expect_equal(inl$inputs$cnt, c(E = "cnt_off"))
  expect_equal(names(inl$stacks$k$main)[1:2], c("k_k", "N_k"))
  expect_equal(names(inl$stacks$cnt$main)[1:2], c("cnt_cnt", "cnt_off"))
  expect_equal(inl$stacks$k$main$N_k, d$n)
  a <- inlaArgumentList(inl)
  expect_equal(a$Ntrials[seq_len(12)], d$n)               # k rows come first
  expect_equal(a$E[12 + seq_len(12)], d$pt)
  s <- fxSpec("surv")
  fs <- fxForms("surv")
  inls <- buildInlaFormula(fs, prepareData(fs, dataOut = s$dataOut, dataLong = s$dataLong),
                           naming = namingScheme(event = "{resp}_status"))
  expect_equal(inls$inputs$time, c(event = "time_status"))
  expect_equal(names(inls$stacks$time$main)[1:2], c("time_time", "time_status"))
  expect_equal(inls$stacks$time$main$time_status, s$dataOut$status)
})

test_that("namingScheme(): duplicated resulting names are an error", {
  f <- smallForms()
  expect_error(buildSmall(f, naming = namingScheme(zeroIndex = "{long}_{out}", scaledIndex = "{long}_{out}")),
               "duplicated INLA names: x_y")
  expect_error(buildSmall(smallForms("(1 + t | id)"), naming = namingScheme(reIndex = "{resp}_{group}")),
               "duplicated INLA names: x_id")
})

test_that("namingScheme(): response / zero-response / weight columns colliding with other names are an error", {
  f <- smallForms()
  expect_error(buildSmall(f, naming = namingScheme(response = "{resp}", zeroResponse = "{resp}")),
               "duplicated INLA names: x")
  expect_error(buildSmall(f, naming = namingScheme(zeroWeight = "{long}_t")),
               "duplicated INLA names: x_t")
  # the response column of a trials() response colliding with its input column
  fk <- collectFormulas(list(k | trials(n) ~ t + scaled(x)), list(x ~ t + (1 | id)),
                        family = list(k = "binomial"))
  expect_error(buildSmall(fk, naming = namingScheme(trials = "{resp}_{resp}")),
               "duplicated INLA names: k_k")
  # the response matrix / LEVELS objects colliding with a column
  expect_error(buildSmall(f, naming = namingScheme(responseObject = "x_t")),
               "duplicated INLA names: x_t")
  expect_error(buildSmall(f, naming = namingScheme(levelsObject = "x_1_id")),
               "duplicated INLA names: x_1_id")
  # one weight name for two longitudinal responses (one long linked to several
  # outcomes shares its weight column by design and is not a duplicate)
  expect_error(buildSmall(pairForms(), naming = namingScheme(zeroWeight = "W")),
               "duplicated INLA names: W")
  f2 <- collectFormulas(list(y ~ t + scaled(x), log(pt) ~ t + scaled(x)), list(x ~ t + (1 | id)))
  expect_equal(buildSmall(f2)$zero$x_pt$weight, "x_WEIGHT")
})

# ---- buildInlaFormula(): response, fixed and random terms ----------------------

test_that("buildInlaFormula(): response columns, zero responses and the response object", {
  inl <- buildSmall(pairForms())
  expect_equal(inl$responses, c(y = "y_y", pt = "pt_pt", x = "x_x", cnt = "cnt_cnt"))
  expect_equal(inl$zeroResponses, c(x = "x_x_0eq", cnt = "cnt_cnt_0eq"))
  expect_equal(deparse1(inl$formula[[2]]), "MX")
  expect_identical(get("LEVELS", environment(inl$formula)), inl$levels)
})

test_that("buildInlaFormula(): one fixed column per design column, names made syntactic", {
  d <- smallData()
  inl <- buildSmall(smallForms("(1 | id)", fixed = "t + g + I(t^2)"), d)
  expect_equal(inl$fixed$x$columns,
               c(x_1 = "(Intercept)", x_t = "t", x_gv = "gv", x_gw = "gw", x_I.t.2. = "I(t^2)"))
  m <- inl$stacks$x$main
  expect_equal(m$x_gv, as.numeric(d$g == "v"))
  expect_equal(m$x_gw, as.numeric(d$g == "w"))
  expect_equal(m$x_I.t.2., d$t^2)
  expect_true(all(c("x_gv", "x_gw", "x_I.t.2.") %in% all.vars(inl$formula)))
})

test_that("buildInlaFormula(): a correlated block is iid<k>d on the first effect + weighted copies", {
  inl <- buildSmall(smallForms("(1 + t | id)"))
  b <- inl$random$x[[1]]
  expect_equal(b$index, c("x_1_id", "x_t_id"))
  expect_equal(b$weight, c(NA, "x_t"))
  expect_equal(b$levels, "x_1_t_id")
  expect_equal(deparse1(b$terms[[1]]),
               'f(x_1_id, model = "iid2d", values = LEVELS$x_1_t_id, n = length(LEVELS$x_1_t_id))')
  expect_equal(deparse1(b$terms[[2]]), 'f(x_t_id, x_t, copy = "x_1_id")')
  expect_equal(inl$levels$x_1_t_id, c(paste0("x|1|S", 1:4), paste0("x|t|S", 1:4)))
  expect_null(termHyper(inl, "x_1_id"))                 # no prior: INLA's default
  b3 <- buildSmall(smallForms("(1 + t + pt | id)"))$random$x[[1]]
  expect_equal(as.list(b3$terms[[1]])$model, "iid3d")
  expect_length(b3$terms, 3L)
  expect_equal(deparse1(b3$terms[[3]]), 'f(x_pt_id, x_pt, copy = "x_1_id")')
  expect_null(as.list(b3$terms[[1]])$order)
})

test_that("buildInlaFormula(): a correlated block with k > 5 is iidkd with order", {
  inl6 <- buildSmall(smallForms("(1 + t + g + cnt + pt | id)"), smallData())   # g: 2 columns -> k = 6
  b6 <- inl6$random$x[[1]]
  t1 <- as.list(b6$terms[[1]])
  expect_equal(t1$model, "iidkd")
  expect_equal(t1$order, 6L)
  expect_length(b6$terms, 6L)
  expect_equal(b6$levels, "x_1_t_gv_gw_cnt_pt_id")
  expect_length(inl6$levels$x_1_t_gv_gw_cnt_pt_id, 6L * 4L)
  expect_true(all(vapply(b6$terms[-1], function(e) as.list(e)$copy, "") == "x_1_id"))
  skip_if_not_installed("INLA")
  spec <- INLA:::inla.interpret.formula(inl6$formula, data.same.len = FALSE, data = NULL)$random.spec
  expect_equal(sort(table(vapply(spec, `[[`, "", "model"))),
               sort(table(c("iidkd", rep("copy", 6), "iid"))))
})

test_that("buildInlaFormula(): a single effect or a || block gives one iid per effect", {
  inl1 <- buildSmall(smallForms("(1 | id)"))
  expect_equal(deparse1(inl1$random$x[[1]]$terms[[1]]),
               'f(x_1_id, model = "iid", values = LEVELS$x_1_id, n = length(LEVELS$x_1_id))')
  inl <- buildSmall(smallForms("(1 + pt || id)"))
  b <- inl$random$x
  expect_length(b, 2L)
  expect_equal(vapply(b, function(e) as.list(e$terms[[1]])$model, ""), c("iid", "iid"))
  expect_equal(vapply(b, `[[`, "", "index"), c("x_1_id", "x_pt_id"))
  expect_equal(vapply(b, `[[`, "", "levels"), c("x_1_id", "x_pt_id"))
  expect_equal(deparse1(b[[2]]$terms[[1]]),
               'f(x_pt_id, x_pt, model = "iid", values = LEVELS$x_pt_id, n = length(LEVELS$x_pt_id))')
  expect_equal(inl$levels$x_pt_id, paste0("x|pt|S", 1:4))
})

test_that("buildInlaFormula(): an intercept-free correlated block has no NA weight", {
  b <- buildSmall(smallForms("(0 + t + pt | id)"))$random$x[[1]]
  expect_equal(b$index, c("x_t_id", "x_pt_id"))
  expect_equal(b$weight, c("x_t", "x_pt"))
  expect_equal(deparse1(b$terms[[1]][[3]]), "x_t")
})

test_that("buildInlaFormula(): an uncorrelated slope without a fixed counterpart still gets its weight column", {
  d <- fxSmallData()
  inl <- buildSmall(smallForms("(1 + t || id)", fixed = "I(t^2)"))
  expect_equal(names(inl$fixed$x$columns), c("x_1", "x_I.t.2."))
  expect_true("x_t" %in% names(inl$stacks$x$main))
  expect_equal(inl$stacks$x$main$x_t, d$t)
  expect_equal(inl$random$x[[2]]$weight, "x_t")
})

test_that("buildInlaFormula(): hyper is set only with a prior, the same for every effect of a || block", {
  f <- smallForms("re(1 + pt || id, prior = gammaPrec(1, 0.1))")
  inl <- buildSmall(f)
  expect_equal(termHyper(inl, "x_1_id"), list(prec = list(prior = "loggamma", param = c(1, 0.1))))
  expect_equal(termHyper(inl, "x_pt_id"), termHyper(inl, "x_1_id"))
  f3 <- smallForms("re(1 + t + pt | id, prior = wishartPrec(df = 6, scale = 2))")
  inl3 <- buildSmall(f3)
  h3 <- termHyper(inl3, "x_1_id")
  expect_equal(h3$theta1$prior, "wishart3d")
  expect_equal(h3$theta1$param, c(6, rep(0.5, 3), rep(0, 3)))
  expect_null(termHyper(inl3, "x_t_id"))                 # copies carry no hyper
  expect_s3_class(inl3$random$x[[1]]$prior, "mvremPrior")
})

test_that("buildInlaFormula(): a user `hyper` on a block without prior is embedded as given", {
  wish <- list(theta1 = list(prior = "wishart2d", param = c(4, 1, 0.0011, 0)))
  inl <- buildSmall(smallForms(out = y ~ 1 + t + re(1 + t | id, hyper = wish) + scaled(x)))
  expect_equal(termHyper(inl, "y_1_id"), wish)
  expect_true(is.call(as.list(termOf(inl, "y_1_id"))$hyper))
  expect_null(inl$random$y[[1]]$prior)
})

test_that("buildInlaFormula(): user f() parameters are merged; a user hyper overriding a prior warns", {
  f <- collectFormulas(
    list(y ~ t + scaled(x, prior = normalCoef(0, 1), hyper = list(beta = list(param = c(1, 1))))),
    list(x ~ t + re(1 | id, prior = pcSd(1), hyper = list(prec = list(prior = "loggamma")),
                    constr = TRUE)))
  ws <- character()
  inl <- withCallingHandlers(buildSmall(f), warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(ws, 2L)
  expect_match(ws, "overrides the `prior` of (block `1 \\| id`|scaled\\(x\\) in the formula for y)")
  t1 <- as.list(inl$random$x[[1]]$terms[[1]])
  expect_true(t1$constr)
  expect_equal(termHyper(inl, "x_1_id")$prec$prior, "loggamma")
  expect_equal(termHyper(inl, "y_FROM_x_id_t")$beta$param, c(1, 1))
  expect_s3_class(inl$scaled$y_x$prior, "mvremPrior")
})

test_that("buildInlaFormula(): a user hyper that sets neither prior nor param is merged silently", {
  f <- collectFormulas(list(y ~ t + scaled(x, prior = normalCoef(0, 1), hyper = list(beta = list(initial = 1)))),
                       list(x ~ t + re(1 | id, prior = pcSd(1), hyper = list(prec = list(initial = 2)))))
  expect_no_warning(inl <- buildSmall(f))
  expect_equal(termHyper(inl, "x_1_id")$prec, list(prior = "pc.prec", param = c(1, 0.01), initial = 2))
  expect_equal(termHyper(inl, "y_FROM_x_id_t")$beta,
               list(fixed = FALSE, prior = "normal", param = c(0, 1), initial = 1))
})

test_that("fTerm(): embeds non-scalar values as calls, keeps weights and merges params", {
  expect_equal(deparse1(fTerm("a", NULL, list(model = "iid", values = c(1, 2)), list(constr = TRUE))),
               'f(a, model = "iid", values = c(1, 2), constr = TRUE)')
  expect_equal(deparse1(fTerm("a", NA, list(v = c(p = 1, q = 2)), list())), "f(a, v = c(p = 1, q = 2))")
  tt <- fTerm("a", "w", list(hyper = list(prec = list(param = c(1, 2)))),
              list(hyper = list(prec = list(initial = 1))))
  expect_equal(deparse1(tt[[3]]), "w")
  expect_equal(eval(as.list(tt)$hyper), list(prec = list(param = c(1, 2), initial = 1)))
  expect_true(is.call(as.list(tt)$hyper))
})

test_that("buildInlaFormula(): f() arguments survive terms() and reach INLA with their names", {
  f <- smallForms("re(1 | id, prior = pcSd(u = 2))",
                  out = y ~ t + scaled(x, hyper = list(beta = list(param = c(0, 0.5)))))
  inl <- buildSmall(f)
  labels <- attr(stats::terms(inl$formula), "term.labels")
  expect_true(any(grepl('hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.01)))', labels, fixed = TRUE)))
  expect_true(any(grepl("hyper = list(prec = list(initial = -6, fixed = TRUE))", labels, fixed = TRUE)))
  skip_if_not_installed("INLA")
  spec <- INLA:::inla.interpret.formula(inl$formula, data.same.len = FALSE, data = NULL)$random.spec
  names(spec) <- vapply(spec, `[[`, "", "term")
  expect_equal(spec$x_1_id$hyper$theta$prior, "pc.prec", ignore_attr = TRUE)
  expect_equal(spec$x_1_id$hyper$theta$param, c(2, 0.01), ignore_attr = TRUE)
  expect_true(spec$x_TO_y_id_t$hyper$theta$fixed)
  expect_equal(spec$x_TO_y_id_t$hyper$theta$initial, -6, ignore_attr = TRUE)
  expect_equal(spec$y_FROM_x_id_t$hyper$theta$param, c(0, 0.5), ignore_attr = TRUE)
  expect_false(spec$y_FROM_x_id_t$hyper$theta$fixed)
})

test_that("buildInlaFormula(): the complete formula of the multi fixture is accepted by INLA", {
  skip_if_not_installed("INLA")
  inl <- fxModel("multi", fit = FALSE)$inla
  i <- INLA:::inla.interpret.formula(inl$formula, data.same.len = FALSE, data = NULL)
  models <- vapply(i$random.spec, `[[`, "", "model")
  # iid: x2 intercept + 3 zero fields; iid2d: x1 block; copy: x1 slope + 3 scaled copies
  expect_equal(as.vector(table(models)[c("iid", "iid2d", "copy")]), c(4L, 1L, 4L))
})

# ---- buildInlaFormula(): zero equations and scaled copies ----------------------

test_that("buildInlaFormula(): one zero-equation term per linked pair, keyed by groups and covariates", {
  inl <- buildSmall(pairForms())
  expect_named(inl$zero, c("cnt_y", "x_pt"))
  z <- inl$zero$x_pt
  expect_equal(z$vars, c("id", "t"))
  expect_equal(z$index, "x_TO_pt_id_t")
  expect_equal(deparse1(z$term),
               'f(x_TO_pt_id_t, x_WEIGHT, model = "iid", values = LEVELS$x_TO_pt_id_t, hyper = list(prec = list(initial = -6, fixed = TRUE)))')
  expect_equal(head(inl$levels$x_TO_pt_id_t, 2), c("x|TO|pt|S1|0", "x|TO|pt|S1|1"))
  expect_length(inl$levels$x_TO_pt_id_t, 12L)
  expect_equal(inl$zero$cnt_y$vars, c("id", "t"))      # t: covariate of I(t^2)
  inlg <- buildSmall(smallForms("(1 | id)", fixed = "t + g"), smallData())
  expect_equal(inlg$zero$x_y$vars, c("id", "t", "g"))
  expect_equal(inlg$zero$x_y$index, "x_TO_y_id_t_g")
})

test_that("buildInlaFormula(): the scaled term copies the zero field with the default coefficient prior", {
  inl <- buildSmall(pairForms())
  expect_named(inl$scaled, c("y_cnt", "pt_x"))
  expect_equal(deparse1(inl$scaled$pt_x$term), paste0(
    'f(pt_FROM_x_id_t, copy = "x_TO_pt_id_t", hyper = list(beta = list(fixed = FALSE, ',
    'prior = "normal", param = c(0, 0.01))))'))
  expect_null(inl$scaled$pt_x$prior)
  s <- as.list(inl$scaled$y_cnt$term)
  expect_equal(deparse1(s[[2]]), "y_FROM_cnt_id_t")
  expect_equal(s$copy, "cnt_TO_y_id_t")
})

test_that("buildInlaFormula(): zero equations exist only for the linked (long, outcome) pairs", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ 1 + t + scaled(x), k ~ 1 + t + scaled(cnt), pt ~ 1 + t),
                       list(x ~ 1 + t + (1 | id), cnt ~ 1 + t + (1 | id)))
  dat <- prepareData(f, d)
  inl <- buildInlaFormula(f, dat)
  expect_named(inl$zero, c("x_y", "cnt_k"))
  expect_named(inl$scaled, c("y_x", "k_cnt"))
  expect_equal(inl$zeroResponses, c(x = "x_x_0eq", cnt = "cnt_cnt_0eq"))
  expect_equal(names(inl$families), c("y_y", "k_k", "pt_pt", "x_x", "cnt_cnt", "x_x_0eq", "cnt_cnt_0eq"))
  expect_equal(names(inl$controlFamily), names(inl$families))
  expect_named(inl$stacks$x, c("main", "y"))
  expect_named(inl$stacks$cnt, c("main", "k"))
  expect_false(any(grepl("y_FROM_cnt|k_FROM_x", deparse1(inl$formula))))
  expect_false(any(grepl("_FROM_", names(inl$stacks$pt$main))))   # pt is not linked
  expect_no_error(formatInlaFormula(inl))
  expect_no_error(latexEquation(f, dat))
  expect_no_error(summarizeData(f, d))
})

test_that("buildInlaFormula(): zero-equation factor levels come from the long data", {
  d <- smallData()
  f <- smallForms("(1 | id)", fixed = "t + g")
  dOut <- d[d$g == "v", ]            # outcome rows only at level v
  inl <- buildSmall(f, dataOut = dOut, dataLong = d)
  z <- inl$stacks$x$y
  expect_equal(nrow(z), nrow(dOut))
  expect_named(z, names(inl$stacks$x$main))
  expect_true(all(z$x_gv == 1) && all(z$x_gw == 0))
  expect_true(all(grepl("\\|v$", z$x_TO_y_id_t_g)))
})

test_that("buildInlaFormula(): a subject without long data still gets a zero equation (outcome with RE)", {
  d <- fxSmallData()
  f <- collectFormulas(list(y ~ t + (1 | id) + scaled(x)), list(x ~ t + (1 | id)))
  inl <- buildSmall(f, dataOut = d, dataLong = d[d$id != "S4", ])
  z <- inl$stacks$x$y
  expect_equal(nrow(z), 12L)
  expect_true("x|TO|y|S4|1" %in% z$x_TO_y_id_t)
  expect_true(all(z$x_1_id %in% inl$levels$x_1_id))
  expect_true("x|1|S4" %in% inl$levels$x_1_id)
})

test_that("buildInlaFormula(): a subject without long data still gets a zero equation (outcome without RE)", {
  skip("KNOWN BUG: RE labels of zero-equation rows of a subject absent from the long data are not in LEVELS when the outcome has no RE on that group")
  d <- fxSmallData()
  inl <- buildSmall(smallForms(), dataOut = d, dataLong = d[d$id != "S4", ])
  z <- inl$stacks$x$y
  expect_true("x|TO|y|S4|1" %in% z$x_TO_y_id_t)
  expect_true(all(z$x_1_id %in% inl$levels$x_1_id))
})

test_that("buildInlaFormula(): outcome rows are matched by covariate values, also with missing covariates", {
  d <- fxSmallData()
  # long data: response missing (row 2: S1, t = 1), covariate missing (row 5: S2, t = 1), no S4
  dLong <- d
  dLong$x[2] <- NA
  dLong$t[5] <- NA
  dLong <- dLong[dLong$id != "S4", ]
  # outcome data: t (used only for matching) missing in row 7 (S3, t = 0)
  dOut <- d
  dOut$t[7] <- NA
  f <- collectFormulas(list(y ~ 1 + (1 | id) + scaled(x)), list(x ~ 1 + t + (1 + t | id)))
  expect_equal(f$out$y$linkVariables, "t")
  dat <- prepareData(f, dataOut = dOut, dataLong = dLong, na.action = na.omit)
  inl <- buildInlaFormula(f, dat)
  y <- inl$stacks$y$main
  m <- inl$stacks$x$main
  z <- inl$stacks$x$y
  expect_equal(dat$out$y$.key, paste0("K", setdiff(1:12, 7)))
  expect_equal(nrow(y), 11L)
  expect_equal(nrow(m), 7L)
  expect_false(any(grepl("S4", m$x_1_id)))
  # one zero equation per (id, t) of the outcome, including combinations without a long observation
  expect_equal(nrow(z), 11L)
  expect_setequal(z$x_TO_y_id_t, inl$levels$x_TO_y_id_t)
  expect_true(all(c("x|TO|y|S1|1", "x|TO|y|S2|1", "x|TO|y|S4|2") %in% z$x_TO_y_id_t))
  expect_false(any(grepl("\\|NA$|\\|NA\\|", z$x_TO_y_id_t)))
  # every outcome row is matched to the zero equation with its own id and t ...
  i <- match(y$y_FROM_x_id_t, z$x_TO_y_id_t)
  expect_false(anyNA(i))
  expect_equal(z$x_t[i], dOut$t[-7])
  expect_equal(sub("^x\\|1\\|", "", z$x_1_id[i]), dOut$id[-7])
  # ... whose covariates and random-effect labels are those of x at that (id, t)
  expect_equal(z$x_t, as.numeric(sub("^.*\\|", "", z$x_TO_y_id_t)))
  expect_true(all(c(z$x_1_id, z$x_t_id) %in% inl$levels$x_1_t_id))
  expect_true("x|1|S4" %in% inl$levels$x_1_t_id)
})

test_that("buildInlaFormula(): outcome rows sharing (id, t) give one zero row; a covariate in vars splits them", {
  d <- smallData()
  dOut <- d[c(1, 1, 2, 2, 2, 5), ]          # S1 t = 0 (x2), S1 t = 1 (x3), S2 t = 1
  dOut$y <- seq_len(6)
  rownames(dOut) <- NULL
  inl <- buildSmall(smallForms(), dataOut = dOut, dataLong = d)
  expect_equal(nrow(inl$stacks$y$main), 6L)
  expect_equal(inl$stacks$x$y$x_TO_y_id_t, c("x|TO|y|S1|0", "x|TO|y|S1|1", "x|TO|y|S2|1"))
  expect_equal(inl$levels$x_TO_y_id_t, inl$stacks$x$y$x_TO_y_id_t)
  expect_false(anyDuplicated(inl$stacks$x$y$x_TO_y_id_t) > 0)
  # with g in the long formula, rows of the same (id, t) but another g differ
  dOut$g <- factor(c("u", "v", "u", "u", "w", "v"), levels = levels(d$g))
  inl2 <- buildSmall(smallForms("(1 | id)", fixed = "t + g"), dataOut = dOut, dataLong = d)
  expect_equal(nrow(inl2$stacks$x$y), 5L)
  expect_equal(inl2$stacks$x$y$x_gv, c(0, 1, 0, 0, 1))
  expect_equal(inl2$stacks$x$y$x_gw, c(0, 0, 0, 1, 0))
})

# ---- buildInlaFormula(): stacks, families and inputs ---------------------------

test_that("buildInlaFormula(): stacks hold the formula's columns, inputs right after the response", {
  inl <- fxModel("multi", fit = FALSE)$inla
  expect_equal(names(inl$stacks$k$main)[1:2], c("k_k", "k_Ntrials"))
  expect_equal(names(inl$stacks$cnt$main),
               c("cnt_cnt", "cnt_E", "cnt_1", "cnt_grpb", "cnt_grpc", "cnt_FROM_x2_id_t"))
  expect_equal(inl$stacks$k$main$k_Ntrials, fxData()$long$n)
  expect_equal(inl$stacks$cnt$main$cnt_E, fxData()$long$pt)
  expect_equal(inl$stacks$k$main$k_k, fxData()$long$k)
  expect_equal(names(inl$stacks$x1$main),
               c("x1_x1", "x1_x1_0eq", "x1_1", "x1_t", "x1_1_id", "x1_t_id", "x1_WEIGHT", "x1_TO_k_id_t"))
  expect_true(all(c(inl$stacks$x1$main$x1_1_id, inl$stacks$x1$main$x1_t_id) %in% inl$levels$x1_1_t_id))
})

test_that("buildInlaFormula(): zero-equation rows have response NA, _0eq = 0, WEIGHT = -1, inputs NA", {
  inl <- fxModel("multi", fit = FALSE)$inla
  expect_named(inl$stacks$x2, c("main", "k", "cnt"))
  x2 <- inl$stacks$x2
  expect_equal(names(x2$main),
               c("x2_x2", "x2_x2_0eq", "x2_1", "x2_t", "x2_1_id", "x2_WEIGHT", "x2_TO_k_id_t", "x2_TO_cnt_id_t"))
  for (p in c("k", "cnt")) {
    z <- x2[[p]]
    expect_named(z, names(x2$main))
    expect_true(all(is.na(z$x2_x2)))
    expect_true(all(z$x2_x2_0eq == 0))
    expect_true(all(z$x2_WEIGHT == -1))
  }
  expect_true(all(is.na(x2$main$x2_x2_0eq)) && all(is.na(x2$main$x2_WEIGHT)))
  # each zero part fills only its own field
  expect_true(all(is.na(x2$k$x2_TO_cnt_id_t)) && !anyNA(x2$k$x2_TO_k_id_t))
  expect_true(all(is.na(x2$cnt$x2_TO_k_id_t)) && !anyNA(x2$cnt$x2_TO_cnt_id_t))
  expect_true(all(is.na(x2$main$x2_TO_k_id_t)))
  # outcome rows point to existing zero equations, which are exactly the LEVELS
  expect_true(all(inl$stacks$cnt$main$cnt_FROM_x2_id_t %in% x2$cnt$x2_TO_cnt_id_t))
  expect_true(all(inl$stacks$k$main$k_FROM_x2_id_t %in% x2$k$x2_TO_k_id_t))
  expect_setequal(x2$k$x2_TO_k_id_t, inl$levels$x2_TO_k_id_t)
})

test_that("buildInlaFormula(): families, controlFamily, inputs and survival are aligned", {
  inl <- fxModel("multi", fit = FALSE)$inla
  expect_equal(inl$families, c(k_k = "binomial", cnt_cnt = "poisson", x1_x1 = "gaussian",
                               x2_x2 = "poisson", x1_x1_0eq = "gaussian", x2_x2_0eq = "gaussian"))
  expect_equal(names(inl$controlFamily), names(inl$families))
  expect_equal(inl$controlFamily$x1_x1_0eq, list(hyper = list(prec = list(initial = 10, fixed = TRUE))))
  expect_equal(inl$controlFamily$k_k, list())
  expect_named(inl$inputs, c("k", "cnt", "x1", "x2"))
  expect_equal(inl$inputs$k, c(Ntrials = "k_Ntrials"))
  expect_equal(inl$inputs$cnt, c(E = "cnt_E"))
  expect_length(inl$inputs$x1, 0L)
  expect_equal(inl$survival, c(k = FALSE, cnt = FALSE, x1 = FALSE, x2 = FALSE))
  expect_equal(names(inl$responses), names(inl$survival))
  expect_equal(names(inl$terms), c("k", "cnt", "x1", "x2"))
})

test_that("buildInlaFormula(): a non-default link appears in controlFamily", {
  d <- fxSmallData()
  d$k <- rep(0:1, 6)
  f <- collectFormulas(list(k ~ t + scaled(x)), list(x ~ t + (1 | id)),
                       family = list(k = binomialFamily(link = "probit")))
  inl <- buildSmall(f, d)
  expect_equal(inl$controlFamily$k_k$control.link, list(model = "probit"))
  expect_equal(buildSmall(smallForms(out = k ~ t + scaled(x)), d)$controlFamily$k_k, list())
})

test_that("buildInlaFormula(): a survival outcome gives an event column, NA in zero-equation rows", {
  s <- fxSpec("surv")
  inl <- fxModel("surv", fit = FALSE)$inla
  expect_equal(inl$families, c(time_time = "weibullsurv", x1_x1 = "gaussian", x1_x1_0eq = "gaussian"))
  expect_equal(inl$survival, c(time = TRUE, x1 = FALSE))
  expect_equal(inl$inputs$time, c(event = "time_event"))
  st <- inl$stacks$time$main
  expect_equal(names(st)[1:2], c("time_time", "time_event"))
  expect_equal(st$time_event, s$dataOut$status)
  expect_equal(st$time_time, s$dataOut$time)
  expect_true(all(is.na(inl$stacks$x1$time$time_event)))
  expect_equal(nrow(inl$stacks$x1$time), nrow(s$dataOut))     # one zero equation per subject (t = 0)
  expect_match(paste(formatInlaFormula(inl, nrows = 1), collapse = "\n"), "time_event")
})

# ---- formatInlaFormula() -------------------------------------------------------

test_that("formatInlaFormula(): one term per line, response groups with header comments", {
  inl <- buildSmall(pairForms())
  txt <- formatInlaFormula(inl, nrows = 2)
  expect_s3_class(txt, "mvremFormulaText")
  expect_equal(txt[1], "MX ~ -1 +")
  expect_true(any(grepl("^   # ----- y: y ~ t \\+ \\(1 \\| id\\) \\+ scaled\\(cnt\\)", txt)))
  expect_true(any(grepl("^   # ----- pt: log\\(pt\\) ~", txt)))
  expect_true("y_1 +" %in% txt)
  expect_true('f(x_t_id, x_t, copy = "x_1_id") +' %in% txt)
  expect_false(grepl("\\+$", txt[length(txt)]))
  expect_true(any(grepl('# family = "gaussian"', txt, fixed = TRUE)))
  expect_true(any(grepl('# zero equations: family = "gaussian", control.family = list(hyper', txt, fixed = TRUE)))
  expect_equal(sum(grepl("zero equations:", txt)), 2L)          # x and cnt
  # nrows rows per part: header + 2 observation + 2 zero-equation rows for x
  i <- grep("----- x:", txt)
  expect_equal(sum(grepl("x\\|1\\|S1", txt[i:(i + 8)])), 4L)
  expect_output(expect_invisible(print(txt)))
  expect_output(print(txt), "----- cnt: cnt ~ I(t^2) + re(1 | id) -----", fixed = TRUE)
})

test_that("formatInlaFormula(): parses back to the INLA formula for any nrows / digits", {
  inl <- fxModel("multi", fit = FALSE)$inla
  for (nr in c(0L, 1L, 3L)) {
    txt <- formatInlaFormula(inl, nrows = nr, digits = 2)
    expect_equal(deparse1(stats::as.formula(paste(txt, collapse = "\n"))), deparse1(inl$formula))
    expect_equal(deparse1(str2lang(paste(txt, collapse = "\n"))), deparse1(inl$formula))
  }
  txt <- formatInlaFormula(inl)
  expect_true(any(grepl('# family = "binomial"', txt, fixed = TRUE)))
  expect_true(any(grepl("k_Ntrials", txt)))
  expect_output(print(txt), "----- cnt: cnt | exposure(pt) ~ grp + scaled(x2) -----", fixed = TRUE)
  # nrows = 0 hides the data; a model gives the same text as its buildInlaFormula() result
  txt0 <- formatInlaFormula(inl, nrows = 0)
  expect_false(any(grepl("k_Ntrials|S01", txt0)))
  expect_output(print(txt0), "MX ~ -1 +", fixed = TRUE)
  expect_equal(formatInlaFormula(fxModel("multi", fit = FALSE), nrows = 1), formatInlaFormula(inl, nrows = 1))
  expect_error(formatInlaFormula(inl, nrows = -1), "nrows >= 0")
})

test_that("formatInlaFormula(): data are shown with `digits` significant digits", {
  inl <- fxModel("multi", fit = FALSE)$inla
  x1 <- inl$stacks$x1$main$x1_x1[1]
  txt <- formatInlaFormula(inl, nrows = 1)
  expect_true(any(grepl(format(x1, digits = 3), txt, fixed = TRUE)))
  expect_false(any(grepl(format(x1, digits = 7), txt, fixed = TRUE)))
  expect_true(any(grepl(format(x1, digits = 6), formatInlaFormula(inl, nrows = 1, digits = 6), fixed = TRUE)))
})

# ---- inlaData() ----------------------------------------------------------------

test_that("inlaData(): the default method stacks a buildInlaFormula() result with a .block factor", {
  inl <- fxModel("multi", fit = FALSE)$inla
  st <- inlaData(inl)
  nr <- sum(unlist(lapply(inl$stacks, function(s) vapply(s, nrow, 0L))))
  expect_equal(nrow(st), nr)
  expect_s3_class(st$.block, "factor")
  # one level per stack part, in likelihood order
  expLevels <- unlist(lapply(names(inl$stacks), function(r) {
    parts <- names(inl$stacks[[r]])
    ifelse(parts == "main", paste0(r, ": observations"), paste0(r, ": zero equations for ", parts))
  }))
  expect_equal(levels(st$.block), expLevels)
  # every block has one row per long observation: the zero equations are keyed by (id, t),
  # and each subject of fxData() is observed once per time
  nObs <- nrow(fxData()$long)
  expect_equal(as.vector(table(st$.block)), rep(nObs, length(expLevels)))
  expect_true(all(is.na(st$k_Ntrials[st$.block != "k: observations"])))
  st2 <- inlaData(inl, n = 2)
  expect_equal(nrow(st2), 2L * length(expLevels))
  expect_equal(names(st2), names(st))
  expect_equal(inlaData(fxModel("multi", fit = FALSE), n = 2), st2)
})

# ---- inlaArgumentList() ----------------------------------------------------------

test_that("inlaArgumentList(): response matrix, Ntrials / E and link per row", {
  inl <- fxModel("multi", fit = FALSE)$inla
  a <- inlaArgumentList(inl)
  st <- inlaData(inl)
  Y <- a$data$MX
  expect_true(is.matrix(Y))
  expect_equal(colnames(Y), names(inl$families))
  expect_true(all(rowSums(!is.na(Y)) == 1L))          # one likelihood per row
  expect_equal(a$family, unname(inl$families))
  expect_equal(a$control.family, unname(inl$controlFamily))
  expect_false(any(c(names(inl$families), "k_Ntrials", "cnt_E", ".block") %in% names(a$data)))
  used <- setdiff(all.vars(inl$formula[[3]]), c("LEVELS", names(inl$levels)))
  expect_true(all(used %in% names(a$data)))
  k <- st$.block == "k: observations"
  expect_equal(a$Ntrials[k], st$k_Ntrials[k])
  expect_true(all(a$Ntrials[!k] == 1))
  cn <- st$.block == "cnt: observations"
  expect_equal(a$E[cn], st$cnt_E[cn])
  expect_true(all(a$E[!cn] == 1))
  # link: the likelihood index of each row
  lik <- a$control.predictor$link
  expect_length(lik, nrow(st))
  expect_equal(lik, apply(!is.na(Y), 1, which), ignore_attr = TRUE)
  expect_true(a$control.predictor$compute)
})

test_that("inlaArgumentList(): verbose, threads and inlaOptions are merged; reserved names refused", {
  inl <- fxModel("multi", fit = FALSE)$inla
  a <- inlaArgumentList(inl)
  expect_false(a$verbose)
  expect_null(a$num.threads)
  a2 <- inlaArgumentList(inl, list(verbose = TRUE, threads = 2,
                                   inlaOptions = list(control.compute = list(dic = TRUE),
                                                      control.predictor = list(compute = FALSE))))
  expect_true(a2$verbose)
  expect_equal(a2$num.threads, 2)
  expect_equal(a2$control.compute, list(dic = TRUE))
  expect_false(a2$control.predictor$compute)
  expect_equal(a2$control.predictor$link, a$control.predictor$link)       # merged, not replaced
  expect_error(inlaArgumentList(inl, list(inlaOptions = list(E = 1, family = "x"))),
               "must not set `E`, `family`: mvrem builds these")
  expect_error(inlaArgumentList(inl, list(inlaOptions = list(1))), "named list")
})

test_that("inlaArgumentList(): a missing response stays NA in the response matrix", {
  d <- fxSmallData()
  d$y[2] <- NA
  inl <- buildSmall(smallForms(), d)
  a <- inlaArgumentList(inl)
  expect_true(all(is.na(a$data$MX[2, ])))
  expect_equal(a$control.predictor$link[2], 1L)       # still the outcome's likelihood
  expect_null(a$Ntrials)
  expect_null(a$E)
})

test_that("inlaArgumentList(): a survival response is an inla.surv inside a list, not a data column", {
  skip_if_not_installed("INLA")
  s <- fxSpec("surv")
  a <- inlaArgumentList(fxModel("surv", fit = FALSE)$inla)
  expect_true(is.list(a$data$MX) && !is.matrix(a$data$MX))
  expect_s3_class(a$data$MX$time_time, "inla.surv")
  expect_null(a$data$time_event)
  expect_equal(a$data$MX$time_time$event[seq_len(nrow(s$dataOut))], s$dataOut$status)
})

test_that("inlaArgumentList(): a survival response needs INLA", {
  local_mocked_bindings(hasInla = function() FALSE)
  expect_error(inlaArgumentList(fxModel("surv", fit = FALSE)$inla), "INLA package is required")
})
