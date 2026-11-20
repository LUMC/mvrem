# Extended tests of the plots: plot.mvremModel() / plotLinesData(), plotModel(), the result
# plots (plotRandomEffects(), plotVarianceComponents(), plotFixedEffects(), plotAssociations(),
# densityPanels(), pairEllipses()), plotResiduals() / plotResidualsData(), the widget specs
# (htmlPlot.R), plot.mvremPredictions() and the widget's JavaScript (run in node).
# Fixtures: helper-fixtures.R; helpers: helper-plots.R.

# ---- useInteractive() / options --------------------------------------------------------------

test_that("useInteractive(): NULL follows htmlwidgets and options(mvrem.interactive)", {
  withr::local_options(mvrem.interactive = NULL)
  expect_identical(useInteractive(NULL), hasWidgets())
  expect_false(useInteractive(FALSE))
  withr::with_options(list(mvrem.interactive = FALSE), expect_false(useInteractive(NULL)))
  withr::with_options(list(mvrem.interactive = "yes"), expect_false(useInteractive(NULL)))  # only TRUE counts
})

test_that("useInteractive(): an explicit TRUE wins over options(mvrem.interactive = FALSE)", {
  skipNoWidgets()
  withr::local_options(mvrem.interactive = FALSE)
  expect_true(useInteractive(TRUE))
})

test_that("useInteractive(): without htmlwidgets NULL / FALSE are base graphics and TRUE an error", {
  withr::local_options(mvrem.interactive = NULL)
  local_mocked_bindings(hasHtmlwidgets = function() FALSE)
  expect_false(useInteractive(NULL))
  expect_false(useInteractive(FALSE))
  expect_error(useInteractive(TRUE), "needs the htmlwidgets package")
})

test_that("plot functions: interactive = TRUE without htmlwidgets is an error, the default is base graphics", {
  m <- fxModel("multi")
  local_mocked_bindings(hasHtmlwidgets = function() FALSE)
  expect_error(plot(m, interactive = TRUE), "htmlwidgets")
  expect_error(plotResiduals(m, interactive = TRUE), "htmlwidgets")
  expect_error(plotRandomEffects(m, interactive = TRUE), "htmlwidgets")
  expect_error(plotVarianceComponents(m, interactive = TRUE), "htmlwidgets")
  expect_error(plotFixedEffects(m, interactive = TRUE), "htmlwidgets")
  expect_error(plotAssociations(m, interactive = TRUE), "htmlwidgets")
  expect_error(htmlPlot(list()), "htmlwidgets")
  withNullDevice({
    expect_identical(plot(m), m)
    expect_identical(plotRandomEffects(m), m)
    expect_identical(plotFixedEffects(m), m)
    expect_identical(plotResiduals(m), m)
  })
})

test_that("plot functions: options(mvrem.interactive = FALSE) draws every plot with base graphics", {
  m <- fxModel("multi")
  withr::local_options(mvrem.interactive = FALSE)
  withNullDevice({
    expect_invisible(plot(m))
    expect_identical(plot(m), m)
    expect_identical(plotResiduals(m), m)
    expect_identical(plotRandomEffects(m), m)
    expect_identical(plotVarianceComponents(m), m)
    expect_identical(plotFixedEffects(m), m)
    expect_identical(plotAssociations(m), m)
  })
})

test_that("plot(): a widget by default when htmlwidgets is installed, base graphics on request", {
  skipNoWidgets()
  withr::local_options(mvrem.interactive = NULL)
  m0 <- fxModel("beta", fit = FALSE)
  expect_s3_class(plot(m0), "mvremPlot")
  expect_s3_class(plotRandomEffects(fxModel("beta")), "mvremPlot")
  withNullDevice(expect_identical(plot(m0, interactive = FALSE), m0))
})

test_that("plot(): n = NULL defaults to 30 subjects for base graphics and 300 for the widget", {
  m <- fxModel("beta", fit = FALSE)
  orig <- plotLinesData
  seen <- NULL
  local_mocked_bindings(plotLinesData = function(x, response, time, group, n, ...) {
    seen <<- n
    orig(x, response, time, group, n, ...)
  })
  withNullDevice(plot(m, interactive = FALSE))
  expect_equal(seen, 30L)
  withNullDevice(plot(m, n = 7, interactive = FALSE))
  expect_equal(seen, 7)
  skipNoWidgets()
  plot(m, interactive = TRUE)
  expect_equal(seen, 300L)
})

# ---- plotLinesData() / plot.mvremModel() -----------------------------------------------------

test_that("plotLinesData(): panels of an unfitted model, in the order passed, finite ordered limits", {
  m <- fxModel("multi", fit = FALSE)
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE)
  expect_false(d$haveFit)
  expect_equal(d$time, "t")
  expect_equal(d$group, "id")
  expect_equal(vapply(d$panels, `[[`, "", "response"), c("k", "cnt", "x1", "x2"))
  titles <- vapply(d$panels, `[[`, "", "title")
  expect_match(titles[1], "^k .*outcome.*binomial")
  expect_match(titles[2], "^cnt .*outcome.*poisson")
  expect_match(titles[3], "^x1 .*longitudinal.*gaussian")
  expect_match(titles[4], "^x2 .*longitudinal.*poisson")
  for (pn in d$panels) {
    expect_false(pn$missing)
    expect_true(all(is.finite(pn$xlim)) && pn$xlim[1] <= pn$xlim[2])
    expect_true(all(is.finite(pn$ylim)) && pn$ylim[1] <= pn$ylim[2])
    expect_false(pn$useFit)
    expect_null(pn$smoothFit)
    expect_length(pn$subjects, 30)
    expect_true(all(vapply(pn$subjects, function(s) is.null(s$fit) && !is.unsorted(s$x), NA)))
    # every drawn point lies inside the limits
    xs <- unlist(lapply(pn$subjects, `[[`, "x")); ys <- unlist(lapply(pn$subjects, `[[`, "y"))
    expect_true(all(xs >= pn$xlim[1] & xs <= pn$xlim[2] & ys >= pn$ylim[1] & ys <= pn$ylim[2]))
    expect_named(pn$smoothObs, c("x", "y"))
  }
  # the same subjects in every panel
  expect_true(all(vapply(d$panels, function(p) identical(subjectIds(p), subjectIds(d$panels[[1]])), NA)))
})

test_that("plotLinesData(): trials() responses as proportions, exposure() responses as rates", {
  d <- plotLinesData(fxModel("multi", fit = FALSE), NULL, NULL, NULL, 30, NULL, TRUE)
  expect_equal(d$panels[[1]]$ylab, "k / n")
  expect_equal(d$panels[[2]]$ylab, "cnt / pt")
  long <- fxData()$long
  s1 <- d$panels[[1]]$subjects[[1]]
  rows <- which(long$id == s1$id)
  expect_equal(s1$y, (long$k / long$n)[rows][order(long$t[rows])])
  expect_true(all(d$panels[[1]]$ylim >= 0 & d$panels[[1]]$ylim <= 1))
  s2 <- d$panels[[2]]$subjects[[1]]
  rows <- which(long$id == s2$id)
  expect_equal(s2$y, (long$cnt / long$pt)[rows][order(long$t[rows])])
})

test_that("plotLinesData(): a response with trials() and exposure() names both denominators (Q27)", {
  skip_if_not_installed("INLA")
  long <- fxData()$long
  m <- jsm(list(k | trials(n) + exposure(pt) ~ 1 + t + scaled(x1), x1 ~ 1 + t + (1 | id)),
           data = long, family = list(k = inlaFamily("binomial")), fit = FALSE)
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE)
  expect_equal(d$panels[[1]]$ylab, "k / n / pt")
  s1 <- d$panels[[1]]$subjects[[1]]
  rows <- which(long$id == s1$id)
  expect_equal(s1$y, (long$k / long$n / long$pt)[rows][order(long$t[rows])])
})

test_that("plotLinesData(): n, ids and group select the subjects", {
  m <- fxModel("beta", fit = FALSE)
  ids <- allIds()
  one <- plotLinesData(m, NULL, NULL, NULL, 1, NULL, TRUE)
  expect_equal(subjectIds(one$panels[[1]]), ids[1])
  two <- plotLinesData(m, NULL, NULL, NULL, 2, NULL, TRUE)
  expect_equal(subjectIds(two$panels[[1]]), ids[c(1, length(ids))])   # evenly spaced: first and last
  big <- plotLinesData(m, NULL, NULL, NULL, 1000, NULL, TRUE)
  expect_equal(subjectIds(big$panels[[1]]), ids)                        # all of them, no duplicates
  for (n in c(3, 7, 30)) {
    p <- plotLinesData(m, NULL, NULL, NULL, n, NULL, TRUE)$panels[[1]]
    expect_length(p$subjects, n)
    expect_false(anyDuplicated(subjectIds(p)) > 0)
  }
  # ids: exactly those subjects (in the order given), overriding n
  sel <- plotLinesData(m, NULL, NULL, NULL, 1, c("S05", "S03", "S20"), TRUE)
  for (pn in sel$panels) expect_equal(subjectIds(pn), c("S05", "S03", "S20"))
  # ids not in the data are dropped silently; none known -> no subjects, limits from all rows
  mixed <- plotLinesData(m, NULL, NULL, NULL, 30, c("S02", "nobody"), TRUE)
  expect_equal(subjectIds(mixed$panels[[1]]), "S02")
  none <- plotLinesData(m, NULL, NULL, NULL, 30, "nobody", TRUE)
  expect_length(none$panels[[1]]$subjects, 0)
  expect_equal(none$panels[[1]]$xlim, range(fxData()$long$t))
})

test_that("plotLinesData(): group = NA pools all rows into one unconnected pseudo-subject", {
  na <- plotLinesData(fxModel("beta", fit = FALSE), NULL, NULL, NA, 30, NULL, TRUE)
  expect_true(is.na(na$group))
  expect_false(na$panels[[1]]$grouped)
  expect_equal(subjectIds(na$panels[[1]]), "all")
  expect_length(na$panels[[1]]$subjects[[1]]$x, nrow(fxData()$long))
})

test_that("plotLinesData(): response filter and explicit time / group, with clear errors", {
  m <- fxModel("multi", fit = FALSE)
  d <- plotLinesData(m, c("x2", "k"), NULL, NULL, 5, NULL, TRUE)
  expect_equal(vapply(d$panels, `[[`, "", "response"), c("k", "x2"))   # in the order of the model
  expect_error(plotLinesData(m, c("k", "zzz"), NULL, NULL, 5, NULL, TRUE), "Unknown response.*zzz")
  # explicit time / group equal to the defaults change nothing
  e <- plotLinesData(m, NULL, "t", "id", 5, NULL, TRUE)
  expect_equal(e$time, "t")
  expect_equal(e$group, "id")
  expect_equal(e$panels, plotLinesData(m, NULL, NULL, NULL, 5, NULL, TRUE)$panels)
  expect_error(plotLinesData(m, NULL, "nope", NULL, 5, NULL, TRUE), "`time`.*'nope'.*not a variable")
  expect_error(plotLinesData(m, NULL, c("t", "n"), NULL, 5, NULL, TRUE), "one variable")
  expect_error(plotLinesData(m, NULL, NULL, "nope", 5, NULL, TRUE), "`group`.*'nope'.*not a variable")
  expect_error(plotLinesData(m, NULL, NULL, c("id", "grp"), 5, NULL, TRUE), "one variable")
})

test_that("plotLinesData(): a time variable present in only some equations gives empty panels", {
  m <- fxModel("multi", fit = FALSE)
  dp <- plotLinesData(m, NULL, "pt", NULL, 5, NULL, TRUE)
  miss <- vapply(dp$panels, function(p) isTRUE(p$missing), NA)
  expect_equal(miss, c(TRUE, FALSE, TRUE, TRUE))          # pt: exposure of cnt only
  expect_match(dp$panels[[1]]$title, "^k:.*`pt`")
  expect_match(dp$panels[[3]]$title, "^x1:.*`pt`")
  expect_equal(dp$panels[[2]]$xlab, "pt")
})

test_that("plotLinesData(): another grouping variable (site) gives one line per level where it is used", {
  m <- fxModel("priors", fit = FALSE)
  d <- plotLinesData(m, NULL, "t", "site", 30, NULL, TRUE)
  expect_equal(d$group, "site")
  sites <- levels(fxData()$long$site)
  yg <- d$panels[[1]]                                  # the outcome has the random effect on site
  expect_equal(yg$response, "yg")
  expect_equal(subjectIds(yg), sites)
  expect_true(all(vapply(yg$subjects, function(s) length(s$x) == sum(fxData()$long$site == s$id), NA)))
  # the longitudinal equations do not use site: their prepared data have no such column, so
  # no subject lines (the smooth mean of all rows is still there)
  for (pn in d$panels[-1]) {
    expect_length(pn$subjects, 0)
    expect_true(pn$grouped)
    expect_named(pn$smoothObs, c("x", "y"))
  }
})

test_that("plotTimeVariable() / plotGroupVariable(): the default choices and their fallbacks", {
  m <- fxModel("priors", fit = FALSE)
  eqs <- c(m$forms$out, m$forms$long)
  prep <- c(m$data$out, m$data$long)
  expect_equal(plotTimeVariable(eqs, prep, NULL), "t")      # the numeric covariate of every equation
  expect_equal(plotGroupVariable(eqs, prep, NULL), "id")    # four blocks on id, one on site
  expect_equal(plotGroupVariable(eqs, prep, "site"), "site")
  expect_equal(plotGroupVariable(eqs, prep, NA), NA_character_)
  # without random effects: the grouping variables of the formulae, else NA
  noRe <- lapply(eqs, function(p) { p$re <- list(); p })
  expect_equal(plotGroupVariable(noRe, prep, NULL), "id")
  noGroup <- lapply(noRe, function(p) { p$variables$group <- character(); p })
  expect_equal(plotGroupVariable(noGroup, prep, NULL), NA_character_)
  # no numeric covariate at all: an error asking for `time`
  noNum <- lapply(eqs, function(p) { p$variables$covariate <- character(); p$linkVariables <- NULL; p })
  expect_error(plotTimeVariable(noNum, prep, NULL), "No numeric covariate")
  # the covariate used by most equations; ties go to the first candidate
  mm <- fxModel("multi", fit = FALSE)
  eqs2 <- c(mm$forms$out, mm$forms$long)
  prep2 <- c(mm$data$out, mm$data$long)
  eqs2[["cnt"]]$variables$covariate <- c("pt", "t")    # pt is a numeric column of cnt's data
  expect_equal(plotTimeVariable(eqs2["cnt"], prep2, NULL), "pt")
  expect_equal(plotTimeVariable(eqs2[c("cnt", "k")], prep2, NULL), "t")
})

test_that("plotLinesData(): a non-numeric time variable gives a clear error", {
  # B10: all plot functions resolve `time` through plotTimeVariable(), which must reject
  # a factor / character variable itself instead of failing later inside range()
  m <- fxModel("beta", fit = FALSE)
  expect_error(plotLinesData(m, NULL, "id", NULL, 5, NULL, TRUE), "`time` = 'id' must be a numeric variable")
  expect_error(plot(m, time = "id", interactive = FALSE), "`time` = 'id' must be a numeric variable")
  expect_error(plot(m, time = "id", interactive = FALSE), "factor")
  expect_error(plot(m, time = "nope", interactive = FALSE), "not a variable of the model's data")
  mf <- fxModel("multi")
  expect_error(plotResiduals(mf, time = "id", interactive = FALSE), "`time` = 'id' must be a numeric variable")
  expect_error(plotResiduals(mf, type = "time", time = "id", interactive = FALSE), "numeric")
  pr <- fxPredict()
  expect_error(plot(pr, time = "id", interactive = FALSE), "`time` = 'id' must be a numeric variable")
})

test_that("plotLinesData(): fitted values of a fitted model, or none with fitted = FALSE", {
  m <- fxModel("multi")
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE)
  expect_true(d$haveFit)
  for (pn in d$panels) {
    expect_true(pn$useFit)
    expect_true(all(vapply(pn$subjects, function(s) length(s$fit) == length(s$x), NA)))
    fits <- unlist(lapply(pn$subjects, `[[`, "fit"))
    expect_true(all(fits >= pn$ylim[1] & fits <= pn$ylim[2]))
    expect_named(pn$smoothFit, c("x", "y"))
  }
  # binomial / poisson fitted means on the proportion / rate scale, as fittedValues() gives them
  fitsK <- unlist(lapply(d$panels[[1]]$subjects, `[[`, "fit"))
  expect_true(all(fitsK > 0 & fitsK < 1))
  expect_true(all(unlist(lapply(d$panels[[2]]$subjects, `[[`, "fit")) > 0))
  fv <- fittedValues(m, "x1")
  s <- d$panels[[3]]$subjects[[1]]
  rows <- which(fv$id == s$id)
  expect_equal(s$fit, fv$mean[rows][order(fv$t[rows])])
  nf <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, FALSE)
  expect_false(nf$haveFit)
  expect_false(any(vapply(nf$panels, `[[`, NA, "useFit")))
})

test_that("plot() of a survival model: the survival outcome is an empty panel", {
  m <- fxModel("surv")
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE)
  expect_equal(vapply(d$panels, `[[`, "", "response"), c("time", "x1"))
  expect_true(d$panels[[1]]$missing)
  expect_match(d$panels[[1]]$title, "^time:.*survival.*not plotted")
  expect_false(d$panels[[2]]$missing)
  expect_true(d$panels[[2]]$useFit)
  only <- plotLinesData(m, "time", NULL, NULL, 30, NULL, TRUE)
  expect_length(only$panels, 1)
  expect_true(only$panels[[1]]$missing)
  withNullDevice({
    expect_identical(plot(m, interactive = FALSE), m)
    expect_identical(plot(m, "time", interactive = FALSE), m)    # only the empty panel, no legend
    expect_identical(plot(fxModel("surv", fit = FALSE), interactive = FALSE), fxModel("surv", fit = FALSE))
  })
  skipNoWidgets()
  w <- expectCleanWidget(plot(m, interactive = TRUE))
  expect_true(w$x$panels[[1]]$empty)
  expect_equal(w$x$panels[[1]]$title, d$panels[[1]]$title)
  expect_true(length(w$x$panels[[2]]$lines) > 0)
})

test_that("plot(): base graphics on unfitted and fitted models, all argument variants", {
  m0 <- fxModel("multi", fit = FALSE)
  m <- fxModel("multi")
  withNullDevice({
    expect_identical(plot(m0, interactive = FALSE), m0)
    expect_identical(plot(m, interactive = FALSE), m)
    expect_identical(plot(m, fitted = FALSE, interactive = FALSE), m)
    expect_identical(plot(m, time = "t", group = "id", interactive = FALSE), m)
    expect_identical(plot(fxModel("priors"), group = "site", interactive = FALSE), fxModel("priors"))
    expect_no_error(plot(m, n = 1, interactive = FALSE))
    expect_no_error(plot(m, n = 1000, interactive = FALSE))
    expect_no_error(plot(m, ids = "nobody", interactive = FALSE))
    expect_no_error(plot(m, ids = c("S01", "S40"), ncol = 4, interactive = FALSE))
    expect_no_error(plot(m, response = "cnt", group = NA, interactive = FALSE))
    expect_no_error(plot(m, time = "pt", interactive = FALSE))   # empty panels for k, x1, x2
    expect_error(plot(m, time = "nope", interactive = FALSE), "not a variable")
    expect_error(plot(m, group = "nope", interactive = FALSE), "not a variable")
    expect_error(plot(m, response = "nope", interactive = FALSE), "Unknown response")
  })
})

test_that("plot() widget: panels, legend, linking and search of a fitted model", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plot(m, interactive = TRUE))
  x <- w$x
  expect_true(x$linked)
  expect_equal(x$search, "id id")              # "<group> id"
  expect_equal(x$ncol, 2L)                     # 4 panels -> 2 columns
  expect_match(panelTitles(w)[1], "^k .*outcome")
  expect_match(panelTitles(w)[2], "^cnt .*outcome")
  expect_equal(legendKinds(w), c("obs", "obsMean", "fit", "fitMean"))
  for (pn in x$panels) {
    expect_setequal(lineIds(pn, "obs"), allIds())          # n = 300 by default: all subjects
    expect_setequal(lineIds(pn, "fit"), allIds())
    expect_true(pn$connect)
    expect_equal(vapply(pn$curves, `[[`, "", "kind"), c("obsMean", "fitMean"))
    expect_true(pn$xlim[1] <= pn$xlim[2] && pn$ylim[1] <= pn$ylim[2])
  }
})

test_that("plot() widget: fitted = FALSE and an unfitted model have no fitted lines or legend entries", {
  skipNoWidgets()
  nf <- plot(fxModel("multi"), fitted = FALSE, interactive = TRUE)
  expect_equal(legendKinds(nf), c("obs", "obsMean"))
  expect_false(any(unlist(lapply(nf$x$panels, lineKinds)) == "fit"))
  expect_true(all(vapply(nf$x$panels, function(p) identical(vapply(p$curves, `[[`, "", "kind"), "obsMean"), NA)))
  w0 <- expectCleanWidget(plot(fxModel("multi", fit = FALSE), interactive = TRUE))
  expect_equal(legendKinds(w0), c("obs", "obsMean"))
  expect_false(any(unlist(lapply(w0$x$panels, lineKinds)) == "fit"))
})

test_that("plot() widget: n, ids, response and ncol select the subjects and panels", {
  skipNoWidgets()
  m <- fxModel("multi")
  expect_length(lineIds(plot(m, n = 1, interactive = TRUE)$x$panels[[1]]), 1)
  sel <- plot(m, ids = c("S07", "S11"), interactive = TRUE)$x
  for (pn in sel$panels) expect_equal(lineIds(pn), c("S07", "S11"))
  none <- expectCleanWidget(plot(m, ids = "nobody", interactive = TRUE))
  expect_length(none$x$panels[[1]]$lines, 0)
  one <- plot(m, response = "x1", interactive = TRUE)$x
  expect_length(one$panels, 1)
  expect_equal(one$ncol, 1L)
  expect_equal(plot(m, ncol = 3, interactive = TRUE)$x$ncol, 3L)
})

test_that("plot() widget: explicit time / group, group = NA and a time variable missing in equations", {
  skipNoWidgets()
  m <- fxModel("multi")
  e <- plot(m, time = "t", group = "id", interactive = TRUE)
  expect_equal(e$x, plot(m, interactive = TRUE)$x)
  ws <- expectCleanWidget(plot(fxModel("priors"), group = "site", interactive = TRUE))
  expect_equal(ws$x$search, "site id")
  expect_equal(lineIds(ws$x$panels[[1]], "obs"), levels(fxData()$long$site))
  ng <- expectCleanWidget(plot(m, group = NA, interactive = TRUE))
  expect_equal(ng$x$search, "subject id")
  expect_false(ng$x$panels[[1]]$connect)
  expect_equal(lineIds(ng$x$panels[[1]]), "all")
  wp <- expectCleanWidget(plot(m, time = "pt", interactive = TRUE))
  expect_true(isTRUE(wp$x$panels[[1]]$empty) && isTRUE(wp$x$panels[[3]]$empty))
  expect_null(wp$x$panels[[2]]$empty)
})

test_that("plot() widget: saved self-contained, without external resources", {
  skipNoWidgets()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  w <- plot(fxModel("beta", fit = FALSE), interactive = TRUE)
  f <- withr::local_tempfile(fileext = ".html")
  htmlwidgets::saveWidget(w, f, selfcontained = TRUE)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "mvremPlot")
  expect_false(grepl("(src|href)=\"https?://", html))
})

# ---- htmlLines / htmlIntervals / htmlPoints --------------------------------------------------

test_that("htmlIntervals() and htmlPoints() build interval and scatter panels", {
  p <- htmlIntervals("T", c("a", "b"), c("A", "B"), c(1, 2), c(0, 1), c(2, 3))
  expect_equal(p$title, "T")
  expect_true(p$yReverse)
  expect_equal(p$vline, 0)
  expect_equal(p$intervals$pos, 1:2)
  expect_equal(p$intervals$tip, c("A\n1 [0, 2]", "B\n2 [1, 3]"))
  expect_null(p$ylabels)
  expect_equal(htmlIntervals("T", "a", "A", 1, 0, 2, ylabels = TRUE)$ylabels, list(pos = 1L, text = "A"))
  q <- htmlPoints("P", c("s1", "s2"), c(1.234567, -1), c(0.5, 2), "u", "v")
  expect_equal(q$points$x, signif(c(1.234567, -1), 5))      # r5(): 5 significant digits
  expect_equal(q$points$tip[1], "s1\nu = 1.235, v = 0.5")
  expect_equal(c(q$vline, q$hline), c(0, 0))
})

test_that("htmlLines(): missing panels, hline, connect and renamed kinds", {
  skipNoWidgets()
  d <- list(haveFit = FALSE, group = "id", panels = list(
    list(missing = TRUE, title = "gone"),
    list(missing = FALSE, title = "p", xlab = "t", ylab = "y", xlim = c(0, 1), ylim = c(-1, 1),
         grouped = TRUE, connect = FALSE, hline = 0,
         subjects = list(list(id = "a", x = c(0, 1), y = c(-1, 1))), smoothObs = list(x = 0:1, y = 0:1))))
  w <- expectCleanWidget(htmlLines(d, 1L, kinds = c(obs = "res", obsMean = "resMean"),
                                   labels = c(obs = "R", obsMean = "RM")))
  expect_equal(w$x$panels[[1]], list(title = "gone", empty = TRUE))
  p <- w$x$panels[[2]]
  expect_false(p$connect)                       # explicit connect wins over grouped
  expect_equal(p$hline, 0)
  expect_equal(p$lines[[1]]$kind, "res")
  expect_equal(p$curves[[1]]$kind, "resMean")
  expect_equal(vapply(w$x$legend, `[[`, "", "label"), c("R", "RM"))
  expect_equal(w$x$search, "id id")
})

# ---- plotModel() -----------------------------------------------------------------------------

test_that("plotModel(): runs on all fixtures, details and fits, returning the text size invisibly", {
  withNullDevice({
    for (nm in c("beta", "multi", "surv", "priors")) {
      m0 <- fxModel(nm, fit = FALSE)
      m <- fxModel(nm)
      expect_invisible(plotModel(m0))
      for (v in list(plotModel(m0, detail = "inla"), plotModel(m), plotModel(m, detail = "inla"),
                     plotModel(m, fitted = FALSE))) {
        expect_type(v, "double")
        expect_length(v, 1)
      }
    }
    expect_error(plotModel(fxModel("beta", fit = FALSE), detail = "other"), "should be one of")
  })
})

test_that("plotModel(): the text size is the one given when the diagram fits, shrunk otherwise", {
  m <- fxModel("multi", fit = FALSE)
  default <- eval(formals(plotModel.mvremModel)$cex)
  big <- plotModelCex(m, 7, 7)
  expect_equal(big, default)                          # fits on a 7 x 7 inch device: as given
  expect_equal(plotModelCex(m, 7, 7, detail = "inla"), default)
  small <- plotModelCex(m, 2, 2)
  expect_lt(small, big)                               # shrunk on a 2 x 2 inch device
  expect_gt(small, 0)
  large <- plotModelCex(m, 7, 7, cex = 2)
  expect_lt(large, 2)                                 # a large cex shrinks as needed ...
  expect_gte(large, big)                              # ... but not below what fits
  tiny <- plotModelCex(m, 1, 1, detail = "inla")       # shrinking stops at a floor: no endless loop
  expect_lte(tiny, small)
  expect_gt(tiny, 0)
  expect_lt(plotModelCex(fxModel("multi"), 2, 2), big)
})

# ---- result plots ----------------------------------------------------------------------------

test_that("result plots: an unfitted model is a clear error", {
  m0 <- fxModel("multi", fit = FALSE)
  withNullDevice({
    expect_error(plotRandomEffects(m0, interactive = FALSE), "not fitted")
    expect_error(plotRandomEffects(m0, type = "pairs", interactive = FALSE), "not fitted")
    expect_error(plotVarianceComponents(m0, interactive = FALSE), "not fitted")
    expect_error(plotFixedEffects(m0, interactive = FALSE), "not fitted")
    expect_error(plotFixedEffects(m0, type = "density", interactive = FALSE), "not fitted")
    expect_error(plotAssociations(m0, interactive = FALSE), "not fitted")
    expect_error(plotAssociations(m0, type = "density", interactive = FALSE), "not fitted")
    expect_error(plotResiduals(m0, interactive = FALSE), "not fitted")
  })
})

test_that("plotRandomEffects(): one panel per response, term and group, linked by subject id", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plotRandomEffects(m, interactive = TRUE))
  expect_equal(panelTitles(w), c("x1: (Intercept) | id", "x1: t | id", "x2: (Intercept) | id"))
  expect_true(w$x$linked)
  expect_equal(w$x$search, "subject id")
  expect_equal(w$x$ncol, 2L)
  re <- randomEffects(m)
  for (pn in w$x$panels) {
    iv <- pn$intervals
    expect_setequal(iv$id, allIds())
    expect_true(all(iv$lo <= iv$mean & iv$mean <= iv$hi))
    r <- sub(":.*", "", pn$title)
    term <- sub("^[^:]+: (.*) \\| id$", "\\1", pn$title)
    d <- re[re$response == r & re$term == term, ]
    expect_equal(iv$mean, signif(d$mean[match(iv$id, d$level)], 5))
  }
  expect_equal(plotRandomEffects(m, ncol = 1, interactive = TRUE)$x$ncol, 1L)
})

test_that("plotRandomEffects(): term selects by formula term, by label or by the alias 1", {
  skipNoWidgets()
  m <- fxModel("multi")
  titles <- function(...) panelTitles(plotRandomEffects(m, ..., interactive = TRUE))
  expect_equal(titles(term = "1"), c("x1: (Intercept) | id", "x2: (Intercept) | id"))
  expect_equal(titles(term = "(Intercept)"), c("x1: (Intercept) | id", "x2: (Intercept) | id"))
  expect_equal(titles("x1", term = "t"), "x1: t | id")
  expect_equal(titles(term = c("t", "1")), c("x1: (Intercept) | id", "x1: t | id", "x2: (Intercept) | id"))
  expect_error(plotRandomEffects(m, term = "zzz", interactive = TRUE), "No random effects")
})

test_that("plotRandomEffects(): response filter and errors", {
  skipNoWidgets()
  m <- fxModel("multi")
  expect_equal(panelTitles(plotRandomEffects(m, "x2", interactive = TRUE)), "x2: (Intercept) | id")
  expect_error(plotRandomEffects(m, "k", interactive = TRUE), "No random effects")
  expect_error(plotRandomEffects(m, "zzz", interactive = TRUE), "Unknown response")
})

test_that("plotRandomEffects(): order = \"mean\" sorts the subjects by posterior mean", {
  skipNoWidgets()
  m <- fxModel("multi")
  srt <- plotRandomEffects(m, "x2", order = "mean", interactive = TRUE)$x$panels[[1]]$intervals
  expect_false(is.unsorted(srt$mean))
  re <- randomEffects(m, "x2")
  expect_equal(srt$id, re$level[order(re$mean)])
  unsorted <- plotRandomEffects(m, "x2", interactive = TRUE)$x$panels[[1]]$intervals
  expect_equal(unsorted$id, re$level)                 # data order by default
})

test_that("plotRandomEffects(): base graphics variants", {
  m <- fxModel("multi")
  withNullDevice({
    expect_identical(plotRandomEffects(m, interactive = FALSE), m)
    expect_identical(plotRandomEffects(m, "x1", term = "t", order = "mean", ncol = 1, interactive = FALSE), m)
    expect_identical(plotRandomEffects(m, type = "pairs", interactive = FALSE), m)
    expect_identical(plotRandomEffects(fxModel("priors"), type = "pairs", interactive = FALSE), fxModel("priors"))
  })
})

test_that("plotRandomEffects(): the outcome's random intercept on site has one interval per level", {
  skipNoWidgets()
  m <- fxModel("priors")
  w <- expectCleanWidget(plotRandomEffects(m, "yg", interactive = TRUE))
  expect_equal(panelTitles(w), "yg: (Intercept) | site")
  expect_equal(w$x$panels[[1]]$intervals$id, levels(fxData()$long$site))
  expect_equal(w$x$ncol, 1L)
  all <- plotRandomEffects(m, interactive = TRUE)
  expect_equal(panelTitles(all)[1], "yg: (Intercept) | site")    # outcome first, then x1, x2, x3
  expect_length(all$x$panels, 7)
  expect_error(plotRandomEffects(m, "yg", type = "pairs", interactive = TRUE), "at least two effects")
})

test_that("plotRandomEffects(): labels use naming$interceptLabel", {
  skipNoWidgets()
  m <- fxModel("beta")                              # the cached fixture is never modified
  m$inla$naming$interceptLabel <- "b0"
  expect_equal(panelTitles(plotRandomEffects(m, interactive = TRUE)), c("x1: b0 | id", "x1: t | id"))
  expect_equal(panelTitles(plotRandomEffects(m, term = "b0", interactive = TRUE)), "x1: b0 | id")
  expect_equal(panelTitles(plotRandomEffects(m, term = "1", interactive = TRUE)), "x1: b0 | id")
  expect_equal(panelTitles(plotRandomEffects(m, term = "(Intercept)", interactive = TRUE)), "x1: b0 | id")
  wp <- plotRandomEffects(m, type = "pairs", interactive = TRUE)
  expect_equal(panelTitles(wp), "x1: b0 + t | id")
  expect_equal(wp$x$panels[[1]]$xlab, "posterior mean of b0 | id")
  # the ellipse is found under the new label too: the same points, curves and note as by default
  default <- plotRandomEffects(fxModel("beta"), type = "pairs", interactive = TRUE)$x$panels[[1]]
  geometry <- function(p) { p$points$tip <- NULL; p[c("points", "curves", "note")] }  # tips carry the label
  expect_equal(geometry(wp$x$panels[[1]]), geometry(default))
  expect_match(wp$x$panels[[1]]$points$tip[1], "posterior mean of b0 | id = ", fixed = TRUE)
})

test_that("plotRandomEffects(type = \"pairs\"): a point per subject and the ellipses of a correlated block", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plotRandomEffects(m, type = "pairs", interactive = TRUE))
  expect_length(w$x$panels, 1)                         # x2 has a single effect: no pair
  pn <- w$x$panels[[1]]
  expect_equal(pn$title, "x1: (Intercept) + t | id")
  expect_equal(pn$xlab, "posterior mean of (Intercept) | id")
  expect_equal(pn$ylab, "posterior mean of t | id")
  expect_setequal(pn$points$id, allIds())
  re <- randomEffects(m, "x1")
  expect_equal(pn$points$x, signif(re$mean[re$term == "(Intercept)"][match(pn$points$id, re$level[re$term == "(Intercept)"])], 5))
  expect_equal(vapply(pn$curves, `[[`, "", "kind"), c("ellipse", "ellipseBound", "ellipseBound"))
  expect_match(pn$note, "^correlation -?[0-9.]+ \\[-?[0-9.]+, -?[0-9.]+\\]$")
  expect_equal(legendKinds(w), c("ellipse", "ellipseBound"))
  expect_true(w$x$linked)
  expect_error(plotRandomEffects(m, "x2", type = "pairs", interactive = TRUE), "at least two effects")
  expect_error(plotRandomEffects(m, "x1", term = "t", type = "pairs", interactive = TRUE), "at least two effects")
  expectCleanWidget(plotRandomEffects(fxModel("surv"), type = "pairs", interactive = TRUE))
})

test_that("plotRandomEffects(type = \"pairs\"): independent blocks and || blocks have no ellipse", {
  skipNoWidgets()
  m <- fxModel("priors")
  w <- expectCleanWidget(plotRandomEffects(m, type = "pairs", interactive = TRUE))
  expect_equal(panelTitles(w), paste0(c("x1", "x2", "x3"), ": (Intercept) + t | id"))
  expect_length(w$x$panels[[1]]$curves, 3)             # x1: correlated block
  for (i in 2:3) {                                     # x2: two blocks on id; x3: `||`
    expect_null(w$x$panels[[i]]$curves)
    expect_null(w$x$panels[[i]]$note)
    expect_setequal(w$x$panels[[i]]$points$id, allIds())
  }
  expect_equal(legendKinds(w), c("ellipse", "ellipseBound"))
  for (r in c("x2", "x3")) {
    one <- plotRandomEffects(m, r, type = "pairs", interactive = TRUE)
    expect_null(one$x$legend)                          # no ellipse anywhere: no legend
    expect_null(one$x$panels[[1]]$curves)
  }
})

test_that("pairEllipses(): symmetric in the two effects, NULL without a table or a correlation", {
  m <- fxModel("multi")
  vc <- varianceComponents(m)$table
  e <- pairEllipses(vc, "x1", "id", "(Intercept)", "t")
  expect_named(e, c("ellipses", "note"))
  expect_named(e$ellipses, c("mean", "lo", "hi"))
  nPoints <- eval(formals(ellipsePoints)$n)
  expect_length(e$ellipses$mean$x, nPoints)
  # the ellipse's extent along x is the 95% chi-square radius times the SD of the first effect
  sd1 <- vc$mean[vc$response == "x1" & vc$parameter == "sd" & vc$effect == "(Intercept)" & vc$group %in% "id"]
  expect_equal(max(abs(e$ellipses$mean$x)), sqrt(qchisq(0.95, 2)) * sd1, tolerance = 1e-3)
  cr <- vc[vc$response == "x1" & vc$parameter == "correlation", ]
  expect_equal(e$note, sprintf("correlation %s [%s, %s]", signif(cr$mean, 2), signif(cr$q0.025, 2), signif(cr$q0.975, 2)))
  # the effects swapped: the same note, the ellipse mirrored (x extent <-> y extent)
  er <- pairEllipses(vc, "x1", "id", "t", "(Intercept)")
  expect_equal(er$note, e$note)
  expect_equal(max(abs(er$ellipses$mean$x)), max(abs(e$ellipses$mean$y)), tolerance = 1e-3)
  expect_equal(max(abs(er$ellipses$mean$y)), max(abs(e$ellipses$mean$x)), tolerance = 1e-3)
  expect_null(pairEllipses(NULL, "x1", "id", "(Intercept)", "t"))
  expect_null(pairEllipses(vc, "x2", "id", "(Intercept)", "t"))
  expect_null(pairEllipses(vc, "x1", "site", "(Intercept)", "t"))
  expect_null(pairEllipses(vc, "x1", "id", "(Intercept)", "zzz"))
  # a renamed intercept is found under its label
  nm <- namingScheme(interceptLabel = "b0")
  expect_equal(pairEllipses(vc, "x1", "id", "b0", "t", nm), e)
  expect_null(pairEllipses(vc, "x1", "id", "(Intercept)", "t", nm))
})

test_that("ellipsePoints(): the 95% ellipse of a bivariate normal", {
  e <- ellipsePoints(2, 1, 0)
  r <- sqrt(qchisq(0.95, 2))
  expect_equal(max(e$x), 2 * r, tolerance = 1e-3)
  expect_equal(max(e$y), r, tolerance = 1e-3)
  expect_equal(e$x[1], e$x[length(e$x)])              # closed curve
  # with correlation 1 the ellipse collapses onto the diagonal
  d <- ellipsePoints(1, 1, 1)
  expect_equal(d$x, d$y)
})

test_that("plotVarianceComponents(): block, correlation and residual panels", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plotVarianceComponents(m, interactive = TRUE))
  titles <- panelTitles(w)
  expect_true(all(c("correlations", "residual") %in% titles))
  expect_false(w$x$linked)
  expect_null(w$x$search)
  for (pn in w$x$panels) {
    expect_equal(pn$ylab, "posterior density")
    expect_equal(pn$xlab, if (pn$title == "correlations") "correlation" else "sd")
    expect_true(all(vapply(pn$curves, function(cv) all(cv$y >= 0) && length(cv$x) == length(cv$y), NA)))
  }
  expect_equal(vapply(w$x$panels[[which(titles == "residual")]]$curves, `[[`, "", "label"), "x1")   # gaussian only
  wv <- expectCleanWidget(plotVarianceComponents(m, scale = "variance", interactive = TRUE))
  expect_equal(wv$x$panels[[1]]$xlab, "variance")
  wr <- expectCleanWidget(plotVarianceComponents(m, response = "x2", interactive = TRUE))
  expect_false("correlations" %in% panelTitles(wr))
  expect_error(plotVarianceComponents(m, response = "zzz", interactive = TRUE), "Unknown response")
  expectCleanWidget(plotVarianceComponents(fxModel("surv"), interactive = TRUE))
  withNullDevice({
    expect_identical(plotVarianceComponents(m, interactive = FALSE), m)
    expect_identical(plotVarianceComponents(m, scale = "variance", ncol = 1, interactive = FALSE), m)
  })
})

test_that("plotVarianceComponents(): every block of the four responses of the priors model", {
  skipNoWidgets()
  m <- fxModel("priors")
  w <- expectCleanWidget(plotVarianceComponents(m, interactive = TRUE))
  titles <- panelTitles(w)
  expect_equal(titles, c("yg | site [(Intercept)]", "x1 | id [(Intercept), t]", "correlations",
                         "x2 | id [(Intercept)]", "x2 | id [t]", "x3 | id [(Intercept)]", "x3 | id [t]", "residual"))
  labels <- lapply(w$x$panels, function(p) vapply(p$curves, `[[`, "", "label"))
  expect_equal(labels[[2]], c("(Intercept)", "t"))
  expect_equal(labels[[3]], "x1 | id [(Intercept), t]: (Intercept), t")   # the only correlation
  expect_equal(labels[[8]], c("yg", "x1", "x3"))                            # x2 (negbin) has no residual
  expect_equal(w$x$ncol, 3L)                                                # 8 panels -> 3 columns
  vc <- varianceComponents(m)
  expect_equal(w$x$panels[[1]]$curves[[1]]$x, signif(vc$marginals[[1]][, 1], 5))
  two <- plotVarianceComponents(m, response = "x2", interactive = TRUE)
  expect_equal(panelTitles(two), c("x2 | id [(Intercept)]", "x2 | id [t]"))
  withNullDevice(expect_identical(plotVarianceComponents(m, interactive = FALSE), m))
})

test_that("plotFixedEffects(): forest panels per response with the fixedEffects() estimates", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plotFixedEffects(m, interactive = TRUE))
  expect_equal(panelTitles(w), paste0(c("k", "cnt", "x1", "x2"), ": fixed effects"))   # the order passed
  expect_false(w$x$linked)
  expect_null(w$x$search)
  fe <- fixedEffects(m)
  cntP <- w$x$panels[[2]]
  expect_equal(cntP$ylabels$text, c("(Intercept)", "grpb", "grpc"))     # factor levels expanded
  expect_equal(cntP$intervals$mean, signif(fe$mean[fe$response == "cnt"], 5))
  expect_equal(cntP$intervals$lo, signif(fe$q0.025[fe$response == "cnt"], 5))
  expect_true(all(cntP$intervals$lo <= cntP$intervals$hi))
  expect_equal(cntP$vline, 0)
  wr <- plotFixedEffects(m, response = c("x2", "k"), interactive = TRUE)
  expect_equal(panelTitles(wr), paste0(c("k", "x2"), ": fixed effects"))  # the order of the model, not of `response`
  expect_error(plotFixedEffects(m, response = "zzz", interactive = TRUE), "Unknown response")
  expectCleanWidget(plotFixedEffects(fxModel("surv"), interactive = TRUE))
  withNullDevice(expect_identical(plotFixedEffects(m, interactive = FALSE), m))
})

test_that("plotFixedEffects(type = \"density\"): one curve per term, the marginal of that term", {
  skipNoWidgets()
  m <- fxModel("multi")
  wd <- expectCleanWidget(plotFixedEffects(m, type = "density", interactive = TRUE))
  expect_equal(panelTitles(wd), paste0(c("k", "cnt", "x1", "x2"), ": fixed effects"))
  expect_equal(vapply(wd$x$panels[[2]]$curves, `[[`, "", "label"), c("(Intercept)", "grpb", "grpc"))
  for (pn in wd$x$panels) {
    expect_equal(pn$vline, 0)
    expect_equal(pn$xlab, "coefficient")
    r <- sub(": fixed effects$", "", pn$title)
    cols <- m$inla$fixed[[r]]$columns                 # INLA column name -> term
    for (cv in pn$curves) {
      marg <- m$fit$marginals.fixed[[names(cols)[cols == cv$label]]]
      expect_equal(cv$x, signif(marg[, 1], 5))
      expect_equal(cv$y, signif(marg[, 2], 5))
    }
  }
  four <- expectCleanWidget(plotFixedEffects(fxModel("priors"), type = "density", interactive = TRUE))
  expect_equal(panelTitles(four), paste0(c("yg", "x1", "x2", "x3"), ": fixed effects"))
  expect_true(all(vapply(four$x$panels, function(p) identical(vapply(p$curves, `[[`, "", "label"), c("(Intercept)", "t")), NA)))
  withNullDevice(expect_identical(plotFixedEffects(m, "cnt", type = "density", ncol = 1, interactive = FALSE), m))
})

test_that("plotAssociations(): one interval / curve per scaled() term, labelled long -> outcome", {
  skipNoWidgets()
  m <- fxModel("multi")
  a <- associations(m)
  w <- expectCleanWidget(plotAssociations(m, interactive = TRUE))
  expect_length(w$x$panels, 1)
  expect_equal(w$x$panels[[1]]$ylabels$text, c("x1 -> k", "x2 -> k", "x2 -> cnt"))
  expect_equal(w$x$panels[[1]]$intervals$mean, signif(a$mean, 5))
  expect_equal(w$x$panels[[1]]$intervals$hi, signif(a$q0.975, 5))
  wd <- expectCleanWidget(plotAssociations(m, type = "density", interactive = TRUE))
  cv <- wd$x$panels[[1]]$curves
  expect_equal(vapply(cv, `[[`, "", "label"), paste0(a$longitudinal, " -> ", a$outcome))
  # the curve of each label is the marginal of that gamma (the copy's `Beta for <out>_FROM_<long>_...`)
  for (z in cv) {
    out <- sub(".* -> ", "", z$label); long <- sub(" -> .*", "", z$label)
    nm <- grep(sprintf("^Beta for %s_FROM_%s_", out, long), names(m$fit$marginals.hyperpar), value = TRUE)
    expect_length(nm, 1)
    expect_equal(z$x, signif(m$fit$marginals.hyperpar[[nm]][, 1], 5))
  }
  expect_equal(wd$x$panels[[1]]$vline, 0)
  expect_equal(wd$x$ncol, 1L)
  expectCleanWidget(plotAssociations(fxModel("surv"), type = "density", interactive = TRUE))
  withNullDevice({
    expect_identical(plotAssociations(m, interactive = FALSE), m)
    expect_identical(plotAssociations(m, type = "density", interactive = FALSE), m)
  })
})

test_that("plotAssociations(): the three gammas of one outcome", {
  skipNoWidgets()
  m <- fxModel("priors")
  a <- associations(m)
  expect_equal(a$longitudinal, c("x1", "x2", "x3"))
  w <- expectCleanWidget(plotAssociations(m, interactive = TRUE))
  expect_equal(w$x$panels[[1]]$ylabels$text, paste0(c("x1", "x2", "x3"), " -> yg"))
  expect_equal(w$x$panels[[1]]$intervals$mean, signif(a$mean, 5))
  wd <- expectCleanWidget(plotAssociations(m, type = "density", interactive = TRUE))
  expect_equal(vapply(wd$x$panels[[1]]$curves, `[[`, "", "label"), paste0(c("x1", "x2", "x3"), " -> yg"))
  withNullDevice(expect_identical(plotAssociations(m, interactive = FALSE), m))
})

test_that("densityPanels(): base and widget, colours cycle, vline", {
  dens <- function(mu) cbind(seq(mu - 3, mu + 3, length.out = 50), stats::dnorm(seq(-3, 3, length.out = 50)))
  panels <- list(list(title = "A", xlab = "a", curves = lapply(1:7, function(i) list(label = paste0("c", i),
                                                                                      density = dens(i)))),
                 list(title = "B", xlab = "b", curves = list(list(label = "one", density = dens(0)))))
  withNullDevice({
    expect_identical(densityPanels(panels, NULL, FALSE, vline = 0, x = "obj"), "obj")
    expect_invisible(densityPanels(panels, 1L, FALSE, x = "obj"))
  })
  skipNoWidgets()
  w <- expectCleanWidget(densityPanels(panels, NULL, TRUE, vline = 0))
  expect_equal(w$x$ncol, 2L)
  expect_false(w$x$linked)
  cols <- vapply(w$x$panels[[1]]$curves, `[[`, "", "color")
  expect_equal(cols[7], cols[1])                         # six colours, then repeated
  expect_equal(length(unique(cols[1:6])), 6)
  expect_equal(w$x$panels[[1]]$curves[[3]]$id, "A c3")
  expect_equal(w$x$panels[[1]]$vline, 0)
  expect_null(densityPanels(panels, 1L, TRUE)$x$panels[[1]]$vline)
})

# ---- plotResiduals() -------------------------------------------------------------------------

test_that("plotResidualsData(): panels, limits containing 0, connect and hline", {
  m <- fxModel("multi")
  d <- plotResidualsData(m, NULL, "time", NULL, NULL, 30, NULL)
  expect_equal(vapply(d$panels, `[[`, "", "response"), c("k", "cnt", "x1", "x2"))
  expect_equal(d$panels[[1]]$ylab, "residual of k / n")
  expect_equal(d$panels[[2]]$ylab, "residual of cnt / pt")
  res <- residuals(m)
  for (pn in d$panels) {
    expect_true(pn$connect)
    expect_equal(pn$hline, 0)
    expect_length(pn$subjects, 30)
    expect_true(all(is.finite(c(pn$xlim, pn$ylim))) && pn$xlim[1] <= pn$xlim[2])
    expect_true(pn$ylim[1] <= 0 && pn$ylim[2] >= 0)
    ys <- unlist(lapply(pn$subjects, `[[`, "y"))
    expect_true(all(ys >= pn$ylim[1] & ys <= pn$ylim[2]))
    s <- pn$subjects[[1]]
    rows <- which(res$response == pn$response & res$id == s$id)
    expect_equal(s$y, res$residual[rows][order(res$t[rows])])
  }
})

test_that("plotResidualsData(): type = \"fitted\" puts the fitted mean on the x axis, unconnected", {
  m <- fxModel("multi")
  f <- plotResidualsData(m, "cnt", "fitted", NULL, NULL, 1, NULL)
  expect_true(is.na(f$time))
  expect_false(f$panels[[1]]$connect)
  expect_equal(f$panels[[1]]$xlab, "fitted mean")
  expect_length(f$panels[[1]]$subjects, 1)
  expect_false(is.unsorted(f$panels[[1]]$subjects[[1]]$x))
  res <- residuals(m, "cnt")
  s <- f$panels[[1]]$subjects[[1]]
  expect_setequal(s$x, res$mean[res$id == s$id])
})

test_that("plotResidualsData(): n, ids, group and explicit time select the subjects and axes", {
  m <- fxModel("multi")
  sel <- plotResidualsData(m, NULL, "time", NULL, NULL, 30, c("S09", "nobody"))
  for (pn in sel$panels) expect_equal(subjectIds(pn), "S09")
  big <- plotResidualsData(m, "x1", "time", NULL, NULL, 1000, NULL)
  expect_setequal(subjectIds(big$panels[[1]]), allIds())
  na <- plotResidualsData(m, "x1", "time", NULL, NA, 30, NULL)
  expect_false(na$panels[[1]]$connect)
  expect_equal(subjectIds(na$panels[[1]]), "all")
  e <- plotResidualsData(m, "x1", "time", "t", "id", 30, NULL)
  expect_equal(e$panels, plotResidualsData(m, "x1", "time", NULL, NULL, 30, NULL)$panels)
  site <- plotResidualsData(fxModel("priors"), c("yg", "x1"), "time", "t", "site", 30, NULL)
  expect_equal(subjectIds(site$panels[[1]]), levels(fxData()$long$site))   # yg: random effect on site
  expect_length(site$panels[[2]]$subjects, 0)                               # x1: no site column
  expect_error(plotResidualsData(m, NULL, "time", "nope", NULL, 30, NULL), "not a variable")
  expect_error(plotResidualsData(m, NULL, "time", NULL, "nope", 30, NULL), "not a variable")
  expect_error(plotResidualsData(m, "zzz", "time", NULL, NULL, 30, NULL), "Unknown response")
})

test_that("plotResidualsData(): a time variable present in only some equations gives empty panels", {
  m <- fxModel("multi")
  pt <- plotResidualsData(m, NULL, "time", "pt", NULL, 30, NULL)         # exposure of cnt only
  expect_equal(vapply(pt$panels, function(p) isTRUE(p$missing), NA), c(TRUE, FALSE, TRUE, TRUE))
  withNullDevice(expect_identical(plotResiduals(m, time = "pt", interactive = FALSE), m))
})

test_that("plotResiduals(): survival responses dropped, a warning when named, an error when alone", {
  m <- fxModel("surv")
  d <- plotResidualsData(m, NULL, "time", NULL, NULL, 30, NULL)
  expect_equal(vapply(d$panels, `[[`, "", "response"), "x1")
  expect_warning(d2 <- plotResidualsData(m, c("time", "x1"), "time", NULL, NULL, 30, NULL), "time; skipped")
  expect_equal(vapply(d2$panels, `[[`, "", "response"), "x1")
  expect_error(plotResidualsData(m, "time", "time", NULL, NULL, 30, NULL), "survival responses")
  withNullDevice({
    expect_identical(plotResiduals(m, interactive = FALSE), m)
    expect_warning(plotResiduals(m, c("time", "x1"), interactive = FALSE), "skipped")
  })
})

test_that("plotResiduals(): base graphics variants", {
  m <- fxModel("multi")
  withNullDevice({
    expect_identical(plotResiduals(m, interactive = FALSE), m)
    expect_identical(plotResiduals(m, time = "t", group = "id", interactive = FALSE), m)
    expect_identical(plotResiduals(fxModel("priors"), group = "site", interactive = FALSE), fxModel("priors"))
    expect_identical(plotResiduals(m, type = "fitted", n = 1, ncol = 4, interactive = FALSE), m)
    expect_identical(plotResiduals(m, "cnt", group = NA, interactive = FALSE), m)
    expect_identical(plotResiduals(m, ids = "nobody", interactive = FALSE), m)
    expect_error(plotResiduals(m, time = "nope", interactive = FALSE), "not a variable")
  })
})

test_that("plotResiduals() widget: residual lines with hline 0, legend without fitted entries", {
  skipNoWidgets()
  m <- fxModel("multi")
  w <- expectCleanWidget(plotResiduals(m, interactive = TRUE))
  expect_equal(w$x$ncol, 2L)
  expect_true(w$x$linked)
  expect_equal(w$x$search, "id id")
  expect_equal(legendKinds(w), c("res", "resMean"))
  for (pn in w$x$panels) {
    expect_equal(pn$hline, 0)
    expect_setequal(lineIds(pn, "res"), allIds())
    expect_equal(vapply(pn$curves, `[[`, "", "kind"), "resMean")
  }
  wf <- expectCleanWidget(plotResiduals(m, "k", type = "fitted", ids = c("S01", "S02"), interactive = TRUE))
  expect_length(wf$x$panels, 1)
  expect_equal(lineIds(wf$x$panels[[1]]), c("S01", "S02"))
  expect_false(wf$x$panels[[1]]$connect)
  ws <- expectCleanWidget(plotResiduals(fxModel("priors"), group = "site", interactive = TRUE))
  expect_equal(ws$x$search, "site id")
  expect_equal(lineIds(ws$x$panels[[1]]), levels(fxData()$long$site))
  expectCleanWidget(plotResiduals(fxModel("surv"), interactive = TRUE))
})

# ---- plot.mvremPredictions() -----------------------------------------------------------------

test_that("plot.mvremPredictions(): base graphics variants and errors", {
  p <- fxPredict()
  expect_equal(unique(p$response), c("k", "cnt", "x1", "x2"))
  withNullDevice({
    expect_identical(plot(p, interactive = FALSE), p)
    expect_identical(plot(p, time = "t", group = "id", interactive = FALSE), p)
    expect_identical(plot(p, data = FALSE, ncol = 1, interactive = FALSE), p)
    expect_identical(plot(p, n = 1, interactive = FALSE), p)
    expect_identical(plot(p, group = NA, interactive = FALSE), p)
    expect_error(plot(p, time = "nope", interactive = FALSE), "not a variable")
    expect_error(plot(p, time = "x1", interactive = FALSE), "not a column of the predictions")
  })
  pl <- p
  attr(pl, "type") <- "linear"          # the linear-predictor branch: no observed data drawn
  withNullDevice(expect_identical(plot(pl, interactive = FALSE), pl))
})

test_that("plot.mvremPredictions() widget: a predicted line and two bounds per subject over the data", {
  skipNoWidgets()
  p <- fxPredict()
  newIds <- unique(fxNewdata$id)                               # three new subjects and S01
  w <- expectCleanWidget(plot(p, interactive = TRUE))
  expect_true(w$x$linked)
  expect_equal(w$x$search, "id id")
  expect_equal(w$x$ncol, 2L)
  expect_match(panelTitles(w)[1], "^k .*outcome")
  expect_match(panelTitles(w)[3], "^x1 .*longitudinal")
  expect_equal(legendKinds(w), c("obs", "pred", "predInt"))
  for (pn in w$x$panels) {
    expect_true(pn$connect)
    expect_equal(lineIds(pn, "pred"), newIds)
    expect_equal(sum(lineKinds(pn) == "pred"), length(newIds))         # one mean line per subject ...
    expect_equal(sum(lineKinds(pn) == "predInt"), 2 * length(newIds))  # ... and two bounds
    expect_setequal(lineIds(pn, "obs"), allIds())            # n = 300: all observed subjects
    r <- sub(" .*", "", pn$title)
    d <- p[p$response == r, ]
    expect_true(pn$xlim[1] <= min(d$t) && pn$xlim[2] >= max(d$t))   # covers the predictions past the data
    tol <- 1e-4 * diff(pn$ylim)                                      # ylim is rounded to 5 digits
    expect_true(pn$ylim[1] <= min(d$q0.025) + tol && pn$ylim[2] >= max(d$q0.975) - tol)
    pi <- Filter(function(l) l$kind == "predInt", pn$lines)
    expect_true(all(vapply(pi, function(l) isFALSE(l$points), NA)))
    s01 <- Filter(function(l) l$kind == "pred" && l$id == "S01", pn$lines)[[1]]
    expect_equal(s01$x, sort(d$t[d$id == "S01"]))
    expect_equal(s01$y, signif(d$mean[d$id == "S01"][order(d$t[d$id == "S01"])], 5))
  }
  expect_equal(plot(p, time = "t", group = "id", interactive = TRUE)$x, w$x)
})

test_that("plot.mvremPredictions() widget: n, data = FALSE and the linear-predictor type", {
  skipNoWidgets()
  p <- fxPredict()
  # at most n observed subjects (the new subjects, without observations, are among the n picked)
  expect_lte(length(lineIds(plot(p, n = 3, interactive = TRUE)$x$panels[[1]], "obs")), 3)
  wn <- plot(p, data = FALSE, interactive = TRUE)
  expect_length(lineIds(wn$x$panels[[1]], "obs"), 0)
  expect_equal(legendKinds(wn), c("pred", "predInt"))
  pl <- p
  attr(pl, "type") <- "linear"
  wl <- expectCleanWidget(plot(pl, interactive = TRUE))
  expect_equal(legendKinds(wl), c("pred", "predInt"))
  expect_equal(wl$x$panels[[3]]$ylab, "linear predictor of x1")
  expect_length(lineIds(wl$x$panels[[1]], "obs"), 0)
})

test_that("plot.mvremPredictions() widget: group = NA pools the predictions into one line", {
  skipNoWidgets()
  wg <- expectCleanWidget(plot(fxPredict(), group = NA, interactive = TRUE))
  expect_equal(wg$x$search, "subject id")
  expect_equal(lineIds(wg$x$panels[[1]], "pred"), "prediction")
  expect_equal(lineIds(wg$x$panels[[1]], "obs"), "all")
})

# ---- the widget's JavaScript (node) ----------------------------------------------------------

test_that("mvremPlot.js: renders every widget kind in node and its interactions work", {
  skipNoWidgets()
  skip_if_not(hasNode(), "node not available")
  m <- fxModel("multi")
  widgets <- list(
    lines = plot(m, interactive = TRUE),
    survival = plot(fxModel("surv"), interactive = TRUE),           # an empty panel
    caterpillar = plotRandomEffects(m, interactive = TRUE),
    pairs = plotRandomEffects(m, type = "pairs", interactive = TRUE),
    densities = plotVarianceComponents(m, interactive = TRUE),
    forest = plotFixedEffects(m, interactive = TRUE),
    residuals = plotResiduals(m, interactive = TRUE),
    predictions = plot(fxPredict(), interactive = TRUE))
  res <- runWidgetJs(widgets)
  for (nm in names(widgets)) {
    r <- res[[nm]]
    x <- widgets[[nm]]$x
    expect_true(r$ok, label = paste0(nm, ": ", r$error, "\n", r$stack))
    if (!isTRUE(r$ok)) next
    expect_equal(r$panels, length(x$panels), label = paste(nm, "panels"))
    expect_equal(r$legend, length(x$legend), label = paste(nm, "legend entries"))
    expect_identical(r$hasSearch, !is.null(x$search), label = paste(nm, "search box"))
    expect_true(r$css)
    expect_true(r$legendOff, label = paste(nm, "legend toggles"))
    expect_true(r$kindRules, label = paste(nm, "hide rules of the legend kinds"))
    # highlightable elements: lines, intervals, points and curves with an id
    nGroups <- sum(vapply(x$panels, function(p)
      length(p$lines) + length(p$intervals$id) + length(p$points$id) +
        sum(vapply(p$curves, function(cv) !is.null(cv$id), NA)), 0))
    expect_equal(r$groups, nGroups, label = paste(nm, "elements"))
    expect_true(r$tipShown, label = paste(nm, "tooltip"))
    expect_match(r$tipText, "click to pin")
    expect_true(r$highlighted, label = paste(nm, "highlight"))
    expect_true(r$tipHiddenAfterLeave, label = paste(nm, "tooltip hidden"))
    expect_true(r$pinnedAfterLeave, label = paste(nm, "pinned"))
    expect_true(r$clearedByEscape, label = paste(nm, "Escape"))
    if (!is.null(x$search)) {
      expect_gte(r$searchMatches, 1)
      expect_true(r$dimAfterSearch, label = paste(nm, "search dims"))
      expect_false(r$dimAfterClear, label = paste(nm, "search cleared"))
    }
  }
  # the tooltip names the panel, the subject and the values
  expect_match(res$lines$tipText, "^k .*\\nS01 observed\\nt = ")
  expect_match(res$pairs$tipText, "^x1: \\(Intercept\\) \\+ t \\| id\\nS01\\nposterior mean of")
})

test_that("mvremPlot.js: a legend kind unknown to the JavaScript can be switched off", {
  skipNoWidgets()
  skip_if_not(hasNode(), "node not available")
  pn <- list(title = "p", xlab = "x", ylab = "y", xlim = c(0, 1), ylim = c(0, 1), connect = TRUE,
             lines = list(list(id = "a", kind = "newKind", x = c(0, 1), y = c(0, 1))))
  w <- htmlPlot(list(pn), 1L, legend = list(list(kind = "newKind", label = "new", color = "#000000", dash = "")))
  res <- runWidgetJs(list(w = w))$w
  expect_true(res$ok, label = paste(res$error, res$stack))
  expect_true(res$legendOff)
  expect_true(res$kindRules)
})

test_that("renderSpecs.js: a broken spec is reported as an error with a non-zero exit status", {
  skip_if_not_installed("jsonlite")
  skip_if_not(hasNode(), "node not available")
  bad <- withr::local_tempfile(fileext = ".json", lines = "{ not json")
  js <- system.file("htmlwidgets", "mvremPlot.js", package = "mvrem")
  out <- suppressWarnings(system2("node", shQuote(c(test_path("js", "renderSpecs.js"), js, bad)),
                                  stdout = TRUE, stderr = TRUE))
  expect_equal(attr(out, "status"), 1L)
  res <- jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE)
  expect_length(res, 1)
  expect_false(res[[1]]$ok)
  expect_match(res[[1]]$error, "JSON")
})

test_that("shareXlim(): the panels of a plot against time share one x range", {
  pn <- list(list(xlim = c(0, 1.5)), list(xlim = c(0.25, 1.75)), list(missing = TRUE, title = "s"),
             list(xlim = c(NA, NA)))
  sh <- shareXlim(pn)
  expect_equal(sh[[1]]$xlim, c(0, 1.75))
  expect_equal(sh[[2]]$xlim, c(0, 1.75))
  expect_identical(sh[3:4], pn[3:4])
  expect_identical(shareXlim(pn[1]), pn[1])
  # the equations of a model: one range over all panels, also for residuals against time
  m <- fxModel("multi", fit = FALSE)
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE)
  lims <- lapply(Filter(function(p) !isTRUE(p$missing), d$panels), `[[`, "xlim")
  expect_true(length(unique(lims)) == 1L)
})

test_that("plotLinesData(): INLA's predictions of missing observations are open circles (E24)", {
  skip_if_not_installed("INLA")
  d <- fxData()$long
  miss <- c(2L, 7L, 11L)
  d$x1[miss] <- NA
  d$y[5L] <- NA
  m <- suppressMessages(jsm(list(y ~ t + scaled(x1), x1 ~ t + (1 + t | id)), data = d,
                            family = list(y = betaFamily())))
  skip_if(is.null(inlaFit(m)), "INLA fit failed")
  pd <- plotLinesData(m, NULL, NULL, NULL, 1000, NULL, TRUE)
  px <- pd$panels[[2]]
  expect_equal(px$response, "x1")
  pred <- px$predicted
  expect_equal(sum(lengths(lapply(pred, `[[`, "x"))), length(miss))
  # the circles are the fitted means of exactly those rows, and not on the subject lines
  fv <- fittedValues(m, "x1")
  for (s in pred) {
    k <- which(d$id == s$id & is.na(d$x1))
    expect_equal(s$x, d$t[k][order(d$t[k])])
    expect_equal(s$fit, fv$mean[k][order(d$t[k])])
    obs <- px$subjects[[which(vapply(px$subjects, `[[`, "", "id") == s$id)]]
    expect_false(any(s$x %in% obs$x))
  }
  expect_length(pd$panels[[1]]$predicted, 1L)
  expect_true(hasPredicted(pd))
  # legend: a circle entry, only when there are such predictions
  lg <- linesLegend(TRUE, predicted = TRUE)
  expect_equal(lg[[5]]$kind, "fitMissing")
  expect_equal(lg[[5]]$shape, "circle")
  expect_length(linesLegend(TRUE), 4L)
  expect_false(hasPredicted(plotLinesData(fxModel("beta"), NULL, NULL, NULL, 30, NULL, TRUE)))
  # base graphics draws them; the widget gets unjoined hollow lines of kind fitMissing
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(m, interactive = FALSE))
  skipNoWidgets()
  w <- plot(m, interactive = TRUE)
  ln <- Filter(function(l) identical(l$kind, "fitMissing"), w$x$panels[[2]]$lines)
  expect_length(ln, length(pred))
  expect_true(all(vapply(ln, function(l) isTRUE(l$hollow), NA)))
  expect_true(any(vapply(w$x$legend, function(e) identical(e$shape, "circle"), NA)))
  skip_if(!hasNode(), "node not available")
  res <- runWidgetJs(list(w = w))$w
  expect_true(res$ok, label = paste(res$error, res$stack))
})

test_that("plotModel(): a fitted model's boxes show the SDs, correlations and residual SD (E6)", {
  m <- fxModel("priors")
  vc <- varianceComponents(m)
  nodes <- diagramNodes(m$forms, vc)
  lines <- function(id) nodes[[which(vapply(nodes, `[[`, "", "id") == id)]]$lines
  t <- vc$table
  num <- function(v) formatC(signif(v, 2), digits = 2, format = "fg", flag = "#")
  x1 <- t[t$response == "x1", ]
  expect_true(paste0("SD | id: (Intercept) ", num(x1$mean[x1$effect == "(Intercept)" & x1$parameter == "sd"]),
                     ", t ", num(x1$mean[x1$effect == "t" & x1$parameter == "sd"]),
                     "; cor ", num(x1$mean[x1$parameter == "correlation"])) %in% lines("x1"))
  expect_true(paste("residual SD", num(x1$mean[x1$block == "residual"])) %in% lines("x1"))
  expect_equal(sum(startsWith(lines("x3"), "SD | id: ")), 2L)        # `||`: one line per effect
  expect_true(any(startsWith(lines("yg"), "SD | site: (Intercept) ")))
  # no estimates without a fit or with fitted = FALSE
  expect_identical(diagramNodes(m$forms), diagramNodes(fxModel("priors", fit = FALSE)$forms))
  expect_false(any(grepl("^SD|residual SD", unlist(lapply(diagramNodes(m$forms), `[[`, "lines")))))
  expect_no_error(withNullDevice(plotModel(m, fitted = FALSE)))
  expect_no_error(withNullDevice(plotModel(m)))
})

test_that("colorBy: candidates, checks, series values, base and widget renderers", {
  m <- fxModel("priors", fit = FALSE)
  prep <- bothParts(m$data)
  expect_equal(colorCandidates(prep), "site")                    # id has too many levels
  expect_error(colorSpec(prep, "id"), "`colorBy` must be one of .*: site\\.")
  expect_error(colorSpec(prep, 1), "`colorBy` must be one of")
  sp <- colorSpec(prep, "site")
  expect_equal(sp$levels$site, c("A", "B", "C", "D"))
  expect_equal(colorOf(c("B", NA, "Z"), sp, "grey"), c(plotColors$groups[2], "grey", "grey"))
  expect_equal(colorOf("B", colorSpec(prep), "grey"), "grey")    # nothing selected
  d <- plotLinesData(m, NULL, NULL, NULL, 30, NULL, TRUE, colorBy = "site")
  expect_equal(d$colors$selected, "site")
  long <- fxData()$long
  s1 <- d$panels[[1]]$subjects[[1]]
  rows <- which(long$id == s1$id)
  expect_equal(s1$cv$site, as.character(long$site[rows][order(long$t[rows])]))
  expect_no_error(withNullDevice(plot(m, colorBy = "site", interactive = FALSE)))
  expect_error(plot(m, colorBy = "nope", interactive = FALSE), "must be one of")
  # no candidates: no dropdown in the widget, plain grey plot
  expect_length(colorSpec(bothParts(fxModel("beta", fit = FALSE)$data))$vars, 0L)
  skipNoWidgets()
  w <- plot(m, colorBy = "site", interactive = TRUE)
  expect_equal(w$x$colorBy$vars, "site")
  expect_equal(w$x$colorBy$selected, "site")
  ln <- Filter(function(l) l$kind == "obs", w$x$panels[[1]]$lines)
  expect_true(all(vapply(ln, function(l) length(l$cv$site) == length(l$x), NA)))
  expect_null(plot(fxModel("beta", fit = FALSE), interactive = TRUE)$x$colorBy)
  skip_if(!hasNode(), "node not available")
  res <- runWidgetJs(list(w = plot(m, interactive = TRUE)))$w
  expect_true(res$ok, label = paste(res$error, res$stack))
  expect_equal(res$colorOptions, 2L)                             # gray + site
  expect_true(res$colored)
  expect_equal(res$colorLevels, 4L)
  expect_true(res$colorReset)
  # residuals: the same option
  f <- fxModel("priors")
  rd <- plotResidualsData(f, NULL, "time", NULL, NULL, 30, NULL, colorBy = "site")
  expect_true(any(vapply(rd$panels[[1]]$subjects, function(s) length(s$cv$site) > 0L, NA)))
  expect_no_error(withNullDevice(plotResiduals(f, colorBy = "site", interactive = FALSE)))
})

test_that("plot.mvremPredictions(): colour dropdown detected automatically, colorBy", {
  p <- fxPredict("multi")
  expect_no_error(withNullDevice(plot(p, colorBy = "grp", interactive = FALSE)))
  expect_error(plot(p, colorBy = "id", interactive = FALSE), "must be one of")
  skipNoWidgets()
  w <- plot(p, interactive = TRUE)
  expect_equal(w$x$colorBy$vars, "grp")                          # no colorBy needed
  expect_null(w$x$colorBy$selected)
  ln <- unlist(lapply(w$x$panels, `[[`, "lines"), recursive = FALSE)
  pred <- Filter(function(l) l$kind %in% c("pred", "predInt"), ln)
  expect_true(all(vapply(pred, function(l) identical(l$cv$grp, rep("a", length(l$x))), NA)))
  expect_equal(plot(p, colorBy = "grp", interactive = TRUE)$x$colorBy$selected, "grp")
  skip_if(!hasNode(), "node not available")
  res <- runWidgetJs(list(w = w))$w
  expect_true(res$ok, label = paste(res$error, res$stack))
  expect_equal(res$colorOptions, 2L)
  expect_true(res$colored)
})
