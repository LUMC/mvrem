# Tests of E18: one list of formulae for both parts (jsm(form) / collectFormulas(form)),
# splitFormulas(), formulaParts() and update() of either style. The pair style
# (formOut / formLong) is tested in the other files; here it is the reference.

# ---- helpers ---------------------------------------------------------------

# jsm() arguments of fixture `name` without its formulae.
fxOtherArgs <- function(s) s[setdiff(names(s), c("formOut", "formLong"))]

# Unfitted jsm() of fixture `name` with the formulae given as `form` (default: the
# combined list, longitudinal first) and the same model in the pair style.
fxBoth <- function(name, form = NULL) {
  s <- fxSpec(name)
  if (is.null(form)) form <- c(s$formLong, s$formOut)
  if (is.function(form)) form <- form(s)
  build <- function(args) suppressWarnings(suppressMessages(do.call(jsm, c(args, fxOtherArgs(s), fit = FALSE))))
  list(comb = build(list(form = form)),
       pair = build(list(formOut = s$formOut, formLong = s$formLong)))
}

# A collectFormulas() result without the order the formulae were passed in (`order`).
unordered <- function(f) f[setdiff(names(f), "order")]

# The data summary with its equations / missing rows in the pair order (outcomes, then
# longitudinal), whatever the order the formulae were passed in.
pairOrderSummary <- function(m) {
  s <- m$dataSummary
  o <- c(names(m$forms$out), names(m$forms$long))
  for (t in c("equations", "missing")) {
    s[[t]] <- s[[t]][order(match(s[[t]]$response, o), seq_len(nrow(s[[t]]))), , drop = FALSE]
    rownames(s[[t]]) <- NULL
  }
  s
}

# Two models are the same up to the order the formulae were passed in (Q3: the results
# follow it, tested below): parsed formulae, data, data summary, INLA specification
# (stacks, families, terms, ...) and the INLA formula (text and LEVELS; its environment is
# a new one per build).
expectSameModel <- function(a, b) {
  expect_identical(unordered(a$forms), unordered(b$forms))
  expect_identical(a$data, b$data)
  expect_identical(pairOrderSummary(a), pairOrderSummary(b))
  ia <- a$inla; ib <- b$inla
  expect_identical(ia[setdiff(names(ia), "formula")], ib[setdiff(names(ib), "formula")])
  expect_identical(ia$stacks, ib$stacks)
  expect_identical(deparse(inlaFormula(a)), deparse(inlaFormula(b)))
  expect_identical(environment(inlaFormula(a))$LEVELS, environment(inlaFormula(b))$LEVELS)
  expect_equal(inlaFormula(a), inlaFormula(b))
  expect_identical(a$options, b$options)
}

# The message of the error raised by `expr` (fails if there is none).
errorText <- function(expr) conditionMessage(expect_error(expr))

# ---- reference models: combined = pair --------------------------------------

for (nm in c("beta", "multi", "surv", "priors")) {
  test_that(paste0("jsm(form): the combined list gives the pair model (fixture ", nm, ")"), {
    m <- fxBoth(nm)
    expectSameModel(m$comb, m$pair)
    expect_identical(unordered(m$comb$forms), unordered(fxForms(nm)))
    # the order passed: longitudinal first here (fxBoth's default), outcomes first for the pair
    f <- fxForms(nm)
    expect_identical(m$comb$forms$order, c(names(f$long), names(f$out)))
    expect_identical(m$pair$forms$order, c(names(f$out), names(f$long)))
  })
}

test_that("jsm(form): any interleaving keeps the order within each part and gives the same model", {
  orders <- list(
    function(s) c(s$formOut, s$formLong),
    function(s) c(s$formLong, s$formOut),
    function(s) list(s$formLong[[1]], s$formOut[[1]], s$formLong[[2]], s$formOut[[2]]),
    function(s) list(s$formOut[[1]], s$formLong[[1]], s$formOut[[2]], s$formLong[[2]]))
  for (o in orders) {
    m <- fxBoth("multi", o)
    expectSameModel(m$comb, m$pair)
  }
  # the order within a part follows the order given
  m <- fxBoth("multi", function(s) rev(c(s$formOut, s$formLong)))
  expect_identical(names(m$comb$forms$out), c("cnt", "k"))
  expect_identical(names(m$comb$forms$long), c("x2", "x1"))
  expect_identical(m$comb$forms$order, c("x2", "x1", "cnt", "k"))
  expect_identical(m$comb$inla$responses, m$pair$inla$responses[c(2, 1, 4, 3)])
})

test_that("jsm(form): names of the list are ignored", {
  m <- fxBoth("beta", function(s) list(a = s$formLong[[1]], b = s$formOut[[1]]))
  expectSameModel(m$comb, m$pair)
})

test_that("jsm(): old positional pair calls still work and are stored as formOut / formLong", {
  s <- fxSpec("beta")
  m1 <- suppressMessages(jsm(s$formOut, s$formLong, s$data, family = s$family, fit = FALSE))
  m2 <- suppressMessages(jsm(form = s$formOut, formLong = s$formLong, data = s$data, family = s$family, fit = FALSE))
  m3 <- suppressMessages(jsm(formLong = s$formLong, formOut = s$formOut, data = s$data, family = s$family, fit = FALSE))
  expectSameModel(m1, fxBoth("beta")$pair)
  expectSameModel(m2, m1)
  expectSameModel(m3, m1)
  for (m in list(m1, m2, m3)) {
    expect_true(all(c("formOut", "formLong") %in% names(m$call)))
    expect_false("form" %in% names(m$call))
  }
  # the combined style is stored as `form`
  mc <- fxBoth("beta")$comb
  expect_true("form" %in% names(mc$call))
  expect_false(any(c("formOut", "formLong") %in% names(mc$call)))
})

test_that("jsm(form): a single formula with scaled() but no longitudinal formula fails as before", {
  expect_error(jsm(y ~ t + scaled(x1), data = fxSmallData(), fit = FALSE),
               "unknown longitudinal response\\(s\\): x1")
})

# ---- single formulae ---------------------------------------------------------

test_that("collectFormulas(form): single formulae", {
  # a single formula as `form` is a one-element list
  expect_identical(splitFormulas(y ~ t + scaled(x)), list(out = list(y ~ t + scaled(x)), long = list()))
  f <- collectFormulas(list(y ~ t + scaled(x), x ~ t + (1 | id)))
  expect_identical(f, collectFormulas(y ~ t + scaled(x), x ~ t + (1 | id)))
  expect_identical(f, collectFormulas(formOut = y ~ t + scaled(x), formLong = x ~ t + (1 | id)))
  expect_identical(f$order, c("y", "x"))
  expect_identical(collectFormulas(list(x ~ t + (1 | id), y ~ t + scaled(x)))$order, c("x", "y"))
  # a single formula without scaled() is longitudinal and unlinked
  expect_error(collectFormulas(x ~ t), "`x ~ t` has no scaled\\(\\) term and is not the target of one")
})

# ---- collectFormulas() directly ----------------------------------------------

test_that("collectFormulas(form): the combined list gives the pair result, families included", {
  for (nm in c("beta", "multi", "surv", "priors")) {
    s <- fxSpec(nm)
    fl <- collectFormulas(c(s$formLong, s$formOut), family = s$family)
    expect_identical(unordered(fl), unordered(fxForms(nm)))
    expect_identical(fl$order, c(names(fl$long), names(fl$out)))
    expect_identical(collectFormulas(form = c(s$formOut, s$formLong), family = s$family), fxForms(nm))
  }
  # familyOut / familyLong key on response names after the split
  f <- collectFormulas(list(x1 ~ t, y ~ scaled(x1)), familyOut = list(y = "beta"),
                       familyLong = list(x1 = "poisson"))
  expect_equal(f$out$y$family$family, "beta")
  expect_equal(f$long$x1$family$family, "poisson")
  expect_error(collectFormulas(list(x1 ~ t, y ~ scaled(x1)), familyOut = list(x1 = "poisson")),
               "familyOut")
})

test_that("splitFormulas(): outcomes are the formulae with scaled(), order kept per part", {
  s <- splitFormulas(list(x1 ~ t + (1 | id), y ~ t + scaled(x1), x2 ~ t, z ~ scaled(x2) + scaled(x1)))
  expect_identical(s$out, list(y ~ t + scaled(x1), z ~ scaled(x2) + scaled(x1)))
  expect_identical(s$long, list(x1 ~ t + (1 | id), x2 ~ t))
  # LHS transformations and family inputs count by response variable
  s2 <- splitFormulas(list(k | trials(n) ~ scaled(x), log(x) ~ t))
  expect_length(s2$out, 1L)
  expect_length(s2$long, 1L)
})

# ---- errors ------------------------------------------------------------------

test_that("errors: a formula without scaled() that is no target is named", {
  msg <- errorText(collectFormulas(list(x1 ~ t, y ~ scaled(x1), z ~ t + (1 | id))))
  expect_match(msg, "`z ~ t + (1 | id)` has no scaled() term and is not the target of one", fixed = TRUE)
  expect_match(msg, "formOut` / `formLong`", fixed = TRUE)
  # all of them are listed
  msg2 <- errorText(collectFormulas(list(x1 ~ t, y ~ scaled(x1), z ~ t, w ~ 1)))
  expect_match(msg2, "`z ~ t`, `w ~ 1` have no scaled() term and are not the target of one", fixed = TRUE)
  # no scaled() at all: every formula is unlinked
  expect_error(collectFormulas(list(y ~ t, x ~ t)), "`y ~ t`, `x ~ t` have no scaled")
  # through jsm()
  expect_error(jsm(list(x ~ t, y ~ scaled(x), z ~ 1), data = fxSmallData(), fit = FALSE),
               "`z ~ 1` has no scaled")
})

test_that("errors: a scaled() target that is an outcome is named", {
  msg <- errorText(collectFormulas(list(x1 ~ t, y ~ scaled(x1), z ~ t + scaled(y))))
  expect_match(msg, "scaled(y) in `z ~ t + scaled(y)` targets `y ~ scaled(x1)`, which is an outcome",
               fixed = TRUE)
  # an outcome naming itself
  expect_error(collectFormulas(list(x ~ t, y ~ scaled(x) + scaled(y))),
               "scaled\\(y\\) in `y ~ scaled\\(x\\) \\+ scaled\\(y\\)` targets `y ~ scaled\\(x\\) \\+ scaled\\(y\\)`")
  expect_error(jsm(list(x ~ t, y ~ scaled(x), k ~ scaled(y)), data = fxSmallData(), fit = FALSE),
               "scaled\\(y\\) in `k ~ scaled\\(y\\)`.*is an outcome")
})

test_that("errors: problems found by the pair validation are reported as before", {
  # unknown target (reported before the then unlinked formula)
  expect_error(collectFormulas(list(y ~ scaled(x9), x1 ~ t)), "unknown longitudinal response\\(s\\): x9")
  # duplicated responses
  expect_error(collectFormulas(list(y ~ scaled(x), x ~ t, log(x) ~ t)), "more than once: x")
  expect_error(collectFormulas(list(y ~ scaled(x), x ~ t, y ~ 1)), "more than once: y")
  # repeated target in one formula, nested scaled(), syntax errors
  expect_error(collectFormulas(list(y ~ scaled(x) + scaled(x), x ~ t)), "more than once for: x")
  expect_error(collectFormulas(list(y ~ log(scaled(x)), x ~ t)), "top-level")
  expect_error(collectFormulas(list(y ~ scaled(x), ~ t)), "two-sided")
  expect_error(collectFormulas(list(y ~ scaled(x), cbind(a, b) ~ t)), "exactly one variable")
})

test_that("errors: formula arguments", {
  expect_error(collectFormulas(list(y ~ scaled(x), x ~ t), formOut = list(y ~ scaled(x))), "not both")
  expect_error(jsm(list(y ~ scaled(x), x ~ t), formOut = list(y ~ scaled(x)), data = fxSmallData(), fit = FALSE),
               "not both")
  expect_error(collectFormulas(formOut = list(y ~ scaled(x))), "`formOut` needs `formLong`")
  expect_error(collectFormulas(), "No model formulae")
  expect_error(jsm(data = fxSmallData(), fit = FALSE), "No model formulae")
  expect_error(collectFormulas(list(y ~ scaled(x), "x ~ t")), "formula or a non-empty list of formulae")
  expect_error(collectFormulas(list()), "formula or a non-empty list of formulae")
  expect_error(splitFormulas("y ~ t"), "formula or a non-empty list of formulae")
})

# ---- the pair style is still the way to an outcome without scaled() ----------

test_that("an outcome without scaled() needs the pair style", {
  f <- collectFormulas(list(y ~ t + (1 | id), z ~ t + scaled(x1)), list(x1 ~ t + (1 | id)))
  expect_identical(names(f$out), c("y", "z"))
  expect_error(collectFormulas(list(y ~ t + (1 | id), z ~ t + scaled(x1), x1 ~ t + (1 | id))),
               "`y ~ t \\+ \\(1 \\| id\\)` has no scaled")
})

# ---- update() ----------------------------------------------------------------

test_that("update(): both styles; `form` / `formOut` replace all formulae", {
  m <- fxBoth("beta")
  s <- fxSpec("beta")
  # other arguments: the style is kept
  cc <- update(m$comb, threads = 2, evaluate = FALSE)
  expect_true("form" %in% names(cc))
  cp <- update(m$pair, threads = 2, evaluate = FALSE)
  expect_true(all(c("formOut", "formLong") %in% names(cp)))
  mc <- suppressWarnings(suppressMessages(update(m$comb, family = list(y = "beta", x1 = studentTFamily()))))
  mp <- suppressWarnings(suppressMessages(update(m$pair, family = list(y = "beta", x1 = studentTFamily()))))
  expectSameModel(mc, mp)
  expect_equal(mc$forms$long$x1$family$family, "t")
  # pair -> combined
  newForm <- list(x1 ~ t + (1 | id), y ~ 1 + scaled(x1))
  m1 <- suppressMessages(update(m$pair, form = newForm))
  expect_false(any(c("formOut", "formLong") %in% names(m1$call)))
  expect_identical(names(m1$forms$out), "y")
  # combined -> pair
  m2 <- suppressMessages(update(m$comb, formOut = list(y ~ 1 + scaled(x1)), formLong = list(x1 ~ t + (1 | id))))
  expect_false("form" %in% names(m2$call))
  expectSameModel(m1, m2)
  # combined -> combined
  m3 <- suppressMessages(update(m$comb, form = rev(newForm)))
  expectSameModel(m3, m1)
  # formLong alone replaces the longitudinal formulae of a pair
  m4 <- suppressMessages(update(m$pair, formLong = list(x1 ~ t + (1 | id))))
  expect_identical(m4$forms$long, m1$forms$long)
  expect_identical(m4$forms$out, m$pair$forms$out)
})
