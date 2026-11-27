# E37: threads policy, availableThreads(), jsmApply()

test_that("checkThreads(): NULL, whole numbers, \"auto\" and \"A:B\" are accepted", {
  expect_null(checkThreads(NULL))
  expect_equal(checkThreads(3), 3L)
  expect_equal(checkThreads("auto"), "auto")
  expect_equal(checkThreads("4:2"), "4:2")
  for (bad in list(0, 1.5, -1, "four", "4", "4:", ":2", "a:b", c(1, 2), NA, TRUE))
    expect_error(checkThreads(bad), "must be NULL, a whole number")
})

test_that("jsm(threads = ): bad values are refused up front", {
  d <- fxSmallData()
  expect_error(jsm(list(y ~ t + scaled(x), x ~ t + (1 | id)), data = d, fit = FALSE, threads = "lots"),
               "threads")
})

test_that("availableThreads(): the smallest limit wins and is named", {
  noCgroup <- tempfile(); dir.create(noCgroup)
  withr::local_options(mvrem.cores = NULL)
  free <- availableThreads(env = c(), cgroup = noCgroup)
  expect_gte(free, 1L)
  expect_equal(attr(availableThreads(env = c(NSLOTS = "1"), cgroup = noCgroup), "source"), "NSLOTS")
  expect_equal(as.integer(availableThreads(env = c(SLURM_CPUS_PER_TASK = "1", NSLOTS = "64"), cgroup = noCgroup)), 1L)
  # R CMD check: at most 2
  chk <- availableThreads(env = c(`_R_CHECK_LIMIT_CORES_` = "TRUE"), cgroup = noCgroup)
  expect_lte(chk, 2L)
  # a limit above the machine's cores does not raise the result
  expect_equal(as.integer(availableThreads(env = c(NSLOTS = "100000"), cgroup = noCgroup)), as.integer(free))
  # junk values are ignored
  expect_equal(as.integer(availableThreads(env = c(NSLOTS = "many"), cgroup = noCgroup)), as.integer(free))
})

test_that("availableThreads(): container CPU quota (cgroups v2 and v1)", {
  withr::local_options(mvrem.cores = NULL)
  d2 <- tempfile(); dir.create(d2)
  writeLines("150000 100000", file.path(d2, "cpu.max"))
  r <- availableThreads(env = c(), cgroup = d2)
  expect_lte(r, 2L)
  if (parallel::detectCores(logical = FALSE) >= 2) expect_equal(attr(r, "source"), "cgroup")
  writeLines("max 100000", file.path(d2, "cpu.max"))
  expect_null(names(cgroupLimit(d2)))
  d1 <- tempfile(); dir.create(file.path(d1, "cpu"), recursive = TRUE)
  writeLines("100000", file.path(d1, "cpu", "cpu.cfs_period_us"))
  writeLines("100000", file.path(d1, "cpu", "cpu.cfs_quota_us"))
  expect_equal(cgroupLimit(d1), c(cgroup = 1))
  writeLines("-1", file.path(d1, "cpu", "cpu.cfs_quota_us"))     # v1: no quota
  expect_length(cgroupLimit(d1), 0L)
})

test_that("availableThreads(): options(mvrem.cores) overrides the detection", {
  withr::local_options(mvrem.cores = 3)
  expect_equal(as.integer(availableThreads()), 3L)
  expect_equal(attr(availableThreads(), "source"), "options(mvrem.cores)")
  withr::local_options(mvrem.cores = 0)
  expect_error(availableThreads(), "mvrem.cores")
})

test_that("availableThreads(): only the performance cores of an Apple Silicon Mac", {
  withr::local_options(mvrem.cores = NULL)
  noCgroup <- tempfile(); dir.create(noCgroup)
  local_mocked_bindings(macPerformanceCores = function() 1L)
  r <- availableThreads(env = c(), cgroup = noCgroup)
  expect_equal(as.integer(r), 1L)
  expect_equal(attr(r, "source"), "performance cores")
  local_mocked_bindings(macPerformanceCores = function() NA_integer_)
  expect_equal(attr(availableThreads(env = c(), cgroup = noCgroup), "source"), "cores")
  # the real sysctl: a count on Apple Silicon, NA elsewhere
  p <- macPerformanceCores()
  expect_true(is.na(p) || (p >= 1L && p <= parallel::detectCores()))
})

test_that("autoBackend(): fork on Linux and in a terminal on macOS, psock on Windows and in a GUI", {
  none <- c(RSTUDIO = "", POSITRON = "")
  expect_equal(autoBackend("windows", "Windows", "Rgui", none), "psock")
  expect_equal(autoBackend("unix", "Linux", "X11", none), "fork")
  expect_equal(autoBackend("unix", "Linux", "RStudio", c(RSTUDIO = "1", POSITRON = "")), "fork")
  expect_equal(autoBackend("unix", "Darwin", "unknown", none), "fork")     # Rscript
  expect_equal(autoBackend("unix", "Darwin", "X11", none), "fork")         # R in a terminal
  expect_equal(autoBackend("unix", "Darwin", "AQUA", none), "psock")       # R.app
  expect_equal(autoBackend("unix", "Darwin", "RStudio", c(RSTUDIO = "1", POSITRON = "")), "psock")
  expect_equal(autoBackend("unix", "Darwin", "unknown", c(RSTUDIO = "", POSITRON = "1")), "psock")
  expect_true(autoBackend() %in% c("fork", "psock"))
})

test_that("resolveThreads(): NULL, number, A:B, auto and the option of the workers", {
  withr::local_options(mvrem.threads = NULL, mvrem.cores = 4)
  expect_null(resolveThreads(NULL))
  expect_equal(resolveThreads(2), 2L)
  expect_equal(resolveThreads("3:2"), "3:2")
  expect_equal(resolveThreads("auto"), "4:1")
  withr::local_options(mvrem.threads = 2)
  expect_equal(resolveThreads(NULL), 2)            # the default inside jsmApply() workers
  expect_equal(resolveThreads(5), 5L)              # an explicit value wins
})

test_that("inlaArgumentList(): threads reach num.threads as given or resolved", {
  withr::local_options(mvrem.threads = NULL, mvrem.cores = 6)
  skip_if_not_installed("INLA")
  m <- fxBuild("beta")
  a <- function(th) inlaArgumentList(m$inla, list(verbose = FALSE, threads = th, inlaOptions = list()))
  expect_null(a(NULL)$num.threads)
  expect_equal(a(2)$num.threads, 2L)
  expect_equal(a("2:3")$num.threads, "2:3")
  expect_equal(a("auto")$num.threads, "6:1")
})

test_that("jsmApply(): results in order with names, budget divided, errors name the element", {
  withr::local_options(mvrem.cores = NULL)
  res <- jsmApply(list(a = 1, b = 2, c = 3), function(v, k) v * k, k = 10, cores = 3,
                  backend = "fork")
  expect_equal(res, list(a = 10, b = 20, c = 30))
  # each worker gets budget / workers threads as the default of its INLA runs
  th <- jsmApply(1:4, function(i) getOption("mvrem.threads"), cores = 8, workers = 2)
  expect_equal(unlist(th), rep(4L, 4))
  th <- jsmApply(1:2, function(i) getOption("mvrem.threads"), cores = 2)
  expect_equal(unlist(th), c(1L, 1L))
  expect_null(getOption("mvrem.threads"))          # the caller's options are untouched
  expect_equal(jsmApply(list(), identity), list())
  expect_equal(jsmApply(1:2, function(i) NULL, cores = 1), list(NULL, NULL))
  err <- tryCatch(jsmApply(list(ok = 1, bad = 2), function(v) if (v == 2) stop("boom") else v,
                           cores = 2), error = conditionMessage)
  expect_match(err, "`bad` failed: boom", fixed = TRUE)
  expect_error(jsmApply(1:2, identity, cores = 0), "cores")
  expect_error(jsmApply(1:2, identity, workers = 0), "workers")
})

test_that("jsmApply(): a worker that dies is reported", {
  skip_on_os("windows")
  suppressWarnings(expect_error(
    jsmApply(1:2, function(i) tools::pskill(Sys.getpid(), 9), cores = 2, backend = "fork"), "died"))
})

test_that("jsmApply(backend = \"psock\"): separate R sessions give the same results", {
  skip_on_cran()
  skip_if_not(requireNamespace("mvrem", quietly = TRUE), "mvrem is not installed")
  res <- jsmApply(list(a = 1, b = 2), function(v, k) v + k, k = 1, cores = 2, backend = "psock")
  expect_equal(res, list(a = 2, b = 3))
})

# FUN as a user writes it at top level: its environment is the global one
inWorkspace <- function(vars, f) {
  for (v in names(vars)) {
    val <- vars[[v]]
    if (is.function(val)) environment(val) <- globalenv()
    assign(v, val, envir = globalenv())
  }
  withr::defer(rm(list = names(vars), envir = globalenv()), envir = parent.frame())
  environment(f) <- globalenv()
  f
}

test_that("workspaceGlobals(): the workspace variables FUN uses, also through its helpers", {
  f <- inWorkspace(list(.jaData = 10, .jaHelper = function(v) v * .jaScale, .jaScale = 2,
                        .jaUnused = 1),
                   function(i) .jaHelper(i) + .jaData + betaFamily()$hyper$phi)
  g <- workspaceGlobals(f)
  expect_setequal(names(g), c(".jaData", ".jaHelper", ".jaScale"))   # not .jaUnused, not betaFamily
  # a closure (local environment) travels with FUN; the workspace variables it uses do not
  # (as in a function called from the console: its environment's parent is the workspace)
  fromConsole <- eval(quote(function() {
    k <- 3
    h <- function(i) i * k + .jaData
    function(i) h(i)
  }), globalenv())
  expect_equal(names(workspaceGlobals(fromConsole())), ".jaData")
  # recursion terminates; package functions are not copied
  rec <- inWorkspace(list(.jaRec = NULL), function(n) if (n > 0) .jaRec(n - 1) else 0)
  assign(".jaRec", rec, envir = globalenv())
  expect_equal(names(workspaceGlobals(rec)), ".jaRec")
  expect_length(workspaceGlobals(jsm), 0L)
})

test_that("jsmApply(backend = \"psock\"): FUN finds mvrem and the workspace as in the caller", {
  skip_on_cran()
  skip_if_not(requireNamespace("mvrem", quietly = TRUE), "mvrem is not installed")
  # betaFamily(): exported by every installed version of mvrem
  f <- inWorkspace(list(.jaOffset = 10, .jaAdd = function(v) v + .jaOffset),
                   function(i) list(.jaAdd(i), inherits(betaFamily(), "mvremFamily")))
  res <- jsmApply(1:2, f, cores = 2, backend = "psock")
  expect_equal(res, list(list(11, TRUE), list(12, TRUE)))
})

test_that("jsmApply(): parallel INLA fits equal the sequential ones", {
  skip_if_not_installed("INLA")
  skip_on_os("windows")
  withr::local_options(mvrem.cores = NULL)
  # a well-identified model (40 subjects): the 4-subject fxSmallData() gave intercepts that
  # differed by > 100% from run to run, serial or not (INLA's optimum is not unique there)
  d <- fxData()$long
  fitOne <- function(i) {
    m <- jsm(list(yg ~ t + scaled(x1), x1 ~ t + (1 + t | id)), data = d, threads = 1)
    fixedEffects(m)[, c("response", "term", "mean")]
  }
  seq <- lapply(1:2, fitOne)
  par <- jsmApply(1:2, fitOne, cores = 2, backend = "fork")
  # INLA's own run-to-run differences are ~4e-4 (relative)
  expect_equal(par[[1]]$mean, seq[[1]]$mean, tolerance = 0.01)
  expect_equal(par[[2]]$mean, par[[1]]$mean, tolerance = 0.01)
  expect_equal(seq[[2]]$mean, seq[[1]]$mean, tolerance = 0.01)
})
