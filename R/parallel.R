# Parallel computation (E37).
#
# Two levels, both with a core budget that is computed once and divided, never left to INLA's
# default (which is all cores, so concurrent fits oversubscribe the machine):
#   1. inside one INLA run: `threads` of jsm() / refit() / predict() = INLA's `num.threads`
#      (outer threads over the hyperparameter integration points : inner solver threads);
#   2. across independent fits: jsmApply().
# INLA is an external program, so several fits can run from parallel R processes.

# `threads` of jsm() / refit() / predict(): NULL, a whole number >= 1 (outer threads, inner 1:
# INLA's reading of a single integer), "auto", or "A:B" (outer : inner, whole numbers >= 0;
# 0 = INLA chooses).
checkThreads <- function(threads, arg = "threads") {
  if (is.null(threads)) return(invisible(NULL))
  if (isString(threads)) {
    if (threads == "auto" || grepl("^[0-9]+:[0-9]+$", threads)) return(invisible(threads))
  } else if (isNumber(threads) && threads == round(threads) && threads >= 1) {
    return(invisible(as.integer(threads)))
  }
  stop("`", arg, "` must be NULL, a whole number >= 1, \"auto\" or a string \"A:B\" ",
       "(outer threads : inner threads).", call. = FALSE)
}

# The `num.threads` value INLA gets: NULL (INLA's default), a number or an "A:B" string.
# NULL falls back to `options(mvrem.threads = )`, which jsmApply() sets in its workers.
resolveThreads <- function(threads) {
  if (is.null(threads)) threads <- getOption("mvrem.threads")
  checkThreads(threads)
  if (is.null(threads)) return(NULL)
  if (identical(threads, "auto")) return(paste0(availableThreads(), ":1"))
  threads
}

#' Number of cores mvrem may use
#'
#' The smallest of the machine's physical cores (on Apple Silicon only the performance cores:
#' macOS counts the efficiency cores as physical cores too, and INLA runs gain little from them)
#' and every limit that is set: the container's
#' CPU quota (Linux cgroups v2 and v1), the scheduler variables `SLURM_CPUS_PER_TASK` and
#' `NSLOTS`, `MC_CORES`, `OMP_THREAD_LIMIT`, and 2 during `R CMD check`
#' (`_R_CHECK_LIMIT_CORES_`). [parallel::detectCores()] alone ignores all of these. A value
#' set with `options(mvrem.cores = )` overrides the detection. This is the budget of
#' `jsm(threads = "auto")` and [jsmApply()].
#'
#' @param env Environment variables, a named character vector (default: the process's).
#' @param cgroup Directory of the cgroup files (default: `/sys/fs/cgroup`).
#' @return A whole number >= 1, with attribute `source` naming the limit that decided.
#' @examples
#' availableThreads()
#' @export
availableThreads <- function(env = Sys.getenv(), cgroup = "/sys/fs/cgroup") {
  set <- getOption("mvrem.cores")
  if (!is.null(set)) {
    checkCount(set, "options(mvrem.cores)", 1L)
    return(structure(as.integer(set), source = "options(mvrem.cores)"))
  }
  env <- stats::setNames(as.character(env), names(env))
  val <- function(v) if (v %in% names(env)) unname(env[[v]]) else NA_character_
  limits <- c(machineCores(), cgroupLimit(cgroup))
  for (v in c("SLURM_CPUS_PER_TASK", "NSLOTS", "MC_CORES", "OMP_THREAD_LIMIT")) {
    n <- suppressWarnings(as.integer(val(v)))
    if (!is.na(n) && n >= 1L) limits[v] <- n
  }
  chk <- val("_R_CHECK_LIMIT_CORES_")
  if (!is.na(chk) && nzchar(chk) && tolower(chk) != "false")
    limits["_R_CHECK_LIMIT_CORES_"] <- 2
  # on a tie the most specific limit is named (scheduler, cgroup, R CMD check over the cores)
  best <- max(which(limits == min(limits)))
  structure(max(1L, as.integer(limits[[best]])), source = names(limits)[best])
}

# The machine's cores as a named number: `performance cores` on Apple Silicon, else `cores` (physical,
# or logical when the physical count is unknown).
machineCores <- function() {
  perf <- macPerformanceCores()
  if (!is.na(perf)) return(c(`performance cores` = perf))
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(cores)) cores <- suppressWarnings(parallel::detectCores())
  c(cores = if (is.na(cores)) 1L else as.integer(cores))
}

# The performance cores of a Mac with two kinds of cores (Apple Silicon: `hw.nperflevels` = 2,
# level 0 = performance); NA elsewhere (Intel Macs, other systems, sysctl not readable).
macPerformanceCores <- function() {
  if (!identical(Sys.info()[["sysname"]], "Darwin")) return(NA_integer_)
  cmd <- if (file.exists("/usr/sbin/sysctl")) "/usr/sbin/sysctl" else "sysctl"
  out <- tryCatch(suppressWarnings(system2(cmd, c("-n", "hw.nperflevels", "hw.perflevel0.physicalcpu"),
                                           stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())
  v <- suppressWarnings(as.integer(out))
  if (length(v) == 2L && !anyNA(v) && v[1] >= 2L && v[2] >= 1L) v[2] else NA_integer_
}

# The CPU quota of a container as a number of cores (cgroups v2 `cpu.max`, v1
# `cpu.cfs_quota_us` / `cpu.cfs_period_us`), rounded up; empty when there is none.
cgroupLimit <- function(dir) {
  read1 <- function(f) tryCatch(readLines(file.path(dir, f), n = 1L, warn = FALSE),
                                error = function(e) character(), warning = function(w) character())
  v2 <- strsplit(read1("cpu.max"), " ", fixed = TRUE)
  if (length(v2) && length(v2[[1]]) == 2L && v2[[1]][1] != "max") {
    q <- suppressWarnings(as.numeric(v2[[1]]))
    if (!anyNA(q) && q[2] > 0) return(c(cgroup = ceiling(q[1] / q[2])))
  }
  q <- suppressWarnings(as.numeric(read1("cpu/cpu.cfs_quota_us")))
  p <- suppressWarnings(as.numeric(read1("cpu/cpu.cfs_period_us")))
  if (length(q) == 1L && length(p) == 1L && !is.na(q) && !is.na(p) && q > 0 && p > 0)
    return(c(cgroup = ceiling(q / p)))
  c()
}

#' Run several INLA computations in parallel
#'
#' Applies `FUN` to every element of `X` in separate R processes, each of which runs its own
#' INLA program: independent fits (refits of changed models, predictions by refitting for
#' several new subjects, cross-validation folds, simulation replicates) are as many times
#' faster as there are workers. The core budget ([availableThreads()], or `cores`) is divided:
#' with `w` workers each INLA run gets `budget / w` threads (through `options(mvrem.threads = )`
#' in the worker, which `jsm()`, `refit()` and `predict()` use when their `threads` is `NULL`),
#' so that the machine is never oversubscribed. Only base R's parallel package is used.
#'
#' Backends: `"fork"` (copies of the running R process: fast to start, everything of the session
#' available; not on Windows, and not safe in a GUI such as RStudio or R.app, see
#' [parallel::mclapply()]) and `"psock"` (separate R sessions). `"auto"` chooses `"fork"` on Linux
#' and on macOS outside a GUI (R in a terminal, `Rscript`), else `"psock"` (Windows, RStudio,
#' Positron, R.app). The psock workers
#' attach the INSTALLED mvrem (a package only `load_all()`ed is not seen by them) with the
#' caller's library paths; `X`, `FUN` and `...` are copied to them, and so are the variables of
#' the workspace that `FUN` uses (also through functions of the workspace it calls), so `FUN`
#' runs as it would in the caller's session. Packages other than mvrem that `FUN` uses must be
#' loaded inside it (or called with `pkg::`).
#'
#' Memory: every worker holds a fit. A fit kept with INLA's configurations
#' (`control.compute = list(config = TRUE)`) can be hundreds of MB; return only what is needed
#' from `FUN` (a table, not the model), and give `cores` fewer workers for large fits.
#'
#' @param X A list (or vector) of inputs.
#' @param FUN A function of one element of `X` (and `...`); it may call `jsm()`, `refit()`,
#'   `predict()`, ... Its results are returned in the order of `X`.
#' @param ... Further arguments of `FUN`.
#' @param cores Total number of cores to use: `"auto"` ([availableThreads()]) or a whole
#'   number.
#' @param workers Number of parallel R processes: default `min(length(X), cores)`; fewer workers
#'   with more threads each suit few, large fits.
#' @param backend `"auto"`, `"fork"` or `"psock"`.
#' @return A list the length (and names) of `X`. An error in `FUN` stops with a message
#'   naming the element.
#' @examples
#' \dontrun{
#' # four candidate models fitted at once, only their DIC kept
#' forms <- list(...)           # a list of lists of formulae
#' dic <- jsmApply(forms, function(f) {
#'   jsm(f, data = d, criteria = TRUE)$fit$dic$dic
#' })
#' }
#' @export
jsmApply <- function(X, FUN, ..., cores = "auto", workers = NULL,
                     backend = c("auto", "fork", "psock")) {
  backend <- match.arg(backend)
  FUN <- match.fun(FUN)
  X <- as.list(X)
  if (identical(cores, "auto")) cores <- as.integer(availableThreads())
  checkCount(cores, "cores", 1L)
  cores <- as.integer(cores)
  checkCount(workers, "workers", 1L, nullOk = TRUE)
  n <- length(X)
  if (!n) return(X)
  w <- min(n, if (is.null(workers)) cores else as.integer(workers))
  per <- max(1L, cores %/% w)
  worker <- makeWorker(X, FUN, per, list(...))
  if (backend == "auto") backend <- autoBackend()
  if (w == 1L) {
    res <- lapply(seq_len(n), worker)
  } else if (backend == "fork") {
    if (.Platform$OS.type == "windows")
      stop("`backend = \"fork\"` is not available on Windows; use \"psock\".", call. = FALSE)
    res <- parallel::mclapply(seq_len(n), worker, mc.cores = w, mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(w)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    startWorkers(cl, FUN)
    res <- parallel::parLapplyLB(cl, seq_len(n), worker)
  }
  # mclapply returns a try-error (or NULL) for a worker that died, e.g. out of memory
  for (i in seq_len(n)) {
    r <- res[[i]]
    if (!is.list(r) || !any(c("value", "error") %in% names(r)))
      stop("jsmApply(): the worker of element ", elementName(X, i),
           " died (out of memory? use fewer `workers`).", call. = FALSE)
    if (!is.null(r$error))
      stop("jsmApply(): element ", elementName(X, i), " failed: ", r$error, call. = FALSE)
  }
  out <- lapply(res, `[[`, "value")
  names(out) <- names(X)
  out
}

# Prepares separate R sessions (psock) to run FUN as the caller would: the caller's library
# paths, mvrem ATTACHED (FUN may call jsm() etc. without `mvrem::`), and the variables FUN uses
# from the caller's workspace copied into each worker's global environment (a function defined
# at top level is serialised without them).
startWorkers <- function(cl, FUN) {
  ok <- parallel::clusterCall(cl, function(paths) {
    .libPaths(paths)
    suppressPackageStartupMessages(require("mvrem", quietly = TRUE, character.only = TRUE))
  }, .libPaths())
  if (!all(vapply(ok, isTRUE, logical(1))))
    stop("jsmApply(): the R sessions of `backend = \"psock\"` cannot load mvrem; install it ",
         "(e.g. devtools::install()): a package only load_all()ed is not seen by them.", call. = FALSE)
  glob <- workspaceGlobals(FUN)
  if (length(glob)) parallel::clusterExport(cl, names(glob), envir = list2env(glob))
  invisible()
}

# The variables of the caller's workspace (global environment) that FUN refers to, directly or
# through functions it calls that are defined there too: a named list. Variables of local
# environments (FUN defined inside a function) travel with FUN anyway; package namespaces are
# loaded by name. The names are those of all.names(), so a name used only as `d$x` brings `x`
# along if the workspace has one (harmless).
workspaceGlobals <- function(FUN) {
  found <- list()
  seen <- list()                     # functions already scanned (recursion, mutual calls)
  todo <- list(FUN)
  while (length(todo)) {
    f <- todo[[1L]]
    todo <- todo[-1L]
    env <- environment(f)
    if (is.null(env) || any(vapply(seen, identical, logical(1), f))) next
    seen <- c(seen, list(f))
    syms <- unique(c(all.names(body(f)), unlist(lapply(formals(f), function(a)
      if (is.language(a)) all.names(a)))))
    for (nm in setdiff(syms, names(found))) {
      e <- userHolder(nm, env)
      if (is.null(e)) next
      val <- get(nm, envir = e, inherits = FALSE)
      if (identical(e, globalenv())) found[[nm]] <- val
      # a function of the workspace or of a local environment may use workspace variables too
      if (is.function(val)) todo <- c(todo, list(val))
    }
  }
  found
}

# The environment in which FUN (environment `env`) finds `nm`, when that is a user environment
# (a local one or the global one); NULL for packages, base, or a name not found.
userHolder <- function(nm, env) {
  while (!identical(env, emptyenv())) {
    if (isNamespace(env) || identical(env, baseenv())) return(NULL)
    if (exists(nm, envir = env, inherits = FALSE)) return(env)
    if (identical(env, globalenv())) return(NULL)
    env <- parent.env(env)
  }
  NULL
}

# The backend of jsmApply(backend = "auto"): psock on Windows and, on macOS, in a GUI or IDE
# session (forking a process with a GUI event loop is unsafe, see ?parallel::mclapply); fork on
# Linux and in a terminal R or Rscript on macOS.
autoBackend <- function(os = .Platform$OS.type, sysname = Sys.info()[["sysname"]],
                        gui = .Platform$GUI, env = Sys.getenv(c("RSTUDIO", "POSITRON"))) {
  if (os == "windows") return("psock")
  if (!identical(sysname, "Darwin")) return("fork")
  ide <- identical(unname(env[["RSTUDIO"]]), "1") || nzchar(unname(env[["POSITRON"]]))
  if (ide || !gui %in% c("X11", "unknown")) "psock" else "fork"
}

elementName <- function(X, i) if (!is.null(names(X)) && nzchar(names(X)[i])) paste0("`", names(X)[i], "`") else i

# The function a worker runs for element i: FUN with the thread budget `per` as the default of
# every INLA run inside. A closure of the package namespace, so that a PSOCK worker gets
# everything it needs by serialisation.
makeWorker <- function(X, FUN, per, dots) {
  force(X); force(FUN); force(per); force(dots)
  function(i) {
    old <- options(mvrem.threads = per)
    on.exit(options(old))
    tryCatch(list(value = do.call(FUN, c(list(X[[i]]), dots))),
             error = function(e) list(error = conditionMessage(e)))
  }
}
