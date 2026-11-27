# E37 benchmark: threads inside one fit, and several fits sequential vs jsmApply().
# Rscript tools/benchParallel.R <label> [nSubjects=300] [reps=2] [nSim=0]   (needs the INSTALLED mvrem)
suppressMessages(library(mvrem))
args <- commandArgs(TRUE)
label <- if (length(args)) args[1] else "machine"
n <- if (length(args) > 1) as.integer(args[2]) else 300L
reps <- if (length(args) > 2) as.integer(args[3]) else 2L
cat("machine:", label, " R", R.version.string, " INLA", as.character(packageVersion("INLA")),
    " detectCores:", parallel::detectCores(), "/", parallel::detectCores(logical = FALSE),
    " availableThreads:", availableThreads(), "(", attr(availableThreads(), "source"), ")\n")

set.seed(1)
nv <- 5
d <- data.frame(id = rep(sprintf("S%04d", 1:n), each = nv), t = rep(seq(0, 2, length.out = nv), n))
b0 <- rnorm(n, 0, 0.7); b1 <- rnorm(n, 0, 0.3)
i <- match(d$id, sprintf("S%04d", 1:n))
d$x <- 1 + 0.5 * d$t + b0[i] + b1[i] * d$t + rnorm(nrow(d), sd = 0.3)
d$y <- plogis(-0.5 + 0.4 * (1 + 0.5 * d$t + b0[i] + b1[i] * d$t) + rnorm(n, 0, 0.3)[i])
forms <- list(y ~ 1 + t + scaled(x), x ~ 1 + t + (1 + t | id))

secs <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }

# 1. threads inside one fit
budget <- as.integer(availableThreads())
ths <- unique(c(list(NULL, 1, 2, 4), if (budget >= 6) list(6), if (budget >= 8) list(8),
                list("auto", "2:2", "4:2", paste0(budget %/% 2, ":2"))))
one <- do.call(rbind, lapply(ths, function(th) {
  s <- sapply(seq_len(reps), function(r) secs(
    jsm(forms, data = d, family = list(y = betaFamily()), threads = th)))
  data.frame(threads = if (is.null(th)) "NULL" else format(th), median_s = round(median(s), 1),
             min_s = round(min(s), 1))
}))
cat("\n## one fit,", n, "subjects,", reps, "reps\n"); print(one, row.names = FALSE)

# 2. four candidate models: sequential vs jsmApply()
candidates <- list(
  "intercept only" = x ~ 1 + t + (1 | id),
  "correlated"     = x ~ 1 + t + (1 + t | id),
  "independent"    = x ~ 1 + t + (1 + t || id),
  "quadratic"      = x ~ 1 + t + I(t^2) + (1 + t | id))
# data and formulae passed through `...` (works with psock before the globals fix)
fitDic <- function(fx, d, th = NULL) {
  m <- mvrem::jsm(list(y ~ 1 + t + scaled(x), fx), data = d, family = list(y = mvrem::betaFamily()),
                  criteria = TRUE, threads = th)
  m$fit$dic$dic
}
rows <- list()
rows[["sequential, threads auto"]] <- secs(lapply(candidates, fitDic, d = d, th = "auto"))
rows[["sequential, threads NULL"]] <- secs(lapply(candidates, fitDic, d = d))
for (b in c(if (.Platform$OS.type != "windows") "fork", "psock")) {
  rows[[paste0("jsmApply ", b, ", 4 workers")]] <- secs(jsmApply(candidates, fitDic, d = d, backend = b))
  rows[[paste0("jsmApply ", b, ", 2 workers")]] <- secs(jsmApply(candidates, fitDic, d = d, backend = b, workers = 2))
}
if (budget > 6 && Sys.info()[["sysname"]] == "Darwin")
  rows[["jsmApply fork, 4 workers, cores = 6 (P cores)"]] <- secs(jsmApply(candidates, fitDic, d = d, backend = "fork", cores = 6))
cat("\n## four models\n")
print(data.frame(run = names(rows), seconds = round(unlist(rows), 1)), row.names = FALSE)
cat("\npsock start-up, 4 workers, trivial FUN:", round(secs(jsmApply(1:4, function(i) i, backend = "psock")), 2), "s\n")

# 3. a simulation study: nSim replicates of the correlated model, sequential vs jsmApply(fork)
nSim <- if (length(args) > 3) as.integer(args[4]) else 0L
if (nSim > 0) {
  simFit <- function(seed, n) {
    set.seed(seed); nv <- 5
    d <- data.frame(id = rep(sprintf("S%04d", 1:n), each = nv), t = rep(seq(0, 2, length.out = nv), n))
    b0 <- rnorm(n, 0, 0.7); b1 <- rnorm(n, 0, 0.3); i <- match(d$id, sprintf("S%04d", 1:n))
    d$x <- 1 + 0.5 * d$t + b0[i] + b1[i] * d$t + rnorm(nrow(d), sd = 0.3)
    d$y <- plogis(-0.5 + 0.4 * (1 + 0.5 * d$t + b0[i] + b1[i] * d$t) + rnorm(n, 0, 0.3)[i])
    m <- suppressMessages(mvrem::jsm(list(y ~ 1 + t + scaled(x), x ~ 1 + t + (1 + t | id)), data = d,
                                     family = list(y = mvrem::betaFamily())))
    mvrem::associations(m)$mean
  }
  sims <- list()
  sims[["sequential, threads NULL"]] <- secs(lapply(1:nSim, simFit, n = n))
  for (w in unique(pmin(nSim, c(2L, 4L, budget %/% 2L, budget))))
    sims[[sprintf("jsmApply fork, %d workers x %d threads", w, max(1L, budget %/% w))]] <-
      secs(jsmApply(1:nSim, simFit, n = n, backend = "fork", workers = w))
  cat("\n## simulation study,", nSim, "replicates of", n, "subjects\n")
  print(data.frame(run = names(sims), seconds = round(unlist(sims), 1)), row.names = FALSE)
}
