# Simulation behind the data-scaled scalar priors of single-effect and `||`
# random-effect blocks (E15, 2026-11-19): pcSd(u = "data"), halfNormalSd(scale =
# "data"), halfCauchySd(scale = "data"), gammaPrec(rate = "data").
#
# Setting "small": 40 subjects x 4 visits (t = 0, 0.5, 1, 1.5),
#   x = 1 + 0.5 t + b0 + b1 t + e,  b0 ~ N(0, 0.2^2), b1 ~ N(0, 0.1^2), e ~ N(0, 0.2^2)
#   y = 0.5 + 0.8 eta_x + e_y,      e_y ~ N(0, 0.3^2)          (gamma = 0.8)
#   model: x ~ 1 + t + re(1 + t || id, prior = P);  y ~ 1 + scaled(x)
# Setting "div10": "small" with x in other units (x / 10: SDs 0.02 / 0.01, residual
#   SD 0.02, gamma 8): fixed priors are informative in the data's units, data-scaled
#   ones are not.
# Setting "x100": x * 100. Every prior fails here (residual SD 7 x too large):
#   the fixed diffuse precision exp(-6) of the zero-equation field u (SD 20) is not
#   diffuse at this scale (see CLAUDE.md, Open issues), so this setting tests mvrem,
#   not the priors.
# Setting "large": SDs 1 (intercept) and 0.5 (slope): the data-scaled rules must
#   not hurt when the SDs are large.
#
# Priors (each effect of the block): INLA's default (loggamma(1, 5e-5)), fixed
# choices, and candidates of the data-scaled rules with s = the lme4 SD of the
# effect, floored at f x residual SD / rms(covariate of the effect):
# pc u = c s with P(sigma > u) = 0.05, half-normal / half-Cauchy scale = c s,
# gamma rate = shape s^2 (E[tau] = 1/s^2). Candidates are numeric per-effect
# priors in two re() blocks (the same model as `||`); the "API" rows use the
# exported "data" priors in one `||` block (they must equal the chosen candidate).
#
# Run (dev container, 8 cores: ~4 min per 8 seeds and setting):
#   Rscript tools/simDataPriors.R [seeds = 5] [setting = small]
# (SIM_NO_API=1 drops the API rows.) Prints per prior: mean posterior mean of each
# SD and of gamma, relative bias (%), and how many seeds' 95% intervals cover the
# true value (for seeds 1-5 and for all seeds).

args <- commandArgs(trailingOnly = TRUE)
nSeeds <- if (length(args) >= 1L) as.integer(args[1L]) else 5L
setting <- if (length(args) >= 2L) args[2L] else "small"
suppressMessages(devtools::load_all(quiet = TRUE))

truth <- switch(setting,
  small = c(sd0 = 0.2, sd1 = 0.1, sdRes = 0.2, gamma = 0.8, unit = 1),
  large = c(sd0 = 1, sd1 = 0.5, sdRes = 0.2, gamma = 0.8, unit = 1),
  div10 = c(sd0 = 0.02, sd1 = 0.01, sdRes = 0.02, gamma = 8, unit = 0.1),
  x100 = c(sd0 = 20, sd1 = 10, sdRes = 20, gamma = 0.008, unit = 100))

simData <- function(seed, n = 40) {
  set.seed(seed)
  ids <- sprintf("S%02d", seq_len(n))
  d <- data.frame(id = rep(ids, each = 4), t = rep(c(0, 0.5, 1, 1.5), n), stringsAsFactors = FALSE)
  i <- match(d$id, ids)
  b0 <- rnorm(n, sd = truth[["sd0"]]); b1 <- rnorm(n, sd = truth[["sd1"]])
  eta <- truth[["unit"]] * (1 + 0.5 * d$t) + b0[i] + b1[i] * d$t
  d$x <- eta + rnorm(nrow(d), sd = truth[["sdRes"]])
  d$y <- 0.5 + truth[["gamma"]] * eta + rnorm(nrow(d), sd = 0.3)
  d
}

# lme4 SDs of the two effects and the residual SD; the floor of an effect is
# f x residual SD / rms(its covariate) (intercept: 1)
lme4Sds <- function(d) {
  f <- suppressMessages(lme4::lmer(x ~ t + (1 + t || id), d))
  vc <- as.data.frame(lme4::VarCorr(f))
  s <- c(vc$sdcor[vc$var1 %in% "(Intercept)"], vc$sdcor[vc$var1 %in% "t"])
  c(s0 = s[1], s1 = s[2], res = stats::sigma(f), rmsT = sqrt(mean(d$t^2)), singular = lme4::isSingular(f))
}
floored <- function(l, f) pmax(l[c("s0", "s1")], f * l[["res"]] / c(1, l[["rmsT"]]))

# name -> a prior for the whole `||` block (NULL = INLA's default), or a
# function(l) of the lme4 results giving one prior per effect (candidates of the
# data-scaled rules; "f" = the floor factor, 0.5 unless given)
cand <- function(make, f) function(l) lapply(floored(l, f), make)
configs <- list(
  "INLA default" = NULL,
  "pcSd(1, 0.01)" = pcSd(1, 0.01),
  "pcSd(0.5, 0.05)" = pcSd(0.5, 0.05),
  "halfNormalSd(1)" = halfNormalSd(1),
  "halfCauchySd(1)" = halfCauchySd(1),
  "gammaPrec(1, 0.01)" = gammaPrec(1, 0.01),
  "gammaPrec(0.001, 0.001)" = gammaPrec(0.001, 0.001),
  "data pc u=3s, f=0.1" = cand(function(z) pcSd(3 * z, 0.05), 0.1),
  "data pc u=6s, f=0.1" = cand(function(z) pcSd(6 * z, 0.05), 0.1),
  "data pc u=3s" = cand(function(z) pcSd(3 * z, 0.05), 0.5),
  "data pc u=6s" = cand(function(z) pcSd(6 * z, 0.05), 0.5),
  "data pc u=10s" = cand(function(z) pcSd(10 * z, 0.05), 0.5),
  "data halfNormal 2s, f=0.1" = cand(function(z) halfNormalSd(2 * z), 0.1),
  "data halfNormal 4s, f=0.1" = cand(function(z) halfNormalSd(4 * z), 0.1),
  "data halfNormal 2s" = cand(function(z) halfNormalSd(2 * z), 0.5),
  "data halfNormal 4s" = cand(function(z) halfNormalSd(4 * z), 0.5),
  "data halfNormal 8s" = cand(function(z) halfNormalSd(8 * z), 0.5),
  "data halfCauchy 2s, f=0.1" = cand(function(z) halfCauchySd(2 * z), 0.1),
  "data halfCauchy 1s" = cand(function(z) halfCauchySd(z), 0.5),
  "data halfCauchy 2s" = cand(function(z) halfCauchySd(2 * z), 0.5),
  "data halfCauchy 4s" = cand(function(z) halfCauchySd(4 * z), 0.5),
  "data gamma shape 0.5, f=0.1" = cand(function(z) gammaPrec(0.5, 0.5 * z^2), 0.1),
  "data gamma shape 0.25" = cand(function(z) gammaPrec(0.25, 0.25 * z^2), 0.5),
  "data gamma shape 0.5" = cand(function(z) gammaPrec(0.5, 0.5 * z^2), 0.5),
  "data gamma shape 1" = cand(function(z) gammaPrec(1, z^2), 0.5),
  "API pcSd(u = \"data\")" = pcSd(u = "data"),
  "API halfNormalSd(\"data\")" = halfNormalSd(scale = "data"),
  "API halfCauchySd(\"data\")" = halfCauchySd(scale = "data"),
  "API gammaPrec(rate = \"data\")" = gammaPrec(rate = "data"))
if (Sys.getenv("SIM_NO_API") != "") configs <- configs[!startsWith(names(configs), "API")]

fitOne <- function(d, cfg, l) {
  if (is.function(cfg)) {
    pr <- cfg(l)
    long <- list(x ~ 1 + t + re(1 | id, prior = pr[[1]]) + re(0 + t | id, prior = pr[[2]]))
  } else {
    pr <- cfg
    long <- if (is.null(pr)) list(x ~ 1 + t + (1 + t || id)) else list(x ~ 1 + t + re(1 + t || id, prior = pr))
  }
  environment(long[[1]]) <- environment()
  m <- suppressWarnings(jsm(list(y ~ 1 + scaled(x)), long, data = d, threads = "1:1"))
  if (is.null(m$fit)) return(NULL)
  vc <- varianceComponents(m, "x")$table
  a <- associations(m)
  row <- function(eff) vc[vc$parameter == "sd" & vc$effect == eff, c("mean", "q0.025", "q0.975")][1, ]
  out <- rbind(sd0 = row("(Intercept)"), sd1 = row("t"), sdRes = row("residual"),
               gamma = a[1, c("mean", "q0.025", "q0.975")])
  out$par <- rownames(out)
  out
}

runSeed <- function(seed) {
  d <- simData(seed)
  l <- lme4Sds(d)
  res <- lapply(names(configs), function(nm) {
    r <- tryCatch(fitOne(d, configs[[nm]], l), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r$config <- nm; r$seed <- seed
    r
  })
  list(fits = do.call(rbind, res), lme4 = c(seed = seed, l))
}

t0 <- Sys.time()
out <- parallel::mclapply(seq_len(nSeeds), runSeed, mc.cores = min(nSeeds, parallel::detectCores()))
res <- do.call(rbind, lapply(out, `[[`, "fits"))
lme <- do.call(rbind, lapply(out, `[[`, "lme4"))
res$true <- truth[res$par]
res$cover <- res$q0.025 <= res$true & res$true <= res$q0.975

summarise <- function(res) {
  summ <- do.call(rbind, lapply(split(res, list(res$config, res$par), drop = TRUE), function(z)
    data.frame(config = z$config[1], par = z$par[1], mean = mean(z$mean),
               bias = 100 * (mean(z$mean) / z$true[1] - 1), cover = sum(z$cover), n = nrow(z))))
  fmt <- function(p) {
    z <- summ[summ$par == p, ]
    stats::setNames(z[match(names(configs), z$config), c("mean", "bias", "cover")],
                    paste0(p, c(".mean", ".bias%", ".cov")))
  }
  data.frame(prior = names(configs), fmt("sd0"), fmt("sd1"), fmt("gamma"), fmt("sdRes"), check.names = FALSE)
}
cat("Setting ", setting, ": true SDs ", truth[["sd0"]], " / ", truth[["sd1"]], ", residual SD ",
    truth[["sdRes"]], ", gamma ", truth[["gamma"]], "; coverage = seeds whose 95% ",
    "interval covers the truth\n\n", sep = "")
cat("lme4 SDs per seed (s0, s1, residual SD):\n")
print(round(lme, 3))
for (n in unique(c(min(5L, nSeeds), nSeeds))) {
  cat("\nSeeds 1 to ", n, ":\n", sep = "")
  print(format(summarise(res[res$seed <= n, ]), digits = 3), row.names = FALSE)
}
cat("\nElapsed:", format(Sys.time() - t0), "\n")
saveRDS(list(res = res, lme4 = lme), file.path(tempdir(), paste0("simDataPriors-", setting, ".rds")))
