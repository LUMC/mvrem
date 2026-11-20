# Fixtures of test-predictSample.R (predict(method = "sample"), E12), cached in
# `.fx` of helper-fixtures.R like the other fitted fixtures (read-only for the tests).
#
#   fxModelConfig(name)          the fixture `name` fitted with control.compute$config = TRUE
#                                (and criteria = TRUE, as fxModel())
#   fxPredictRefit(name, type)   predict(method = "refit") of fxModelConfig(name) on newdata
#   fxPredictSample(name, type)  predict(method = "sample", nsamples = fxNsamples, seed = 1)
#   fxPredictNewdata(name)       newdata of the comparisons: fxNewdata for "multi" (new subjects
#                                N1-N3, S01 at new times); for "priors" the same subjects with
#                                the site of each subject
#
# The refit of "multi" with type = "mean" is fxPredict() (same model and data: the
# configurations do not change the fit).

fxNsamples <- 4000L

fxModelConfig <- function(name = "multi") {
  key <- paste0(name, "_fit_config")
  testthat::skip_if_not_installed("INLA")
  testthat::skip_if_not_installed("Matrix")
  if (is.null(.fx[[key]])) {
    s <- fxSpec(name)
    m <- suppressWarnings(suppressMessages(do.call(jsm, c(s, list(
      criteria = TRUE, inlaOptions = list(control.compute = list(config = TRUE)))))))
    .fx[[key]] <- if (is.null(inlaFit(m))) "failed" else m
  }
  if (identical(.fx[[key]], "failed")) testthat::skip(paste("INLA fit of fixture", name, "with config failed"))
  .fx[[key]]
}

fxPredictNewdata <- function(name = "multi") {
  nd <- fxNewdata
  if (name == "priors") nd$site <- c("A", "B", "C", "A", "A")   # S01 is at site A
  nd
}

fxPredictRefit <- function(name = "multi", type = "mean") {
  if (name == "multi" && type == "mean") return(fxPredict())
  key <- paste0(name, "_predict_refit_", type)
  if (is.null(.fx[[key]])) {
    t0 <- proc.time()
    p <- suppressMessages(predict(fxModelConfig(name), fxPredictNewdata(name), type = type, method = "refit"))
    attr(p, "elapsed") <- (proc.time() - t0)[["elapsed"]]
    .fx[[key]] <- p
  }
  .fx[[key]]
}

fxPredictSample <- function(name = "multi", type = "mean") {
  key <- paste0(name, "_predict_sample_", type)
  if (is.null(.fx[[key]])) {
    m <- fxModelConfig(name)
    t0 <- proc.time()
    p <- predict(m, fxPredictNewdata(name), type = type, method = "sample", nsamples = fxNsamples, seed = 1)
    attr(p, "elapsed") <- (proc.time() - t0)[["elapsed"]]
    .fx[[key]] <- p
  }
  .fx[[key]]
}

# Standardized differences of the sample method from the refit method: the difference of
# the means in Monte-Carlo standard errors (refit SD / sqrt(nsamples)), the SD ratio and
# the differences of the 2.5% / 97.5% quantiles in refit SDs.
fxPredictDiff <- function(refit, sample, n = fxNsamples) {
  data.frame(refit[c("response", "id", "t")],
             z = (sample$mean - refit$mean) / (refit$sd / sqrt(n)),
             sdRatio = sample$sd / refit$sd,
             dLo = (sample$q0.025 - refit$q0.025) / refit$sd,
             dHi = (sample$q0.975 - refit$q0.975) / refit$sd)
}
