#' Random effects and their variance components
#'
#' Posterior summaries of the subject-level random effects and of the
#' covariance of each random-effect block, read from a fitted model with the
#' names of the model (responses, groups, terms) instead of INLA's row numbers.
#'
#' `varianceComponents()` returns standard deviations (or variances) and
#' correlations with 95% intervals computed from the posterior marginals of the
#' hyperparameters (INLA's precisions are transformed, so the intervals are not
#' those of the precisions inverted at the mean). The matrices `cov` and `cor`
#' are plug-ins built from the posterior means of the standard deviations and
#' correlations. For Gaussian (and other precision-parametrised) responses the
#' residual standard deviation is included, as block `"residual"`.
#'
#' @param x A fitted `mvremModel`, e.g. from [jsm()].
#' @param response Names of the responses to show (default: all).
#' @param scale `"sd"` (default) or `"variance"`: the scale of the table and of
#'   the marginals in `varianceComponents()`.
#' @param ... Further arguments of the methods.
#' @return
#' * `randomEffects()`: a data frame with one row per subject (group level) and
#'   effect: `response`, `group`, `term`, `level`, and the posterior `mean`,
#'   `sd`, `q0.025`, `q0.5`, `q0.975`, `mode`.
#' * `varianceComponents()`: an object of class `mvremVarCov`, a list with
#'   `table` (`response`, `group`, `block`, `parameter` = `sd` / `variance` /
#'   `correlation`, `effect`, `effect2`, `mean`, `sd`, `q0.025`, `q0.5`,
#'   `q0.975`), `cov` and `cor` (lists of named matrices, one per random-effect
#'   block), `marginals` (list parallel to the rows of `table`: posterior
#'   density of each parameter as a two-column matrix `x`, `y`) and `scale`.
#' @examples
#' \dontrun{
#' m <- jsm(list(y ~ t + scaled(x1)), list(x1 ~ t + (1 + t | id)), data = d)
#' randomEffects(m)
#' vc <- varianceComponents(m)
#' vc$cor
#' plotRandomEffects(m)
#' }
#' @name randomEffects
NULL

#' @rdname randomEffects
#' @export
randomEffects <- function(x, ...) UseMethod("randomEffects")
#' @rdname randomEffects
#' @export
randomEffects.mvremModel <- function(x, response = NULL, ...) {
  needFit(x)
  resp <- checkResponses(x, response)
  sr <- x$fit$summary.random
  nm <- x$inla$naming
  rows <- list()
  for (r in intersect(resp, names(x$inla$random))) for (b in x$inla$random[[r]]) {
    lev <- x$inla$groups[[b$group]]
    # a correlated block is ONE INLA field holding all its effects (its other indices
    # are copies); a single effect or a `||` block has one field per effect
    fields <- if (length(b$effects) > 1L) rep(b$index[1L], length(b$effects)) else b$index
    for (e in seq_along(b$effects)) {
      s <- sr[[fields[e]]]
      # the rows of this effect, in the order of the group levels: matched by the labels
      # the model gave INLA (`"x1|t|S01"`), whatever order INLA or the data use
      labels <- fillValue(nm$reValue, resp = r, term = termLabel(b$effects[e], nm), level = lev)
      k <- match(labels, as.character(s$ID))
      if (anyNA(k))
        stop("Random-effect labels of `", fields[e], "` not found in the INLA field: ",
             paste(labels[is.na(k)][1:min(3L, sum(is.na(k)))], collapse = ", "), call. = FALSE)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(response = r, group = b$group, term = b$effects[e], level = lev,
                   stringsAsFactors = FALSE),
        statTable(s[k, setdiff(names(s), "ID"), drop = FALSE]))
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(response = character(), group = character(), term = character(),
               level = character())
  rownames(out) <- NULL
  out
}

# The responses selected by a `response` argument (NULL = all), checked, in the order of
# the model (the order in which the formulae were passed), whatever the order given.
checkResponses <- function(x, response) {
  all <- modelResponses(x)
  if (is.null(response)) return(all)
  bad <- setdiff(response, all)
  if (length(bad))
    stop("Unknown response(s): ", paste(bad, collapse = ", "), ". Responses of the model: ",
         paste(all, collapse = ", "), ".", call. = FALSE)
  all[all %in% response]
}

#' @rdname randomEffects
#' @export
varianceComponents <- function(x, ...) UseMethod("varianceComponents")
#' @rdname randomEffects
#' @export
varianceComponents.mvremModel <- function(x, response = NULL, scale = c("sd", "variance"), ...) {
  needFit(x)
  requireInla()
  scale <- match.arg(scale)
  resp <- checkResponses(x, response)
  tr <- if (scale == "sd") function(p) 1 / sqrt(p) else function(p) 1 / p
  sh <- x$fit$summary.hyperpar
  mh <- x$fit$marginals.hyperpar
  sr <- x$fit$summary.random

  # posterior summary (mean, sd, quantiles) of a transformed hyperparameter + its density
  transformed <- function(name, f) {
    if (is.null(m <- mh[[name]])) return(NULL)
    tm <- tryCatch(suppressWarnings(INLA::inla.tmarginal(f, m)), error = function(e) NULL)
    if (is.null(tm) || !all(is.finite(tm[, 1L])) || !all(is.finite(tm[, 2L]))) {
      # a precision marginal spanning orders of magnitude (an SD near 0) is smoothed into
      # negative precisions by inla.tmarginal(): transform INLA's internal log precision
      lm <- x$fit$internal.marginals.hyperpar[[sub("^Precision", "Log precision", name)]]
      if (is.null(lm)) stop("Cannot transform the marginal of '", name, "'.", call. = FALSE)
      tm <- INLA::inla.tmarginal(function(th) f(exp(th)), lm)
    }
    z <- INLA::inla.zmarginal(tm, silent = TRUE)
    list(stats = c(mean = z$mean, sd = z$sd, q0.025 = z$quant0.025, q0.5 = z$quant0.5,
                   q0.975 = z$quant0.975), density = tm)
  }
  plain <- function(name) {
    if (!name %in% rownames(sh) || is.null(mh[[name]])) return(NULL)
    s <- statTable(sh[name, , drop = FALSE])
    list(stats = unlist(s[c("mean", "sd", "q0.025", "q0.5", "q0.975")]), density = mh[[name]])
  }
  rows <- list(); dens <- list()
  add <- function(r, group, block, parameter, effect, effect2, est) {
    if (is.null(est)) return(invisible())
    rows[[length(rows) + 1L]] <<- data.frame(
      response = r, group = group, block = block, parameter = parameter, effect = effect,
      effect2 = effect2, as.list(est$stats), stringsAsFactors = FALSE)
    dens[[length(dens) + 1L]] <<- est$density
  }
  covs <- list(); cors <- list()

  for (r in intersect(resp, names(x$inla$random))) for (b in x$inla$random[[r]]) {
    k <- length(b$effects)
    block <- sprintf("%s | %s [%s]", r, b$group, paste(b$effects, collapse = ", "))
    correlated <- k > 1L   # one INLA field for all effects of a correlated block
    precNames <- if (correlated) paste0("Precision for ", b$index[1L], " (component ", seq_len(k), ")") else
      paste0("Precision for ", b$index)
    sdv <- rep(NA_real_, k)
    for (j in seq_len(k)) {
      add(r, b$group, block, scale, b$effects[j], NA_character_, transformed(precNames[j], tr))
      s1 <- transformed(precNames[j], function(p) 1 / sqrt(p))
      if (!is.null(s1)) sdv[j] <- s1$stats[["mean"]]
    }
    cm <- diag(k)
    if (correlated) for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
      est <- plain(paste0("Rho", i, ":", j, " for ", b$index[1L]))
      add(r, b$group, block, "correlation", b$effects[i], b$effects[j], est)
      if (!is.null(est)) cm[i, j] <- cm[j, i] <- est$stats[["mean"]]
    }
    dimnames(cm) <- list(b$effects, b$effects)
    covs[[block]] <- outer(sdv, sdv) * cm
    dimnames(covs[[block]]) <- dimnames(cm)
    cors[[block]] <- cm
  }

  # residual / observation precision of the likelihoods
  ht <- hyperTable(x)
  for (j in which(ht$type %in% "family" & ht$parameter %in% c("residual precision", "precision") &
                    ht$response %in% resp))
    add(ht$response[j], NA_character_, "residual", scale, "residual", NA_character_,
        transformed(rownames(sh)[j], tr))

  table <- if (length(rows)) do.call(rbind, rows) else
    data.frame(response = character(), group = character(), block = character(),
               parameter = character(), effect = character(), effect2 = character())
  # the random-effect blocks, then the residual SDs, each by response in the order of the model
  if (nrow(table)) {
    k <- order(table$block == "residual", match(table$response, resp), seq_len(nrow(table)))
    table <- table[k, , drop = FALSE]
    dens <- dens[k]
  }
  rownames(table) <- NULL
  structure(list(table = table, cov = covs, cor = cors, marginals = dens, scale = scale),
            class = "mvremVarCov")
}

# varianceComponents() for the displays that show it as one part among others
# (estimatesTable(), the LaTeX lookup, the pairs plot): if it fails (INLA not installed,
# a marginal that cannot be transformed) they go on without it, with a warning that says
# what is missing instead of silently dropping the rows.
varianceComponentsOrNull <- function(x, what) {
  tryCatch(varianceComponents(x), error = function(e) {
    warning("Variance components left out of ", what, ": ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

#' @export
print.mvremVarCov <- function(x, digits = 3, ...) {
  cat("Variance components (posterior summaries; scale = ", x$scale, ")\n\n", sep = "")
  t <- x$table
  num <- vapply(t, is.numeric, NA)
  t[num] <- lapply(t[num], signif, digits)
  print(t[setdiff(names(t), c("group", "response"))], row.names = FALSE)
  for (b in names(x$cov)) {
    if (nrow(x$cov[[b]]) < 2L) next
    cat("\n", b, ": covariance (from posterior mean SDs and correlations)\n", sep = "")
    print(signif(x$cov[[b]], digits))
  }
  invisible(x)
}
