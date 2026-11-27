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
  sr <- x$fit$summary.random
  resp <- selectResponses(x, response)
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
      labels <- fillValue(nm$reValue, resp = r, term = inlaTermName(b$effects[e], nm), level = lev)
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
  bindRows(rows, data.frame(response = character(), group = character(), term = character(),
                            level = character()))
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
  resp <- selectResponses(x, response)
  tr <- if (scale == "sd") function(p) 1 / sqrt(p) else function(p) 1 / p
  sh <- x$fit$summary.hyperpar
  mh <- x$fit$marginals.hyperpar

  # posterior summary (mean, sd, quantiles) of a transformed precision + its density;
  # `internal` = the name of INLA's internal (log precision) marginal, the fallback
  transformed <- function(name, f, internal = NA) {
    if (is.null(m <- mh[[name]])) return(NULL)
    tm <- tryCatch(suppressWarnings(INLA::inla.tmarginal(f, m)), error = function(e) NULL)
    if (is.null(tm) || !all(is.finite(tm[, 1L])) || !all(is.finite(tm[, 2L]))) {
      # a precision marginal spanning orders of magnitude (an SD near 0) is smoothed into
      # negative precisions by inla.tmarginal(): transform INLA's internal log precision
      lm <- if (!is.na(internal)) x$fit$internal.marginals.hyperpar[[internal]]
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
  # `blockId` = the position of the block in `inl$random[[response]]` (NA for the residual
  # rows): the structural key of a block; `block` is its display label (D4)
  add <- function(r, group, blockId, block, parameter, effect, effect2, est) {
    if (is.null(est)) return(invisible())
    rows[[length(rows) + 1L]] <<- data.frame(
      response = r, group = group, blockId = blockId, block = block, parameter = parameter,
      effect = effect, effect2 = effect2, as.list(est$stats), stringsAsFactors = FALSE)
    dens[[length(dens) + 1L]] <<- est$density
  }
  covs <- list(); cors <- list()
  known <- latentHyperRows(x$inla)

  for (r in intersect(resp, names(x$inla$random))) for (i in seq_along(x$inla$random[[r]])) {
    b <- x$inla$random[[r]][[i]]
    k <- length(b$effects)
    block <- blockLabel(r, b$group, b$effects)
    hr <- blockHyperRows(known, r, i)
    prec <- hr[hr$short == "prec", , drop = FALSE]           # one row per effect, in order
    sdv <- rep(NA_real_, k)
    for (j in seq_len(k)) {
      add(r, b$group, i, block, scale, b$effects[j], NA_character_,
          transformed(prec$.inla[j], tr, prec$internal[j]))
      s1 <- transformed(prec$.inla[j], function(p) 1 / sqrt(p), prec$internal[j])
      if (!is.null(s1)) sdv[j] <- s1$stats[["mean"]]
    }
    cm <- diag(k)
    for (j in which(hr$short == "rho")) {                     # only a correlated block has them
      est <- plain(hr$.inla[j])
      add(r, b$group, i, block, "correlation", hr$effect[j], hr$effect2[j], est)
      e1 <- match(hr$effect[j], b$effects); e2 <- match(hr$effect2[j], b$effects)
      if (!is.null(est)) cm[e1, e2] <- cm[e2, e1] <- est$stats[["mean"]]
    }
    dimnames(cm) <- list(b$effects, b$effects)
    covs[[block]] <- outer(sdv, sdv) * cm
    dimnames(covs[[block]]) <- dimnames(cm)
    cors[[block]] <- cm
  }

  # residual / observation precision of the likelihoods
  ht <- hyperTable(x)
  for (j in which(ht$type %in% "family" & ht$short %in% "prec" & ht$response %in% resp))
    add(ht$response[j], NA_character_, NA_integer_, "residual", scale, "residual", NA_character_,
        transformed(ht$.inla[j], tr))

  table <- bindRows(rows, data.frame(response = character(), group = character(), blockId = integer(),
                                     block = character(), parameter = character(), effect = character(),
                                     effect2 = character()))
  # the random-effect blocks, then the residual SDs, each by response in the order of the model
  if (nrow(table)) {
    k <- responseOrder(x, table, first = table$block == "residual")
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

# The display label of a random-effect block in varianceComponents()$table / $cov / $cor.
blockLabel <- function(r, group, effects) sprintf("%s | %s [%s]", r, group, paste(effects, collapse = ", "))

#' @export
print.mvremVarCov <- function(x, digits = 3, ...) {
  cat("Variance components (posterior summaries; scale = ", x$scale, ")\n\n", sep = "")
  t <- x$table
  num <- vapply(t, is.numeric, NA)
  t[num] <- lapply(t[num], signif, digits)
  print(t[setdiff(names(t), c("group", "response", "blockId"))], row.names = FALSE)
  for (b in names(x$cov)) {
    if (nrow(x$cov[[b]]) < 2L) next
    cat("\n", b, ": covariance (from posterior mean SDs and correlations)\n", sep = "")
    print(signif(x$cov[[b]], digits))
  }
  invisible(x)
}
