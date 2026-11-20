#' Plots of the estimates of a fitted model
#'
#' Base-graphics summaries of the posterior results, one panel per equation or
#' random-effect term (same style as [plot.mvremModel()]):
#'
#' * `plotRandomEffects()`: caterpillar plot of the random effects: posterior
#'   mean and 95% interval of every subject, one panel per response, term and
#'   group; shows the spread of the subjects and outliers. With
#'   `type = "pairs"`, the posterior means of two effects of a block against each
#'   other (e.g. intercept and slope), one point per subject, with the estimated
#'   correlation of a correlated block: the 95% ellipse of the bivariate normal
#'   with the posterior-mean SDs and correlation (solid; wider than the cloud of
#'   posterior means, which are shrunk towards 0), the ellipses with the
#'   correlation at the bounds of its 95% interval (dashed) and the correlation
#'   with its interval as text (no ellipse for effects without an estimated
#'   correlation: `||` or separate blocks).
#' * `plotVarianceComponents()`: posterior densities of the standard deviations
#'   (or variances) of the random effects and the residuals, one panel per
#'   random-effect block, and of the correlations of the blocks with several
#'   effects (one more panel).
#' * `plotFixedEffects()`: posterior mean and 95% interval of the fixed effects,
#'   one panel per response (`type = "forest"`), or their posterior densities,
#'   one curve per term (`type = "density"`).
#' * `plotAssociations()`: the same for the association coefficients
#'   \eqn{\gamma} of the `scaled()` terms.
#'
#' By default each returns an HTML widget (when htmlwidgets is installed; see
#' [htmlPlot()]): the subject, term or curve under the mouse is identified and
#' highlighted, and the same subject is highlighted in all random-effect panels.
#'
#' @inheritParams randomEffects
#' @param term Terms to show (e.g. `"t"`, `"1"` or `"(Intercept)"`); default: all.
#' @param order `"none"` (default, data order) or `"mean"` (subjects sorted by
#'   posterior mean).
#' @param type `plotRandomEffects()`: `"caterpillar"` (default) or `"pairs"`;
#'   `plotFixedEffects()` / `plotAssociations()`: `"forest"` (default) or
#'   `"density"` (see above).
#' @param ncol Number of panel columns (default: automatic).
#' @param interactive `NULL` (default): an HTML widget with mouse-over
#'   identification if the htmlwidgets package is installed, else base graphics;
#'   `TRUE` / `FALSE` force one of them (see [htmlPlot()]).
#' @return The widget, or `x` invisibly for base graphics.
#' @examples
#' \dontrun{
#' plotRandomEffects(m, "x1", order = "mean")
#' plotRandomEffects(m, "x1", type = "pairs")
#' plotRandomEffects(m, "x1", interactive = FALSE)   # base graphics
#' plotVarianceComponents(m)
#' plotFixedEffects(m)
#' plotFixedEffects(m, type = "density")
#' plotAssociations(m)
#' }
#' @name plotResults
NULL

#' @rdname plotResults
#' @export
plotRandomEffects <- function(x, ...) UseMethod("plotRandomEffects")
#' @rdname plotResults
#' @export
plotRandomEffects.mvremModel <- function(x, response = NULL, term = NULL,
                                         order = c("none", "mean"), type = c("caterpillar", "pairs"),
                                         ncol = NULL, interactive = NULL, ...) {
  order <- match.arg(order)
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  re <- randomEffects(x, response)
  re$label <- effectLabel(re$term, x$inla$naming)
  if (!is.null(term))   # `term`: as in the formula, as labelled, or "1" for the intercept
    re <- re[re$term %in% term | re$label %in% term | (re$term == "(Intercept)" & "1" %in% term), , drop = FALSE]
  if (!nrow(re)) stop("No random effects to plot.", call. = FALSE)
  if (type == "pairs") return(plotRandomPairs(x, re, ncol, interactive))

  re$panel <- paste0(re$response, ": ", re$label, " | ", re$group)
  panels <- unique(re$panel)
  rows <- lapply(panels, function(p) {          # the rows of each panel, in the order drawn
    d <- re[re$panel == p, ]
    if (order == "mean") d <- d[base::order(d$mean), ]
    d
  })
  if (isTRUE(interactive)) {
    ps <- Map(function(p, d) htmlIntervals(p, d$level, d$level, d$mean, d$q0.025, d$q0.975, ylab = "subject"),
              panels, rows)
    return(htmlPlot(unname(ps), autoNcol(length(ps), ncol), linked = TRUE, search = "subject id"))
  }
  op <- panelLayout(length(panels), ncol)
  on.exit(graphics::par(op), add = TRUE)
  for (i in seq_along(panels)) {
    p <- panels[i]
    d <- rows[[i]]
    n <- nrow(d)
    graphics::plot(d$mean, seq_len(n), type = "n", ylim = c(n, 1),
                   xlim = range(c(d$q0.025, d$q0.975), finite = TRUE),
                   xlab = "posterior mean and 95% interval", ylab = "subject", main = p, cex.main = 0.95)
    graphics::abline(v = 0, lty = 3, col = plotColors$base$reference)
    graphics::segments(d$q0.025, seq_len(n), d$q0.975, seq_len(n), col = plotColors$base$interval)
    graphics::points(d$mean, seq_len(n), pch = 16, cex = 0.45, col = plotColors$accent)
  }
  invisible(x)
}

# Posterior means of two effects of one block against each other, a point per subject, with
# the 95% ellipse of the estimated distribution of the effects (posterior-mean SDs and
# correlation) and the ellipses with the correlation at the bounds of its 95% interval.
plotRandomPairs <- function(x, re, ncol, interactive) {
  vc <- varianceComponentsOrNull(x, "the pairs plot (no ellipses)")$table
  pairs <- list()
  # one key per (response, group); "\r" does not occur in variable names, so the pairs stay apart
  for (key in unique(paste(re$response, re$group, sep = "\r"))) {
    d <- re[paste(re$response, re$group, sep = "\r") == key, ]
    terms <- unique(d$label)
    if (length(terms) < 2L) next
    for (i in seq_len(length(terms) - 1L)) for (j in (i + 1L):length(terms)) {
      a <- d[d$label == terms[i], ]
      b <- d[d$label == terms[j], ]
      lev <- intersect(a$level, b$level)
      pairs[[length(pairs) + 1L]] <- c(list(
        title = paste0(d$response[1L], ": ", terms[i], " + ", terms[j], " | ", d$group[1L]), level = lev,
        x = a$mean[match(lev, a$level)], y = b$mean[match(lev, b$level)],
        xlab = paste0("posterior mean of ", terms[i], " | ", d$group[1L]),
        ylab = paste0("posterior mean of ", terms[j], " | ", d$group[1L])),
        pairEllipses(vc, d$response[1L], d$group[1L], terms[i], terms[j], x$inla$naming))
    }
  }
  if (!length(pairs)) stop("`type = \"pairs\"` needs a block with at least two effects.", call. = FALSE)
  ellipseCol <- plotColors$ellipse
  haveEllipse <- any(vapply(pairs, function(p) !is.null(p$ellipses), logical(1)))
  if (isTRUE(interactive)) {
    ps <- lapply(pairs, function(p) {
      pn <- htmlPoints(p$title, p$level, p$x, p$y, p$xlab, p$ylab)
      if (!is.null(p$ellipses)) {
        e <- p$ellipses
        pn$curves <- list(
          list(kind = "ellipse", color = ellipseCol, dash = "", width = 1.8, x = r5(e$mean$x), y = r5(e$mean$y)),
          list(kind = "ellipseBound", color = ellipseCol, dash = "4 3", width = 1.2, x = r5(e$lo$x), y = r5(e$lo$y)),
          list(kind = "ellipseBound", color = ellipseCol, dash = "4 3", width = 1.2, x = r5(e$hi$x), y = r5(e$hi$y)))
        pn$note <- p$note
      }
      pn
    })
    legend <- if (haveEllipse)
      list(list(kind = "ellipse", label = "95% ellipse (posterior-mean SDs and correlation)",
                color = ellipseCol, dash = ""),
           list(kind = "ellipseBound", label = "correlation at its 95% bounds", color = ellipseCol, dash = "4 3"))
    return(htmlPlot(ps, autoNcol(length(ps), ncol), legend, linked = TRUE, search = "subject id"))
  }
  op <- panelLayout(length(pairs), ncol)
  on.exit(graphics::par(op), add = TRUE)
  for (p in pairs) {
    e <- p$ellipses
    xlim <- range(c(p$x, e$mean$x, e$lo$x, e$hi$x), finite = TRUE)
    ylim <- range(c(p$y, e$mean$y, e$lo$y, e$hi$y), finite = TRUE)
    graphics::plot(p$x, p$y, pch = 16, cex = 0.6, col = plotColors$accent, xlim = xlim, ylim = ylim,
                   xlab = p$xlab, ylab = p$ylab, main = p$title, cex.main = 0.95)
    graphics::abline(h = 0, v = 0, lty = 3, col = plotColors$base$reference)
    if (!is.null(e)) {
      graphics::lines(e$lo$x, e$lo$y, col = ellipseCol, lty = 2)
      graphics::lines(e$hi$x, e$hi$y, col = ellipseCol, lty = 2)
      graphics::lines(e$mean$x, e$mean$y, col = ellipseCol, lwd = 1.8)
      graphics::legend("topleft", legend = p$note, bty = "n", cex = 0.85, text.col = "grey20")
    }
  }
  invisible(x)
}

# The ellipses of a pair of effects (labels as in randomEffects(): "1", "t") of a correlated
# block, from the varianceComponents() table: NULL when the block has no estimated
# correlation of the two (separate blocks, `||`).
pairEllipses <- function(vc, response, group, ti, tj, naming = namingScheme()) {
  if (is.null(vc)) return(NULL)
  lab <- effectLabel(vc$effect, naming)
  lab2 <- effectLabel(vc$effect2, naming)
  inBlock <- vc$response == response & vc$group %in% group
  sd1 <- vc[inBlock & vc$parameter == "sd" & lab == ti, ]
  sd2 <- vc[inBlock & vc$parameter == "sd" & lab == tj, ]
  cr <- vc[inBlock & vc$parameter == "correlation" &
             ((lab == ti & lab2 == tj) | (lab == tj & lab2 == ti)), ]
  if (nrow(sd1) != 1L || nrow(sd2) != 1L || nrow(cr) != 1L) return(NULL)
  ell <- function(rho) ellipsePoints(sd1$mean, sd2$mean, rho)
  list(ellipses = list(mean = ell(cr$mean), lo = ell(cr$q0.025), hi = ell(cr$q0.975)),
       note = sprintf("correlation %s [%s, %s]", signif(cr$mean, 2), signif(cr$q0.025, 2), signif(cr$q0.975, 2)))
}

# 95% probability ellipse of a bivariate normal with SDs s1, s2 and correlation rho (centred at 0)
ellipsePoints <- function(s1, s2, rho, n = 120L) {
  th <- seq(0, 2 * pi, length.out = n)
  r <- sqrt(stats::qchisq(0.95, 2))
  list(x = r * s1 * cos(th), y = r * s2 * (rho * cos(th) + sqrt(max(0, 1 - rho^2)) * sin(th)))
}

#' @rdname plotResults
#' @export
plotVarianceComponents <- function(x, ...) UseMethod("plotVarianceComponents")
#' @rdname plotResults
#' @export
plotVarianceComponents.mvremModel <- function(x, response = NULL, scale = c("sd", "variance"),
                                              ncol = NULL, interactive = NULL, ...) {
  interactive <- useInteractive(interactive)
  vc <- varianceComponents(x, response, scale)
  t <- vc$table
  if (!nrow(t)) stop("No variance components to plot.", call. = FALSE)
  t$panel <- ifelse(t$parameter == "correlation", "correlations",
                    ifelse(t$block == "residual", "residual", t$block))
  # the blocks (responses in the order of the model; correlations after the first correlated
  # block), then the residual SDs
  panels <- unique(c(setdiff(t$panel, "residual"), intersect("residual", t$panel)))
  label <- function(p, k) {
    if (p == "correlations") paste0(t$block[k], ": ", t$effect[k], ", ", t$effect2[k]) else
      if (p == "residual") t$response[k] else t$effect[k]
  }
  ps <- lapply(panels, function(p) {
    k <- which(t$panel == p)
    lab <- label(p, k)
    list(title = p, xlab = if (p == "correlations") "correlation" else vc$scale,
         curves = Map(function(l, m) list(label = l, density = m), lab, vc$marginals[k]))
  })
  densityPanels(ps, ncol, interactive, vline = NULL, x)
}

# Panels of posterior densities, interactive or base graphics: `panels` = list of
# list(title, xlab, curves = list(list(label, density = two-column matrix))); `vline` = a
# reference line (e.g. 0); returns the widget or `x` invisibly.
densityPanels <- function(panels, ncol, interactive, vline = NULL, x = NULL) {
  cols <- plotColors$densities
  colsOf <- function(p) cols[(seq_along(p$curves) - 1L) %% length(cols) + 1L]
  if (isTRUE(interactive)) {
    ps <- lapply(panels, function(p) {
      cc <- colsOf(p)
      c(list(title = p$title, xlab = p$xlab, ylab = "posterior density",
             curves = lapply(seq_along(p$curves), function(i) {
               cv <- p$curves[[i]]
               list(id = paste(p$title, cv$label), label = cv$label, color = cc[i], kind = "curve",
                    width = 2, x = r5(cv$density[, 1]), y = r5(cv$density[, 2]))
             })),
        if (!is.null(vline)) list(vline = vline))
    })
    return(htmlPlot(ps, autoNcol(length(ps), ncol), linked = FALSE, search = NULL))
  }
  op <- panelLayout(length(panels), ncol)
  on.exit(graphics::par(op), add = TRUE)
  for (p in panels) {
    ds <- lapply(p$curves, `[[`, "density")
    graphics::plot(range(c(vline, unlist(lapply(ds, function(m) m[, 1]))), finite = TRUE),
                   c(0, max(unlist(lapply(ds, function(m) m[, 2])), finite = TRUE)),
                   type = "n", xlab = p$xlab, ylab = "posterior density", main = p$title, cex.main = 0.95)
    if (!is.null(vline)) graphics::abline(v = vline, lty = 3, col = plotColors$base$reference)
    cc <- colsOf(p)
    for (i in seq_along(ds)) graphics::lines(ds[[i]], col = cc[i], lwd = 2)
    graphics::legend("topright", legend = vapply(p$curves, `[[`, "", "label"), col = cc, lwd = 2,
                     bty = "n", cex = 0.75)
  }
  invisible(x)
}

# Forest plot: estimates with 95% intervals, one row per label. The left margin grows with
# the longest label (term names), so it overrides the margin of panelLayout() for its panel.
forestPlot <- function(tab, labels, main) {
  n <- nrow(tab)
  graphics::par(mar = c(4, max(6, 0.6 * max(nchar(labels))), 2.5, 1))
  graphics::plot(tab$mean, seq_len(n), type = "n", ylim = c(n + 0.5, 0.5), yaxt = "n", ylab = "",
                 xlim = range(c(tab$q0.025, tab$q0.975), finite = TRUE),
                 xlab = "posterior mean and 95% interval", main = main, cex.main = 0.95)
  graphics::abline(v = 0, lty = 3, col = plotColors$base$reference)
  graphics::axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.85)
  graphics::segments(tab$q0.025, seq_len(n), tab$q0.975, seq_len(n), col = plotColors$base$forestInterval, lwd = 2)
  graphics::points(tab$mean, seq_len(n), pch = 16, col = plotColors$accent)
}

#' @rdname plotResults
#' @export
plotFixedEffects <- function(x, ...) UseMethod("plotFixedEffects")
#' @rdname plotResults
#' @export
plotFixedEffects.mvremModel <- function(x, response = NULL, type = c("forest", "density"), ncol = NULL,
                                        interactive = NULL, ...) {
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  fe <- fixedEffects(x)
  fe <- fe[fe$response %in% checkResponses(x, response), , drop = FALSE]
  resp <- orderResponses(x, fe$response)
  if (type == "density") {
    ps <- lapply(resp, function(r) {
      cols <- x$inla$fixed[[r]]$columns          # INLA names -> terms
      list(title = paste0(r, ": fixed effects"), xlab = "coefficient",
           curves = Map(function(l, nm) list(label = l, density = x$fit$marginals.fixed[[nm]]),
                        unname(cols), names(cols)))
    })
    return(densityPanels(ps, ncol, interactive, vline = 0, x))
  }
  if (isTRUE(interactive)) {
    ps <- lapply(resp, function(r) {
      d <- fe[fe$response == r, ]
      htmlIntervals(paste0(r, ": fixed effects"), paste(r, d$term), d$term, d$mean, d$q0.025, d$q0.975,
                    ylabels = TRUE)
    })
    return(htmlPlot(ps, autoNcol(length(ps), ncol), linked = FALSE, search = NULL))
  }
  op <- panelLayout(length(resp), ncol)
  on.exit(graphics::par(op), add = TRUE)
  for (r in resp) {
    d <- fe[fe$response == r, ]
    forestPlot(d, d$term, paste0(r, ": fixed effects"))
  }
  invisible(x)
}

#' @rdname plotResults
#' @export
plotAssociations <- function(x, ...) UseMethod("plotAssociations")
#' @rdname plotResults
#' @export
plotAssociations.mvremModel <- function(x, type = c("forest", "density"), interactive = NULL, ...) {
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  a <- associations(x)
  lab <- paste0(a$longitudinal, " -> ", a$outcome)
  if (type == "density") {
    nm <- associationRows(x)$.inla               # parallel to the rows of associations()
    p <- list(title = "associations (gamma)", xlab = "gamma",
              curves = Map(function(l, n) list(label = l, density = x$fit$marginals.hyperpar[[n]]), lab, nm))
    return(densityPanels(list(p), 1L, interactive, vline = 0, x))
  }
  if (isTRUE(interactive))
    return(htmlPlot(list(htmlIntervals("associations (gamma)", lab, lab, a$mean, a$q0.025, a$q0.975,
                                       ylabels = TRUE)), 1L, linked = FALSE, search = NULL))
  op <- graphics::par(no.readonly = TRUE)   # forestPlot() changes the margins
  on.exit(graphics::par(op), add = TRUE)
  forestPlot(a, lab, "associations (gamma)")
  invisible(x)
}
