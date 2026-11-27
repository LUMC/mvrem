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
#' @param colorBy,shapeBy `plotRandomEffects(type = "pairs")`: a categorical variable of the data
#'   colouring / shaping the subjects' points, as in [plot.mvremModel()] (a subject takes its
#'   value when it is the same in all its rows); the widget offers the variables in dropdowns and
#'   hides / shows a level by a click on its legend entry. Ignored by the other plots.
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
                                         ncol = NULL, interactive = NULL, colorBy = NULL,
                                         shapeBy = NULL, ...) {
  order <- match.arg(order)
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  re <- randomEffects(x, response)
  re$label <- modelEffectLabel(x, re$term)
  if (!is.null(term))   # `term`: as in the formula, as labelled, or "1" for the intercept
    re <- re[re$term %in% term | re$label %in% term | (re$term == "(Intercept)" & "1" %in% term), , drop = FALSE]
  if (!nrow(re)) stop("No random effects to plot.", call. = FALSE)
  if (type == "pairs") return(plotRandomPairs(x, re, ncol, interactive, colorBy, shapeBy))

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
    return(htmlPlot(unname(ps), autoNcol(length(ps), ncol), linked = TRUE, search = "subject"))
  }
  withPanels(length(panels), ncol, function(i) forestPlot(rows[[i]], NULL, panels[i], ylab = "subject"))
  invisible(x)
}

# Posterior means of two effects of one block against each other, a point per subject, with
# the 95% ellipse of the estimated distribution of the effects (posterior-mean SDs and
# correlation) and the ellipses with the correlation at the bounds of its 95% interval.
plotRandomPairs <- function(x, re, ncol, interactive, colorBy = NULL, shapeBy = NULL) {
  vc <- varianceComponentsOrNull(x, "the pairs plot (no ellipses)")$table
  # colorBy / shapeBy: a point is a subject (a level of the block's group), so it takes the
  # subject's value of the variable (NA when the value varies within the subject)
  prep <- plotData(x)
  groups <- groupSpec(prep, colorBy, shapeBy)
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
        cv = subjectGroupValues(prep, d$group[1L], lev, groups$vars),
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
      pn <- htmlPoints(p$title, p$level, p$x, p$y, p$xlab, p$ylab, p$cv)
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
    return(htmlPlot(ps, autoNcol(length(ps), ncol), legend, linked = TRUE, search = "subject",
                    groupBy = groups))
  }
  withPanels(length(pairs), ncol, function(k) {
    p <- pairs[[k]]
    e <- p$ellipses
    xlim <- range(c(p$x, e$mean$x, e$lo$x, e$hi$x), finite = TRUE)
    ylim <- range(c(p$y, e$mean$y, e$lo$y, e$hi$y), finite = TRUE)
    graphics::plot(p$x, p$y, cex = 0.6, xlim = xlim, ylim = ylim,
                   pch = pchOf(if (!is.null(groups$shape)) p$cv[[groups$shape]], groups, length(p$x)),
                   col = colorOf(if (!is.null(groups$color)) p$cv[[groups$color]], groups, plotColors$accent),
                   xlab = p$xlab, ylab = p$ylab, main = p$title, cex.main = 0.95)
    graphics::abline(h = 0, v = 0, lty = 3, col = plotColors$base$reference)
    if (!is.null(e)) {
      graphics::lines(e$lo$x, e$lo$y, col = ellipseCol, lty = 2)
      graphics::lines(e$hi$x, e$hi$y, col = ellipseCol, lty = 2)
      graphics::lines(e$mean$x, e$mean$y, col = ellipseCol, lwd = 1.8)
      graphics::legend("topleft", legend = p$note, bty = "n", cex = 0.85, text.col = "grey20")
    }
    if (k == 1L) groupLegends(groups)
  })
  invisible(x)
}

# The value of each variable `vars` per subject `levels` of grouping `group`, over the data of all
# responses: the subject's value when it is the same in all its rows, else NA (named list of
# character vectors parallel to `levels`; variables absent from the data with the group dropped).
subjectGroupValues <- function(prep, group, levels, vars) {
  out <- lapply(vars, function(v) {
    rows <- lapply(prep, function(df) if (!is.null(df[[group]]) && !is.null(df[[v]]))
      data.frame(g = as.character(df[[group]]), v = as.character(df[[v]])))
    rows <- do.call(rbind, rows)
    if (is.null(rows)) return(NULL)
    rows <- unique(rows[!is.na(rows$v), , drop = FALSE])
    one <- vapply(split(rows$v, rows$g), function(z) if (length(z) == 1L) z else NA_character_, "")
    unname(one[as.character(levels)])
  })
  Filter(Negate(is.null), stats::setNames(out, vars))
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
  withPanels(length(panels), ncol, function(k) {
    p <- panels[[k]]
    ds <- lapply(p$curves, `[[`, "density")
    panelFrame(list(title = p$title, xlab = p$xlab, ylab = "posterior density"), vline = vline,
               xlim = range(c(vline, unlist(lapply(ds, function(m) m[, 1]))), finite = TRUE),
               ylim = c(0, max(unlist(lapply(ds, function(m) m[, 2])), finite = TRUE)))
    cc <- colsOf(p)
    for (i in seq_along(ds)) graphics::lines(ds[[i]], col = cc[i], lwd = 2)
    graphics::legend("topright", legend = vapply(p$curves, `[[`, "", "label"), col = cc, lwd = 2,
                     bty = "n", cex = 0.75)
  })
  invisible(x)
}

# Forest plot: estimates with 95% intervals, one row per label, the first on top. With `labels`
# (term names) the left margin grows with the longest one, overriding panelLayout()'s margin for
# its panel; without them (the caterpillar of random effects: many thin rows) the rows are
# unlabelled, with `ylab`.
forestPlot <- function(tab, labels, main, ylab = "") {
  n <- nrow(tab)
  labelled <- !is.null(labels)
  if (labelled) graphics::par(mar = c(4, max(6, 0.6 * max(nchar(labels))), 2.5, 1))
  panelFrame(list(title = main, xlab = "posterior mean and 95% interval", ylab = ylab), vline = 0,
             xlim = range(c(tab$q0.025, tab$q0.975), finite = TRUE),
             ylim = if (labelled) c(n + 0.5, 0.5) else c(n, 1), yaxt = if (labelled) "n" else "s")
  if (labelled) graphics::axis(2, at = seq_len(n), labels = labels, las = 1, cex.axis = 0.85)
  b <- plotColors$base
  graphics::segments(tab$q0.025, seq_len(n), tab$q0.975, seq_len(n),
                     col = if (labelled) b$forestInterval else b$interval, lwd = if (labelled) 2 else 1)
  graphics::points(tab$mean, seq_len(n), pch = 16, cex = if (labelled) 1 else 0.45, col = plotColors$accent)
}

# The forest / density / interval-widget rendering shared by plotFixedEffects() and
# plotAssociations(): `groups` = one entry per panel, list(title, tab (mean, q0.025, q0.975),
# labels, ids (for the widget), marginals (named list of densities, for type "density")).
intervalPlot <- function(groups, type, ncol, interactive, x) {
  if (type == "density") {
    ps <- lapply(groups, function(g)
      list(title = g$title, xlab = g$xlab,
           curves = Map(function(l, m) list(label = l, density = m), g$labels, g$marginals)))
    return(densityPanels(ps, ncol, interactive, vline = 0, x))
  }
  if (isTRUE(interactive)) {
    ps <- lapply(groups, function(g)
      htmlIntervals(g$title, g$ids, g$labels, g$tab$mean, g$tab$q0.025, g$tab$q0.975, ylabels = TRUE))
    return(htmlPlot(ps, autoNcol(length(ps), ncol), linked = FALSE, search = NULL))
  }
  withPanels(length(groups), ncol, function(i) forestPlot(groups[[i]]$tab, groups[[i]]$labels, groups[[i]]$title))
  invisible(x)
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
  fe <- fe[fe$response %in% selectResponses(x, response), , drop = FALSE]
  groups <- lapply(orderResponses(x, fe$response), function(r) {
    d <- fe[fe$response == r, ]
    cols <- x$inla$fixed[[r]]$columns          # INLA names -> terms, parallel to the rows of d
    list(title = paste0(r, ": fixed effects"), xlab = "coefficient", tab = d, labels = d$term,
         ids = paste(r, d$term), marginals = x$fit$marginals.fixed[names(cols)])
  })
  intervalPlot(groups, type, ncol, interactive, x)
}

#' @rdname plotResults
#' @export
plotAssociations <- function(x, ...) UseMethod("plotAssociations")
#' @rdname plotResults
#' @export
plotAssociations.mvremModel <- function(x, response = NULL, type = c("forest", "density"),
                                        interactive = NULL, ...) {
  type <- match.arg(type)
  interactive <- useInteractive(interactive)
  a <- associations(x)
  nm <- associationRows(x)$.inla               # parallel to the rows of associations()
  keep <- a$outcome %in% selectResponses(x, response)
  if (!any(keep)) stop("No association for the response(s) given.", call. = FALSE)
  a <- a[keep, , drop = FALSE]; nm <- nm[keep]
  lab <- linkLabel(a$longitudinal, a$outcome)
  intervalPlot(list(list(title = "associations (gamma)", xlab = "gamma", tab = a, labels = lab, ids = lab,
                         marginals = x$fit$marginals.hyperpar[nm])), type, 1L, interactive, x)
}
