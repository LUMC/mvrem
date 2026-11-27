#' Diagram of a model
#'
#' Draws the structure of a model as a path diagram: one box per equation (the
#' longitudinal responses on the left, the outcomes on the right) with its
#' response, family (and link, `trials()` / `exposure()` / `event()` / `weights()` input) and
#' random-effect blocks, and one arrow per association from a longitudinal
#' response to an outcome, labelled with its parameter `gamma`.
#'
#' For a fitted model the arrows carry the posterior mean and the 95% credible
#' interval of `gamma`, and each box adds the posterior means of the standard
#' deviations (and correlations) of its random-effect blocks and of the residual
#' standard deviation (see [varianceComponents()]); for a model that was only constructed (`fit = FALSE`)
#' they carry its symbol, and `gamma` with a `prior` or other parameters in
#' `scaled()` is marked with `*`.
#'
#' @param x A model (class `mvremModel`), e.g. from [jsm()].
#' @param detail `"model"` (default): the equations and associations. `"inla"`:
#'   also the zero equations, as INLA sees them: each association passes through
#'   the field `u` of its zero equation (`0 = eta - u`, then the outcome uses
#'   `gamma * u`).
#' @param fitted If `TRUE` (default) and the model is fitted, show the estimates
#'   of `gamma` and the standard deviations.
#' @param cex Text size (shrunk automatically when the diagram does not fit).
#' @param ... Ignored.
#' @return The text size used (`cex`, shrunk if the diagram did not fit the
#'   device at the size given), invisibly.
#' @examples
#' \dontrun{
#' m <- jsm(list(y ~ 1 + t + scaled(x1)), list(x1 ~ 1 + t + (1 + t | id)),
#'          data = d, family = list(y = betaFamily()), fit = FALSE)
#' plotModel(m)
#' plotModel(m, detail = "inla")
#' }
#' @export
plotModel <- function(x, ...) UseMethod("plotModel")

#' @rdname plotModel
#' @export
plotModel.mvremModel <- function(x, detail = c("model", "inla"), fitted = TRUE, cex = 0.9, ...) {
  detail <- match.arg(detail)
  showFit <- isTRUE(fitted) && !is.null(x$fit)
  vc <- if (showFit) varianceComponentsOrNull(x, "the model diagram")
  nodes <- diagramNodes(x$forms, vc, x$inla$naming)
  links <- diagramLinks(c(x$forms$out, latentParts(x$forms)))   # E38: latent links too
  est <- if (showFit) associations(x)

  op <- graphics::par(mar = c(0.5, 0.5, 0.5, 0.5), xaxs = "i", yaxs = "i")
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new()
  g <- diagramGeometry(nodes, detail, cex)
  graphics::plot.window(c(0, 1), c(0, 1))

  drawDiagramBoxes(nodes, g)
  segs <- drawDiagramArrows(links, g, detail)
  # the labels after all arrows, so that they lie on top
  labels <- lapply(seq_len(nrow(links)), function(i) gammaLabel(links[i, ], est))
  drawArrowLabels(segs, labels, g)
  if (detail == "inla" && nrow(links)) drawInlaHeader(g)
  invisible(g$cex)
}

# Layout constants of the diagram. Sizes in inches or in text lines (one line = the device's
# character height x `lineHeight`); `user` = units of the plot window, 0..1 in x and y.
diagramStyle <- list(
  col = "#003c7b", fillLong = "#dcecf7", fillOut = "#fdebd2",   # text and borders, box fills
  fillLatent = "#e3f1e0",  # box fill of a latent equation (middle column, E38)
  lineHeight = 1.25,       # line spacing (character heights)
  boxPadIn = 0.2,          # box width beyond its widest line (inches)
  columnGapIn = c(model = 1.5, inla = 2),   # space between the two columns of boxes (inches)
  boxLinesFit = 0.9,       # extra lines per box when checking that a column fits the height
  boxPadLines = 0.5,       # box height beyond its lines (lines)
  boxGapLines = 0.4,       # vertical gap between boxes (lines)
  textTopLines = 0.25,     # extra space above the first text line in a box (lines)
  shrink = 0.92,           # cex factor per step while the diagram does not fit
  minCex = 0.3,            # stop shrinking below this
  boxLwd = 1.5,
  arrowEndSpread = 0.7,    # share of a box's height over which several arrow ends spread
  arrow = list(model = c(length = 0.1, lwd = 1.5), inla = c(length = 0.08, lwd = 1.2)),
  uHalfWidth = 0.06,       # detail = "inla": half width of the box of u (user)
  uHalfHeightLines = 0.7,  # ... and its half height (lines)
  uCex = 0.9,              # ... text size of u relative to cex
  labelPadW = 1.1,         # gamma label background: width and height relative to the text
  labelPadH = 1.6,
  # positions tried for a gamma label, as fractions along its arrow, in order of preference
  labelFracs = c(0.5, 0.35, 0.65, 0.25, 0.75, 0.15, 0.85),
  crossSamples = 60,       # points per arrow when counting arrows that cross a label
  headerLines = 1.2,       # detail = "inla": header above the top box (lines)
  headerShift = 0.03)      # ... moved outwards from the arrow midpoints (user)

# The boxes of the diagram: one per equation, longitudinal (side "L") before outcomes ("R");
# lines: response, family / link / inputs, random-effect blocks; with variance components
# `vc` (varianceComponents() of a fit) also their estimates (sdLines()).
diagramNodes <- function(forms, vc = NULL, naming = namingScheme()) {
  boxLines <- function(p) {
    fam <- p$family
    link <- fam$link
    fl <- paste0(fam$family, if (!is.null(link)) paste0(", ", link, " link"))
    aux <- vapply(names(p$aux), function(a) paste0(a, "(", deparse1(p$aux[[a]]), ")"), "")
    re <- vapply(p$re, function(b) paste0("(", deparse1(b$term), ")"), "")
    c(p$lhs, paste0(fl, if (length(aux)) paste0("; ", paste(aux, collapse = ", "))),
      if (length(re)) paste("re:", re) else "no random effects",
      sdLines(vc, p$response, naming))
  }
  latentLines <- function(p) {
    re <- vapply(p$re, function(b) paste0("(", deparse1(b$term), ")"), "")
    c(paste0("latent(", p$response, ")"), "no observations",
      if (length(re)) paste("re:", re) else "no random effects",
      sdLines(vc, p$response, naming))
  }
  c(lapply(forms$long, function(p) list(id = p$response, lines = boxLines(p), side = "L")),
    lapply(latentParts(forms), function(p) list(id = p$response, lines = latentLines(p), side = "M")),
    lapply(forms$out, function(p) list(id = p$response, lines = boxLines(p), side = "R")))
}

# The arrows: one row per scaled() term (from = longitudinal, to = outcome); `prior` marks a
# gamma with a prior or other parameters given in scaled().
diagramLinks <- function(outs) {
  rows <- lapply(Filter(function(o) length(o$scaled) > 0L, outs), function(o)
    data.frame(from = scaledTargets(o), to = o$response,
               prior = vapply(o$scaled, function(s) !is.null(s$prior) || length(s$params) > 0L, NA),
               coef = vapply(o$scaled, function(s) if (is.null(s$coef)) NA_real_ else s$coef, 0),
               stringsAsFactors = FALSE))
  bindRows(rows, data.frame(from = character(), to = character(), prior = logical(), coef = numeric(),
                            stringsAsFactors = FALSE))
}

# The geometry of the diagram on the current device (after plot.new()): the text size `cex`
# (shrunk until both columns of boxes and the gap between them fit), the box width `wx`, the
# line height `lh` (inches), user units per inch `yu`, the box centres `xL` / `xR`, the arrow
# ends `x0` / `x1` and the boxes `pos` (per id: centre `y`, height `h`).
diagramGeometry <- function(nodes, detail, cex) {
  st <- diagramStyle
  pin <- graphics::par("pin")
  lineH <- function(cex) graphics::par("cin")[2L] * cex * st$lineHeight
  colLines <- function(side) sum(vapply(Filter(function(n) n$side == side, nodes),
                                       function(n) length(n$lines) + st$boxLinesFit, 0))
  nCol <- if (any(vapply(nodes, function(n) n$side == "M", NA))) 3L else 2L   # E38: latent column
  repeat {
    wIn <- max(unlist(lapply(nodes, function(n) graphics::strwidth(n$lines, "inches", cex)))) + st$boxPadIn
    okW <- nCol * wIn + (nCol - 1L) * st$columnGapIn[[detail]] <= pin[1L]
    okH <- max(colLines("L"), colLines("M"), colLines("R")) * lineH(cex) <= pin[2L]
    if ((okW && okH) || cex < st$minCex) break
    cex <- cex * st$shrink
  }
  lh <- lineH(cex)
  xu <- 1 / pin[1L]; yu <- 1 / pin[2L]          # user units per inch (x, y in 0..1)
  wx <- wIn * xu
  xL <- wx / 2; xR <- 1 - wx / 2                 # box centres
  list(cex = cex, lh = lh, yu = yu, wx = wx, xL = xL, xR = xR, x0 = xL + wx / 2, x1 = xR - wx / 2,
       xc = c(L = xL, M = 0.5, R = xR),
       side = stats::setNames(vapply(nodes, `[[`, "", "side"), vapply(nodes, `[[`, "", "id")),
       pos = c(stackBoxes(nodes, "L", lh, yu), stackBoxes(nodes, "M", lh, yu), stackBoxes(nodes, "R", lh, yu)))
}

# The boxes of one side, stacked and centred vertically: per id c(y = centre, h = height).
stackBoxes <- function(nodes, side, lh, yu) {
  st <- diagramStyle
  ns <- Filter(function(n) n$side == side, nodes)
  h <- vapply(ns, function(n) (length(n$lines) + st$boxPadLines) * lh * yu, 0)
  gap <- st$boxGapLines * lh * yu
  tot <- sum(h) + gap * (length(ns) - 1L)
  top <- 0.5 + tot / 2
  y <- numeric(length(ns)); cur <- top
  for (k in seq_along(ns)) { y[k] <- cur - h[k] / 2; cur <- cur - h[k] - gap }
  stats::setNames(Map(function(yy, hh) c(y = yy, h = hh), y, h), vapply(ns, `[[`, "", "id"))
}

# The boxes of the diagram: one per node at its computed position, the longitudinal ones on
# the left, with the node's text lines (response, family, RE blocks, fitted SDs).
drawDiagramBoxes <- function(nodes, g) {
  st <- diagramStyle
  for (id in names(g$pos)) {
    n <- Filter(function(z) z$id == id, nodes)[[1L]]
    p <- g$pos[[id]]
    xc <- g$xc[[n$side]]
    fill <- switch(n$side, L = st$fillLong, M = st$fillLatent, R = st$fillOut)
    graphics::rect(xc - g$wx / 2, p["y"] - p["h"] / 2, xc + g$wx / 2, p["y"] + p["h"] / 2,
                   col = fill, border = st$col, lwd = st$boxLwd)
    ys <- p["y"] + p["h"] / 2 - (seq_along(n$lines) - 0.5 + st$textTopLines) * g$lh * g$yu
    graphics::text(xc, ys, n$lines, cex = g$cex, col = st$col, font = c(2L, rep(1L, length(n$lines) - 1L)))
  }
}

# The heights where the arrows leave / reach the boxes: several arrows at one box are spread
# over its edge, so that they do not start / end at one point.
arrowEnds <- function(links, pos, column) {
  out <- rep(NA_real_, nrow(links))
  for (id in unique(links[[column]])) {
    k <- which(links[[column]] == id)
    h <- pos[[id]]["h"] * diagramStyle$arrowEndSpread
    out[k] <- pos[[id]]["y"] + if (length(k) == 1L) 0 else seq(h / 2, -h / 2, length.out = length(k))
  }
  out
}

# Draws the arrows (detail = "inla": through a box of the zero-equation field u) and returns
# per arrow the segment its label goes on (x0, y0, x1, y1): the whole arrow, or the part after u.
drawDiagramArrows <- function(links, g, detail) {
  st <- diagramStyle
  yFrom <- arrowEnds(links, g$pos, "from")
  yTo <- arrowEnds(links, g$pos, "to")
  a <- st$arrow[[detail]]
  lapply(seq_len(nrow(links)), function(i) {
    # from the right edge of the source's column to the left edge of the target's
    x0 <- g$xc[[g$side[[links$from[i]]]]] + g$wx / 2
    x1 <- g$xc[[g$side[[links$to[i]]]]] - g$wx / 2
    if (detail == "model") {
      graphics::arrows(x0, yFrom[i], x1, yTo[i], length = a[["length"]], col = st$col, lwd = a[["lwd"]])
      return(c(x0, yFrom[i], x1, yTo[i]))
    }
    uw <- st$uHalfWidth
    xm <- (x0 + x1) / 2; ym <- (yFrom[i] + yTo[i]) / 2
    graphics::arrows(x0, yFrom[i], xm - uw, ym, length = a[["length"]], col = st$col, lwd = a[["lwd"]])
    graphics::arrows(xm + uw, ym, x1, yTo[i], length = a[["length"]], col = st$col, lwd = a[["lwd"]])
    graphics::rect(xm - uw, ym - g$lh * g$yu * st$uHalfHeightLines, xm + uw, ym + g$lh * g$yu * st$uHalfHeightLines,
                   col = "white", border = st$col)
    graphics::text(xm, ym, bquote(u^(.(linkLabel(links$from[i], links$to[i])))),
                   cex = g$cex * st$uCex, col = st$col)
    c(xm + uw, ym, x1, yTo[i])
  })
}

# The label of an arrow (one row of diagramLinks()): the estimate of gamma with its 95%
# interval when `est` (associations()) has it, else the symbol (`*`: prior / parameters set).
gammaLabel <- function(link, est) {
  # scaled(x, coef = c): fixed, no estimate (E38)
  if (!is.null(link$coef) && !is.na(link$coef))
    return(bquote(gamma == .(format(link$coef)) ~ "(fixed)"))
  symbol <- if (link$prior) quote(gamma * "*") else quote(gamma)
  if (is.null(est)) return(symbol)
  e <- est[est$outcome == link$to & est$longitudinal == link$from, ]
  if (nrow(e) != 1L) return(symbol)
  bquote(gamma == .(sprintf("%.2f", e$mean)) ~ (.(sprintf("%.2f", e[["q0.025"]])) *
                                               ", " * .(sprintf("%.2f", e[["q0.975"]]))))
}

# Draws the labels of the arrows on white backgrounds. Each label goes to the first position
# along its arrow (`segs`) that does not overlap a label placed before, choosing among those
# the one crossed by the fewest other arrows. The preferred positions of successive arrows
# differ (0.3, 0.5, 0.7 of the way) so that the labels of parallel arrows do not line up.
drawArrowLabels <- function(segs, labels, g) {
  st <- diagramStyle
  placed <- list()
  overlaps <- function(r) any(vapply(placed, function(q)
    r[1L] < q[3L] && r[3L] > q[1L] && r[2L] < q[4L] && r[4L] > q[2L], NA))
  crossed <- function(r, seg) sum(vapply(segs, function(q) {
    if (identical(q, seg)) return(0)
    t <- seq(0, 1, length.out = st$crossSamples)
    px <- q[1L] + t * (q[3L] - q[1L]); py <- q[2L] + t * (q[4L] - q[2L])
    any(px > r[1L] & px < r[3L] & py > r[2L] & py < r[4L])
  }, 0))
  for (i in seq_along(segs)) {
    seg <- segs[[i]]
    lab <- labels[[i]]
    pref <- 0.3 + 0.4 * ((i - 1L) %% 3L) / 2
    w <- graphics::strwidth(lab, cex = g$cex) * st$labelPadW
    h <- graphics::strheight(lab, cex = g$cex) * st$labelPadH
    cand <- lapply(unique(c(pref, st$labelFracs)), function(f) {   # box (x0, y0, x1, y1) + centre
      xm <- seg[1L] + f * (seg[3L] - seg[1L]); ym <- seg[2L] + f * (seg[4L] - seg[2L])
      c(xm - w / 2, ym - h / 2, xm + w / 2, ym + h / 2, xm, ym)
    })
    free <- Filter(function(r) !overlaps(r[1:4]), cand)
    if (!length(free)) free <- cand                       # nothing free: all positions
    r <- free[[which.min(vapply(free, crossed, 0, seg = seg))]]   # fewest crossing arrows, then first
    placed[[length(placed) + 1L]] <- r
    graphics::rect(r[1L], r[2L], r[3L], r[4L], col = "white", border = NA)
    graphics::text(r[5L], r[6L], lab, cex = g$cex, col = st$col)
  }
}

# detail = "inla": what each half of the arrows means, above the boxes.
drawInlaHeader <- function(g) {
  st <- diagramStyle
  x0 <- g$x0; x1 <- g$x1
  top <- max(vapply(g$pos, function(p) p["y"] + p["h"] / 2, 0)) + st$headerLines * g$lh * g$yu
  graphics::text((x0 + (x0 + x1) / 2) / 2 - st$headerShift, top, "0 = eta - u", cex = g$cex, col = st$col)
  graphics::text(((x0 + x1) / 2 + x1) / 2 + st$headerShift, top, "outcome uses gamma * u", cex = g$cex, col = st$col)
}

# E6: the estimates of the variance components of response `r` in its box (posterior means,
# 2 significant digits): per random-effect block "SD | id: (Intercept) 0.46, t 0.18; cor 0.57",
# then "residual SD 0.22". Empty without `vc`.
sdLines <- function(vc, r, naming) {
  if (is.null(vc)) return(character())
  t <- vc$table[vc$table$response == r, , drop = FALSE]
  num <- function(v) formatC(signif(v, 2), digits = 2, format = "fg", flag = "#")
  out <- character()
  for (b in unique(t$block[t$block != "residual"])) {
    tb <- t[t$block == b, , drop = FALSE]
    sd <- tb[tb$parameter != "correlation", , drop = FALSE]
    co <- tb[tb$parameter == "correlation", , drop = FALSE]
    out <- c(out, paste0(
      if (vc$scale == "sd") "SD" else "variance", " | ", sd$group[1L], ": ",
      paste(effectLabel(sd$effect, naming), num(sd$mean), collapse = ", "),
      if (nrow(co)) paste0("; cor ", paste(
        if (nrow(co) > 1L) paste0(effectLabel(co$effect, naming), "~", effectLabel(co$effect2, naming), " ") else "",
        num(co$mean), sep = "", collapse = ", "))))
  }
  res <- t[t$block == "residual", , drop = FALSE]
  if (nrow(res)) out <- c(out, paste("residual SD", num(res$mean[1L])))
  out
}
