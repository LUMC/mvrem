#' Prior distributions for random-effect blocks
#'
#' Constructors of prior distributions for the hyperparameters of a
#' random-effect block, given as `re(..., prior = <constructor>)`. Each name
#' tells the distribution and the quantity it is placed on: the precision
#' \eqn{\tau}, the standard deviation \eqn{\sigma = \tau^{-1/2}}, or the
#' precision matrix \eqn{\mathbf{Q}} of a correlated block.
#'
#' * `gammaPrec()`: \eqn{\tau \sim \mathrm{Gamma}(shape, rate)} (INLA `loggamma`;
#'   the defaults are INLA's default for `iid` effects).
#' * `pcSd()`: penalised-complexity prior with \eqn{P(\sigma > u) = \alpha},
#'   i.e. \eqn{\sigma \sim \mathrm{Exp}(-\log(\alpha) / u)} (INLA `pc.prec`).
#' * `halfNormalSd()`: \eqn{\sigma \sim \mathcal{N}^{+}(0, scale^2)} (INLA
#'   `logtnormal` with mean 0 and precision \eqn{1/scale^2}).
#' * `halfCauchySd()`: \eqn{\sigma \sim \mathrm{C}^{+}(0, scale)}, the
#'   half-Cauchy prior of Gelman (2006) (INLA `expression:` prior on
#'   \eqn{\log \tau}).
#' * `wishartPrec()`: \eqn{\mathbf{Q} \sim \mathcal{W}_k(df, \mathbf{V})} with
#'   \eqn{\mathrm{E}[\mathbf{Q}] = df\,\mathbf{V}} (INLA `wishart<k>d` /
#'   `wishartkd`, whose parameter \eqn{\mathbf{R} = \mathbf{V}^{-1}}).
#'
#' Scalar priors are for single effects and uncorrelated (`||`) blocks; the same
#' prior is used for every effect of the block (split the block into several
#' `re()` terms for different priors). A correlated block with \eqn{k > 1}
#' effects accepts only `wishartPrec()`. The dimension \eqn{k} is known once the
#' data is (factor levels expand), so `wishartPrec()` is checked against the
#' block by [buildInlaFormula()].
#'
#' @param shape,rate Shape and rate of the gamma distribution.
#' @param u,alpha Upper bound and tail probability: \eqn{P(\sigma > u) = \alpha}.
#' @param scale For `halfNormalSd()` / `halfCauchySd()`: scale of the
#'   distribution of \eqn{\sigma}. For `wishartPrec()`: the scale matrix
#'   \eqn{\mathbf{V}}, a symmetric positive-definite \eqn{k \times k} matrix or
#'   a positive number \eqn{s} meaning \eqn{s\,\mathbf{I}_k}.
#' @param df Degrees of freedom of the Wishart distribution, \eqn{df > k + 1};
#'   `NULL` means \eqn{k + 2}.
#' @return An object of class `mvremPrior`.
#' @examples
#' pcSd(u = 1)
#' wishartPrec(df = 4, scale = diag(c(1, 0.1)))
#' @name priors
NULL

#' @rdname priors
#' @export
gammaPrec <- function(shape = 1, rate = 5e-5) {
  checkPositive(shape, "shape")
  checkPositive(rate, "rate")
  newPrior("gamma", "prec", list(shape = shape, rate = rate))
}

#' @rdname priors
#' @export
pcSd <- function(u, alpha = 0.01) {
  checkPositive(u, "u")
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  newPrior("pc", "sd", list(u = u, alpha = alpha))
}

#' @rdname priors
#' @export
halfNormalSd <- function(scale) {
  checkPositive(scale, "scale")
  newPrior("halfNormal", "sd", list(scale = scale))
}

#' @rdname priors
#' @export
halfCauchySd <- function(scale) {
  checkPositive(scale, "scale")
  newPrior("halfCauchy", "sd", list(scale = scale))
}

#' @rdname priors
#' @export
wishartPrec <- function(df = NULL, scale = 1) {
  if (!is.null(df)) checkPositive(df, "df")
  if (is.matrix(scale)) {
    if (!is.numeric(scale) || anyNA(scale) || nrow(scale) != ncol(scale) ||
        !isSymmetric(unname(scale)) ||
        any(eigen(scale, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("`scale` must be a symmetric positive-definite matrix.", call. = FALSE)
    scale <- unname(scale)
  } else {
    checkPositive(scale, "scale")
  }
  newPrior("wishart", "precMatrix", list(df = df, scale = scale))
}

newPrior <- function(dist, target, args) {
  structure(list(dist = dist, target = target, args = args), class = "mvremPrior")
}

checkPositive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0)
    stop("`", name, "` must be a single positive number.", call. = FALSE)
}

#' @export
format.mvremPrior <- function(x, ...) {
  a <- x$args
  num <- function(v) vapply(v, function(z) format(signif(z, 4)), "")
  txt <- switch(x$dist,
    gamma = paste0("Gamma prior on precision: tau ~ Gamma(shape = ", num(a$shape),
                   ", rate = ", num(a$rate), ")   [INLA: loggamma, param = c(",
                   num(a$shape), ", ", num(a$rate), ")]"),
    pc = paste0("PC prior on SD: P(sigma > ", num(a$u), ") = ", num(a$alpha),
                ", i.e. sigma ~ Exp(", num(-log(a$alpha) / a$u), ")   [INLA: pc.prec, param = c(",
                num(a$u), ", ", num(a$alpha), ")]"),
    halfNormal = paste0("Half-normal prior on SD: sigma ~ N+(0, sd = ", num(a$scale),
                        ")   [INLA: logtnormal, param = c(0, ", num(1 / a$scale^2), ")]"),
    halfCauchy = paste0("Half-Cauchy prior on SD: sigma ~ C+(0, ", num(a$scale),
                        ")   [INLA: expression: prior on log precision]"),
    wishart = paste0("Wishart prior on precision matrix: Q ~ W(df = ",
                     if (is.null(a$df)) "k + 2" else num(a$df), ", scale = ",
                     scaleText(a$scale, num), ")   [INLA: wishart<k>d / wishartkd, R = scale^-1]")
  )
  paste("<mvrem prior>", txt)
}

#' @export
print.mvremPrior <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

scaleText <- function(s, num) {
  if (!is.matrix(s)) return(paste(num(s), "* I"))
  if (all(s[row(s) != col(s)] == 0)) return(paste0("diag(", paste(num(diag(s)), collapse = ", "), ")"))
  paste0("rbind(", paste(apply(s, 1L, function(r) paste0("c(", paste(num(r), collapse = ", "), ")")),
                         collapse = ", "), ")")
}

# Wishart settings resolved for a block of k effects: df, k x k scale matrix.
wishartArgs <- function(prior, k) {
  a <- prior$args
  df <- if (is.null(a$df)) k + 2 else a$df
  V <- if (is.matrix(a$scale)) a$scale else diag(a$scale, k)
  list(df = df, scale = V)
}

# Check that `prior` suits a block of k effects (correlated: `|` block).
checkBlockPrior <- function(prior, k, correlated, term) {
  if (is.null(prior)) return(invisible())
  blk <- paste0("`", deparse1(term), "`")
  if (correlated && k > 1L) {
    if (prior$dist != "wishart")
      stop("Correlated block ", blk, " (k = ", k, ") needs wishartPrec(); use `||` or ",
           "separate re() blocks for scalar priors.", call. = FALSE)
    if (is.matrix(prior$args$scale) && nrow(prior$args$scale) != k)
      stop("wishartPrec() scale of block ", blk, " must be a ", k, " x ", k, " matrix.",
           call. = FALSE)
    df <- wishartArgs(prior, k)$df
    if (df <= k + 1)
      stop("wishartPrec() of block ", blk, ": df must be > k + 1 = ", k + 1L, ".",
           call. = FALSE)
  } else if (prior$dist == "wishart") {
    stop("wishartPrec() needs a correlated block with >= 2 effects (block ", blk,
         "); use gammaPrec(), pcSd(), halfNormalSd() or halfCauchySd().", call. = FALSE)
  }
  invisible()
}

# INLA `hyper` list of a prior for a block of k effects.
priorHyper <- function(prior, k) {
  a <- prior$args
  prec <- function(p, param) list(prec = c(list(prior = p), if (!is.null(param)) list(param = param)))
  switch(prior$dist,
    gamma = prec("loggamma", c(a$shape, a$rate)),
    pc = prec("pc.prec", c(a$u, a$alpha)),
    halfNormal = prec("logtnormal", c(0, 1 / a$scale^2)),
    halfCauchy = prec(halfCauchyExpression(a$scale), NULL),
    wishart = {
      w <- wishartArgs(prior, k)
      R <- solve(w$scale)
      list(theta1 = list(prior = if (k <= 5L) paste0("wishart", k, "d") else "wishartkd",
                         param = c(w$df, diag(R), R[lower.tri(R)])))
    }
  )
}

# Density of sigma = exp(-theta / 2) ~ C+(0, s), on theta = log precision
# (Jacobian |d sigma / d theta| = sigma / 2).
halfCauchyExpression <- function(s) {
  paste0("expression: sigma = exp(-theta/2); s = ", format(s, digits = 15),
         "; logdens = log(2) - log(pi) - log(s) - log(1 + (sigma/s)^2);",
         " return(logdens + log(sigma) - log(2));")
}

# LaTeX lines of the prior of an RE block: tau = precision symbol(s) (scalar
# priors) or Q (wishart); idx = effect indices for the "k = ..." suffix.
priorLatex <- function(prior, k, tau, Q, idx) {
  if (is.null(prior)) return(character())
  a <- prior$args
  suffix <- if (length(idx) > 1L) paste0(", \\quad k = ", paste(idx, collapse = ", ")) else ""
  sd <- paste0("\\left(", tau, "\\right)^{-1/2}")
  switch(prior$dist,
    gamma = paste0(tau, " &\\sim \\mathrm{Gamma}\\left(", numLatex(a$shape), ", ",
                   numLatex(a$rate), "\\right)", suffix),
    pc = paste0(sd, " &\\sim \\mathrm{Exp}\\left(", numLatex(-log(a$alpha) / a$u),
                "\\right), \\quad P\\left(", sd, " > ", numLatex(a$u), "\\right) = ",
                numLatex(a$alpha), suffix),
    halfNormal = paste0(sd, " &\\sim \\mathcal{N}^{+}\\left(0, \\left(",
                        numLatex(1 / a$scale^2), "\\right)^{-1}\\right)", suffix),
    halfCauchy = paste0(sd, " &\\sim \\mathrm{C}^{+}\\left(0, ", numLatex(a$scale), "\\right)",
                        suffix),
    wishart = {
      kk <- if (is.matrix(a$scale)) nrow(a$scale) else k
      w <- wishartArgs(prior, kk)
      paste0(Q, " &\\sim \\mathcal{W}_{", kk, "}\\left(", numLatex(w$df), ", ",
             matrixLatex(a$scale, kk), "\\right)")
    }
  )
}

# Number with 3 significant digits; 5e-05 -> 5 \cdot 10^{-5}.
numLatex <- function(x) {
  s <- formatC(signif(x, 3), digits = 3, format = "g")
  s <- trimws(s)
  if (grepl("e", s)) {
    parts <- strsplit(s, "e", fixed = TRUE)[[1L]]
    s <- paste0(if (parts[1L] != "1") paste0(parts[1L], " \\cdot "), "10^{",
                as.integer(parts[2L]), "}")
  }
  s
}

matrixLatex <- function(s, k) {
  if (!is.matrix(s)) {
    I <- paste0("\\mathbf{I}_{", k, "}")
    return(if (s == 1) I else paste0(numLatex(s), "\\,", I))
  }
  rows <- apply(s, 1L, function(r) paste(vapply(r, numLatex, ""), collapse = " & "))
  paste0("\\begin{pmatrix}", paste(rows, collapse = " \\\\ "), "\\end{pmatrix}")
}
