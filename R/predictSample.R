# predict(method = "sample"): predictions for new rows from joint posterior draws
# of a fit, without a new INLA run (see ?predict.mvremModel).

# Error unless the fit keeps INLA's integration points (`config = TRUE`).
needConfig <- function(x) {
  if (!hasConfig(x))
    stop("`method = \"sample\"` needs a fit with `control.compute = list(config = TRUE)`. ",
         "Fit the model with `jsm(..., inlaOptions = list(control.compute = list(config = TRUE)))` ",
         "(or `m <- update(m, inlaOptions = list(control.compute = list(config = TRUE)))`), ",
         "or use `method = \"refit\"`.", call. = FALSE)
  invisible(x)
}

# TRUE if the fit keeps INLA's integration points (`config = TRUE`): predict(method =
# "auto") then samples.
hasConfig <- function(x) !is.null(x$fit$misc$configs)

# Error of a model `method = "sample"` cannot handle (class `mvremSampleUnsupported`):
# predict(method = "auto") falls back to the refit.
sampleUnsupported <- function(...)
  stop(structure(class = c("mvremSampleUnsupported", "error", "condition"),
                 list(message = paste0(...), call = NULL)))

# Joint posterior draws of a fit with config = TRUE: an integration point per draw,
# chosen with probability proportional to its posterior density (as
# INLA::inla.posterior.sample()), then the latent field from INLA's Gaussian
# approximation at that point (improved mean, precision Q; linear constraints
# applied by conditioning). Returns list(latent = matrix (latent x draws, without
# the linear predictor), hyper = matrix (hyperparameters on the user scale, named as
# summary.hyperpar x draws), rows = function(tag) -> rows of `latent`).
posteriorDraws <- function(fit, n) {
  cs <- fit$misc$configs
  ct <- cs$contents
  # in INLA's compact mode (.preopt) the configurations exclude the linear predictor;
  # otherwise it is the first part of the field: it is sampled and not used
  if (isTRUE(cs$.preopt)) {
    keep <- !ct$tag %in% c("APredictor", "Predictor")
    ct <- list(tag = ct$tag[keep], start = ct$start[keep] - ct$start[keep][1L] + 1L, length = ct$length[keep])
  }
  nx <- length(cs$config[[1L]]$mean)
  hyperNames <- rownames(fit$summary.hyperpar)
  ld <- vapply(cs$config, function(cf) cf$log.posterior, 0)
  idx <- sample.int(length(ld), n, replace = TRUE, prob = exp(ld - max(ld)))
  latent <- matrix(NA_real_, nx, n)
  hyper <- matrix(NA_real_, length(hyperNames), n, dimnames = list(hyperNames, NULL))
  constr <- cs$constr
  for (k in sort(unique(idx))) {
    cols <- which(idx == k)
    cf <- cs$config[[k]]
    Q <- symmetricPrecision(cf$Q)
    ch <- Matrix::Cholesky(Q, LDL = FALSE, perm = TRUE)   # P Q P' = L L'
    z <- matrix(stats::rnorm(nx * length(cols)), nx)
    x <- as.matrix(Matrix::solve(ch, Matrix::solve(ch, z, system = "Lt"), system = "Pt"))
    if (!is.null(constr) && length(constr$e)) {
      A <- matrix(constr$A, ncol = nx)
      QiAt <- as.matrix(Matrix::solve(ch, t(A), system = "A"))
      x <- x - QiAt %*% solve(A %*% QiAt, A %*% x - constr$e)
    }
    mu <- if (!is.null(cf$improved.mean)) cf$improved.mean else cf$mean
    latent[, cols] <- x + mu
    if (length(hyperNames)) {
      th <- vapply(seq_along(hyperNames), function(j) fit$misc$from.theta[[j]](cf$theta[[j]]), 0)
      hyper[, cols] <- th
    }
  }
  rows <- function(tag) {
    i <- match(tag, ct$tag)
    if (is.na(i)) stop("Internal: no latent component '", tag, "' in the fit.", call. = FALSE)
    seq.int(ct$start[i], length.out = ct$length[i])
  }
  list(latent = latent, hyper = hyper, rows = rows)
}

# INLA stores one triangle of the symmetric precision matrix of a configuration
# (a full matrix is taken as symmetric).
symmetricPrecision <- function(Q) {
  S <- Matrix::summary(Q)
  uplo <- if (all(S$i >= S$j)) "L" else "U"
  Matrix::forceSymmetric(Q, uplo = uplo)
}

# Inverse link of a family (INLA's link functions; links with parameters are not
# supported).
inverseLink <- function(fam) {
  link <- fam$link
  if (is.null(link) || identical(link, "default")) {
    links <- familyLinks(fam$family)
    link <- if (length(links) > 1L) links[2L] else "identity"
  }
  f <- tryCatch(getExportedValue("INLA", paste0("inla.link.inv", link)), error = function(e) NULL)
  if (is.null(f) || link %in% c("robit", "sn", "gevit", "cgevit", "quantile", "pquantile", "powerlogit",
                                 "logitoffset", "sslogit", "loga"))
    sampleUnsupported("`method = \"sample\"` does not support the link '", link, "'; use `method = \"refit\"`.")
  function(eta) { out <- f(as.vector(eta)); dim(out) <- dim(eta); out }
}

# Summaries of draws (rows = predictions, columns = draws), as the INLA tables.
drawSummary <- function(s) {
  q <- t(apply(s, 1L, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE))
  mode <- apply(s, 1L, function(v) {
    if (!all(is.finite(v)) || stats::sd(v) == 0) return(v[1L])
    d <- stats::density(v)
    d$x[which.max(d$y)]
  })
  data.frame(mean = rowMeans(s), sd = apply(s, 1L, stats::sd), q0.025 = q[, 1L], q0.5 = q[, 2L],
             q0.975 = q[, 3L], mode = mode)
}

# Predictions by sampling: `inlNew` is buildInlaFormula() of the data with the new
# rows appended (the last rows of each predicted response's observation stack).
predictSampled <- function(object, inlNew, newdata, response, type, nsamples, seed) {
  requireInla()
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("`method = \"sample\"` needs the Matrix package.", call. = FALSE)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  fit <- object$fit
  inlOld <- object$inla
  eqs <- bothParts(object$forms)
  links <- lapply(stats::setNames(response, response), function(r)
    if (type == "mean") inverseLink(eqs[[r]]$family))   # checked before sampling
  dr <- posteriorDraws(fit, nsamples)
  hyperRow <- function(name, what) {
    if (!name %in% rownames(dr$hyper))
      sampleUnsupported("`method = \"sample\"`: the fit has no estimated ", what, " ('", name, "'); ",
                        "use `method = \"refit\"`.")
    dr$hyper[name, ]
  }

  # random effects of a block at the labels of the rebuilt model (rows) x draws:
  # labels of the fit -> their sampled values; new labels (new groups) -> draws from
  # N(0, Q^-1) with the sampled hyperparameters, one per group level and draw
  blockCache <- new.env(parent = emptyenv())
  blockDraws <- function(r, i) {
    key <- paste(r, i)
    if (!is.null(blockCache[[key]])) return(blockCache[[key]])
    b <- inlNew$random[[r]][[i]]
    field <- b$index[1L]
    valsNew <- inlNew$levels[[b$levels]]
    valsOld <- inlOld$levels[[inlOld$random[[r]][[i]]$levels]]
    k <- length(b$effects)
    B <- matrix(NA_real_, length(valsNew), nsamples)
    pos <- match(valsNew, valsOld)
    have <- !is.na(pos)
    B[have, ] <- dr$latent[dr$rows(field)[pos[have]], , drop = FALSE]
    if (any(!have)) {
      nlev <- length(valsNew) %/% k
      miss <- which(!have)
      lev <- sort(unique((miss - 1L) %% nlev + 1L))   # new group levels (positions)
      if (k == 1L) {
        tau <- hyperRow(paste0("Precision for ", field), "precision")
        B[lev, ] <- matrix(stats::rnorm(length(lev) * nsamples), length(lev)) *
          rep(1 / sqrt(tau), each = length(lev))
      } else {
        if (k > 5L)
          sampleUnsupported("`method = \"sample\"` does not support correlated blocks of more than five ",
                            "effects; use `method = \"refit\"`.")
        tau <- vapply(seq_len(k), function(e)
          hyperRow(sprintf("Precision for %s (component %d)", field, e), "precision"), numeric(nsamples))
        rho <- list()
        for (e1 in seq_len(k - 1L)) for (e2 in (e1 + 1L):k)
          rho[[paste(e1, e2)]] <- hyperRow(sprintf("Rho%d:%d for %s", e1, e2, field), "correlation")
        z <- array(stats::rnorm(length(lev) * k * nsamples), c(length(lev), k, nsamples))
        for (s in seq_len(nsamples)) {
          R <- diag(k)
          for (e1 in seq_len(k - 1L)) for (e2 in (e1 + 1L):k)
            R[e1, e2] <- R[e2, e1] <- rho[[paste(e1, e2)]][s]
          sdv <- 1 / sqrt(tau[s, ])
          U <- chol(R * outer(sdv, sdv))           # t(U) %*% U = covariance
          v <- z[, , s, drop = FALSE]
          dim(v) <- c(length(lev), k)
          b1 <- v %*% U                            # rows: levels, columns: effects
          for (e in seq_len(k)) B[(e - 1L) * nlev + lev, s] <- b1[, e]
        }
      }
    }
    blockCache[[key]] <- B
    B
  }

  # linear predictor (rows of `st` x draws) of response r at stack rows `st`
  etaRows <- function(r, st) {
    n <- nrow(st)
    eta <- matrix(0, n, nsamples)
    for (col in names(inlNew$fixed[[r]]$columns)) {
      if (!col %in% names(inlOld$fixed[[r]]$columns))
        stop("Internal: fixed-effect column '", col, "' is not in the fit.", call. = FALSE)
      eta <- eta + outer(st[[col]], dr$latent[dr$rows(col), ])
    }
    for (i in seq_along(inlNew$random[[r]])) {
      b <- inlNew$random[[r]][[i]]
      B <- blockDraws(r, i)
      vals <- inlNew$levels[[b$levels]]
      for (e in seq_along(b$index)) {
        w <- if (is.na(b$weight[e])) 1 else st[[b$weight[e]]]
        eta <- eta + w * B[match(st[[b$index[e]]], vals), , drop = FALSE]
      }
    }
    for (sc in Filter(function(s) s$out == r, inlNew$scaled)) {
      z <- inlNew$zero[[paste0(sc$long, "_", r)]]
      zst <- inlNew$stacks[[sc$long]][[r]]
      zrows <- zst[match(st[[sc$index]], zst[[z$index]]), , drop = FALSE]
      gamma <- hyperRow(paste0("Beta for ", sc$index), "association")
      eta <- eta + etaRows(sc$long, zrows) * rep(gamma, each = n)
    }
    eta
  }

  lapply(stats::setNames(response, response), function(r) {
    st <- inlNew$stacks[[r]]$main
    n <- nrow(newdata[[r]])
    st <- st[seq.int(nrow(st) - n + 1L, nrow(st)), , drop = FALSE]
    eta <- etaRows(r, st)
    if (type == "mean") eta <- links[[r]](eta)
    drawSummary(eta)
  })
}
