# Assembly and execution of the INLA call, and friendly summaries of its
# results. Nothing here is exported: users go through jsm() and the accessors
# of R/mvremModel.R, and never need INLA's interface.

# Arguments INLA::inla() must receive from mvrem (not overridable).
reservedInlaArguments <- c("formula", "family", "data", "control.family", "Ntrials", "E", "scale")

# Stacked data of all responses as one data frame (rows in the order of the INLA
# likelihoods; columns absent from a part are NA), with the block of each row.
stackRows <- function(inl, n = Inf) {
  parts <- list()
  blocks <- character()
  for (r in names(inl$stacks)) for (p in names(inl$stacks[[r]])) {
    s <- utils::head(inl$stacks[[r]][[p]], n)
    block <- stackBlockLabel(r, p)
    s$.block <- rep_len(block, nrow(s))   # rep_len: a zero-row part (n = 0) takes no length-1 value
    parts[[length(parts) + 1L]] <- s
    blocks <- c(blocks, block)
  }
  out <- rbindFill(parts)
  out$.block <- factor(out$.block, blocks)
  rownames(out) <- NULL
  out
}

# The survival responses of a buildInlaFormula() result (names of `inl$survival` that are TRUE).
survivalResponses <- function(inl) names(inl$survival)[vapply(inl$survival, isTRUE, NA)]

# The `.block` label of a stack part of response r: "y: observations" (part "main"),
# "x1: zero equations for y" (the zero part for outcome y) or "eta: latent equations".
stackBlockLabel <- function(r, part = "main") {
  if (part == "main") return(paste0(r, ": observations"))
  if (part == latentStackPart) return(paste0(r, ": latent equations"))
  paste0(r, ": zero equations for ", part)
}

# The name of a latent equation's only stack part (its zero rows; E38).
latentStackPart <- ".latent"

# Positions of the observation rows of response r among the rows of stackRows() (= INLA's
# fitted values / linear predictor), given the `.block` column as character.
observationRows <- function(block, r) which(block == stackBlockLabel(r))

# rbind() of data frames with different columns: the union of the columns in the
# order of their first appearance (those in `last` moved to the end), NA where a
# part lacks a column.
rbindFill <- function(parts, last = NULL) {
  cols <- unique(unlist(lapply(parts, names)))
  cols <- c(setdiff(cols, last), intersect(last, cols))
  do.call(rbind, unname(lapply(parts, function(o) {
    for (cn in setdiff(cols, names(o))) o[[cn]] <- rep_len(NA, nrow(o))   # rep_len: zero-row parts
    o[cols]
  })))
}

# The arguments of INLA::inla(): formula, one likelihood per response (and zero
# equation), control.family, the data list (response matrix + all columns of the
# formula), Ntrials / E / scale (one value per row: the input of the likelihood the row
# belongs to; NA -> 1, e.g. zero equations and responses without the input) and the
# user's options.
inlaArgumentList <- function(inl, options = list()) {
  st <- stackRows(inl)
  resp <- names(inl$families)
  inputs <- unlist(unname(lapply(inl$inputs, as.list)), recursive = FALSE)  # list of Ntrials = col, ...
  data <- as.list(st[setdiff(names(st), c(resp, ".block", unlist(inputs, use.names = FALSE)))])
  surv <- survivalResponses(inl)
  if (length(surv)) {
    # a survival response is an inla.surv object (time, event); the response
    # object is then a list of the response columns (one likelihood each)
    requireInla()
    Y <- as.list(st[resp])
    for (r in surv) {
      col <- inl$responses[[r]]
      Y[[col]] <- INLA::inla.surv(time = Y[[col]], event = st[[inl$inputs[[r]][["event"]]]])
    }
    data[[inl$naming$responseObject]] <- Y
  } else {
    data[[inl$naming$responseObject]] <- as.matrix(st[resp])
  }

  args <- list(formula = inl$formula, family = unname(inl$families), data = data,
               control.family = unname(inl$controlFamily))
  for (a in c("Ntrials", "E", "scale")) {
    cols <- unlist(inputs[names(inputs) == a], use.names = FALSE)
    if (!length(cols)) next
    v <- rep(NA_real_, nrow(st))
    for (cn in cols) v[!is.na(st[[cn]])] <- st[[cn]][!is.na(st[[cn]])]
    v[is.na(v)] <- 1
    args[[a]] <- v
  }

  # fitted values of rows with a missing response: use the link of the row's likelihood
  # (INLA's default would be the identity)
  lik <- unlist(lapply(names(inl$stacks), function(r) {
    s <- inl$stacks[[r]]
    vapply(names(s), function(p) match(if (p == "main") inl$responses[[r]] else
      inl$zeroResponses[[r]], names(inl$families)), 0L) |> rep(vapply(s, nrow, 0L))
  }), use.names = FALSE)
  args$control.predictor <- list(compute = TRUE, link = lik)

  args$verbose <- isTRUE(options$verbose)
  checkInlaOptions(options$inlaOptions, options$threads)
  threads <- resolveThreads(options$threads)
  if (!is.null(threads)) args$num.threads <- threads
  if (length(options$inlaOptions))
    args <- utils::modifyList(args, options$inlaOptions, keep.null = TRUE)
  args
}

# `inlaOptions` must be a named list without the arguments mvrem builds from the model;
# `threads` see checkThreads(). Checked by jsm(), refit() and
# predict() up front, and again when the inla() arguments are assembled.
checkInlaOptions <- function(inlaOptions, threads = NULL) {
  if (is.null(inlaOptions)) inlaOptions <- list()
  if (!isNamedList(inlaOptions))
    stop("`inlaOptions` must be a named list.", call. = FALSE)
  bad <- intersect(names(inlaOptions), reservedInlaArguments)
  if (length(bad))
    stop("`inlaOptions` must not set ", paste0("`", bad, "`", collapse = ", "),
         ": mvrem builds ", if (length(bad) > 1L) "these" else "this", " from the model.",
         call. = FALSE)
  checkThreads(threads)
  invisible()
}

requireInla <- function() {
  if (!hasInla())
    stop("The INLA package is required to fit the model. Install it (see the README) or ",
         "use `fit = FALSE` to build the model without fitting it.", call. = FALSE)
}

# Run INLA for a model object; returns the fit, or the error condition.
runInla <- function(x) {
  requireInla()
  args <- inlaArgumentList(x$inla, x$options)
  tryCatch(do.call(INLA::inla, args), error = function(e) e)
}

# ---------------------------------------------------------------------------
# Results with readable names

firstLine <- function(cond) sub("\n.*", "", conditionMessage(cond))

needFit <- function(x) {
  if (is.null(x$fit))
    stop("The model is not fitted", if (!is.null(x$inlaError)) paste0(" (INLA failed: ",
         firstLine(x$inlaError), ")"), ". Use refit() to fit it.", call. = FALSE)
  invisible(x)
}

# INLA summary table -> data frame with readable column names.
statTable <- function(s) {
  s <- as.data.frame(s, check.names = FALSE)
  s <- s[setdiff(names(s), "kld")]
  names(s) <- ifelse(grepl("quant$", names(s)), paste0("q", sub("quant$", "", names(s))),
                     names(s))
  rownames(s) <- NULL
  s
}

# fixedEffects(): one row per fixed-effect column of every response, in the model's order,
# with INLA's summary matched by the column's INLA name (`.inla`).
fixedEffectsTable <- function(x) {
  needFit(x)
  sf <- x$fit$summary.fixed
  rows <- list()
  for (r in orderResponses(x, names(x$inla$fixed))) {
    cols <- x$inla$fixed[[r]]$columns
    rows[[r]] <- data.frame(response = r, term = unname(cols), .inla = names(cols),
                            stringsAsFactors = FALSE)
  }
  key <- do.call(rbind, rows)
  tab <- statTable(sf[key$.inla, , drop = FALSE])
  out <- cbind(key[c("response", "term")], tab)
  rownames(out) <- NULL
  out
}

# The latent-field hyperparameters the model gives INLA, one row per INLA row name
# (`.inla`) with mvrem's labels: per random block either `Precision for <index>` (a single
# effect; each effect of a `||` block is its own block) or, for a correlated block of k
# effects (one iid<k>d / iidkd field on the first index), `Precision for <index> (component
# j)` and `Rho i:j for <index>`; per scaled() term `Beta for <scaled index>`. This is the ONLY
# place that spells INLA's row names (D1): consumers look rows up by `response`, `block` (the
# position of the block in `inl$random[[response]]`), `effect` / `effect2` or `link`, and read
# `.inla` (the row of `summary.hyperpar` / `marginals.hyperpar`), `internal` (the name of the
# internal marginal, log precision) and `short` ("prec", "rho", "beta").
latentHyperRows <- function(inl) {
  effLabel <- function(e) effectLabel(e, inl$naming)
  rows <- list()
  add <- function(.inla, response, parameter, type, short, link = NA, group = NA, effect = NA,
                  effect2 = NA, block = NA, internal = NA)
    rows[[length(rows) + 1L]] <<- data.frame(
      .inla = .inla, response = response, parameter = parameter, type = type, short = short,
      link = link, group = group, effect = effect, effect2 = effect2, block = block,
      internal = internal, stringsAsFactors = FALSE)
  for (r in names(inl$random)) for (i in seq_along(inl$random[[r]])) {
    b <- inl$random[[r]][[i]]
    k <- length(b$effects)
    if (k == 1L) {
      add(paste0("Precision for ", b$index), r, paste0("precision of ", effectGroupLabel(b$effects, b$group, inl$naming)),
          "random", "prec", group = b$group, effect = b$effects, block = i,
          internal = paste0("Log precision for ", b$index))
      next
    }
    for (j in seq_len(k))
      add(paste0("Precision for ", b$index[1L], " (component ", j, ")"), r,
          paste0("precision of ", effectGroupLabel(b$effects[j], b$group, inl$naming)), "random", "prec",
          group = b$group, effect = b$effects[j], block = i,
          internal = paste0("Log precision for ", b$index[1L], " (component ", j, ")"))
    for (i2 in seq_len(k - 1L)) for (j in (i2 + 1L):k)
      add(paste0("Rho", i2, ":", j, " for ", b$index[1L]), r,
          paste0("correlation of ", effectGroupLabel(b$effects[c(i2, j)], b$group, inl$naming)),
          "random", "rho", group = b$group, effect = b$effects[i2], effect2 = b$effects[j], block = i)
  }
  for (sc in inl$scaled)
    add(paste0("Beta for ", sc$index), sc$out, paste0("gamma (", linkLabel(sc$long, sc$out), ")"),
        "association", "beta", link = sc$long)
  bindRows(rows, data.frame(.inla = character(), response = character(), parameter = character(),
                            type = character(), short = character(), link = character(),
                            group = character(), effect = character(), effect2 = character(),
                            block = integer(), internal = character()))
}

# The rows of latentHyperRows() of one random-effect block (response r, position i in
# `inl$random[[r]]`): its precisions in effect order (`short == "prec"`), then the correlations.
blockHyperRows <- function(known, r, i) known[known$response %in% r & known$block %in% i, , drop = FALSE]

# The value of a FIXED likelihood hyperparameter of response r (INLA's `short.name`), which
# `summary.hyperpar` does not list: from.theta(initial) of the fit's `all.hyper`; NULL if the
# hyperparameter is estimated or does not exist.
fixedFamilyHyper <- function(x, r, short) {
  i <- match(x$inla$responses[[r]], names(x$inla$families))
  for (h in x$fit$all.hyper$family[[i]]$hyper)
    if (isTRUE(h$fixed) && isTRUE(h$short.name == short)) return(as.vector(h$from.theta(h$initial)))
  NULL
}

# Columns of hyperTable() that hyperparameters() drops: the family's / latent term's short
# name (`prec`, `phi`, `rho`, `beta`, ...), the linked longitudinal response of an association,
# the group and effect(s) of a random-effect precision / correlation, INLA's row name.
hyperInternalColumns <- c("short", "link", "group", "effect", "effect2", ".inla")

# Hyperparameters in the order of INLA's summary: likelihoods (non-fixed ones),
# then the latent terms in formula order. Rows are labelled from the model;
# anything not recognised keeps INLA's own name (`type` NA).
hyperTable <- function(x) {
  needFit(x)
  inl <- x$inla
  sh <- x$fit$summary.hyperpar
  n <- if (is.null(sh)) 0L else nrow(sh)
  lab <- data.frame(response = rep(NA_character_, n), parameter = rownames(sh), type = NA_character_,
                    short = NA_character_, link = NA_character_, group = NA_character_,
                    effect = NA_character_, effect2 = NA_character_,
                    .inla = if (n) rownames(sh) else character(), stringsAsFactors = FALSE)
  if (n == 0L) return(cbind(lab, statTable(data.frame(mean = numeric(), sd = numeric()))))
  rows <- rownames(sh)

  # likelihood hyperparameters: one row per non-fixed theta, likelihoods in order. They are
  # matched by POSITION: INLA lists them first, in likelihood order, under descriptive names
  # ("Precision for the Gaussian observations[2]") that do not name the response.
  fam <- x$fit$all.hyper$family
  respOfColumn <- c(stats::setNames(names(inl$responses), inl$responses))
  k <- 0L
  for (i in seq_along(fam)) {
    col <- names(inl$families)[i]
    for (h in fam[[i]]$hyper) {
      if (isTRUE(h$fixed)) next
      # survival likelihoods list ten cure-model coefficients (beta1, ...) that are only
      # estimated with INLA's `cure` covariates, which mvrem does not use
      r <- unname(respOfColumn[col])
      if (isTRUE(inl$survival[r]) && !is.null(h$short.name) && grepl("^beta[0-9]+$", h$short.name)) next
      k <- k + 1L
      if (k > n) break
      short <- if (is.null(h$short.name)) "hyper" else h$short.name
      lab$response[k] <- r
      lab$type[k] <- "family"
      lab$short[k] <- short
      lab$parameter[k] <- if (short == "prec") {
        if (identical(unname(inl$families[col]), "gaussian")) "residual precision" else "precision"
      } else short
    }
  }

  # latent-field hyperparameters: INLA's row names are known from the model (latentHyperRows())
  known <- latentHyperRows(inl)
  j <- which(is.na(lab$type) & rows %in% known$.inla)
  k <- match(rows[j], known$.inla)
  cols <- c("response", "parameter", "type", "short", "link", "group", "effect", "effect2")
  lab[j, cols] <- known[k, cols]
  out <- cbind(lab, statTable(sh))
  rownames(out) <- NULL
  out
}
