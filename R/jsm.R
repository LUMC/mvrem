#' @import dplyr
#' @import tidyr
#' @import rlang
#' @import lme4
#' @import R6
#' @import INLA
#' @importFrom purrr imap imap_chr list_flatten list_transpose map map_chr map_dfr discard
#' @importFrom stats setNames
#' @importFrom utils capture.output

#'
#' @title Ensure Named Formulas
#' @description Ensures that a list of formulas has names corresponding to their left-hand side variables.
#' @param forms A formula or a list of formulas.
#' @return A named list of formulas, where the names correspond to the left-hand side variables.
#'
.forms2namedForms <- function(forms) {
  if (rlang::is_formula(forms)) forms <- list(forms)

  stopifnot(is.list(forms))
  stopifnot(all(sapply(forms, inherits, "formula")))

  lhsNames <- sapply(forms, \(f) all.vars(f[[2L]]))
  if (!is.null(names(forms))) {
    missingNames <- names(forms) == ""
    names(forms)[missingNames] <- lhsNames[missingNames]
  } else {
    names(forms) <- lhsNames
  }
  forms
}

#'
#' @title Ensure Named Data
#' @description Ensures that a list of data frames has names corresponding to the formulas.
#' @details
#' If a single data frame is provided, it will be used for all formulas.
#' If a named list of data frames is provided, it must have the same length as the list of formulas, and the names must correspond to the formulas.
#' If a non-named list of data frames is provided, the names will be assigned based on the formulas keeping the order consistent.
#'
#' @param data A data frame or a list of data frames.
#' @param namedForms A named list of formulas.
#' @return A named list of data frames, where the names correspond to the formulas.
#'
.data2namedData <- function(data, namedForms) {
  if (is.data.frame(data)) {
    # single data frame: use the same data frame for all formulas
    return (namedForms |> lapply(\(f) data))
  } else if (is.list(data)) {
    # list of data frames: use the corresponding data frame for each formula
    stopifnot(length(data) == length(namedForms))
    stopifnot(all(sapply(data, is.data.frame)))

    if (!is.null(names(data))) {
      stopifnot(all(names(data) %in% names(namedForms)))
      return (data[names(namedForms)])
    } else {
      return (setNames(data, names(namedForms)))
    }
  } else {
    stop("data must be a data frame or a list of data frames")
  }
}

.families2named_families <- function(families, namedForms, defaultFamily) {
  if (is.null(families)) families <- rep(defaultFamily, length(namedForms))
  stopifnot(is.character(families))

  lhsNames <- setNames(nm = names(namedForms))
  if (!is.null(names(families))) {
    stopifnot(all(names(families) %in% lhsNames))
    families <- sapply(lhsNames, \(n) {
      ifelse(n %in% names(families), families[n], defaultFamily)
    })
  } else {
    stopifnot(length(families) == length(lhsNames))
    names(families) <- lhsNames
  }
  families
}

# ------------------------------------------------------------------------

# ----- vars -----

.isIntercept <- \(lmeVars) {
  lmeVars == "(Intercept)"
}

.renameIntercept <- \(lmeVars, pars) {
  purrr::map_chr(lmeVars, \(lmeVar) {
    if (lmeVar == "(Intercept)") pars$interceptId else lmeVar
  })
}

.mergeLmeVars <- \(lmeVars, pars) {
  paste(lmeVars, collapse = pars$sep)
}

.lmeVar2zeroEqLmeVar <- \(lmeVar, pars) {
  paste(lmeVar, pars$zeroBlockSuffix, sep = pars$sep)
}

.lmeVar2inlaVar <- \(lmeVar, prefixes, pars) {
  paste(c(prefixes, lmeVar), collapse = pars$sep)
}

# ----- cols -----

.mergeCols <- \(fr, prefixes, pars) {
  col <- do.call(paste, c(fr, sep = pars$valSep))
  if (!is.null(prefixes)) {
    col <- sapply(col, \(v) paste(c(prefixes, v), collapse = pars$valSep))
  }
  col
}

.getAllBlockColGrps <- \(blocks) {
  blocks |>
    purrr::map(\(block) block$blockColGrp) |>
    dplyr::bind_rows()
}

# ----- levels -----

.addLevels <- \(block, keyInlaVar, colInlaVar, levels, pars) {
  names(levels) <- NULL
  levels <- c(block$levels[[keyInlaVar]], levels)
  levels <- levels[!duplicated(levels)]
  block$levels[[keyInlaVar]] <- levels
  block
}

.getAllLevels <- \(blocks) {
  blocks |>
    purrr::map(\(block) block$levels) |>
    purrr::list_transpose() |>
    purrr::map(\(x) unlist(x)) |>
    purrr::map(\(x) setNames(x[!duplicated(x)], nm = NULL))
}

.getLevelsExpr <- \(keyInlaVar, pars) {
  paste0(pars$levelsInlaVar, "$", keyInlaVar)
}

# ----- single block operations -----

.initBlocks <- function(
  namedForms, namedData, namedFamilies, isOut, na.action
) {
  purrr::imap(namedForms, \(f, n) {
    d <- as.data.frame(namedData[[n]])
    stopifnot(all(all.vars(f) %in% names(d)))

    if (is.null(rownames(d))) {
      rownames(d) <- sprintf("%s_%d", n, seq_len(nrow(d)))
    } else {
      stopifnot(!any(duplicated(rownames(d))))
    }

    list(
      formula = deparse(f),
      data = d,
      family = namedFamilies[[n]],
      isOut = isOut,
      parsedFormula = lme4::glFormula(f, data = d, na.action = na.action),
      outVar = as.character(f[[2L]]),
      fr = list(),
      colGrp2cols = list(
        "L" = character(), "L0" = character(),
        "R" = character(), "R0" = character()
      ),
      terms = list(),
      levels = list(),
      familiesFArgs = list(),
      blockColGrp = list()
    )
  })
}

.addCol <- \(block, inlaVar, col, colGrps) {
  names(col) <- NULL
  block$fr[[inlaVar]] <- col
  for (cg in colGrps) {
    if (!(inlaVar %in% block$colGrp2cols[[cg]])) {
      block$colGrp2cols[[cg]] <- c(block$colGrp2cols[[cg]], inlaVar)
    }
  }
  block
}

.addFamilyFArgs <- \(block, inlaVar, family, fArgs) {
  block$familiesFArgs <- block$familiesFArgs |>
    c(list(
      list(
        inlaVar = inlaVar, family = family, fArgs = fArgs
      )
    ))
  block
}

.addBlockColGrp <- \(block, blockId, colGrp, inlaVar) {
  block$blockColGrp <- block$blockColGrp |> 
    c(
      list(
        list(.blockId = blockId, .colGrp = colGrp, .inlaVar = inlaVar)
      )
    )
  block
}

.addTerm <- \(block, inlaVar) {
  block$terms <- c(block$terms, inlaVar)
  block
}

.addFTerm <- \(block, grpInlaVar, inlaVar, fArgs) {
  args <- c(
    grpInlaVar, inlaVar,
    purrr::imap_chr(
      fArgs, \(v, k) paste0(
        k, " = ", deparse(v, width.cutoff = 500L),
        collapse = ""
      )
    )
  )
  .fOf <- \(args) paste0("f(", paste(args, collapse = ", "), ")")
  .addTerm(block, .fOf(args))
}

.buildFArgs <- function(defaultArgs, inlaVar, userPars) {
  args <- defaultArgs
  if (!is.null(userPars) && (inlaVar %in% names(userPars))) {
    args <- utils::modifyList(args %||% list(), userPars[[inlaVar]])
  }

  args
}

# ----- multiple blocks operations -----

.initL <- function(
  outForms, longForms, outData, longData, outFamily, longFamily,
  pars
) {
  if (is.null(longData)) {
    longData <- outData
  }

  l <- list()
  l$pars <- pars

  outNamedForms <- .forms2namedForms(outForms)
  outNamedData <- .data2namedData(outData, outNamedForms)
  outNamedFamilies <- .families2named_families(
    outFamily, outNamedForms, defaultFamily = "gaussian"
  )
  outBlocks <- .initBlocks(
    namedForms = outNamedForms, namedData = outNamedData,
    namedFamilies = outNamedFamilies,
    isOut = TRUE, na.action = pars$na.action
  )

  longNamedForms <- .forms2namedForms(longForms)
  longNamedData <- .data2namedData(longData, longNamedForms)
  longNamedFamilies <- .families2named_families(
    longFamily, longNamedForms, defaultFamily = "gaussian"
  )
  longBlocks <- .initBlocks(
    namedForms = longNamedForms, namedData = longNamedData,
    namedFamilies = longNamedFamilies,
    isOut = FALSE, na.action = pars$na.action
  )

  stopifnot(all(!duplicated(c(names(outBlocks), names(longBlocks)))))

  l$blocks <- c(outBlocks, longBlocks)
  l$longBlocksIds <- names(longBlocks)
  l$outBlocksIds <- names(outBlocks)
  l
}

.addFixedEffects <- function(l, pars) {
  l$blocks <- purrr::imap(l$blocks, \(block, blockId) {
    pf <- block$parsedFormula
    cnms <- pf$reTrms$cnms
    fr <- pf$fr
    X <- as.data.frame(pf$X)
    Xfr <- X |> cbind(fr[setdiff(colnames(fr), colnames(X))])
    grpsR <- if (block$isOut) c("R") else c("R", "R0")

    #block$lmeFr <- Xfr

    # ----- outcome variable (LHS) -----
    outInlaVar <- .lmeVar2inlaVar(block$outVar, blockId, pars)
    block <- .addCol(block, outInlaVar, fr[[block$outVar]], c("L"))
    block <- .addBlockColGrp(block, blockId, "R", outInlaVar)
    block <- .addFamilyFArgs(
      block, outInlaVar, block$family, list()
    )

    if (!block$isOut) {
      outZeroEqInlaVar <-
        .lmeVar2inlaVar(.lmeVar2zeroEqLmeVar(block$outVar, pars), blockId, pars)
      block <- .addCol(block, outZeroEqInlaVar, fr[[block$outVar]] * 0, c("L0"))
      block <- .addBlockColGrp(block, blockId, "R0", outZeroEqInlaVar)
      block <- .addFamilyFArgs(
        block, outZeroEqInlaVar, pars$zeroBlockFamily, pars$zeroBlockFArgs
      )
    }

    # ----- fixed effect variables (RHS) -----
    fixEfVars <- colnames(X)
    for (v in fixEfVars) {
      inlaV <- .lmeVar2inlaVar(.renameIntercept(v, pars), blockId, pars)
      block <- .addCol(block, inlaV, Xfr[[v]], grpsR)
      block <- .addTerm(block, inlaV)
    }
    for (v in sort(setdiff(unique(unlist(cnms)), fixEfVars))) {
      inlaV <- .lmeVar2inlaVar(.renameIntercept(v, pars), blockId, pars)
      col <- if (.isIntercept(v)) rep(1.0, nrow(Xfr)) else Xfr[[v]]
      block <- .addCol(block, inlaV, col, grpsR)
    }

    block
  })
  l
}

.addRandomEffects <- function(l, pars, userPars) {
  l$blocks <- purrr::imap(l$blocks, \(block, blockId) {
    pf <- block$parsedFormula
    cnms <- pf$reTrms$cnms
    fr <- pf$fr
    #X <- as.data.frame(pf$X)
    #Xfr <- X |> cbind(fr[setdiff(colnames(fr), colnames(X))])
    grpsR <- if (block$isOut) c("R") else c("R", "R0")

    # ----- random effect variables (RHS) -----
    for (i in seq_along(cnms)) {
      rndEfVars <- cnms[[i]]
      rndEfId <- names(cnms)[[i]]
      rndEfGrpVars <- sort(strsplit(rndEfId, ":")[[1]])

      grpInlaVar <- .lmeVar2inlaVar(
        .mergeLmeVars(rndEfGrpVars, pars), blockId, pars
      )
      if (length(rndEfVars) == 1L) {
        if (.isIntercept(rndEfVars)) {
          # ----- add intercept col -----
          grpInla1Var <- .lmeVar2inlaVar(
            .mergeLmeVars(rndEfGrpVars, pars),
            c(blockId, pars$interceptId), pars
          )
          grpInla1Col <- .mergeCols(
            fr[rndEfGrpVars], c(blockId, pars$interceptId), pars
          )
          block <- .addCol(block, grpInla1Var, grpInla1Col, grpsR)

          # ----- add levels -----
          lvs <- grpInla1Col
          lvs <- lvs[!duplicated(lvs)]
          lvsInlaVar <- .lmeVar2inlaVar(
            .mergeLmeVars(rndEfGrpVars, pars),
            c(blockId, pars$interceptId), pars
          )
          block <- .addLevels(block, lvsInlaVar, lvsInlaVar, lvs, pars)

          # ----- add intercept term -----
          args <- .buildFArgs(pars$rndEfInterceptFArgs, grpInla1Var, userPars)
          args$values <- str2lang(.getLevelsExpr(lvsInlaVar, pars))
          block <- .addFTerm(block, grpInla1Var, NULL, args)
        } else {
          # ----- add slope col -----
          rndEfTVar <- rndEfVars[[1L]]
          grpInlaTVar <- .lmeVar2inlaVar(
            .mergeLmeVars(rndEfGrpVars, pars), c(blockId, rndEfTVar), pars
          )
          grpInlaTCol <- .mergeCols(
            fr[rndEfGrpVars], c(blockId, rndEfTVar), pars
          )
          block <- .addCol(block, grpInlaTVar, grpInlaTCol, grpsR)

          # ----- add levels -----
          lvs <- grpInlaTCol
          lvs <- lvs[!duplicated(lvs)]
          lvsInlaVar <- .lmeVar2inlaVar(
            .mergeLmeVars(rndEfGrpVars, pars), c(blockId, rndEfTVar), pars
          )
          block <- .addLevels(block, lvsInlaVar, lvsInlaVar, lvs, pars)

          # ----- add slope term -----
          args <- .buildFArgs(pars$rndEfSlopeFArgs, grpInlaTVar, userPars)
          args$values <- str2lang(.getLevelsExpr(lvsInlaVar, pars))
          slopeInlaVar <- .lmeVar2inlaVar(rndEfTVar, blockId, pars)
          block <- .addFTerm(block, grpInlaTVar, slopeInlaVar, args)
        }
      } else if (length(rndEfVars) == 2L && .isIntercept(rndEfVars[[1L]])) {
        # ----- add intercept col -----
        grpInla1Var <- .lmeVar2inlaVar(
          .mergeLmeVars(rndEfGrpVars, pars), c(blockId, pars$interceptId), pars
        )
        grpInla1Col <- .mergeCols(
          fr[rndEfGrpVars], c(blockId, pars$interceptId), pars
        )
        block <- .addCol(block, grpInla1Var, grpInla1Col, grpsR)

        # ----- add slope col -----
        rndEfTVar <- rndEfVars[[2L]]
        grpInlaTVar <- .lmeVar2inlaVar(
          .mergeLmeVars(rndEfGrpVars, pars), c(blockId, rndEfTVar), pars
        )
        grpInlaTCol <- .mergeCols(
          fr[rndEfGrpVars], c(blockId, rndEfTVar), pars
        )
        block <- .addCol(block, grpInlaTVar, grpInlaTCol, grpsR)

        # ----- add levels -----
        lvs <- c(grpInla1Col, grpInlaTCol)
        lvs <- lvs[!duplicated(lvs)]
        lvsInlaVar <- .lmeVar2inlaVar(
          .mergeLmeVars(rndEfGrpVars, pars),
          c(blockId, pars$interceptId, rndEfTVar),
          pars
        )
        block <- .addLevels(block, lvsInlaVar, lvsInlaVar, lvs, pars)

        # ----- add intercept term -----
        args <- .buildFArgs(pars$rndEfInterceptSlopeFArgs, grpInla1Var, userPars)
        args$values <- .getLevelsExpr(lvsInlaVar, pars) |>
          str2lang()
        args$n <- paste0("length(", .getLevelsExpr(lvsInlaVar, pars), ")") |>
          str2lang()
        block <- .addFTerm(block, grpInla1Var, NULL, args)

        # ----- add slope term -----
        args <- .buildFArgs(pars$rndEfInterceptSlopeCopyFArgs, grpInlaTVar, userPars)
        args$copy <- grpInla1Var
        slopeInlaVar <- .lmeVar2inlaVar(rndEfTVar, blockId, pars)
        block <- .addFTerm(block, grpInlaTVar, slopeInlaVar, args)
      } else {
        stop("Unsupported random effect structure")
      }
    }
    block
  })
  l
}

.addCopyVariables <- \(l, blockId2forcedLmeVars, pars, userPars) {
  # ----- identify lme variables for each block -----
  l$blockId2lmeVars <- purrr::imap(l$blocks, \(block, blockId) {
    lmeVars <- c(
      colnames(block$parsedFormula$X) |> purrr::discard(.isIntercept),
      colnames(block$parsedFormula$fr)
    )
    lmeVars[!duplicated(lmeVars)]
  })

  # ----- helper fun: calc shared key lme vars between two blocks -----
  .calcKeyLmeVars <- \(longBlockId, outBlockId) {
    sort(intersect(
      c(
        l$blockId2lmeVars[[longBlockId]],
        blockId2forcedLmeVars[[longBlockId]]
      ),
      c(
        l$blockId2lmeVars[[outBlockId]],
        blockId2forcedLmeVars[[outBlockId]]
      )
    ))
  }

  # ----- helper fun: build frame with all key lme vars -----
  .getExFrX <- \(block, keyLmeVars) {
    exFrX <- data.frame(block$parsedFormula$fr, block$parsedFormula$X)
    missingKeyLmeVars <- keyLmeVars[!(keyLmeVars %in% colnames(exFrX))]
    if (length(missingKeyLmeVars) > 0L) {
      exFrX <- cbind(
        exFrX,
        block$data[
          rownames(block$parsedFormula$fr), missingKeyLmeVars,
          drop = FALSE
        ]
      )
    }
    exFrX
  }

  # ----- to long blocks add sources -----
  l$blocks[l$longBlocksIds] <- l$blocks[l$longBlocksIds] |>
    purrr::imap(\(longBlock, longBlockId) {
      frX <- data.frame(longBlock$parsedFormula$fr, longBlock$parsedFormula$X)
      wgtInlaVar <- .lmeVar2inlaVar(
        c(),
        c(longBlockId, pars$weightColMark),
        pars
      )
      wgtInlaCol <- frX[[longBlock$outVar]] * 0. - 1.
      longBlock <- .addCol(longBlock, wgtInlaVar, wgtInlaCol, c("R0"))

      for (outBlockId in l$outBlocksIds) {
        keyLmeVars <- .calcKeyLmeVars(longBlockId, outBlockId)
        frX <- .getExFrX(longBlock, keyLmeVars)
        inlaVar <- .lmeVar2inlaVar(
          .mergeLmeVars(keyLmeVars, pars),
          c(longBlockId, pars$toColMark, outBlockId),
          pars
        )
        inlaCol <- .mergeCols(
          frX[keyLmeVars],
          c(longBlockId, pars$toColMark, outBlockId),
          pars
        )
        longBlock <- .addCol(longBlock, inlaVar, inlaCol, c("R0"))
        longBlock <- .addLevels(longBlock, inlaVar, inlaVar, inlaCol, pars)
        args <- .buildFArgs(pars$fromFArgs, inlaVar, userPars)
        args$values <- str2lang(.getLevelsExpr(inlaVar, pars))
        longBlock <- .addFTerm(longBlock, inlaVar, wgtInlaVar, args)
      }

      longBlock
    })

  # ----- to out blocks add destinations -----
  l$blocks[l$outBlocksIds] <- l$blocks[l$outBlocksIds] |>
    purrr::imap(\(outBlock, outBlockId) {
      for (longBlockId in l$longBlocksIds) {
        keyLmeVars <- .calcKeyLmeVars(longBlockId, outBlockId)
        frX <- .getExFrX(outBlock, keyLmeVars)
        fromInlaVar <- .lmeVar2inlaVar(
          .mergeLmeVars(keyLmeVars, pars),
          c(longBlockId, pars$toColMark, outBlockId),
          pars
        )
        toInlaVar <- .lmeVar2inlaVar(
          .mergeLmeVars(keyLmeVars, pars),
          c(outBlockId, pars$fromColMark, longBlockId),
          pars
        )
        inlaCol <- .mergeCols(
          frX[keyLmeVars],
          c(longBlockId, pars$toColMark, outBlockId),
          pars
        )
        outBlock <- .addCol(outBlock, toInlaVar, inlaCol, c("R"))
        outBlock <- .addLevels(outBlock, fromInlaVar, toInlaVar, inlaCol, pars)
        args <- .buildFArgs(pars$toFArgs, toInlaVar, userPars)
        args$copy <- fromInlaVar
        outBlock <- .addFTerm(outBlock, toInlaVar, NULL, args)
      }
      outBlock
    })

  # ----- build matching stats -----
  l$copyVarStats <- purrr::map(l$longBlocksIds, \(longBlockId) {
    purrr::map(l$outBlocksIds, \(outBlockId) {
      keyLmeVars <- .calcKeyLmeVars(longBlockId, outBlockId)
      srcInlaVar <- .lmeVar2inlaVar(
        .mergeLmeVars(keyLmeVars, pars),
        c(longBlockId, pars$toColMark, outBlockId),
        pars
      )
      dstInlaVar <- .lmeVar2inlaVar(
        .mergeLmeVars(keyLmeVars, pars),
        c(outBlockId, pars$fromColMark, longBlockId),
        pars
      )
      srcCol <- l$blocks[[longBlockId]]$fr[[srcInlaVar]]
      dstCol <- l$blocks[[outBlockId]]$fr[[dstInlaVar]]

      .cntClass <- \(col, one, many) {
        dplyr::tibble(id = col) |>
          dplyr::count(id, name = "cnt") |>
          dplyr::mutate(
            cnt = dplyr::if_else(cnt > 1L, many, one)
          )
      }

      dplyr::full_join(
        .cntClass(srcCol, one = "1", many = "N") |> dplyr::rename(src = cnt),
        .cntClass(dstCol, one = "1", many = "M") |> dplyr::rename(dst = cnt),
        by = c("id")
      ) |>
        tidyr::replace_na(list(src = "0", dst = "0")) |>
        dplyr::count(src, dst, name = "cnt") |>
        dplyr::mutate(
          longBlockId = longBlockId,
          outBlockId = outBlockId,
          relation = paste0(src, "-", dst),
          matchedVars = paste(keyLmeVars, collapse = ",")
        ) |>
        dplyr::select(longBlockId, outBlockId, relation, cnt, matchedVars)
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()

  if (any(l$copyVarStats$relation %in% c("1-M", "N-1", "N-M"))) {
    txt <- paste(capture.output(print(l$copyVarStats)), collapse = "\n")
    warning(paste(
      c(
        "Some data rows are not uniquely identified by the variables in their formula,",
        "leading to multimapping relations (1-M, N-1, N-M) between long and out rows.",
        "Use e.g.: blockId2forcedLmeVars = list(x1 = c('t'), ...) to force usage of additional match variables",
        "for specific equation blocks (here, forces usage of variable 't' in block for 'x1').",
        "",
        "Count summary:",
        txt
      ),
      collapse = "\n"
    ))
  }

  l
}

# ------------------------------------------------------------------------

#'
#' @title Default Parameters for JSM unless Overridden by User
#'
#' @description
#' This list contains the default parameters used by the JSM package,
#' unless the user overriddes (some) elements
#' in the `pars` argument of the `jsm` function.
#'
#' @export
#'
#' @seealso \code{\link{jsm}}
defaultPars <- list(
  na.action = "na.pass",

  # ----- utility variables and composition parts -----
  #'
  #' @param lhsInlaVar The name used for the INLA LHS variable
  #' (the matrix with the outcome variables).
  lhsInlaVar = "MX",
  #'
  #' @param outcomeIdInlaVar The name used for the variable holding the
  #' vector indicating which outcome each row corresponds to
  #' (see INLA `control.predictor` documentation).
  #'
  outcomeIdInlaVar = "OUTCOME",
  #'
  #' @param levelsInlaVar The name used for the variable holding
  #' a list with the levels generated for each random effect in the model.
  #' Used in the `f(...)` terms of the INLA formula.
  #'
  levelsInlaVar = "LEVELS",
  #'
  #' @param interceptId The string used to construct parts of identifiers
  #' representing various intercept terms in the model.
  #'
  interceptId = "1",
  #'
  #' @param slopeId The string used to construct parts of identifiers
  #' representing various slope terms in the model.
  #'
  slopeId = "S",

  # ----- separators -----
  #'
  #' @param sep The separator used in constructing variable identifiers.
  #'
  sep = "_",
  #'
  #' @param valSep The separator used for values (levels)
  #' of various variables.
  #'
  valSep = "|",

  # ----- zero block -----
  #'
  #' @param zeroBlockSuffix The suffix added to the name of a variable
  #' used in the "zero block" in the model.
  #'
  zeroBlockSuffix = "0eq",
  #'
  #' @param zeroBlockFamily The INLA likelihood distribution family
  #' used for the "zero block" in the model. See INLA `family`
  #' documentation.
  #'
  zeroBlockFamily = "gaussian",
  #'
  #' @param zeroBlockFArgs The list of arguments specifying additional
  #' prior hypermarameters for the likelihood in the INLA.
  #' See INLA `control.family` documentation.
  #'
  zeroBlockFArgs = list(
    hyper = list(
      prec = list(initial = 1, fixed = TRUE)
    )
  ),

  # ----- random effects arguments -----
  #'
  #' @param rndEfInterceptFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` terms
  #' for the intercept-only random effects.
  #'
  rndEfInterceptFArgs = list( # intercept only
    model = "iid",
    values = "???",
    hyper = list()
  ),
  #'
  #' @param rndEfSlopeFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` terms
  #' for the slope-only random effects.
  #'
  rndEfSlopeFArgs = list( # slope only
    model = "iid",
    values = "???",
    hyper = list()
  ),
  #'
  #' @param rndEfInterceptSlopeFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` intercept term
  #' for the intercept-and-slope random effects.
  #'
  rndEfInterceptSlopeFArgs = list( # intercept and slope, intercept part
    model = "iid2d",
    values = "???",
    n = NA_integer_,
    hyper = list(
      theta1 = list(
        prior = "wishart2d", param = c(3, 1.000000, 0.004000, 0)
      )
    )
  ),
  #'
  #' @param rndEfInterceptSlopeCopyFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` slope term
  #' for the intercept-and-slope random effects.
  #'
  rndEfInterceptSlopeCopyFArgs = list( # intercept and slope, slope part
    copy = "???"
  ),

  # ----- copy vars from long to out blocks -----
  #'
  #' @param fromColMark The string used to construct parts of identifiers
  #' representing the "scaled" variable in the longitudinal "zero block".
  #'
  fromColMark = "FROM",
  #'
  #' @param toColMark The string used to construct parts of identifiers
  #' representing the "scaled" variable in the outcome terms of the model.
  #'
  toColMark = "TO",
  #'
  #' @param weightColMark The string used to construct parts of identifiers
  #' representing the column of `1s` weights.
  #'
  weightColMark = "WEIGHT",
  #'
  #' @param fromFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` terms
  #' for the longitudinal "zero block" terms of the
  #' model "scaled" variable.
  #'
  fromFArgs = list(
    model = "iid",
    values = "???",
    hyper = list()
  ),
  #'
  #' @param toFArgs The list specifying arguments of
  #' prior hypermarameters in INLA formula `f(...)` terms
  #' for the outcome terms of the model "scaled" variable.
  #'
  toFArgs = list(
    copy = "???",
    hyper = list(
      beta = list(
        fixed = FALSE,
        prior = "normal",
        initial = 0,
        param = c(0, 1 / (5^2)) # mean, precision = 1/(sd^2)
      )
    )
  )
)

# ------------------------------------------------------------------------

#'
#' A class with methods for translation of Joint Spatial Model (JSM)
#' longitudinal and outcome formulas/data into a format suitable for INLA,
#' as well as managing the model fitting and retrieval of INLA results.
#'
#' Use the \code{\link{jsm}} function to create a new instance of this class.
#'
JSMMdl <- R6::R6Class(
  classname = "JSMMdl",
  private = list(
    fit = TRUE,
    l = list(),
    inlaRes = NULL
  ),
  public = list(
    #'
    #' @description
    #' Create a new JSM model object, used by the `jsm` function.
    #' @returns A new `JSMMdl` object.
    #'
    initialize = function(
      ..., fit = TRUE, blockId2forcedLmeVars = list(),
      pars = NULL,
      userPars = list()
    ) {
      if (is.null(pars)) pars <- list()
      stopifnot(is.list(pars))
      pars <- utils::modifyList(defaultPars, pars)

      private$fit <- fit
      private$inlaRes <- NULL

      l <- .initL(..., pars = pars)
      l <- .addFixedEffects(l, pars = pars)
      l <- .addRandomEffects(l, pars = pars, userPars = userPars)
      l <- .addCopyVariables(
        l, blockId2forcedLmeVars = blockId2forcedLmeVars,
        pars = pars, userPars = userPars
      )
      private$l <- l

      if (private$fit) {
        warning("Auto fitting is disabled. Use callInla() method.")
      }
    },
    #'
    #' @description
    #' Get the internal list representing the state model. Not for general use.
    #' It will be disabled in future versions of the package.
    #' @returns The internal list `l` containing the model components.
    #'
    getL = function() {
      private$l
    },
    #'
    #' @description
    #' Get the families for each INLA outcome variable.
    #' @returns A named character vector mapping INLA outcome variables
    #'   to their families (e.g. `c(y_y = "gaussian", x_x = "binomial")`).
    #'   Note: order is not guaranteed to correspond
    #'   to the order of INLA outcome variables.
    #'
    getFamilies = function(defaultFamily = "gaussian") {
      l <- self$getL()
      familesFArgs <- l$blocks |>
        purrr::map(\(bl) bl$familiesFArgs) |>
        purrr::list_flatten()
      inlaVars <- familesFArgs |>
        purrr::map_chr(\(ff) ff$inlaVar)
      families <- familesFArgs |>
        purrr::map_chr(\(ff) ff$family %||% defaultFamily)
      setNames(families, nm = inlaVars)
    },
    #'
    #' @description
    #' Get the INLA `control.family` list for each INLA outcome variable.
    #' @returns A named list mapping INLA outcome variables to their
    #'   control family lists. Note: order is not guaranteed to correspond
    #'   to the order of INLA outcome variables.
    #'
    getControlFamilies = function(defaultFArgs = list()) {
      l <- self$getL()
      familesFArgs <- l$blocks |>
        purrr::map(\(bl) bl$familiesFArgs) |>
        purrr::list_flatten()
      inlaVars <- familesFArgs |>
        purrr::map_chr(\(ff) ff$inlaVar)
      controlFamilies <- familesFArgs |>
        purrr::map(\(ff) ff$fArgs %||% defaultFArgs)
      setNames(controlFamilies, nm = inlaVars)
    },
    getFrList = function(blockIds = NULL, colGrps = NULL) {
      l <- self$getL()
      blocks <- l$blocks
      if (!is.null(blockIds)) blocks <- blocks[blockIds]

      purrr::imap(blocks, \(block, blockId) {
        if (is.null(colGrps)) {
          cols <- names(block$fr)
        } else {
          cols <- unlist(block$colGrp2cols[colGrps], use.names = FALSE)
        }
        stopifnot(all(cols %in% names(block$fr)))
        if (length(cols) == 0) {
          return(NULL)
        } else {
          block$fr[cols]
        }
      })
    },
    getFrs = function(colGrps, blockIds = NULL, n = +Inf, addBlockId = FALSE) {
      l <- self$getL()
      if (is.null(blockIds)) blockIds <- names(l$blocks)
      d <- purrr::imap(l$blocks[blockIds], \(block, blockId) {
        purrr::map_dfr(colGrps, \(colGrp) {
          ll <- self$getFrList(blockId, colGrps = colGrp)
          names(ll) <- NULL
          dd <- ll |> as.data.frame() |> head(n)
          if (addBlockId && nrow(dd) > 0L) {
            dd$.blockId <- blockId
            dd$.colGrp <- colGrp
          }
          dd
        }) |>
          dplyr::bind_rows()
      }) |>
        dplyr::bind_rows()
      rownames(d) <- NULL
      d
    },
    getLevels = function() {
      l <- self$getL()
      purrr::map(l$blocks, \(bl) bl$levels)
    },
    #'
    #' @description
    #' Get the INLA formula text for the model, optionally including comments and a preview of the data.
    #' @param addComments Logical, whether to include comments in the formula text.
    #' @param n Number of rows of data to include in the comments.
    #' @returns A character string representing the formula text.
    #'
    getFormulaText = function(addComments = TRUE, n = 3L) {
      l <- self$getL()
      allTerms <- c(
        list(list(.minusOne = list(TRUE, paste0(l$pars$lhsInlaVar, " ~ -1")))),
        purrr::imap(l$blocks, \(bl, blockId) {
          c(
            if (addComments) {
              ll <- list(
                list(FALSE, paste0("   #")),
                list(
                  FALSE,
                  paste0("   # ----- ", blockId, ": ", bl$formula, " -----")
                ),
                list(FALSE, paste0("   #"))
              )
              if (n > 0L) {
                d <- dplyr::bind_cols(
                  self$getFrs(
                    c("L", "L0"), blockIds = blockId, n = n
                  ),
                  self$getFrs(
                    c("R", "R0"), blockIds = blockId, n = n
                  )
                )
                dTxt <- d |>
                  print(row.names = FALSE) |>
                  capture.output() |>
                  lapply(\(x) list(FALSE, paste0("   # ", x)))
                ll <- c(ll, dTxt)
                ll <- c(ll, list(list(FALSE, paste0("   #"))))
              }
              ll
            },
            purrr::map(bl$terms, \(t) list(TRUE, t))
          )
        })
      ) |>
        unlist(recursive = FALSE)
      addPlusIdxs <- allTerms |> sapply(\(x) x[[1]]) |> which() |> head(-1L)
      for (i in addPlusIdxs) {
        allTerms[[i]][[2]] <- paste0(allTerms[[i]][[2]], " +")
      }
      paste0(sapply(allTerms, \(x) x[[2]]), collapse = "\n")
    },
    #'
    #' @description
    #' Get the INLA formula for the model as an R formula object.
    #'   This formula is used when `callInla()` is called.
    #'
    getFormula = function(...) {
      as.formula(self$getFormulaText(...))
    },
    #'
    #' @description
    #' Get the data list for the INLA model, including the left-hand side,
    #'   right-hand side, and levels variables.
    #'
    getData = function() {
      l <- self$getL()
      lhs <- self$getFrs(c("L", "L0"), addBlockId = FALSE)

      allBlockColGrps <- .getAllBlockColGrps(l$blocks)
      rhs <- self$getFrs(c("R", "R0"), addBlockId = TRUE) |>
        dplyr::left_join(allBlockColGrps, by = c(".blockId", ".colGrp"))
      rhs[[l$pars$outcomeIdInlaVar]] <- rhs[[".inlaVar"]]
      rhs <- rhs |>
        dplyr::select(-.blockId, -.colGrp, -.inlaVar)

      ret <- list()
      ret[[l$pars$lhsInlaVar]] <- lhs
      ret <- c(ret, rhs)

      ret[[l$pars$levelsInlaVar]] <- .getAllLevels(l$blocks)

      ret
    },
    #'
    #' Call INLA to fit the model.
    #'
    #' @description
    #' Fit the Joint Spatial Model (JSM) using INLA with the formulas
    #' obtained with `getFormula()`, and data provided by `getData()`.
    #' Moreover, it uses other methods of the class to get the necessary
    #' components for the INLA call, such as families and control settings.
    #'
    #' @return
    #' Currently the function returns the result of the INLA call.
    #'
    callInla = function(...) {
      l <- self$getL()

      data <- self$getData()
      lhsInlaVars <- colnames(data[[l$pars$lhsInlaVar]])
      linkCol <- match(data[[l$pars$outcomeIdInlaVar]], lhsInlaVars)

      families <- self$getFamilies()[lhsInlaVars]
      controlFamilies <- self$getControlFamilies()[lhsInlaVars]

      private$inlaRes <- INLA::inla(
        formula = self$getFormula(),
        data = data,
        control.compute = list(
          dic = TRUE, cpo = TRUE, waic = TRUE
        ),
        family = families,
        control.family = setNames(controlFamilies, NULL), # INLA bug
        # control.inla = list(control.vb = list(emergency = 30)), # NOTE: unknown parameter
        control.predictor = list(compute = TRUE, link = linkCol),
        ...
      )
      private$inlaRes
    },
    plotData2 = function(
    ) {
      # self <- jsmMdl; private <- self$.__enclos_env__$private
      l <- self$getL()
      d <- self$getData()
      d$.variable <- d[[l$pars$lhsInlaVar]]
      d$.value <- d[[l$pars$outcomeIdInlaVar]] |>
        purrr::imap_dbl(\(v, idx) d$MX[idx, v])
    },
    #'
    #' Provides a visualization of the data used in the model,
    #' with options to specify the x-axis, color, and shape variables,
    #' as well as which blocks to include.
    #'
    plotData = function(
      xVar = NULL, colorVar = NULL, shapeVar = NULL,
      blockIds = NULL, addFamily = TRUE
    ) {
      l <- self$getL()
      blocks <- l$blocks
      if (!is.null(blockIds)) {
        blocks <- blocks[blockIds]
      }
      d <- purrr::imap(blocks, \(block, blockId) {
        fr <- block$parsedFormula$fr
        dd <- tidyr::pivot_longer(
          fr, cols = block$outVar, names_to = ".variable", values_to = ".value"
        )
        if (addFamily) {
          dd <- dd |>
            dplyr::mutate(
              .variable = paste0(.variable, " (", block$family, ")")
            )
        }
        dd
      }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(id = factor(id, levels = id[!duplicated(id)]))

      userVars <- colnames(d) |>
        purrr::keep(\(n) (n != ".value" && n != ".variable"))
      contVars <- userVars |>
        purrr::keep(\(n) is.numeric(d[[n]]))
      catVars <- userVars |>
        purrr::keep(\(n) !is.numeric(d[[n]]))

      if (is.null(xVar)) {
        stopifnot(length(contVars) > 0L)
        xVar <- contVars[[1L]]
      }
      if (is.null(colorVar)) {
        if (length(catVars) > 0L) {
          colorVar <- catVars[[1L]]
        }
      }
      if (is.null(shapeVar)) {
        if (length(catVars) > 1L) {
          shapeVar <- catVars[[2L]]
        }
      }

      p <- ggplot2::ggplot(d) +
        ggplot2::aes(x = !!rlang::sym(xVar), y = .value)
      if (!is.null(colorVar)) {
        p <- p + ggplot2::aes(color = !!rlang::sym(colorVar))
      }
      if (!is.null(shapeVar)) {
        p <- p + ggplot2::aes(shape = !!rlang::sym(shapeVar))
      }
      p <- p +
        ggplot2::geom_line(linewidth = 0.25, alpha = 0.25) +
        ggplot2::geom_point(size = 0.25) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~.variable, scales = "free_y") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab(NULL)
      p
    }
  )
)

# ------------------------------------------------------------------------
# jsm - top level function
# ------------------------------------------------------------------------

#'
#' Create a Joint Spatial Model (JSM) for INLA using the specified long and outcome formulas and data.
#'
#' @param formLong,formOut Formulas for the longitudinal and outcome model parts.
#'
#'  Named lists of formulas for the longitudinal/outcome parts of the model.
#'  If a name is missing, the corresponding formula outcome variable is used.
#'
#' @param dataLong,dataOut Data for each of the longitudinal and outcome formulas.
#'
#'  Named lists of data frames corresponding to each formula.
#'  If names are missing, matching to the formulas is done by order.
#'  When instead of a list, a single data frame is provided, it is used for all formulas.
#'
#' @param familyLong,familyOut Distribution families for each formula.
#'
#'  Named vector of families for each formula (`gaussian` used by default).
#'
#' @param fit Logical, whether to auto-fit the model (deprecated, use callInla instead).
#'
#' @param ... Additional arguments passed to JSMMdl$new().
#'
#' @param pars The parameters for the model construction. These parameters override
#' the default parameters described in \code{\link{defaultPars}}.
#'
#' @return An object of class \code{\link{JSMMdl}}.
#'
#' @examples
#' \dontrun{
#' d <- read.csv(...)
#' mdl <- jsm(
#'   fit = FALSE,
#'   # ----- longitudinal model -----
#'   formLong = list(x ~ 1 + t + (1 + t|id)),
#'   dataLong = d, familyLong = "gaussian",
#'   # ----- outcome model -----
#'   formOut = list(y ~ 1 + t + (1 + t|id)),
#'   dataOut = d, familyOut = "beta"
#' )
#' mdl$plotData()
#'
#' formulaText <- mdl$getFormulaText(addComments = TRUE, n = 4L)
#' cat(formulaText, "\n")
#'
#' mdl$callInla()
#' }
#'
#' @seealso \code{\link{defaultPars}}
#' @export
#'
jsm <- function(
  formLong, dataLong,
  formOut, dataOut = dataLong,
  familyLong = NULL, familyOut = NULL,
  fit = TRUE,
  ...,
  pars = NULL,
  userPars = list()
) {
  JSMMdl$new(
    outForms = formOut,
    longForms = formLong,
    outData = dataOut,
    longData = dataLong,
    outFamily = familyOut,
    longFamily = familyLong,
    ...,
    fit = fit,
    pars = pars,
    userPars = userPars
  )
}

# ------------------------------------------------------------------------
