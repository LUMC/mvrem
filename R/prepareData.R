# lme4 control that switches off lme4's sanity checks of the design (observations vs
# levels / random effects, groups with one level): mvrem uses lme4 only to build grouping
# factors (prepareData()) or a rough covariance for a data-scaled prior (lme4BlockFit()),
# where such designs are legitimate (e.g. a fixture with one row per subject). `singular`:
# also silence the singular-fit check (a singular fit is handled by the caller).
lme4Control <- function(fun = c("glmer", "lmer"), singular = TRUE) {
  fun <- match.arg(fun)
  args <- list(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", check.nlev.gtr.1 = "ignore")
  if (singular) args$check.conv.singular <- "ignore"
  do.call(if (fun == "lmer") lme4::lmerControl else lme4::glmerControl, args)
}

# Name of the internal row-key column added to every data frame.
keyColumn <- ".key"

#' Prepare the data of a joint model
#'
#' Matches the outcome and longitudinal data to the parsed formulae (see
#' [collectFormulas()]), checks that all variables used by a formula are
#' present and adds a row-key column `.key`. [lme4::glFormula()] is then called
#' on the simplified `lme4Formula` to select the variables, build the grouping
#' factors (e.g. `site:id`) and handle rows with missing values.
#'
#' @param forms Result of [collectFormulas()].
#' @param data Data for all formulae (outcome and longitudinal): a single data
#'   frame, or a list of data frames named by response (any response of the
#'   model). Used for each part for which `dataOut` / `dataLong` is not given.
#' @param dataOut,dataLong Data for the outcome / longitudinal formulae only,
#'   overriding `data` for that part: a single data frame used for all
#'   formulae of the part, or a list of data frames named by the responses of
#'   the part (as in `collectFormulas()$out` / `$long`). At least one of `data`
#'   and the part-specific argument is required for each part.
#' @param na.action Strategy for rows with missing values in the variables of a
#'   formula: [stats::na.pass()] (default: keep all rows), [stats::na.omit()]
#'   (drop the rows) or [stats::na.fail()] (stop with an error), another
#'   function, or its name as a string. Kept rows may have a missing response
#'   (INLA treats it as a value to predict); a missing covariate, grouping
#'   variable, `trials()` / `exposure()` / `weights()` value or link variable in a kept row is
#'   an error (use `na.omit` to drop such rows).
#' @return A list with elements `out` and `long`, each a list of data frames
#'   named by response. Each data frame contains
#'   * `.key`: character row key, `"K"` followed by the row names of the user's
#'     data frame if it has any (they must be unique), otherwise by the row
#'     number. Rows coming from the same user data frame (or subsets of it) thus
#'     share keys across formulae,
#'   * the variables used by the formula (response, `trials()` / `exposure()` / `event()` / `weights()`
#'     variables, covariates, grouping variables), untransformed; outcome data also contain the `linkVariables`
#'     (variables of all longitudinal formulae, needed for the zero equations),
#'   * one factor column per grouping term, named as the term (e.g. `id`,
#'     `site:id`); a grouping variable is thereby replaced by its factor version.
#'
#'   Only rows kept by `na.action` are returned (all rows by default). The responses (and
#'   `trials()` / `exposure()` values, `weights()` positive and finite) are checked against the support of the
#'   formula's family (see [families]), e.g. a beta response must be in (0, 1).
#' @export
prepareData <- function(forms, data = NULL, dataOut = NULL, dataLong = NULL,
                        na.action = stats::na.pass) {
  na.action <- match.fun(na.action)
  list(
    out = assignData(resolveData(forms, "out", data, dataOut), forms$out, na.action),
    long = assignData(resolveData(forms, "long", data, dataLong), forms$long, na.action)
  )
}

assignData <- function(data, parsed, na.action) {
  argName <- attr(data, "argName")
  res <- lapply(names(parsed), function(r) {
    df <- addKey(data[[r]], argName)
    res <- selectData(parsed[[r]], df, na.action)
    checkResponse(parsed[[r]], res)
    res
  })
  stats::setNames(res, names(parsed))
}

# Data of one part: the part-specific argument if given, else `data`, matched
# to the part's formulae (see matchData()). A list `data` may be named by the
# responses of both parts; the part's elements are used. The result has the
# name of the argument used as attribute `argName` (for messages).
resolveData <- function(forms, part, data, dataPart) {
  argName <- if (part == "out") "dataOut" else "dataLong"
  parsed <- forms[[part]]
  if (is.null(dataPart)) {
    if (is.null(data))
      stop("No data for the ", if (part == "out") "outcome" else "longitudinal",
           " formulae: give `data` or `", argName, "`.", call. = FALSE)
    if (is.data.frame(data) || !is.list(data)) {
      argName <- "data"
    } else {
      nms <- names(data)
      if (is.null(nms) || any(nms == "") || anyDuplicated(nms))
        stop("`data` must be a data frame or a list of data frames with ",
             "unique names matching the responses.", call. = FALSE)
      allResp <- c(names(forms$out), names(forms$long))
      extra <- setdiff(nms, allResp)
      if (length(extra))
        stop("Data frame(s) in `data` do not correspond to any formula: ",
             paste(extra, collapse = ", "), call. = FALSE)
      argName <- "data"
      data <- data[nms %in% names(parsed)]
    }
    dataPart <- data
  }
  res <- matchData(dataPart, parsed, argName)
  attr(res, "argName") <- argName
  res
}

# Data frame per response (a single data frame is shared), with the
# formula's variables checked to be present.
matchData <- function(data, parsed, argName) {
  resp <- names(parsed)

  if (is.data.frame(data)) {
    data <- stats::setNames(rep(list(data), length(resp)), resp)
  } else if (is.list(data)) {
    nms <- names(data)
    if (is.null(nms) || any(nms == "") || anyDuplicated(nms))
      stop("`", argName, "` must be a data frame or a list of data frames with ",
           "unique names matching the responses.", call. = FALSE)
    missingResp <- setdiff(resp, nms)
    if (length(missingResp))
      stop("No data frame in `", argName, "` for response(s): ",
           paste(missingResp, collapse = ", "), call. = FALSE)
    extra <- setdiff(nms, resp)
    if (length(extra))
      stop("Data frame(s) in `", argName, "` do not correspond to any formula: ",
           paste(extra, collapse = ", "), call. = FALSE)
    notDf <- nms[!vapply(data, is.data.frame, TRUE)]
    if (length(notDf))
      stop("Element(s) of `", argName, "` are not data frames: ",
           paste(notDf, collapse = ", "), call. = FALSE)
    data <- data[resp]
  } else {
    stop("`", argName, "` must be a data frame or a named list of data frames.",
         call. = FALSE)
  }

  data <- lapply(data, plainDataFrame)
  for (r in resp) checkVariables(data[[r]], parsed[[r]], argName)
  data
}

# A plain data frame: tibbles and other subclasses (no row names, list-like `[`)
# break the row-name based matching of prepareData().
plainDataFrame <- function(d) if (identical(class(d), "data.frame")) d else as.data.frame(d)

# Select the formula's variables and the rows kept by na.action, and add the
# grouping factors, using lme4::glFormula(); rows are identified by .key.
# lme4's sample-size checks are switched off: they target lme4 fits, not INLA.
selectData <- function(p, df, na.action) {
  formula <- p$lme4Formula
  rownames(df) <- df[[keyColumn]]
  # Only the response may stay missing (INLA treats a missing response as a
  # prediction); covariates, grouping variables, trials / exposure and link
  # variables must be complete in the rows that are kept.
  pred <- unique(c(p$variables$covariate, p$variables$group, p$variables$aux, p$linkVariables))
  if (length(pred)) {
    kept <- na.action(df[pred])
    bad <- pred[vapply(pred, function(v) anyNA(kept[[v]]), NA)]
    if (length(bad))
      stop("Missing values in ", paste0("`", bad, "`", collapse = ", "), " (formula of '", p$response,
           "'): only the response may be missing with this `na.action`. Use ",
           "`na.action = na.omit` to drop such rows.", call. = FALSE)
  }
  if (!length(p$re)) {
    frame <- stats::model.frame(formula, df, na.action = na.action)
    groups <- list()
  } else {
    ctrl <- lme4Control("glmer", singular = FALSE)
    glf <- lme4::glFormula(formula, data = df, family = stats::gaussian,
                           na.action = na.action, control = ctrl)
    frame <- glf$fr
    groups <- glf$reTrms$flist
  }
  vars <- unique(unlist(p$variables, use.names = FALSE))
  res <- df[rownames(frame), unique(c(keyColumn, vars, p$linkVariables)), drop = FALSE]
  # the model frame above covers only the variables of the lme4 formula; trials / exposure /
  # weights (aux) and the link variables are not in it, so na.action is applied once more
  # to drop (na.omit) or reject (na.fail) rows missing one of them
  if (length(p$variables$aux) || length(p$linkVariables)) {
    res <- na.action(res)
    attr(res, "na.action") <- NULL
  }
  for (g in names(groups)) res[[g]] <- groups[[g]][rownames(frame) %in% rownames(res)]
  rownames(res) <- NULL
  res
}

checkVariables <- function(df, p, argName) {
  vars <- unique(c(unlist(p$variables, use.names = FALSE), p$linkVariables))
  miss <- setdiff(vars, names(df))
  if (length(miss))
    stop("Variable(s) used in formula for '", p$response, "' not found in `",
         argName, "`: ", paste(miss, collapse = ", "), call. = FALSE)
  df
}

# Key = "K" + user-provided row names if present (must be unique), else "K" + row number.
addKey <- function(df, argName) {
  if (keyColumn %in% names(df))
    stop("Column name '", keyColumn, "' is reserved for internal use.", call. = FALSE)
  hasRowNames <- .row_names_info(df) > 0L
  ids <- if (hasRowNames) rownames(df) else seq_len(nrow(df))
  if (anyNA(ids) || anyDuplicated(ids))
    stop("Row names of a data frame in `", argName, "` must be unique and not missing.",
         call. = FALSE)
  df[[keyColumn]] <- paste0("K", ids)
  df
}
