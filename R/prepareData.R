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
#' @param dataOut,dataLong Either a single data frame used for all formulae of
#'   that part, or a list of data frames named by response (as in
#'   `collectFormulas()$out` / `$long`).
#' @param forms Result of [collectFormulas()].
#' @param na.action Strategy for rows with missing values in the variables of a
#'   formula: a function such as [stats::na.omit()] (default, drop the rows) or
#'   [stats::na.fail()] (stop with an error), or its name as a string.
#' @return A list with elements `out` and `long`, each a list of data frames
#'   named by response. Each data frame contains
#'   * `.key`: character row key, `"K"` followed by the row names of the user's
#'     data frame if it has any (they must be unique), otherwise by the row
#'     number. Rows coming from the same user data frame (or subsets of it) thus
#'     share keys across formulae,
#'   * the variables used by the formula (response, covariates, grouping
#'     variables), untransformed; outcome data also contain the `linkVariables`
#'     (variables of the `scaled()` targets, needed for the zero equations),
#'   * one factor column per grouping term, named as the term (e.g. `id`,
#'     `site:id`); a grouping variable is thereby replaced by its factor version.
#'
#'   Only rows kept by `na.action` are returned.
#' @export
prepareData <- function(dataOut, dataLong, forms, na.action = stats::na.omit) {
  na.action <- match.fun(na.action)
  list(
    out = assignData(dataOut, forms$out, "dataOut", na.action),
    long = assignData(dataLong, forms$long, "dataLong", na.action)
  )
}

assignData <- function(data, parsed, argName, na.action) {
  data <- matchData(data, parsed, argName)
  res <- lapply(names(parsed), function(r) {
    df <- addKey(data[[r]], argName)
    selectData(parsed[[r]], df, na.action)
  })
  stats::setNames(res, names(parsed))
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

  for (r in resp) checkVariables(data[[r]], parsed[[r]], argName)
  data
}

# Select the formula's variables and the rows kept by na.action, and add the
# grouping factors, using lme4::glFormula(); rows are identified by .key.
# lme4's sample-size checks are switched off: they target lme4 fits, not INLA.
selectData <- function(p, df, na.action) {
  formula <- p$lme4Formula
  rownames(df) <- df[[keyColumn]]
  if (!length(p$re)) {
    frame <- stats::model.frame(formula, df, na.action = na.action)
    groups <- list()
  } else {
    ctrl <- lme4::glmerControl(check.nobs.vs.nlev = "ignore",
                               check.nobs.vs.nRE = "ignore",
                               check.nlev.gtr.1 = "ignore")
    glf <- lme4::glFormula(formula, data = df, family = stats::gaussian,
                           na.action = na.action, control = ctrl)
    frame <- glf$fr
    groups <- glf$reTrms$flist
  }
  vars <- unique(unlist(p$variables, use.names = FALSE))
  res <- df[rownames(frame), c(keyColumn, vars, p$linkVariables), drop = FALSE]
  if (length(p$linkVariables)) {
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
