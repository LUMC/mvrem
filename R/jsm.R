#' Fit a multivariate joint model
#'
#' Wraps the whole pipeline: [collectFormulas()] (formulae and families),
#' [prepareData()] (data), [summarizeData()], [buildInlaFormula()] (INLA
#' formula and stacked data) and, unless `fit = FALSE`, the call of INLA. The
#' user never has to deal with INLA's interface: priors, families and
#' associations are given with the functions of this package and the results are
#' read with [fixedEffects()], [hyperparameters()], [associations()] or
#' `summary()`.
#'
#' See [collectFormulas()] for the formula syntax (`re()`, `scaled()`), the
#' families and priors, and [prepareData()] for the data arguments. Family
#' inputs follow the response on the left-hand side, as in brms:
#' `k | trials(n) ~ ...`, `cnt | exposure(pt) ~ ...`, `time | event(status) ~ ...`
#' and observation weights `x | weights(w) ~ ...` (gaussian, Student-t, gamma,
#' beta: the precision of row i is multiplied by `w[i]`; see [families]).
#'
#' The formulae are given as one list (`form`; outcomes are the formulae with
#' `scaled()` terms, see [splitFormulas()]) or, as before, as the pair
#' `formOut` / `formLong` (positionally `jsm(outcomes, longitudinal, data)`),
#' which is needed only for an outcome without a `scaled()` term. A pair-style
#' model stores its call with `formOut` / `formLong` whichever way it was
#' written.
#'
#' @inheritParams collectFormulas
#' @inheritParams prepareData
#' @param time Name of the time variable (optional): a numeric covariate of at
#'   least one formula. It is the horizontal axis of the plots and is reported by
#'   [timeVariable()]; the model's equations do not depend on it. Default: the
#'   numeric covariate used by most formulae (a message names the other
#'   candidates when there are several).
#' @param naming Naming scheme of the INLA objects, see [namingScheme()].
#' @param fit If `FALSE`, the model is only constructed (formulae, data, INLA
#'   formula and stacked data, LaTeX equations are all available); INLA is not
#'   needed. Fit it later with [refit()].
#' @param verbose Show INLA's progress output.
#' @param threads Number of threads for INLA (default: INLA's choice).
#' @param inlaOptions Advanced: named list of further arguments of
#'   `INLA::inla()`, e.g. `list(control.compute = list(dic = TRUE, waic = TRUE))`.
#'   The formula, families and data are built by this package and cannot be set.
#' @param criteria If `TRUE`, INLA also computes DIC, WAIC and CPO, so that
#'   `modelFit()` works (slower). Merged under `inlaOptions$control.compute`.
#' @return An object of class `c("jsm", "mvremModel")`. Use the accessors
#'   ([modelFormulas()], [modelData()], [dataSummary()], [inlaFormula()],
#'   [inlaData()], [inlaFit()], [fixedEffects()], [hyperparameters()],
#'   [associations()]), [latexEquation()], [formatInlaFormula()], `summary()`,
#'   `coef()`, [refit()] and `update()`.
#'
#'   If INLA fails, the model is returned unfitted with a warning (and
#'   [refit()] can be tried again).
#' @examples
#' \dontrun{
#' m <- jsm(
#'   list(y ~ 1 + t + re(1 | id) + scaled(x1),
#'        x1 ~ 1 + t + re(1 + t | id, prior = wishartPrec())),
#'   data   = d,
#'   family = list(y = betaFamily())
#' )
#' # the same model in the pair style
#' m2 <- jsm(
#'   formOut  = list(y ~ 1 + t + re(1 | id) + scaled(x1)),
#'   formLong = list(x1 ~ 1 + t + re(1 + t | id, prior = wishartPrec())),
#'   data     = d,
#'   family   = list(y = betaFamily())
#' )
#' summary(m)
#' associations(m)
#' latexEquation(m)
#' }
#' @export
jsm <- function(form = NULL, formLong = NULL, data = NULL, dataOut = NULL, dataLong = NULL,
                family = NULL, familyOut = NULL, familyLong = NULL, time = NULL,
                na.action = stats::na.pass,
                naming = getOption("mvrem.naming", namingScheme()),
                fit = TRUE, verbose = FALSE, threads = NULL, inlaOptions = list(),
                criteria = FALSE, formOut = NULL) {
  call <- match.call()
  # the pair style is stored as formOut / formLong whatever the argument name used, so that
  # update(m, formOut = ...) works for every pair model (as before E18)
  if (!is.null(call$form) && !is.null(call$formLong) && is.null(call$formOut))
    names(call)[names(call) == "form"] <- "formOut"
  if (isTRUE(fit)) requireInla()
  if (isTRUE(criteria))
    inlaOptions <- utils::modifyList(
      list(control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)), inlaOptions)
  forms <- collectFormulas(form, formLong, family, familyOut, familyLong, formOut)
  dat <- prepareData(forms, data, dataOut, dataLong, na.action)
  summ <- summarizeData(forms, data, dataOut, dataLong)
  inl <- buildInlaFormula(forms, dat, naming)
  tv <- resolveTimeVariable(forms, dat, time)

  x <- structure(
    list(call = call, forms = forms, data = dat, dataSummary = summ, inla = inl, time = tv,
         options = list(verbose = verbose, threads = threads, inlaOptions = inlaOptions),
         fit = NULL, inlaError = NULL),
    class = c("jsm", "mvremModel"))
  if (isTRUE(fit)) x <- refit(x)
  x
}
