# Helpers of the build tests (test-build.R): small models on fxSmallData()
# built without INLA.
#
#   smallData()                fxSmallData() with the factor `g` not collinear with `t`
#                              (lme4 reports a rank-deficient design for `t + g` otherwise)
#   smallForms(re, out, fixed) collectFormulas() of `out` (default y ~ t + scaled(x)) and
#                              x ~ <fixed> + <re> (strings, default "t" and "(1 | id)")
#   buildSmall(forms, data, ..., naming)
#                              buildInlaFormula(forms, prepareData(forms, data, ...), naming);
#                              data = fxSmallData() unless given
#   termOf(inl, index)         the f() call of the field `index` in the INLA formula
#   termHyper(inl, index)      its evaluated `hyper` argument (NULL if none)

smallData <- function() fxSmallData()   # kept as the name used throughout this file (helpers are sourced alphabetically: lazy)

smallForms <- function(re = "(1 | id)", out = y ~ t + scaled(x), fixed = "t") {
  collectFormulas(list(out), list(stats::as.formula(paste("x ~", fixed, "+", re))))
}

buildSmall <- function(forms, data = fxSmallData(), ...,
                       naming = getOption("mvrem.naming", namingScheme())) {
  buildInlaFormula(forms, prepareData(forms, data, ...), naming = naming)
}

termOf <- function(inl, index) {
  calls <- unlist(inl$terms, recursive = FALSE)
  for (cl in calls)
    if (is.call(cl) && identical(deparse1(cl[[2L]]), index)) return(cl)
  stop("no f() term for index ", index)
}

termHyper <- function(inl, index) eval(as.list(termOf(inl, index))$hyper)
