# mvrem 0.1.0

First public release. The interface may still change.

* `jsm()` builds and fits multivariate joint models with INLA from lme4-style formulae: several
  longitudinal markers and outcomes, each with its own family, fixed and random effects; an outcome
  is linked to a marker by `scaled(x)` (the marker's linear predictor times an estimated
  coefficient), also when the two are measured at different times.
* Families: Gaussian, Student t, lognormal, gamma, beta, Poisson, negative binomial, binomial,
  beta-binomial, zero-inflated Poisson and negative binomial, Weibull and exponential survival
  (marker evaluated at one time per subject), and any other INLA likelihood via `inlaFamily()`;
  inputs `trials()`, `exposure()`, `event()` and observation `weights()`.
* Priors for random effects (`gammaPrec()`, `pcSd()`, `halfNormalSd()`, `halfCauchySd()`,
  `wishartPrec()`), optionally scaled from an lme4 fit (`"data"`), and `normalCoef()` for
  associations.
* Results with readable names: `summary()`, `estimatesTable()`, `fixedEffects()`, `associations()`,
  `hyperparameters()`, `randomEffects()`, `varianceComponents()`, `fittedValues()`, `residuals()`,
  `modelFit()`; predictions for new subjects and times (`predict()`, by refitting or by sampling the
  posterior).
* Plots, interactive when htmlwidgets is installed: data and fitted values, residuals, random
  effects, variance components, fixed effects, associations, predictions; a model diagram
  (`plotModel()`) and the model equations in LaTeX (`latexEquation()`).
* Not yet: a time-varying ("current value") survival association. The default priors are not final.
