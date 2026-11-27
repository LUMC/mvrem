# mvrem 0.2.0

* Latent equations (E38): `latent(eta) ~ 1 + t + (1 | id) + scaled(x1) + scaled(x2)` is an
  unobserved process combining fixed and random effects with scaled copies of longitudinal
  responses; outcomes link to it with `scaled(eta)` (markers -> latent -> outcomes). Its options
  go in the family list, `family = list(eta = latentSpec(center, centerTargets, logPrec))`. The
  latent is evaluated at the values of its variables in the data of its outcomes, through one
  field shared by all of them. `fittedValues(m, "eta")` gives the latent process at these points;
  `associations()` gets a column `level`. Its scale and sign are identified only by the priors
  unless a coefficient is fixed: `scaled(x, coef = 1)` (new, for any `scaled()`). Checked against
  a hand-written INLA model of the same structure (identical results at the same
  hyperparameters). Not yet: `predict(method = "sample")` (refits instead), latents of latents.
* `jsm(zeroLogPrec = )` sets the fixed log precision of the zero-equation fields (default -6, as
  before; also `options(mvrem.zeroLogPrec = )`); it is kept by `update()` and `predict()`.
* `modelFit()` has a row `observations` (the criteria over the observation rows of all responses)
  before `overall` (INLA's totals, which also count the zero-equation rows: these grow with the
  number of links, so compare models with `observations`).

* Parallel computation (E37): `threads` of `jsm()` / `refit()` / `predict()` also takes `"auto"` (the
  cores `availableThreads()` finds: container CPU quota, SLURM / SGE variables, 2 under `R CMD check`,
  only the performance cores on Apple Silicon) and `"A:B"` (outer : inner INLA threads); `jsmApply()`
  runs independent fits in parallel R processes (forked on Linux and on macOS in a terminal or
  `Rscript`, separate sessions on Windows and in RStudio / Positron / R.app) and divides the core budget.
  The separate sessions attach mvrem and get the workspace variables that `FUN` uses, so `FUN`
  runs as in the caller's session. Tested on Linux and macOS (Apple Silicon); Windows is untested.
* Colour and shape (E30): `plot()`, `plotResiduals()`, `plot()` of predictions and
  `plotRandomEffects(type = "pairs")` take `colorBy` and `shapeBy` (a categorical variable of the
  data, including columns no formula uses); the widget offers all such variables in dropdowns,
  and clicking a level in a legend hides / shows that subset.
* The widget's search box lists the matching subjects in a dropdown and also finds them by a
  subject-level variable (e.g. a clinic).
* `survivalCurves()` arguments renamed before their first release: `colorBy` (was `by`), `nTimes`
  (was `n`), `seed = NULL` (was 1), as in `plot()` / `predict()`; a missing value of `colorBy` is the
  level `"(missing)"`.
* `plotAssociations()` takes `response`; `plot()` of predictions takes `response` and `ids`.
* Arguments are validated up front: `fit`, `verbose`, `criteria` must be TRUE / FALSE, `threads` a
  positive whole number, `"auto"` or `"A:B"`, `inlaOptions` a named list (also with `fit = FALSE`), `n` / `ncol` /
  `nsamples` / `nTimes` whole numbers.
* Bugs fixed: `offset()` in a formula was dropped silently (now an error: give the exposure on the
  left-hand side); a formula without fixed effects (`x ~ 0 + (1 | id)`) failed; a data frame without
  rows, or `na.action = na.omit` dropping every row of a response, failed with an internal error;
  a link variable of different type in the outcome and the longitudinal data failed with an
  internal error; a categorical covariate with one value gave R's contrasts error; a numeric grouping
  variable was offered as a time variable; `survivalCurves()` failed when `colorBy` had missing
  values; the model survival curve of a Weibull outcome with a fixed shape was all `NA`; a subject
  with outcome data but no marker data got marker random-effect labels outside the model's levels
  when the outcome (or latent equation) had no random effect on that group.
* `varianceComponents()$table` has a column `blockId` (the position of the block in the
  formula) next to the display label `block`.

* `plot()` draws a survival outcome as its Kaplan-Meier curve (censored times as open circles,
  linked to the subjects in the widget) with the model's survival curve averaged over the subjects;
  a 95% band when the fit keeps INLA's configurations (`config = TRUE`); one curve per level with
  `colorBy`. The curves are also available as a data frame from the new `survivalCurves()`.
* `plot()` of predictions with `colorBy`: the observations behind the predictions are faint, the
  predicted means and bounds opaque.

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
