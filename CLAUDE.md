# mvrem

R package under development: tools for constructing and fitting Multivariate Joint Models
with the Integrated Nested Laplace Approximation (INLA). Models have several longitudinal and
outcome processes with fixed and random effects; every outcome is associated with every
longitudinal process through a scaled copy of its linear predictor.

## User-facing API (target)

```r
jsmMdl <- jsm(
  formOut = list(
    z ~ 1 + t + (1 + t | id),            # plain lme4 notation -> default priors
    y ~ 1 + t + re(1 + t | id, ...)      # correlated intercept/slope; ... = prior / f() args
              + scaled(x1, ...)          # parameters of the association with x1 (optional)
  ),                                     # z and y are linked to x1, x2, x3 in any case
  formLong = list(
    x1 ~ 1 + t + re(0 + t | id, ...) + re(1 | id, ...),  # two independent RE blocks
    x2 ~ 1 + t + re(1 + t | id, ...),
    x3 ~ 1 + t + re(1 + t || id, ...)                    # `||`: uncorrelated (lme4 semantics)
  ),
  dataLong = d,
  dataOut = d,
  familyOut = list(y = betaFamily(phi = gammaPrec(1, 0.1))),  # unlisted -> gaussianFamily()
  familyLong = list(x3 = poissonFamily()),                    # INLA names ("poisson") also ok
  fit = FALSE                            # FALSE: build the model object only, don't call INLA
)
```

### Syntax rules

- **Left-hand side**: an expression of exactly one variable (`log(y)` ok, `cbind(a, b)` not);
  that variable is the response name used everywhere (list names, `scaled()` targets, data
  list names). Family inputs follow brms-style after `|`: `k | trials(n) ~ ...` (binomial-type,
  INLA `Ntrials`), `cnt | exposure(pt) ~ ...` (count families, INLA `E`); any expression inside,
  one unnamed argument each, only these two.
- **Random effects**: `(... | g)` (lme4 notation) = `re(... | g)` with default priors; `|`
  correlated, `||` uncorrelated. `re(..., prior = <constructor>)` sets the block's prior (see
  Priors); it is not an `f()` argument.
- **Associations**: every outcome is ALWAYS linked to EVERY longitudinal response (zero
  equation + scaled copy with an estimated scaling coefficient). `scaled(x, ...)` (bare name,
  no `~`) only sets non-default parameters of the link to `x`. Only in `formOut`; `x` must be a
  response in `formLong`; at most one `scaled()` per target in a formula.
- **Extra arguments** of `re()` / `scaled()` must be named; evaluated in the formula's
  environment; merged into the default `f()` arguments (`modifyList`, nested lists such as
  `hyper` merged). User parameters are given ONLY this way, no `userPars` argument keyed by
  INLA names.
- `re()` / `scaled()` must be top-level additive terms (not nested in other calls).
- **Families**: `familyOut` / `familyLong` are named lists (names = responses of that part) of
  family objects or INLA likelihood names (see Families).

## Status

Pipeline: `collectFormulas()` (step 1) -> `prepareData()` (step 2) -> `buildInlaFormula()`
(step 3); around it `summarizeData()`, `formatInlaFormula()`, `latexEquation()`.

Not done yet: the `jsm()` wrapper; the INLA call (stacks combined into one data list,
`family = inl$families`, `control.family = inl$controlFamily`, `Ntrials` / `E` from the
`inl$inputs` columns with `NA` -> 1); fitting; survival responses (`inla.surv`). A
hand-assembled `inla()` call (beta, probit binomial with trials, Poisson with exposure,
Gaussian long) ran and recovered simulated values with the zero-equation field's fixed
diffuse precision (see Open issues).

## Step 1: formula collection — `R/parseFormula.R`

`parseJsmFormula(formula, part = "out" / "long")` returns:
- `response` (the LHS variable), `lhs` (response expression as written, without `| ...`),
  `responseExpr` (its language object), `aux` (named list of `trials` / `exposure`
  expressions), `part`, `formula` (as given)
- `fixedFormula`: fixed part only (no REs, no `scaled()`); `lme4Formula`: `re(bar, ...)` ->
  `(bar)`, `scaled()` removed, plain bar terms kept
- `variables`: `response` / `covariate` / `group` / `aux` via `all.vars()` (`log(y)` -> `y`,
  `I(t^2)` -> `t`). Covariates: fixed terms + left of `|`; groups: right of `|` (`site:id` ->
  `site`, `id`); aux: variables of `trials()` / `exposure()`. `scaled()` targets and
  parameters are not variables.
- `re`: one entry per RE block (`re()` and plain `( | )` alike): `term`, `group`,
  `correlated` (`FALSE` for `||`), `params`, `prior` (`mvremPrior` or `NULL`; not in `params`)
- `scaled` (outcome only): the formula's `scaled()` terms as written (`target`, `params`)

`collectFormulas(formOut, formLong, familyOut = NULL, familyLong = NULL)` returns
`list(out, long, variables)`; `out` / `long` named by response, `variables` = data frame
(`variable`, `role` (`response` / `covariate` / `group` / `aux`), `part`, `formula`). It adds:
- `scaled` (outcome): replaced by one entry per long response (`formLong` order), `params`
  from the matching `scaled()` term or `list()`
- `linkVariables` (outcome): groups + covariates of ALL long formulae not already in the
  outcome formula; the outcome data must contain them (zero equations live at outcome rows)
- `family`: `mvremFamily` (see Families). Elements of `familyOut` / `familyLong` may be family
  objects or INLA names (named character vectors too): a name -> its constructor with
  defaults (`familyFromName()`), else `inlaFamily(name)`; unlisted -> `gaussianFamily()`.
  Names validated only if INLA is installed (else a warning).

Errors: one-sided formulae; LHS without exactly one variable; unsupported, repeated or
named-argument terms after `|` on the LHS; `trials()` / `exposure()` not used by the family
(`checkAux()`); duplicate responses (also `y` and `log(y)`); non-bar first argument of `re()`;
non-name first argument of `scaled()`; unnamed extra arguments; nested `re()` / `scaled()`;
`scaled()` in a long formula; unknown or repeated `scaled()` target; `prior` in `re()` not an
`mvremPrior`; `prior` in `scaled()` (not supported yet); family lists not a named list of
families / names, unnamed or duplicated names, names not responses of that part (message says
if it belongs to the other part), unknown INLA name, survival likelihood.

## Step 2: data preparation — `R/prepareData.R`

`prepareData(dataOut, dataLong, forms, na.action = na.omit)` returns `list(out, long)` with one
data frame per response. `dataOut` / `dataLong`: one data frame shared by that part, or a list
named by response variable (`y` for `log(y) ~ ...`) matching the responses exactly.

Each data frame contains:
- `.key`: `"K"` + the user's explicit row names (must be unique) or `"K"` + row number; the
  same row gets the same key across formulae (also for subsets of one data frame)
- the formula's variables, untransformed (incl. aux), + `linkVariables` for outcomes
- one factor column per grouping term, named as the term (`id`, `site:id`; replaces a
  same-named variable)

`lme4::glFormula()` on `lme4Formula` only selects rows and builds the grouping factors (no
design matrices kept; lme4's nobs / nlev checks disabled); formulae without REs use
`model.frame()`. `na.action` (function or name, e.g. `na.omit` / `na.fail`) applies per formula
to its variables, aux and link variables; other columns are ignored. `checkResponse()` then
checks the (transformed) response against the family's `support`: `real`, `positive`, `unit`
(0, 1), `count`, `binomial` (0..trials; 0/1 without `trials()`), `any` (`inlaFamily()`);
trials must be positive integers, exposure positive.

Errors: missing / extra / unnamed list elements, non-data-frame elements, formula or link
variables missing in the data, user data containing `.key`, non-unique row names, responses /
trials / exposure outside the family's support.

`summarizeData(forms, dataOut, dataLong)` (`R/summary.R`; raw data as for `prepareData()`,
shares `matchData()`) -> class `mvremDataSummary`: `equations` (part, response, INLA family,
rows, complete, dropped, groups = levels per grouping term among complete rows; long first),
`missing` (per equation and variable, incl. aux and link variables), `zero` (per (long, out)
pair: matching variables, number of unique combinations in complete outcome rows).
`complete` = rows kept by `prepareData(na.action = na.omit)`.

## Step 3: INLA formula — `R/inlaFormula.R`

`buildInlaFormula(forms, data, naming = getOption("mvrem.naming", namingScheme()))`.

### Naming

All names come from `namingScheme(...)` templates (placeholders `{resp}`, `{term}`, `{terms}`,
`{group}`, `{level}`, `{long}`, `{out}`, `{vars}`, `{values}`; names made syntactic; unknown
templates and duplicated resulting names -> error). Default scheme (user-provided reference,
must be kept):
- response matrix `MX`; columns `y_y`, `x1_x1`, zero equations `x1_x1_0eq`
- family inputs `<resp>_Ntrials` (`trials`), `<resp>_E` (`exposure`)
- fixed `<resp>_<term>`, intercept `1` (`y_1`, `y_t`); RE slopes use this column as weight
- RE index `<resp>_<term>_<group>` (`y_1_id`, `y_t_id`), `values = LEVELS$<resp>_<terms>_<group>`
  (`y_1_t_id`), labels `"y|1|S001"`, `n = length(LEVELS$...)`
- zero-equation field `<long>_TO_<out>_<vars>` (`x1_TO_y_id_t`), weight `<long>_WEIGHT`,
  labels `"x1|TO|y|S001|0.5"`; `{vars}` = groups + covariates of the long formula
- scaled copy `<out>_FROM_<long>_<vars>` (`y_FROM_x1_id_t`)

### Terms

`f()` arguments are embedded as the calls that build them (`fTerm()`: `str2lang(deparse1())`
for non-scalar values): `terms()` deparses embedded list objects WITHOUT names
(`hyper = list(list("pc.prec", ...))`), and INLA then silently ignored every `hyper` (RE
priors, scaled `beta`, zero field). Tested via `INLA:::inla.interpret.formula()`; tests read
`hyper` with `eval()`.

Grouped per response, outcomes first, then longitudinal:
- fixed: one column per design-matrix column (factor levels expand)
- random: `|` with k > 1 -> `iid<k>d` (`iidkd` + `order` for k > 5) on the first effect, `copy`
  for the others; a single effect or `||` -> one `iid` per effect. Without `prior` no `hyper`
  (INLA default); with `prior`, `hyper = priorHyper(prior, k)` after `n`; then `params` merged
  (warning if a user `hyper` sets `prior` / `param` of it)
- zero equation, per (long, out) pair (all pairs): `f(x1_TO_y_id_t, x1_WEIGHT, model = "iid",
  values = LEVELS$x1_TO_y_id_t, hyper = list(prec = list(initial = -6, fixed = TRUE)))`
  (`zeroFieldHyper`: fixed diffuse precision of `u`), i.e. `0 = eta_x1 - u`
- scaled copy, per pair: `f(y_FROM_x1_id_t, copy = "x1_TO_y_id_t", hyper = list(beta =
  list(fixed = FALSE, prior = "normal", param = c(0, 0.01))))`, then the pair's `scaled()`
  `params` merged

### Result

- `responses`, `zeroResponses` (column names)
- `families`: INLA likelihood per response column (responses, then zero responses =
  `"gaussian"`), in the order of the INLA likelihoods
- `controlFamily`: INLA `control.family`, same names / order as `families`:
  `controlFamily(family)` per response; zero equations `zeroControlFamily` (placeholder)
- `inputs`: per response the family-input columns, `c(Ntrials = "y_Ntrials")` /
  `c(E = "z_E")` (empty if none)
- `fixed`, `random` (blocks with `group`, `lhs`, `effects`, `index`, `weight`, `levels`,
  `terms`), `zero`, `scaled`, `levels` (= `LEVELS`), `naming`, `terms` (per response),
  `formulas` (the user's), `formula` (its environment holds `LEVELS`)
- `stacks`: per response the stacked data with the formula's column names (input columns
  right after the response). `main` = observation rows; long responses also have one
  zero-equation part per outcome: response `NA`, `_0eq` = 0, WEIGHT = -1, inputs `NA`,
  covariates and RE labels from the outcome rows, factor levels from the long data. One
  zero-equation row per UNIQUE combination of `{vars}` values in the outcome data
  (= `LEVELS$<zeroIndex>`); outcome rows are matched by these values (not `.key`), so
  combinations without a long observation (missing values, subjects absent from the long
  data) still get a zero equation (tested with artificially missing covariates).

## Families — `R/families.R`

Class `mvremFamily`: `family` (INLA name), `kind`, `link`, `linkDefault`, `hyper` (INLA hyper
short name -> prior), `specs` (LaTeX symbol, accepted prior constructors), `control`, `aux`
(allowed LHS inputs), `support`; `format` / `print` methods.

| Constructor | INLA likelihood | Hyper | LHS input | Support |
|---|---|---|---|---|
| `gaussianFamily(prec, link)` | `gaussian` | `prec` | | real |
| `studentTFamily(prec, dof)` | `t` | `prec`, `dof` | | real |
| `lognormalFamily(prec)` | `lognormal` | `prec` | | positive |
| `gammaFamily(prec, link)` | `gamma` | `prec` | | positive |
| `betaFamily(phi, link)` | `beta` | `phi` | | unit |
| `poissonFamily(link)` | `poisson` | | exposure | count |
| `negBinomialFamily(size, link)` | `nbinomial` | `size` | exposure | count |
| `binomialFamily(link)` | `binomial` | | trials | binomial |
| `betaBinomialFamily(rho, link)` | `betabinomial` | `rho` | trials | binomial |
| `zeroInflatedPoissonFamily(type, prob)` | `zeroinflatedpoisson0/1` | `prob` | exposure | count |
| `zeroInflatedNegBinomialFamily(type, size, prob)` | `zeroinflatednbinomial0/1` | `size`, `prob` | exposure | count |
| `inlaFamily(family, link, hyper)` | any non-survival | any | both | any |

- All take `control` (other `control.family` entries, merged last; warning if its `hyper`
  overrides a prior argument).
- Hyper arguments: `NULL` (INLA default); an `mvremPrior` (`prec` / `phi`: `gammaPrec()`,
  `pcSd()`, `halfNormalSd()`, `halfCauchySd()`; `size`: `gammaPrec()` only; converted with
  `priorHyper(prior, 1)`); or a raw INLA hyper list (`prior`, `param`, `initial`, `fixed`,
  `to.theta`, `from.theta`).
- Links validated against INLA's list (if installed). `controlFamily(fam)` (exported) ->
  `control.link = list(model = link)` only if the link differs from INLA's default, `hyper`,
  then `control`.
- Survival likelihoods (INLA `survival = TRUE`) are rejected: deferred (need `inla.surv`, one
  row per subject, association at event times).

## Priors — `R/priors.R`

Constructors (class `mvremPrior`: `dist`, `target` = `"prec"` / `"sd"` / `"precMatrix"`,
`args`; `format` / `print`), arguments validated on construction:
- `gammaPrec(shape = 1, rate = 5e-5)` -> `prec = list(prior = "loggamma", param = c(shape, rate))`
- `pcSd(u, alpha = 0.01)` -> `"pc.prec"`, `c(u, alpha)` (sigma ~ Exp(-log(alpha) / u))
- `halfNormalSd(scale)` -> `"logtnormal"`, `c(0, 1 / scale^2)` (INLA: sd ~ N(mean, 1/prec),
  truncated)
- `halfCauchySd(scale)` -> `prior = "expression: ..."` on theta = log tau (Jacobian
  log(sigma / 2); verified against logtnormal written the same way; INLA marks `expression:` as
  to be removed)
- `wishartPrec(df = NULL, scale = 1)`: Q ~ W_k(df, V), E[Q] = df V; `df = NULL` -> k + 2, scalar
  scale -> s I -> `theta1 = list(prior = "wishart<k>d" / "wishartkd" (k > 5), param = c(df,
  diag(R), R[lower.tri(R)]))`, R = V^-1 (order R12, R13, ..., R23, ... per INLA docs for
  k <= 5; iidkd documented only for k = 3, assumed the same; INLA accepts the k = 6 vector)

`checkBlockPrior()` (in `buildInlaFormula()`, and `latexEquation()` with data, as k depends on
factor expansion): a `|` block with k > 1 needs `wishartPrec()` (scale k x k, df > k + 1 as in
the INLA docs); single / `||` blocks need a scalar prior (the same for each effect).

## Display — `R/formatFormula.R`

`formatInlaFormula(inl, nrows = 4, digits = 3)` -> class `mvremFormulaText` (printed with
`cat`): one term per line; each response group preceded by `#` comments with the user
formula, `family = ..., control.family = ...` (long: also a `zero equations:` line) and the
first `nrows` rows of each stack part (`nrows = 0` hides data; numbers with `digits`
significant digits). Parses back to the same formula.

## LaTeX — `R/latex.R`

`latexEquation(x, data = NULL)`: `x` = a parsed formula or the `collectFormulas()` result;
`data` = the matching step-2 data frame / `prepareData()` result (expands factors into
`\mathbb{I}[g_i = v]`; without data a factor is one term). Returns class `mvremLatex`, the
body of an `aligned` environment; `format(x, delim = "$$")` wraps it. In the Rmd, rendered
math uses `results = "asis"` + `cat(format(..., delim = "$$"))`.

- Names: observation index `i`; `t` -> `t_i`, `x1` -> `x_{1,i}`, other -> `\mathit{name}_i`.
- Effects: `\beta^{(r)}_k` (k = 0 intercept); `b^{(r)}_{k,g[i]}` shares k with the fixed effect
  of the same column, else the next free k.
- Families (`familyLatex()`): link of the mean (`\operatorname{logit}`, `\log`, `\Phi^{-1}`) +
  a distribution line, e.g. `y_i \sim \mathrm{Beta}(\mu_i \phi, (1 - \mu_i)\phi)`,
  `\mathrm{Binomial}(n_i, p_i)` / `\mathrm{Bernoulli}(p_i)`, `\mathrm{Poisson}(E_i \lambda_i)`,
  `\mathrm{NegBin}(E_i \mu_i, \kappa)`, `\mathrm{Gamma}(\phi, \phi / \mu_i)`, zero-inflated
  `\pi\,\delta_0 + (1 - \pi)\,\mathrm{Poisson}(...)` (type 0: `\mathrm{Poisson}^{+}`);
  gaussian (identity) / lognormal (`\log y`) / t: residual `\varepsilon^{(r)}_i` with
  `\mathcal{N}` / `t_{\nu}`; `inlaFamily()`: generic `g_r(E[lhs])` + `lhs \sim \mathrm{name}`.
  Hyper symbols `\tau^{(r)}_\varepsilon`, `\phi^{(r)}`, `\kappa^{(r)}`, `\nu^{(r)}`,
  `\rho^{(r)}`, `\pi^{(r)}`, `\theta^{(r)}_{\mathrm{name}}`. Standalone parsed formula without
  family: long = gaussian, outcome = generic link only.
- Associations: each `scaled` entry (standalone parsed outcome: only its written `scaled()`
  terms) -> `\gamma_{y,x}\,u^{(x \to y)}_{\mathrm{id}[i],t_i}`; outcome and long are matched
  by the VALUES of the long formula's groups + covariates (not by row). A whole model adds per
  pair the zero equation `0 = <fixed + RE terms of x, no index> - u^{(x \to y)}_{\mathrm{id},t},
  (id, t) \in D_y` and `u \sim N(0, \tau^{-1})`.
- Priors (`priorLatex()`): one line after the distribution of the block / family; tau symbol
  `\tau^{(r)}_{k,g}` (`, \quad k = 0, 1` for several effects), SD as `\left(\tau\right)^{-1/2}`:
  `\mathrm{Gamma}(a, b)`, `\mathrm{Exp}(\lambda), \quad P(... > u) = \alpha`,
  `\mathcal{N}^{+}(0, (prec)^{-1})`, `\mathrm{C}^{+}(0, s)`,
  `\mathbf{Q}^{(r)}_{g} \sim \mathcal{W}_{k}(df, pmatrix | s\,\mathbf{I}_k)`; numbers with 3
  significant digits (`numLatex()`, `5 \cdot 10^{-5}`). No prior given: no line.
- Line order: long equations, zero equations, outcome equations, then distributions (family /
  residual / RE, long before outcome; `|` block -> MVN `N(0, Q^{-1})`, else univariate
  normals), `u` lines last. All normals use the precision parametrisation (as INLA):
  `N(0, (\tau)^{-1})`, residual `\tau_\varepsilon`; no `\sigma` / `\Sigma`.

## Open issues and placeholders

Defaults that are placeholders — don't build on them or tune them without the user:
- RE blocks without `prior`: no `hyper` (INLA default)
- scaled copy: `beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.01))`
- zero-equation field: `zeroFieldHyper` = `list(prec = list(initial = -6, fixed = TRUE))`
  (fixed diffuse precision of `u`), not user-settable. INLA's default (estimated precision)
  shrinks `u` and the long predictor to 0 (degenerate fits); the fixed diffuse value recovered
  simulated values in the hand-assembled fit
- zero-equation likelihood: `zeroControlFamily` = `list(hyper = list(prec = list(initial = 10,
  fixed = TRUE)))`, not user-settable

Still open:
- priors in `scaled()` (reuse the prior constructors)
- survival families; Gaussian observation weights (INLA `scale`)
- the pass-through of named `re()` / `scaled()` arguments into `f()` stays as is

## Development

- Tests: `Rscript -e 'devtools::test()'`; docs: `Rscript -e 'roxygen2::roxygenise()'`.
- Code style: camelCase names, roxygen2 docs with markdown.
- Dependencies: Imports `stats`, `utils`, `lme4`; Suggests `testthat`, `INLA`
  (`Additional_repositories` set). Avoid `lme4::findbars()` (moved to reformulas, warns); use
  the parsed `re` entries instead.
- `README.md`: short info, installation (INLA repository, then
  `remotes::install_github("LUMC/mvrem", ref = "devel-ai")` or a local clone), one example
  whose printed output is pasted in — re-run it and update the output after API changes.
- Examples: `CODE_EXAMPLES.Rmd`, rendered to `CODE_EXAMPLES.html` with
  `Rscript -e 'rmarkdown::render("CODE_EXAMPLES.Rmd")'`; update and re-render after API
  changes. Chunks printing INLA formulas use `class.output = "r", comment = ""`. The HTML is
  fully self-contained (no CDN): math via vendored KaTeX 0.16.11 in `tools/katex/` (js, css
  with woff2-only fonts, `katex-header.html` render script; MIT `LICENSE`), included with
  `pandoc_args: --katex` + `includes: in_header`. Don't switch back to `math_method:
  katex/mathjax` — rmarkdown marks those as external (CDN), giving tofu glyphs offline.
- Presentation: `presentation/mvrem-demo.qmd` (Quarto reveal.js; simple model -> INLA formula
  -> stacked-data tile plot -> LaTeX -> families/priors -> hand-assembled fit (200 subjects;
  with 40, gamma is noticeably attenuated for some seeds) -> status -> Extras incl. the dev
  container). Render with `quarto render presentation/mvrem-demo.qmd`; fully self-contained
  HTML (`embed-resources: true`). Math: `html-math-method: katex` only makes pandoc emit raw
  TeX in `span.math` (`plain` converts TeX to HTML — math is lost); Quarto's own KaTeX loader
  loads js/css from the `url` at run time (not embedded), so `presentation/katex-header.html`
  embeds the vendored KaTeX (js, css, woff2 fonts) and renders the spans itself, removing the
  `math` class so Quarto's loader finds nothing. Every slide heading has a short explicit
  id (`{#step3}`, used to refer to slides; new slides need one); `slide-ids.html` shows it as
  `[id]` left of the slide number (title slide: `[title]`). Re-render after API changes.
