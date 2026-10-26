# mvrem

R package under development.

## Description

Provides tools for constructing and fitting Multivariate Joint Models using the
Integrated Nested Laplace Approximation (INLA) approach. The package facilitates
the specification of complex hierarchical models with multiple longitudinal and
outcome processes, handling both fixed and random effects.

## User-facing API (target)

```r
jsmMdl <- jsm(
  formOut = list(
    z ~ 1 + t + (1 + t|id),        # plain lme4 notation -> default priors
    y ~ 1 + t + re(1 + t|id, ...)  # correlated intercept/slope RE; ... = distribution/hyperpriors
              + scaled(x1, ...)    # association with longitudinal process x1; ... = extra params
              + scaled(x2, ...)    # scaled() allowed ONLY in outcome formulae
  ),
  formLong = list(
    x1 ~ 1 + t + re(0 + t|id, ...) + re(1|id, ...),  # two independent RE blocks
    x2 ~ 1 + t + re(1 + t|id, ...),
    x3 ~ 1 + t + re(1 + t||id, ...)                  # `||` = uncorrelated REs (lme4 semantics)
  ),
  dataLong = d,
  dataOut = d,
  familyOut = c(y = "beta"),   # named; INLA likelihood names; unspecified -> "gaussian"
  familyLong = c(x3 = "poisson"),
  fit = FALSE                  # FALSE: build model object only, don't call INLA
)
```

Syntax rules:
- `(... | id)` is lme4 notation, equivalent to `re(... | id)` with default priors.
- `re(... | g, ...)` random effects; `|` correlated, `||` uncorrelated.
- `scaled(x, ...)` (bare name, no `~`) links an outcome to longitudinal submodel `x` via an
  estimated scaling coefficient. Only valid in `formOut`; target must be a response in `formLong`.
- Extra arguments in `re()` / `scaled()` must be named; evaluated in the formula's environment.
  They are merged into the default `f()` arguments (`modifyList`, nested lists like `hyper`
  merged). User parameters are given ONLY this way — no `userPars` argument keyed by INLA names.
- `re()` / `scaled()` must be top-level additive terms (not nested in other calls).
- LHS must contain exactly one variable (`log(y)` ok, `cbind(a, b)` not); that variable is the
  response name used everywhere (list names, `scaled()` targets, data list names).

## Pipeline

`collectFormulas()` (step 1) -> `prepareData()` (step 2) -> `buildInlaFormula()` (step 3)
-> `formatInlaFormula()` (display); `latexEquation()`; `summarizeData()`. Not yet done: `jsm()` wrapper, INLA
call (`family = inl$families`, `control.family`, zero-eq precision), fitting, multi-column
responses (e.g. survival `inla.surv`, binomial `Ntrials`).

### Step 1: formula collection — `R/parseFormula.R`

`collectFormulas(formOut, formLong, familyOut = NULL, familyLong = NULL)` returns
`list(out, long, variables)`; `out`/`long` named
by response; `variables` = data frame (`variable`, `role`, `part`, `formula`) over all formulae.
Per formula (`parseJsmFormula(formula, part)`):
- `response` (single LHS variable), `lhs` (LHS as written), `part` (`"out"` / `"long"`)
- `formula`: as given by the user
- `fixedFormula`: fixed part only (no REs, no `scaled()`)
- `lme4Formula`: `re(bar, ...)` -> `(bar)`, `scaled()` removed, plain bar terms kept
- `variables`: `response` / `covariate` / `group` from the AST via `all.vars()` (`log(y)` -> `y`,
  `I(t^2)` -> `t`). Covariates: fixed terms + left of `|`; groups: right of `|`
  (`site:id` -> `site`, `id`). `scaled()` targets and params are not variables.
- `re`: one entry per RE block (`re()` and plain `( | )` alike): `term`, `group`,
  `correlated` (`FALSE` for `||`), `params`
- `scaled` (outcome only): `target`, `params`
- `family` (set by `collectFormulas`): from `familyOut` / `familyLong` — named character
  vectors (names = response variables of that part, values = INLA likelihood names,
  `inlaFamilies()` = `names(INLA::inla.models()$likelihood)`); unlisted -> `"gaussian"`.
  Validated only if INLA installed (else warning). Errors: unknown family, names not
  responses of that part (message says if it belongs to the other part), unnamed/duplicated.
- `linkVariables` (outcome only, set by `collectFormulas`): groups + covariates of its
  `scaled()` targets not already in the outcome formula (needed for zero-equation rows).

Errors: `scaled()` in long formula; unknown `scaled()` target; duplicate responses (also
`y` and `log(y)`); LHS without exactly one variable; non-bar first arg of `re()`; non-name
first arg of `scaled()`; unnamed extra args; nested `re()`/`scaled()`; one-sided formulae.

### Step 2: data preparation — `R/prepareData.R`

`prepareData(dataOut, dataLong, forms, na.action = na.omit)` returns `list(out, long)`: per
response one data frame. `dataOut` / `dataLong`: one data frame (shared by that part) or a
list named by response variable (`y` for `log(y) ~ ...`); list must match responses exactly.

Each returned data frame contains:
- `.key`: character, `"K"` + user's explicit row names (checked unique) or `"K"` + row number;
  same row -> same key across formulae (also for subsets of one data frame)
- the formula's variables (untransformed) + `linkVariables` for outcomes
- one factor column per grouping term, named as the term (`id`, `site:id`; replaces a
  same-named variable)

`lme4::glFormula()` on `lme4Formula` only selects rows and builds grouping factors (no design
matrices kept); lme4 nobs/nlev checks disabled; formulae without REs use `model.frame`.
`na.action` (function or name, e.g. `na.omit` / `na.fail`) applies per formula to its
variables (and outcome link variables); other columns ignored.

Errors: missing/extra/unnamed list elements, non-data-frame elements, formula variables
missing in data, user data containing `.key`, non-unique row names.

### Step 3: INLA formula — `R/inlaFormula.R`

`buildInlaFormula(forms, data, naming = getOption("mvrem.naming", namingScheme()))`.

Naming is configurable: all names come from `namingScheme(...)` templates (placeholders
`{resp}`, `{term}`, `{terms}`, `{group}`, `{level}`, `{long}`, `{out}`, `{vars}`, `{values}`;
names made syntactic; unknown templates and duplicated resulting names -> error). Default
scheme (user-provided reference, must be kept):
- response matrix `MX`; columns `y_y`, `x1_x1`, zero-eq `x1_x1_0eq`
- fixed `<resp>_<term>`, intercept `1` (`y_1`, `y_t`); RE slopes use this column as weight
- RE index `<resp>_<term>_<group>` (`y_1_id`, `y_t_id`), `values = LEVELS$<resp>_<terms>_<group>`
  (`y_1_t_id`), labels `"y|1|S001"`, `n = length(LEVELS$...)`
- zero-eq field `<long>_TO_<out>_<vars>` (`x1_TO_y_id_t`), weight `<long>_WEIGHT`, labels
  `"x1|TO|y|S001|0.5"`; `{vars}` = groups + covariates of the long formula
- scaled `<out>_FROM_<long>_<vars>` (`y_FROM_x1_id_t`)

Terms (formula grouped per response: outcomes first, then longitudinal):
- fixed: one column per design-matrix column (factor levels expand)
- random: `|` with k > 1 -> `iid<k>d` (`iidkd` + `order` for k > 5) on first effect, `copy`
  for the others; single effect or `||` -> one `iid` per effect. Default prior: INLA's.
- zero equation, per (long, out) pair from `scaled()`: `f(x1_TO_y_id_t, x1_WEIGHT,
  model = "iid", values = LEVELS$x1_TO_y_id_t, hyper = list())`, i.e. `0 = eta_x1 - u`
- scaled: `f(y_FROM_x1_id_t, copy = "x1_TO_y_id_t", hyper = list(beta = list(fixed = FALSE,
  prior = "normal", param = c(0, 0.01))))`

Returns `responses`, `zeroResponses`, `families` (family per response column: responses,
then zero responses = `"gaussian"`; order of INLA likelihoods), `fixed`, `random` (blocks with `group`, `lhs`,
`effects`, `index`, `weight`, `levels`, `terms`), `zero`, `scaled`, `levels` (= `LEVELS`),
`naming`, `terms` (per response), `formulas` (user's), `stacks`, `formula` (env holds `LEVELS`).

`stacks`: per response the stacked data with the formula's column names. `main` = observation
rows; longitudinal responses also have one zero-equation part per outcome using them
(response `NA`, `_0eq` = 0, WEIGHT = -1, covariates and RE labels from the outcome rows,
factor levels from the long data). Zero-eq rows: one per UNIQUE combination of `{vars}`
values in the outcome data (= `LEVELS$<zeroIndex>`); outcome rows are matched to them via these
values (not `.key`), so combinations without a long observation (missing values, subjects
absent from long data) still get a zero equation. Tested with artificially missing covariates.

### Data summary — `R/summary.R`

`summarizeData(forms, dataOut, dataLong)` (raw data, same input forms as `prepareData`) ->
class `mvremDataSummary` with `equations` (part, response, family, rows, complete, dropped,
groups = levels per grouping term among complete rows; long first), `missing` (per equation
and variable incl. link variables), `zero` (per `scaled()` pair: matching vars, number of
unique combos in complete outcome rows). `complete` = rows kept by `prepareData(na.omit)`.
Data matching shared with `prepareData` via `matchData()`.

### Display — `R/formatFormula.R`

`formatInlaFormula(inl, nrows = 4, digits = 3)` -> class `mvremFormulaText` (printed with `cat`): formula
one term per line, each response group preceded by `#` comments with the user formula and the
first `nrows` rows of each stack part (`nrows = 0` hides data; numbers shown with `digits`
significant digits). Parses back to the same formula.

### LaTeX — `R/latex.R`

`latexEquation(x, data = NULL)`: `x` = parsed formula (step 1 element) or `collectFormulas()`
result; `data` = matching step-2 data frame / `prepareData()` result (expands factors into
`\mathbb{I}[g_i = v]`; without data a factor is one term). Returns class `mvremLatex`: body of
an `aligned` environment; `format(x, delim = "$$")` wraps it. Notation: observation index `i`;
`\beta^{(r)}_k` (k = 0 intercept); `b^{(r)}_{k,g[i]}` shares k with the fixed effect of the
same column, else next free k; gaussian family: identity + `\varepsilon^{(r)}_i` (+ its
normal line); other families: generic link `g_r(E[lhs])` + line `lhs \sim \mathrm{family}`;
standalone parsed formula without family: long = gaussian, outcome = generic link only;
`scaled(x)` -> `\gamma_{y,x}\,u^{(x \to y)}_{\mathrm{id}[i],t_i}` — outcome and long are matched
by the VALUES of the long formula's groups + covariates (not by row index); whole model adds
a zero equation per pair `0 = <fixed + RE terms of x, no index> - u^{(x \to y)}_{\mathrm{id},t},
(id, t) \in D_y` and
`u ~ N(0, tau^{-1})`; line order: long equations, zero eqs, outcome equations, then
distributions (family/residual/RE; long before outcome), `u` lines last; RE
distributions after equations (`|` block -> MVN `N(0, Q^{-1})`, else univariate normals). All normals use the precision
parametrisation (as INLA): `N(0, (\tau)^{-1})`, residual `\tau_\varepsilon`; no `\sigma` / `\Sigma`.
Names: `t` -> `t_i`, `x1` -> `x_{1,i}`, other -> `\mathit{name}_i`. In the Rmd, rendered
math uses `results = "asis"` + `cat(format(..., delim = "$$"))`.

## Deferred: prior specification

Priors are not designed yet and will be developed in a future step. Until then:
- current prior-related defaults are placeholders only — don't build on them or tune them:
  RE blocks get no default `hyper` (INLA default); scaled terms use
  `beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.01))`; zero-eq field
  `hyper = list()`
- the reference output's `wishart2d` default (`param = c(3, 1, 0.004, 0)`) is not implemented
- zero-equation likelihood precision (`control.family`) is not set anywhere yet
- the pass-through of named `re()` / `scaled()` arguments into `f()` stays as is

## Development

- Tests: `Rscript -e 'devtools::test()'`; docs: `Rscript -e 'roxygen2::roxygenise()'`.
- Code style: camelCase names, roxygen2 docs with markdown.
- Dependencies: Imports `stats`, `utils`, `lme4`; Suggests `testthat`, `INLA`
  (`Additional_repositories` set). Avoid `lme4::findbars()` (moved to reformulas, warns); use the
  parsed `re` entries instead.
- `README.md`: short info, installation (devtools from source + INLA repo), one example whose
  printed output is pasted in — re-run it and update the output after API changes.
- Examples: `CODE_EXAMPLES.Rmd`, rendered to `CODE_EXAMPLES.html` via
  `Rscript -e 'rmarkdown::render("CODE_EXAMPLES.Rmd")'`; update and re-render after API
  changes. Chunks printing INLA formulas use `class.output = "r", comment = ""`.
  The HTML is fully self-contained (no CDN): math via vendored KaTeX 0.16.11 in `tools/katex/`
  (js, css with woff2-only fonts, `katex-header.html` render script; MIT `LICENSE`), included
  with `pandoc_args: --katex` + `includes: in_header`. Don't switch back to `math_method:
  katex/mathjax` — rmarkdown marks those as external (CDN), giving tofu glyphs offline.
