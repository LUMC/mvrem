# mvrem

Multivariate Joint Models with INLA.

`mvrem` provides tools for constructing and fitting Multivariate Joint Models using the
Integrated Nested Laplace Approximation (INLA) approach. It facilitates the specification of
complex hierarchical models with multiple longitudinal and outcome processes, handling both
fixed and random effects.

Models are written in familiar lme4-style formulae, extended with:

- `re(1 + t | id, ...)` — random effects with extra (prior) parameters; `|` correlated,
  `||` uncorrelated; plain `(1 + t | id)` is also accepted;
- `scaled(x1, ...)` — in an outcome formula, the (scaled) linear predictor of the
  longitudinal submodel `x1`.

In both, the other named `...` arguments are passed to INLA's `f()`, e.g. `scaled(x1, hyper = list(beta =
list(prior = "normal", param = c(0, 0.01))))`.

> **Status:** under development. Model specification, data preparation and construction of
> the INLA formula are implemented; fitting the model (`jsm()`) is not available yet, and the
> default priors are not finalised.

## Installation

The package is developed on GitHub: <https://github.com/LUMC/mvrem> (current development
branch: [`devel-ai`](https://github.com/LUMC/mvrem/tree/devel-ai)).

INLA is not on CRAN; it is used to validate likelihood families and (in the future) to fit
models. Install it first from the INLA repository:

```r
install.packages("INLA",
                 repos = c(getOption("repos"),
                           INLA = "https://inla.r-inla-download.org/R/stable"),
                 dep = TRUE)
```

Then install `mvrem` from the `devel-ai` branch:

```r
# install.packages("remotes")
remotes::install_github("LUMC/mvrem", ref = "devel-ai")
```

or from a local clone:

```sh
git clone -b devel-ai https://github.com/LUMC/mvrem.git
```

```r
devtools::install("mvrem")
```

### Development container (no local INLA installation)

Instead of installing INLA locally, you can use the development container in
[`.devcontainer/`](.devcontainer/). Its image (based on `rocker/geospatial`, R 4.5) comes
with INLA, `lme4` and `devtools` already installed:

1. clone the repository (see above) and open it in VS Code with the *Dev Containers*
   extension;
2. run **Dev Containers: Reopen in Container**. The image is built on first use, and the
   repository is mounted at `/home/rstudio/project`. RStudio Server is also available on
   port 8787;
3. inside the container, load the package with `devtools::load_all()` or install it with
   `devtools::install()`.

## Example

A beta-distributed outcome `y` associated with the trajectory of a longitudinal marker `x1`,
with priors on the random effects:

```r
library(mvrem)

# simulated data: 50 subjects, 4 visits each
set.seed(1)
n <- 50
d <- data.frame(id = sprintf("S%03d", rep(seq_len(n), each = 4)),
                t  = rep(c(0, 0.5, 1, 1.5), n))
b <- rnorm(n, sd = 0.5)[as.integer(factor(d$id))]
d$x1 <- 1 + 0.5 * d$t + b + rnorm(nrow(d), sd = 0.2)        # longitudinal marker
d$y  <- plogis(-1 + 0.3 * d$t + 0.8 * (d$x1 - 1) +
               rnorm(nrow(d), sd = 0.3))                     # outcome in (0, 1)

# 1. parse the formulae
forms <- collectFormulas(
  formOut   = list(y ~ 1 + t + re(1 | id, prior = pcSd(u = 1)) + scaled(x1)),
  formLong  = list(x1 ~ 1 + t + re(1 + t | id,
                                  prior = wishartPrec(df = 4, scale = diag(c(1, 0.1))))),
  familyOut = c(y = "beta")
)

# 2. assign the data to the formulae (variables selected, incomplete rows dropped)
dat <- prepareData(dataOut = d, dataLong = d, forms = forms)

# 3. build the INLA formula and the stacked data
inl <- buildInlaFormula(forms, dat)
formatInlaFormula(inl, nrows = 2)
```

```r
MX ~ -1 +
   #
   # ----- y: y ~ 1 + t + re(1 | id, prior = pcSd(u = 1)) + scaled(x1) -----
   #
   #    y_y y_1 y_t   y_1_id   y_FROM_x1_id_t
   #  0.241   1 0.0 y|1|S001   x1|TO|y|S001|0
   #  0.294   1 0.5 y|1|S001 x1|TO|y|S001|0.5
   #
y_1 +
y_t +
f(y_1_id, model = "iid", values = LEVELS$y_1_id, n = length(LEVELS$y_1_id), hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
f(y_FROM_x1_id_t, copy = "x1_TO_y_id_t", hyper = list(beta = list(fixed = FALSE, prior = "normal", param = c(0, 0.01)))) +
   #
   # ----- x1: x1 ~ 1 + t + re(1 + t | id, prior = wishartPrec(df = 4, scale = diag(c(1, 0.1)))) -----
   #
   #  x1_x1 x1_x1_0eq x1_1 x1_t   x1_1_id   x1_t_id x1_WEIGHT     x1_TO_y_id_t
   #  0.766        NA    1  0.0 x1|1|S001 x1|t|S001        NA             <NA>
   #  0.814        NA    1  0.5 x1|1|S001 x1|t|S001        NA             <NA>
   #     NA         0    1  0.0 x1|1|S001 x1|t|S001        -1   x1|TO|y|S001|0
   #     NA         0    1  0.5 x1|1|S001 x1|t|S001        -1 x1|TO|y|S001|0.5
   #
x1_1 +
x1_t +
f(x1_1_id, model = "iid2d", values = LEVELS$x1_1_t_id, n = length(LEVELS$x1_1_t_id), hyper = list(theta1 = list(prior = "wishart2d", param = c(4, 1, 10, 0)))) +
f(x1_t_id, x1_t, copy = "x1_1_id") +
f(x1_TO_y_id_t, x1_WEIGHT, model = "iid", values = LEVELS$x1_TO_y_id_t, hyper = list())
```

The model equations as LaTeX (body of an `aligned` environment):

```r
latexEquation(forms, dat)
```

$$
\begin{aligned}
x_{1,i} &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t_{i} + b^{(x_{1})}_{0,\mathrm{id}[i]} + b^{(x_{1})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{1})}_{i} \\
0 &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t + b^{(x_{1})}_{0,\mathrm{id}} + b^{(x_{1})}_{1,\mathrm{id}}\, t - u^{(x_{1} \to y)}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{y} \\
g_{y}\left(\mathrm{E}\left[y_{i}\right]\right) &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + b^{(y)}_{0,\mathrm{id}[i]} + \gamma_{y,x_{1}}\,u^{(x_{1} \to y)}_{\mathrm{id}[i],t_{i}} \\
\varepsilon^{(x_{1})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{1})}_{\varepsilon}\right)^{-1}\right) \\
\left(b^{(x_{1})}_{0,l}, b^{(x_{1})}_{1,l}\right)^\top &\sim \mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(x_{1})}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
\mathbf{Q}^{(x_{1})}_{\mathrm{id}} &\sim \mathcal{W}_{2}\left(4, \begin{pmatrix}1 & 0 \\ 0 & 0.1\end{pmatrix}\right) \\
y_{i} &\sim \mathrm{beta} \\
b^{(y)}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(y)}_{0,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
\left(\tau^{(y)}_{0,\mathrm{id}}\right)^{-1/2} &\sim \mathrm{Exp}\left(4.61\right), \quad P\left(\left(\tau^{(y)}_{0,\mathrm{id}}\right)^{-1/2} > 1\right) = 0.01 \\
u^{(x_{1} \to y)}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{1} \to y}\right)^{-1}\right)
\end{aligned}
$$

More examples (data lists, missing values, naming schemes, validation errors) are in
`CODE_EXAMPLES.Rmd` / `CODE_EXAMPLES.html`.

## Priors of random effects
The random effects are assumed to follow a Gaussian distribution $mathbf{b}_i \sim N(0, $\methbf{D}$)$, where $\mathbf{D}$ denotes the variance-covariance matrix of the random effects. In the case of a single random effect $k = 1$ (e.g., a random intercept or a random slope only), $\mathbf{D}$ reduces to a scalar variance. 
The prior of a random-effect variance-covariance matrix $\mathbf{D}$ is given in `re()` with `prior =` and one of these
functions. Each name says the distribution and the quantity it applies to (precision
$\tau$, standard deviation $\sigma = \tau^{-1/2}$, or the precision matrix $\mathbf{Q}$):

| Function | Prior | Use for | Default parameters
|---|---|---|---|
| `gammaPrec(shape, rate, fixed, initial)` | $\tau \sim \mathrm{Gamma}(shape, rate)$ | single / `\|\|` effects | shape = 1; rate = 5e-5; fixed = FALSE; initial = 4
| `pcSd(u)` | PC prior, $P(\sigma > u) = \alpha$ | single / `\|\|` effects | alpha = 0.01
| `halfNormalSd(scale)` | $\sigma \sim \mathcal{N}^{+}(0, scale^2)$ | single / `\|\|` effects |
| `halfCauchySd(scale)` | $\sigma \sim \mathrm{C}^{+}(0, scale)$ | single / `\|\|` effects |
| `wishartPrec(df, scale matrix)` | $\mathbf{Q} \sim \mathcal{W}_k(df, \mathbf{V})$ | correlated (`\|`) blocks | df = k + 2; scale = diag(1, 2)

The option `fixed` determines if a fixed value for the parameter is assumed. If `fixed = TRUE`, then the parameter is fixed to the `initial` value.

```r
x1 ~ 1 + t + re(1 | id, prior = pcSd(u = 1))                  # P(sigma > 1) = 0.01
x1 ~ 1 + t + re(1 + t || id, prior = halfNormalSd(scale = 2))  # same prior for both effects
x1 ~ 1 + t + re(1 + t | id, prior = wishartPrec(df = 4, scale = diag(c(1, 0.1))))
```

Blocks without `prior` (including plain `(1 | id)`) use INLA's default prior.

## Prior of the scaling parameter
In INLA formulation, the scaling parameter is also a random effects for which its assumed a Gaussian distribution. By defult $\gamma \sim N(0, 0.0001)$. The code formulation can be `normal(mean, precision, fixed, initial)`.

