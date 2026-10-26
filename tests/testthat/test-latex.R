forms <- collectFormulas(
  formOut = list(z ~ 1 + t + (1 + t | id),
                 log(y) ~ 1 + t + re(1 + t | id) + scaled(x1) + scaled(x2)),
  formLong = list(x1 ~ 1 + t + re(0 + t | id) + re(1 | id),
                  x2 ~ t - 1 + re(1 | id),
                  x3 ~ 1 + I(t^2) + g + age:t + re(1 + t || site:id)),
  familyOut = c(y = "beta")
)
d <- data.frame(id = rep(1:3, each = 2), site = rep(c("a", "b"), 3), t = rep(0:1, 3),
                g = factor(rep(c("u", "v", "w"), 2)), age = 1:6,
                x1 = 1, x2 = 1, x3 = 1, z = 1, y = 0.5)
dat <- prepareData(d, d, forms)

lines <- function(x) strsplit(unclass(x), " \\\\\n", fixed = TRUE)[[1]]

test_that("single longitudinal formula", {
  tex <- latexEquation(forms$long$x2)
  expect_s3_class(tex, "mvremLatex")
  expect_match(tex, "^\\\\begin\\{aligned\\}\n")
  expect_match(tex, "\n\\\\end\\{aligned\\}$")
  expect_match(tex, "x_{2,i} &= \\beta^{(x_{2})}_{1}\\, t_{i} + b^{(x_{2})}_{2,\\mathrm{id}[i]}",
               fixed = TRUE)                      # no intercept: numbering starts at 1
  expect_match(tex, "\\varepsilon^{(x_{2})}_{i}", fixed = TRUE)
  expect_match(tex, "\\mathcal{N}\\left(0, \\left(\\tau^{(x_{2})}_{2,\\mathrm{id}}\\right)^{-1}\\right)",
               fixed = TRUE)                      # precision parametrisation
  expect_match(tex, "\\left(\\tau^{(x_{2})}_{\\varepsilon}\\right)^{-1}", fixed = TRUE)
  expect_false(grepl("sigma", tex))
})

test_that("random effects share the index of the fixed effect of the same column", {
  tex <- latexEquation(forms$long$x1)
  expect_match(tex, "b^{(x_{1})}_{1,\\mathrm{id}[i]}\\, t_{i} + b^{(x_{1})}_{0,\\mathrm{id}[i]}",
               fixed = TRUE)
  # residual + two uncorrelated blocks -> three univariate normals
  expect_length(lines(tex), 4)
})

test_that("correlated block gives a multivariate normal", {
  tex <- latexEquation(forms$out$z)
  expect_match(tex, "\\left(b^{(z)}_{0,l}, b^{(z)}_{1,l}\\right)^\\top &\\sim", fixed = TRUE)
  expect_match(tex, "\\mathcal{N}\\left(\\mathbf{0}, \\left(\\mathbf{Q}^{(z)}_{\\mathrm{id}}\\right)^{-1}\\right)",
               fixed = TRUE)
})

test_that("outcome: link, transformed response and scaled terms", {
  tex <- latexEquation(forms$out$y)
  expect_match(tex, "g_{y}\\left(\\mathrm{E}\\left[\\log\\left(y_{i}\\right)\\right]\\right) &=",
               fixed = TRUE)
  expect_match(tex, "\\gamma_{y,x_{1}}\\,u^{(x_{1} \\to y)}_{i} + \\gamma_{y,x_{2}}\\,u^{(x_{2} \\to y)}_{i}",
               fixed = TRUE)
  expect_false(grepl("varepsilon", tex))
  expect_match(tex, "\\log\\left(y_{i}\\right) &\\sim \\mathrm{beta}", fixed = TRUE)
})

test_that("gaussian family: identity and residual; standalone formulae", {
  tex <- latexEquation(forms$out$z)                 # gaussian by default
  expect_match(tex, "^\\\\begin\\{aligned\\}\nz_\\{i\\} &= ")
  expect_match(tex, "\\varepsilon^{(z)}_{i} &\\sim \\mathcal{N}", fixed = TRUE)

  p <- parseJsmFormula(y ~ t, part = "out")         # no family known: generic link
  expect_match(latexEquation(p), "g_{y}", fixed = TRUE)
  expect_false(grepl("\\sim", latexEquation(p)))
  p <- parseJsmFormula(x ~ t, part = "long")        # longitudinal: gaussian
  expect_match(latexEquation(p), "\\varepsilon", fixed = TRUE)

  fp <- collectFormulas(y ~ t, n ~ t, familyLong = c(n = "poisson"))
  tex <- latexEquation(fp$long$n)
  expect_match(tex, "g_{n}\\left(\\mathrm{E}\\left[n_{i}\\right]\\right)", fixed = TRUE)
  expect_match(tex, "n_{i} &\\sim \\mathrm{poisson}", fixed = TRUE)
})

test_that("zero equations: linear predictor of the long submodel minus u", {
  l <- lines(latexEquation(forms, dat))
  z1 <- l[grepl("^0 &= ", l)][1]
  # indexed by the values of the matching variables (id, t), not by observation
  expect_equal(z1, paste0(
    "0 &= \\beta^{(x_{1})}_{0} + \\beta^{(x_{1})}_{1}\\, t + b^{(x_{1})}_{1,\\mathrm{id}}\\, t",
    " + b^{(x_{1})}_{0,\\mathrm{id}} - u^{(x_{1} \\to y)}_{\\mathrm{id},t},",
    " \\quad \\left(\\mathrm{id}, t\\right) \\in \\mathcal{D}_{y}"))
  # the outcome uses u at its own values of id and t
  expect_match(l[grepl("^g_\\{y\\}", l)],
               "\\gamma_{y,x_{1}}\\,u^{(x_{1} \\to y)}_{\\mathrm{id}[i],t_{i}}", fixed = TRUE)
  expect_false(grepl("varepsilon", z1))
  # no zero equations without scaled() / for a single formula
  f0 <- collectFormulas(y ~ t, x ~ t)
  expect_false(any(grepl("^0 &=", lines(latexEquation(f0)))))
  expect_false(grepl("0 &=", latexEquation(forms$out$y), fixed = TRUE))
})

test_that("covariate expressions and data-based factor expansion", {
  tex <- latexEquation(forms$long$x3)
  expect_match(tex, "\\beta^{(x_{3})}_{1}\\, {t_{i}}^{2}", fixed = TRUE)
  expect_match(tex, "\\beta^{(x_{3})}_{2}\\, g_{i}", fixed = TRUE)
  expect_match(tex, "\\mathit{age}_{i} \\cdot t_{i}", fixed = TRUE)
  expect_match(tex, "\\mathrm{site:id}[i]", fixed = TRUE)

  texD <- latexEquation(forms$long$x3, dat$long$x3)
  expect_match(texD, "\\beta^{(x_{3})}_{2}\\, \\mathbb{I}\\left[g_{i} = \\text{v}\\right]",
               fixed = TRUE)
  expect_match(texD, "\\beta^{(x_{3})}_{3}\\, \\mathbb{I}\\left[g_{i} = \\text{w}\\right]",
               fixed = TRUE)
  # slope t has no fixed counterpart -> next free index
  expect_match(texD, "b^{(x_{3})}_{5,\\mathrm{site:id}[i]}\\, t_{i}", fixed = TRUE)
})

test_that("whole model: equations first, then distributions", {
  tex <- latexEquation(forms, dat)
  l <- lines(tex)
  expect_length(l, 5 + 2 + 1 + 4 + 7 + 2)  # equations, zero eqs, family, residuals, REs, u
  # longitudinal equations, zero equations, outcome equations, distributions
  expect_match(l[1], "\nx_\\{1,i\\} &=")
  expect_match(l[3], "^x_\\{3,i\\} &=")
  expect_match(l[4:5], "^0 &= ")
  expect_match(l[6], "z_{i} &=", fixed = TRUE)
  expect_match(l[7], "^g_\\{y\\}")
  expect_true(all(grepl("\\\\sim", l[8:21])))
  expect_match(l[8], "^\\\\varepsilon\\^\\{\\(x_\\{1\\}\\)\\}")   # long distributions first
  expect_match(l[20:21], "\\\\left\\(\\\\tau_\\{x_\\{[12]\\} \\\\to y\\}\\\\right\\)\\^\\{-1\\}")
  expect_match(l[20:21], "^u\\^\\{\\(x_\\{[12]\\} \\\\to y\\)\\}_\\{\\\\mathrm\\{id\\},t\\} &\\\\sim")
  expect_equal(latexEquation(forms), latexEquation(forms, NULL))
})

test_that("format wraps in delimiters; invalid input errors", {
  tex <- latexEquation(forms$out$z)
  expect_match(format(tex, delim = "$$"), "^\\$\\$\\\\begin.*\\\\end\\{aligned\\}\\$\\$$")
  expect_output(print(tex), "begin{aligned}", fixed = TRUE)
  expect_error(latexEquation(list(a = 1)), "parsed formula")
  expect_error(latexEquation(forms$out$z, dat), "data frame")
  expect_error(latexEquation(forms, d), "prepareData")
})

test_that("variable names and expressions", {
  expect_equal(latexName("t", "i"), "t_{i}")
  expect_equal(latexName("x12", "i"), "x_{12,i}")
  expect_equal(latexName("my_var", "i"), "\\mathit{my\\_var}_{i}")
  expect_equal(exprLatex(quote(sqrt(t) / exp(u))), "\\frac{\\sqrt{t_{i}}}{\\exp\\left(u_{i}\\right)}")
  expect_equal(exprLatex(quote(ns(t, 3))), "\\operatorname{ns}\\left(t_{i}, 3\\right)")
})
