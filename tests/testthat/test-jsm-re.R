test_that("re captures tail arguments without evaluating them", {
  translated <- mvrem:::.translateRe(
    y ~ re(1 + t | id, par = 1L, custom = unresolved_value)
  )

  expect_true(is.list(translated$params))
  expect_identical(
    translated$params$id$par,
    1L
  )
  expect_identical(
    translated$params$id$custom,
    quote(unresolved_value)
  )
})