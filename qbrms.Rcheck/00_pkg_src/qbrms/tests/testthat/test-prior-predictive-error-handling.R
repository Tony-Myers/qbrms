library(testthat)

test_that("prior predictive sampling fails gracefully when variables missing", {
  expect_error(
    qbrms::qbrms(y ~ x, data = data.frame(a=1:5, b=6:10), family = gaussian()),
    "Variables not found in data: y, x"
  )
})

test_that("prior predictive sampling errors on unknown family", {
  expect_error(
    qbrms::qbrms(y ~ x, data = data.frame(y=rnorm(10), x=rnorm(10)), family = "fakefamily"),
    "Family 'fakefamily' not yet supported in qbrms"
  )
})
