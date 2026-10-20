test_that("format_report_pvals keeps small p-values", {
  x <- data.frame(
    method = c("a", "b", "c"),
    b = c(0.5, 0.1, 0.2),
    pval = c(1e-50, 0.0123, NA),
    Q_pval = c(0.5, 3.2e-12, 1)
  )
  out <- format_report_pvals(x)
  expect_equal(out$pval, c("1e-50", "0.0123", NA))
  expect_equal(out$Q_pval, c("0.5", "3.2e-12", "1"))
  expect_equal(out$b, x$b)
  expect_equal(out$method, x$method)
  expect_null(format_report_pvals(NULL))
  expect_false(grepl("\\| *0 *\\|", paste(knitr::kable(out), collapse = "\n")))
})
