# Tests for mr_rucker_bootstrap().
#
# 1. mr_rucker() returns a list with one element per exposure-outcome combination,
#    but mr_rucker_bootstrap() accesses the per-combination result directly
#    ($rucker/$Q/$res/$selected). It must unwrap [[1]], otherwise it errors with
#    "values must be length 1, but FUN(X[[1]]) result is length 0".
#
# 2. Regression test for issue #684: the bootstrap draws are pre-generated as a
#    matrix and were filled in a way that scrambled which SNP each draw came from,
#    inflating the bootstrapped standard errors. The effect is worst when the
#    number of SNPs divides nboot; with nsnp = 50 and nboot = 1000 the Rucker
#    mean/median SE was inflated from ~0.06-0.09 to ~0.50.

load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("mr_rucker_bootstrap() runs and returns the expected structure", {
  set.seed(1234)
  rb <- suppressMessages(mr_rucker_bootstrap(dat))

  expect_named(
    rb,
    c("rucker", "res", "bootstrap_estimates", "boostrap_q", "q_plot", "e_plot")
  )
  # $rucker is the unwrapped single-combination result, not the combo list.
  expect_true(all(c("rucker", "Q", "res", "selected") %in% names(rb$rucker)))
  expect_equal(nrow(rb$bootstrap_estimates), 1000L)
  expect_true(all(c("Rucker mean", "Rucker median") %in% rb$res$Method))
})

test_that("mr_rucker_bootstrap() SEs are not inflated when nsnp divides nboot (issue #684)", {
  # nsnp = 50 divides the default nboot = 1000 - the worst case for the bug.
  dat50 <- dat[1:50, ]
  set.seed(1234)
  rb <- suppressMessages(mr_rucker_bootstrap(dat50))

  mean_se <- rb$res$SE[rb$res$Method == "Rucker mean"]
  median_se <- rb$res$SE[rb$res$Method == "Rucker median"]

  # Correct SEs are ~0.06-0.09; the bug inflated them to ~0.50.
  expect_lt(mean_se, 0.25)
  expect_lt(median_se, 0.25)
})

test_that("mr_rucker_cooksdistance() runs and returns a well-formed result", {
  # Like mr_rucker_bootstrap(), this accessed the per-combination result
  # directly ($cooksdistance/$selected/$rucker) while mr_rucker() returns a
  # combo list, so the Cook's distance loop silently never ran and a malformed
  # object was returned. It must unwrap [[1]].
  set.seed(1234)
  cd <- suppressMessages(mr_rucker_cooksdistance(dat))

  expect_true(is.data.frame(cd$selected))
  expect_true(is.data.frame(cd$rucker))
  expect_equal(cd$selected$Method, "Rucker (CD)")
})
