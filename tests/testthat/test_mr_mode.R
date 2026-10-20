# Regression tests for issue #684.
#
# In v0.6.30-v0.7.8 the mr_mode() bootstrap pre-generated its random draws with
#   matrix(rnorm(nboot * n, mean = BetaIV, sd = seBetaIV), nrow = nboot, ncol = n)
# rnorm() recycles the length-n mean/sd vectors element-wise, but matrix() fills
# column-by-column, so each column (SNP) did not keep its own mean/SE - the SNPs'
# parameters were scrambled across the matrix. This left the point estimates
# untouched but inflated the bootstrap standard errors (and hence p-values). The
# fix lays the parameters out with rep(..., each = nboot) so column j draws from
# SNP j. Correct SEs for this data are ~0.1; the bug inflated them ~3x (and much
# more when the number of SNPs divides nboot).

load(system.file("extdata", "test_commondata.RData", package = "TwoSampleMR"))

test_that("mr_mode() weighted/penalised bootstrap SEs are not inflated (issue #684)", {
  set.seed(1234)
  m <- mr_mode(dat)

  # Point estimates are unaffected by the bug.
  expect_equal(m$b[m$method == "Simple mode"], 0.340, tolerance = 1e-2)
  expect_equal(m$b[m$method == "Weighted mode"], 0.379, tolerance = 1e-2)

  # Weighted and penalised modes are corrupted at any number of SNPs because the
  # bootstrapped betas get mispaired with the fixed-order weights.
  expect_equal(m$se[m$method == "Weighted mode"], 0.105, tolerance = 1e-1)
  expect_equal(m$se[m$method == "Penalised mode"], 0.109, tolerance = 1e-1)

  # Anti-inflation ceiling: the buggy SEs were ~0.30+, the correct ones ~0.10.
  expect_lt(m$se[m$method == "Weighted mode"], 0.2)
  expect_lt(m$se[m$method == "Penalised mode"], 0.2)
})

test_that("mr_mode() simple-mode SE is not inflated when nsnp divides nboot (issue #684)", {
  # The simple mode is only corrupted when gcd(nsnp, nboot) > 1. Use nsnp = 50,
  # which divides the default nboot = 1000 - the worst case, where the buggy
  # fill made every bootstrap draw come from a single SNP.
  dat50 <- dat[1:50, ]
  set.seed(1234)
  m <- mr_mode(dat50)

  # Correct simple-mode SE is ~0.17; the bug inflated it to ~0.57.
  expect_lt(m$se[m$method == "Simple mode"], 0.35)
})

test_that("mr_mode() removes SNPs with beta.exposure of 0 (issue #570)", {
  d <- dat
  d$beta.exposure[1:2] <- 0
  set.seed(1234)
  expect_warning(m <- mr_mode(d), "Removing 2 SNP\\(s\\) with beta.exposure of 0")
  set.seed(1234)
  m0 <- mr_mode(dat[-(1:2), ])
  expect_equal(m, m0)
  expect_true(all(m$nsnp == 77))

  expect_warning(
    w <- mr_weighted_mode(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome),
    "Removing 2 SNP\\(s\\)"
  )
  expect_true(is.finite(w$b))
  expect_equal(w$nsnp, 77)
})

test_that("mr_mode() returns NA when fewer than 3 SNPs have non-zero beta.exposure (issue #570)", {
  d <- dat[1:4, ]
  d$beta.exposure[1:2] <- 0
  warnings <- character()
  w <- withCallingHandlers(
    mr_simple_mode(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome),
    warning = function(cnd) {
      warnings <<- c(warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(w, list(b = NA, se = NA, pval = NA, nsnp = NA))
  expect_match(warnings[1], "Removing 2 SNP\\(s\\)")
  expect_match(warnings[2], "Need at least 3 SNPs")
})
