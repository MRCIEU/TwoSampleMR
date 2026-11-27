# The fixed effects standard errors in mr_rucker() should equal the
# closed-form fixed effects values, (X' W X)^-1 with W = 1 / se.outcome^2,
# whatever the residual standard error. These data have a residual standard
# error of about 0.70 in both the IVW and the Egger model (under-dispersion),
# where the fixed effects SEs used to be too small by that factor (#630).

bx <- c(0.07, 0.119, 0.142, 0.078, 0.06, 0.12, 0.103, 0.131, 0.146, 0.061)
by <- c(0.0159, 0.0382, 0.0441, 0.0255, 0.0268, 0.0726, 0.0252, 0.0308, 0.0552, 0.003)
se_out <- rep(0.02, 10)
underdispersed_dat <- data.frame(
  SNP = paste0("rs", 1:10),
  exposure = "X",
  outcome = "Y",
  id.exposure = "X",
  id.outcome = "Y",
  beta.exposure = bx,
  se.exposure = 0.01,
  beta.outcome = by,
  se.outcome = se_out,
  mr_keep = TRUE
)

test_that("mr_rucker() fixed effects SEs match the closed form when the residual SE is below 1", {
  r <- suppressMessages(mr_rucker(underdispersed_dat))[[1]]
  expect_lt(summary(r$lmod_ivw)$sigma, 1)
  expect_lt(summary(r$lmod_egger)$sigma, 1)

  w <- 1 / se_out^2
  se_ivw_fe <- 1 / sqrt(sum(w * bx^2))
  X <- cbind(1, bx)
  V_egger <- solve(t(X) %*% (w * X))

  rk <- r$rucker
  expect_equal(rk$SE[rk$Method == "IVW fixed effects"], se_ivw_fe)
  expect_equal(rk$SE[rk$Method == "Egger fixed effects"], sqrt(V_egger[2, 2]))
  expect_equal(r$intercept$SE[r$intercept$Method == "Egger fixed effects"], sqrt(V_egger[1, 1]))

  # Same as the package's own fixed effects IVW
  fe <- mr_ivw_fe(bx, by, underdispersed_dat$se.exposure, se_out)
  expect_equal(rk$SE[rk$Method == "IVW fixed effects"], fe$se)
})
