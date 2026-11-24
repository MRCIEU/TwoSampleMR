context("LDSC")

test_that("the LDSC functions are placeholders until opengwasr can serve them (opengwas/opengwasr#26)", {
  expect_error(ldsc_h2("ukb-b-19953"), class = "twosamplemr_not_yet_available")
  expect_error(
    ldsc_rg("ukb-b-19953", "ebi-a-GCST003116"),
    "opengwas/opengwasr#26",
    class = "twosamplemr_not_yet_available"
  )
})
