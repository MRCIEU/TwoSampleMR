test_that("format_d() removes duplicated proxy SNPs within each outcome without warning", {
  d <- data.frame(
    rsid = c("rs1", "rs2", "rs1", "rs2"),
    chr = "1",
    position = c(100, 200, 100, 200),
    beta = 0.1,
    se = 0.01,
    n = 1000,
    p = 1e-8,
    eaf = 0.3,
    ea = "A",
    nea = "G",
    trait = c("Trait A", "Trait A", "Trait B", "Trait B"),
    id = c("id-a", "id-a", "id-b", "id-b"),
    proxy = TRUE,
    target_snp = c("rs1", "rs2", "rs1", "rs2"),
    proxy_snp = "rs3",
    target_a1 = "A",
    target_a2 = "G",
    proxy_a1 = "A",
    proxy_a2 = "G"
  )
  expect_no_warning(out <- format_d(d))
  expect_equal(nrow(out), 2)
  expect_setequal(out$id.outcome, c("id-a", "id-b"))
})
