test_that("BLUPF90 preparation supports an arbitrary trait vector", {
  for (k in c(1L, 2L, 3L, 6L, 10L)) {
    set.seed(30 + k)
    p <- data.frame(ID = paste0("a", 1:30), Sire = NA, Dam = NA)
    d <- data.frame(ID = p$ID, Sex = rep(c("F", "M"), 15))
    traits <- paste0("trait", seq_len(k))
    for (t in traits) d[[t]] <- rnorm(30)
    job <- prepare_blupf90(p, d, traits, fixed = ~Sex)
    expect_equal(job$n_traits, k)
    expect_equal(dim(job$G_start), c(k, k))
    expect_equal(nrow(job$parameter_map), k * (k + 1L))
    expect_equal(sum(grepl("OPTION se_covar_function", job$parameter_lines)), k + k * (k - 1L) / 2L)
    expect_true(any(job$parameter_lines == "OPTION method VCE"))
    expect_true(any(job$parameter_lines == "OPTION sol se"))
    expect_false(any(job$data == job$missing_code))
  }
})
test_that("sparse inbreeding agrees with tabular relationships", {
  p <- data.frame(ID = letters[1:5], Sire = c(NA, NA, "a", "a", "c"), Dam = c(NA, NA, "b", "b", "d"))
  sparse <- pedigree_inbreeding(p)
  dense <- pedigree_relationship(p)
  expect_equal(sparse, dense$inbreeding, ignore_attr = TRUE)
})
test_that("v2 requires no sommer installation for BLUPF90 preparation", {
  p <- data.frame(ID = letters[1:5], Sire = NA, Dam = NA)
  d <- data.frame(ID = p$ID, w = 1:5)
  expect_s3_class(prepare_blupf90(p, d, "w"), "blupf90_job")
  expect_error(run_genetic_analysis(p, d, "w",
    backend = "blupf90",
    blupf90_control = list(executable = "/does/not/exist")
  ), "executable")
})
