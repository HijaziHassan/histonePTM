# tests/testthat/test-misc_shiftNumbers.R

test_that("misc_shiftNumbers shifts K, S, T positions correctly", {
  expect_equal(
    misc_shiftNumbers("tmaNt-K19ac-K24tma", by = -1),
    "tmaNt-K18ac-K23tma"
  )
  expect_equal(
    misc_shiftNumbers("S10ph-S15ph", by = -1),
    "S9ph-S14ph"
  )
  expect_equal(
    misc_shiftNumbers("T3ph-T6ph", by = -1),
    "T2ph-T5ph"
  )
  expect_equal(
    misc_shiftNumbers("K27me3-S10ph-T3ph", by = -1),
    "K26me3-S9ph-T2ph"
  )
  expect_equal(
    misc_shiftNumbers("Phospho (S10); Acetyl (K27); Phospho (T3)", by = -1),
    "Phospho (S9); Acetyl (K26); Phospho (T2)"
  )
})

test_that("misc_shiftNumbers handles multiple matches in one string", {
  expect_equal(
    misc_shiftNumbers("K1-K2-K3", by = 1),
    "K2-K3-K4"
  )
  expect_equal(
    misc_shiftNumbers("S1-S2-S3", by = 1),
    "S2-S3-S4"
  )
})

test_that("misc_shiftNumbers handles labels with no matches", {
  expect_equal(
    misc_shiftNumbers("No numbers here"),
    "No numbers here"
  )
  expect_equal(
    misc_shiftNumbers(""),
    ""
  )
})

test_that("misc_shiftNumbers respects custom patterns", {
  expect_equal(
    misc_shiftNumbers("R2me1-R5me2", pattern = "(?<=R)[0-9]+", by = -1),
    "R1me1-R4me2"
  )
})
