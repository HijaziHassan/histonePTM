test_that("multiplication works", {

  renamed_ptm <- ptm_rename(PTM = "TMAyl_correct (Any N-term); Butyryl (K28); Trimethyl (K37); Propionyl (K38)",
             lookup = histptm_lookup,
             software = "Proline")


  expect_equal(renamed_ptm, "tmaNt-K28bu-K37me3-K38pr")
})
