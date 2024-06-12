test_that("PTM conversion Sklyine works", {
  expect_equal(ptm_beautify('K[+72.021129]SAPSTGGVK[+98.1]K[+56]PHR',
                            lookup =shorthistptm_mass,
                            software = 'Skyline',
                            residue = 'remove')
  ,


  "KlaSAPSTGGVKme2KprPHR")
})

test_that("PTM conversion from Sklyine works", {
  expect_equal(ptm_beautify('K[+72.021129]SAPSTGGVK[+98.1]K[+56]PHR',
                            lookup =shorthistptm_mass,
                            software = 'Skyline',
                            residue = 'keep')
               ,


               "la-me2-pr")
})



test_that("PTM conversion from Proline works", {
  expect_equal(ptm_beautify('Propionyl (Any N-term); Crotonyl (K1); Propionyl (K6)',
                            lookup =histptm_lookup,
                            software = 'Proline',
                            residue = 'keep')
               ,


               "prNt-K1cr-K6pr")
})


test_that("PTM conversion from Proline works", {
  expect_equal(ptm_beautify('Propionyl (Any N-term); Crotonyl (K1); Propionyl (K6)',
                            lookup = histptm_lookup,
                            software = 'Proline',
                            residue = 'remove')
               ,


               "prNt-cr-pr")
})


test_that("PTM conversion from Proline works", {
  expect_equal(ptm_beautify('K[+112.1]SAPSIGGVK[+28]K[+56]PHR',
                            lookup = shorthistptm_mass,
                            software = 'Skyline',
                            residue = 'remove')
               ,


               "prNt-pr-me2-pr")
})



