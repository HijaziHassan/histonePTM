# tests/testthat/test-misc_wrapLabels.R

test_that("misc_wrapLabels handles semicolon lists correctly", {
  # Simple semicolon list
  expect_equal(
    misc_wrapLabels("A; B; C"),
    "A\nB\nC"
  )

  # Long semicolon segments that exceed max_width
  long_label <- "VeryLongModificationName (K1); AnotherLong (K2); Short (K3)"
  result <- misc_wrapLabels(long_label, max_width = 20)
  # Expected: first line "VeryLongModificationName (K1);", second "AnotherLong (K2);", third "Short (K3)"
  expect_true(stringr::str_count(result, "\n") == 2)
  # Check each line length <=20 (approx)
  lines <- strsplit(result, "\n")[[1]]
  expect_true(all(nchar(lines) <= 20))
})

test_that("misc_wrapLabels handles hyphen chains correctly", {
  # Fits within default width
  expect_equal(
    misc_wrapLabels("tmaNt-K19ac-K24tma"),
    "tmaNt-K19ac-K24tma"
  )

  # Breaks with small max_width
  expect_equal(
    misc_wrapLabels("tmaNt-K19ac-K24tma", max_width = 10),
    "tmaNt-\nK19ac-\nK24tma"
  )

  # Complex chain with many parts
  long_chain <- "K27me3-K36me2-K79me1-K80me2"
  result <- misc_wrapLabels(long_chain, max_width = 12)
  expected <- "K27me3-\nK36me2-\nK79me1-\nK80me2"
  expect_equal(result, expected)
})

test_that("misc_wrapLabels handles plain text (space-separated)", {
  # Short text stays as is
  expect_equal(
    misc_wrapLabels("short text"),
    "short text"
  )

  # Long text gets wrapped on spaces
  long_text <- "This is a very long label with many words"
  result <- misc_wrapLabels(long_text, max_width = 15)
  # Check that lines are not longer than 15 characters
  lines <- strsplit(result, "\n")[[1]]
  expect_true(all(nchar(lines) <= 15))
  # Ensure no hyphen was introduced
  expect_false(grepl("-", result))
})

test_that("misc_wrapLabels handles mixed input vector", {
  labels <- c(
    "A; B; C",
    "K27me3-K36me2-K79me1",
    "plain long text for wrapping"
  )
  result <- misc_wrapLabels(labels, max_width = 10)
  expect_length(result, 3)
  # Check that each result has newlines where appropriate
  expect_true(grepl("\n", result[1]))
  expect_true(grepl("\n", result[2]))
  expect_true(grepl("\n", result[3]))
})

test_that("misc_wrapLabels handles edge cases", {
  # Empty vector
  expect_equal(misc_wrapLabels(character(0)), character(0))

  # NA values
  expect_equal(misc_wrapLabels(NA_character_), NA_character_)

  # Empty string
  expect_equal(misc_wrapLabels(""), "")

  # String with only spaces
  expect_equal(misc_wrapLabels("   "), "")

  # Unicode minus (U+2212) recognized as hyphen
  result <- misc_wrapLabels("K27me3−K36me2", max_width = 10)
  expect_equal(result, "K27me3−\nK36me2")
})

test_that("misc_wrapLabels prioritizes semicolon over hyphen", {
  # Mixed: semicolon present, should use semicolon logic, not hyphen
  mixed <- "A; B-C"
  result <- misc_wrapLabels(mixed)
  expect_true(grepl("\n", result))  # semicolon triggers newline
  expect_false(grepl("-\n", result))  # hyphen not used as line break
})

test_that("misc_wrapLabels respects different max_width values", {
  label <- "K27me3-K36me2-K79me1"
  # Width 20 -> one line
  expect_equal(misc_wrapLabels(label, max_width = 20), label)
  # Width 12 -> breaks at each hyphen
  expect_true(grepl("\n", misc_wrapLabels(label, max_width = 12)))
  # Check that line lengths <= max_width
  lines <- strsplit(misc_wrapLabels(label, max_width = 12), "\n")[[1]]
  expect_true(all(nchar(lines) <= 12))
})
