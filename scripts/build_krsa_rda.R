#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# ==============================================================================
# Configuration
# ==============================================================================

input_file <- "results/mappings/krsa_family_mapping.tsv"
output_dir <- "results/krsa_rda"

# Preserve the existing KRSA naming convention.
coverage_object_name <- "KRSA_coverage_STK_PamChip_87202_v1"
mapping_object_name  <- "KRSA_Mapping_STK_PamChip_87202_v1"

coverage_file <- file.path(
  output_dir,
  paste0(coverage_object_name, ".rda")
)

mapping_file <- file.path(
  output_dir,
  paste0(mapping_object_name, ".rda")
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================================
# Load New Mapping
# ==============================================================================

x <- read_tsv(
  input_file,
  show_col_types = FALSE
)

required <- c(
  "Peptide",
  "Kinase"
)

missing <- setdiff(
  required,
  names(x)
)

if (length(missing) > 0L) {
  stop(
    "Input mapping is missing required column(s): ",
    paste(missing, collapse = ", ")
  )
}

x <- x %>%
  transmute(
    Peptide = trimws(as.character(Peptide)),
    Kinase  = trimws(as.character(Kinase))
  ) %>%
  filter(
    !is.na(Peptide),
    nzchar(Peptide),
    !is.na(Kinase),
    nzchar(Kinase)
  ) %>%
  distinct(
    Peptide,
    Kinase
  ) %>%
  arrange(
    Peptide,
    Kinase
  )

cat(
  "Loaded",
  nrow(x),
  "unique peptide-family relationships.\n"
)

cat(
  "Unique peptides:",
  n_distinct(x$Peptide),
  "\n"
)

cat(
  "Unique kinase families:",
  n_distinct(x$Kinase),
  "\n\n"
)

# ==============================================================================
# 1. KRSA Coverage Object
# ==============================================================================
#
# Existing structure:
#
#     Kin   Substrates
# 1 BARK1 ACM1_421_433
# 2 BARK1 ACM1_444_456
# ...
#
# One row per kinase-family / substrate relationship.
# ==============================================================================

KRSA_coverage_STK_PamChip_87102_v2 <- x %>%
  transmute(
    Kin = Kinase,
    Substrates = Peptide
  ) %>%
  distinct() %>%
  arrange(
    Kin,
    Substrates
  ) %>%
  as.data.frame()

# ==============================================================================
# 2. KRSA Mapping Object
# ==============================================================================
#
# Existing structure:
#
# # A tibble:
#   Substrates       Kinases
#   <chr>            <chr>
#   ACM1_421_433     BARK1 PKCA ...
#
# One row per substrate, with kinase families collapsed into a
# space-delimited string.
# ==============================================================================

KRSA_Mapping_STK_PamChip_87102_v1 <- x %>%
  group_by(Peptide) %>%
  summarise(
    Kinases = paste(
      sort(unique(Kinase)),
      collapse = " "
    ),
    .groups = "drop"
  ) %>%
  transmute(
    Substrates = Peptide,
    Kinases = Kinases
  ) %>%
  arrange(
    Substrates
  )

# ==============================================================================
# Validation
# ==============================================================================

expected_pairs <- nrow(x)

coverage_pairs <- nrow(
  KRSA_coverage_STK_PamChip_87102_v2
)

mapping_pairs <- sum(
  lengths(
    strsplit(
      KRSA_Mapping_STK_PamChip_87102_v1$Kinases,
      " ",
      fixed = TRUE
    )
  )
)

if (coverage_pairs != expected_pairs) {
  stop(
    "Coverage validation failed: expected ",
    expected_pairs,
    " relationships but generated ",
    coverage_pairs,
    "."
  )
}

if (mapping_pairs != expected_pairs) {
  stop(
    "Mapping validation failed: expected ",
    expected_pairs,
    " relationships but recovered ",
    mapping_pairs,
    " from collapsed mapping."
  )
}

if (anyDuplicated(
  KRSA_coverage_STK_PamChip_87102_v2[
    c("Kin", "Substrates")
  ]
)) {
  stop(
    "Duplicate Kin/Substrates relationships remain in coverage object."
  )
}

if (anyDuplicated(
  KRSA_Mapping_STK_PamChip_87102_v1$Substrates
)) {
  stop(
    "Duplicate Substrates remain in mapping object."
  )
}

# ==============================================================================
# Save .rda Files
# ==============================================================================

save(
  KRSA_coverage_STK_PamChip_87102_v2,
  file = coverage_file,
  compress = "xz"
)

save(
  KRSA_Mapping_STK_PamChip_87102_v1,
  file = mapping_file,
  compress = "xz"
)

# ==============================================================================
# Report
# ==============================================================================

cat("Created:\n")
cat("  ", coverage_file, "\n", sep = "")
cat("  ", mapping_file, "\n\n", sep = "")

cat("Coverage preview:\n")
print(
  head(
    KRSA_coverage_STK_PamChip_87102_v2
  )
)

cat("\nMapping preview:\n")
print(
  head(
    KRSA_Mapping_STK_PamChip_87102_v1
  )
)

cat("\nValidation passed.\n")
cat(
  "Total peptide-family relationships:",
  expected_pairs,
  "\n"
)