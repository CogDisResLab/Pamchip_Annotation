#' Calculate Mapping Coverage Summary
#'
#' Summarize peptide- and family-level coverage across the integrated
#' kinase-evidence pipeline.
#'
#' @param chip_sites Validated chip-site annotations.
#' @param master_evidence Harmonized source-native evidence table after
#'   peptide/site linkage.
#' @param family_evidence Family-level collapsed evidence table.
#'
#' @return A named list containing overall coverage metrics, source-specific
#'   coverage, concordance metrics, and family-level counts.
calculate_mapping_coverage <- function(
  chip_sites,
  master_evidence,
  family_evidence
) {

  if (!is.data.frame(chip_sites)) {
    stop(
      "`chip_sites` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  if (!is.data.frame(master_evidence)) {
    stop(
      "`master_evidence` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  if (!is.data.frame(family_evidence)) {
    stop(
      "`family_evidence` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  required_chip <- c(
    "peptide_id"
  )

  required_master <- c(
    "peptide_id",
    "evidence_source"
  )

  required_family <- c(
    "peptide_id",
    "kinase_family",
    "signor_supported",
    "kinase_library_supported",
    "gps6_supported",
    "predictive_concordant",
    "concordant"
  )

  missing_chip <- setdiff(
    required_chip,
    names(chip_sites)
  )

  missing_master <- setdiff(
    required_master,
    names(master_evidence)
  )

  missing_family <- setdiff(
    required_family,
    names(family_evidence)
  )

  if (length(missing_chip) > 0L) {
    stop(
      "`chip_sites` is missing required column(s): ",
      paste(missing_chip, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_master) > 0L) {
    stop(
      "`master_evidence` is missing required column(s): ",
      paste(missing_master, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_family) > 0L) {
    stop(
      "`family_evidence` is missing required column(s): ",
      paste(missing_family, collapse = ", "),
      call. = FALSE
    )
  }

  valid_peptide <- function(x) {
    !is.na(x) &
      nzchar(
        trimws(
          as.character(x)
        )
      )
  }

  unique_chip_peptides <- unique(
    as.character(
      chip_sites$peptide_id[
        valid_peptide(
          chip_sites$peptide_id
        )
      ]
    )
  )

  total_chip_peptides <- length(
    unique_chip_peptides
  )

  mapped_family_peptides <- unique(
    as.character(
      family_evidence$peptide_id[
        valid_peptide(
          family_evidence$peptide_id
        ) &
          !is.na(
            family_evidence$kinase_family
          ) &
          nzchar(
            trimws(
              as.character(
                family_evidence$kinase_family
              )
            )
          )
      ]
    )
  )

  n_mapped_family_peptides <- length(
    mapped_family_peptides
  )

  n_unmapped_peptides <- total_chip_peptides -
    n_mapped_family_peptides

  pct_mapped <- if (total_chip_peptides > 0L) {
    100 *
      n_mapped_family_peptides /
      total_chip_peptides
  } else {
    NA_real_
  }

  pct_unmapped <- if (total_chip_peptides > 0L) {
    100 *
      n_unmapped_peptides /
      total_chip_peptides
  } else {
    NA_real_
  }

  source_peptide_count <- function(flag_column) {

    x <- family_evidence

    supported <- !is.na(
      x[[flag_column]]
    ) &
      x[[flag_column]] %in% TRUE

    length(
      unique(
        as.character(
          x$peptide_id[
            supported &
              valid_peptide(
                x$peptide_id
              )
          ]
        )
      )
    )
  }

  source_pair_count <- function(flag_column) {

    supported <- !is.na(
      family_evidence[[flag_column]]
    ) &
      family_evidence[[flag_column]] %in% TRUE

    sum(
      supported,
      na.rm = TRUE
    )
  }

  n_signor_peptides <- source_peptide_count(
    "signor_supported"
  )

  n_kl_peptides <- source_peptide_count(
    "kinase_library_supported"
  )

  n_gps6_peptides <- source_peptide_count(
    "gps6_supported"
  )

  n_predictive_concordant_peptides <- source_peptide_count(
    "predictive_concordant"
  )

  n_concordant_peptides <- source_peptide_count(
    "concordant"
  )

  source_summary <- tibble::tibble(
    evidence_source = c(
      "PhosphoSIGNOR",
      "Kinase Library",
      "GPS6",
      "Predictive concordant",
      "Multi-source concordant"
    ),

    peptide_count = c(
      n_signor_peptides,
      n_kl_peptides,
      n_gps6_peptides,
      n_predictive_concordant_peptides,
      n_concordant_peptides
    ),

    peptide_percent = if (total_chip_peptides > 0L) {
      100 *
        c(
          n_signor_peptides,
          n_kl_peptides,
          n_gps6_peptides,
          n_predictive_concordant_peptides,
          n_concordant_peptides
        ) /
        total_chip_peptides
    } else {
      rep(
        NA_real_,
        5L
      )
    },

    peptide_family_pairs = c(
      source_pair_count(
        "signor_supported"
      ),
      source_pair_count(
        "kinase_library_supported"
      ),
      source_pair_count(
        "gps6_supported"
      ),
      source_pair_count(
        "predictive_concordant"
      ),
      source_pair_count(
        "concordant"
      )
    )
  )

  family_counts <- family_evidence |>
    dplyr::filter(
      !is.na(.data$kinase_family),
      nzchar(
        trimws(
          as.character(
            .data$kinase_family
          )
        )
      )
    ) |>
    dplyr::count(
      .data$kinase_family,
      name = "peptide_family_rows",
      sort = TRUE
    )

  family_source_counts <- family_evidence |>
    dplyr::filter(
      !is.na(.data$kinase_family),
      nzchar(
        trimws(
          as.character(
            .data$kinase_family
          )
        )
      )
    ) |>
    dplyr::group_by(
      .data$kinase_family
    ) |>
    dplyr::summarise(
      signor_rows =
        sum(
          .data$signor_supported %in% TRUE,
          na.rm = TRUE
        ),

      kinase_library_rows =
        sum(
          .data$kinase_library_supported %in% TRUE,
          na.rm = TRUE
        ),

      gps6_rows =
        sum(
          .data$gps6_supported %in% TRUE,
          na.rm = TRUE
        ),

      predictive_concordant_rows =
        sum(
          .data$predictive_concordant %in% TRUE,
          na.rm = TRUE
        ),

      concordant_rows =
        sum(
          .data$concordant %in% TRUE,
          na.rm = TRUE
        ),

      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$concordant_rows
      ),
      dplyr::desc(
        .data$gps6_rows
      ),
      .data$kinase_family
    )

  linked_master <- master_evidence |>
    dplyr::filter(
      valid_peptide(
        .data$peptide_id
      )
    )

  unlinked_master <- master_evidence |>
    dplyr::filter(
      !valid_peptide(
        .data$peptide_id
      )
    )

  master_linkage <- tibble::tibble(
    metric = c(
      "Total master evidence rows",
      "Linked master evidence rows",
      "Unlinked master evidence rows"
    ),

    value = c(
      nrow(master_evidence),
      nrow(linked_master),
      nrow(unlinked_master)
    )
  )

  overall <- tibble::tibble(
    metric = c(
      "Total chip peptides",
      "Peptides with >=1 kinase-family mapping",
      "Peptides without kinase-family mapping",
      "Percent chip peptides mapped",
      "Percent chip peptides unmapped",
      "Total family-evidence rows",
      "Unique kinase families"
    ),

    value = c(
      total_chip_peptides,
      n_mapped_family_peptides,
      n_unmapped_peptides,
      pct_mapped,
      pct_unmapped,
      nrow(family_evidence),
      dplyr::n_distinct(
        family_evidence$kinase_family,
        na.rm = TRUE
      )
    )
  )

  list(
    overall = overall,
    source_coverage = source_summary,
    family_counts = family_counts,
    family_source_counts = family_source_counts,
    master_linkage = master_linkage
  )
}

#' Perform Mapping Quality Control
#'
#' Run structural and consistency checks across the integrated kinase mapping
#' pipeline. The function is designed for both pipeline gating and report
#' generation: hard schema failures stop immediately, whereas expected biological
#' exclusions (for example, curated SIGNOR phosphosites that are not present on
#' the current chip) are summarized rather than treated as errors.
#'
#' @param chip_sites Validated chip-site annotations.
#' @param master_evidence Harmonized evidence table after peptide/site linkage.
#' @param family_evidence Family-level collapsed evidence table.
#' @param family_mapping Primary family-level mapping table.
#' @param experimental_mapping Experimental/PhosphoSIGNOR mapping table.
#' @param kinase_library_mapping Kinase Library mapping table.
#' @param gps6_mapping GPS6 mapping table.
#' @param predictive_concordant_mapping GPS6 + Kinase Library concordant mapping.
#' @param concordant_mapping Multi-source concordant mapping.
#'
#' @return Named list containing QC status, summary metrics, duplicate checks,
#'   source consistency checks, linkage diagnostics, and issue tables.
perform_mapping_qc <- function(
  chip_sites,
  master_evidence,
  family_evidence,
  family_mapping,
  experimental_mapping = NULL,
  kinase_library_mapping = NULL,
  gps6_mapping = NULL,
  predictive_concordant_mapping = NULL,
  concordant_mapping = NULL
) {

  # --------------------------------------------------------------------------
  # Validation helpers
  # --------------------------------------------------------------------------

  require_df <- function(x, name) {
    if (!is.data.frame(x)) {
      stop(
        "`", name, "` must be a data.frame or tibble.",
        call. = FALSE
      )
    }
  }

  require_columns <- function(x, columns, name) {

    missing <- setdiff(
      columns,
      names(x)
    )

    if (length(missing) > 0L) {
      stop(
        "`",
        name,
        "` is missing required column(s): ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }

  valid_chr <- function(x) {
    !is.na(x) &
      nzchar(
        trimws(
          as.character(x)
        )
      )
  }

  count_unique_pairs <- function(x) {

    if (
      nrow(x) == 0L ||
      !"peptide_id" %in% names(x) ||
      !"kinase_family" %in% names(x)
    ) {
      return(0L)
    }

    x |>
      dplyr::filter(
        valid_chr(.data$peptide_id),
        valid_chr(.data$kinase_family)
      ) |>
      dplyr::distinct(
        .data$peptide_id,
        .data$kinase_family
      ) |>
      nrow()
  }

  duplicate_pairs <- function(x, object_name) {

    if (
      nrow(x) == 0L ||
      !"peptide_id" %in% names(x) ||
      !"kinase_family" %in% names(x)
    ) {
      return(
        tibble::tibble(
          object = character(),
          peptide_id = character(),
          kinase_family = character(),
          n = integer()
        )
      )
    }

    x |>
      dplyr::filter(
        valid_chr(.data$peptide_id),
        valid_chr(.data$kinase_family)
      ) |>
      dplyr::count(
        .data$peptide_id,
        .data$kinase_family,
        name = "n"
      ) |>
      dplyr::filter(
        .data$n > 1L
      ) |>
      dplyr::mutate(
        object = object_name,
        .before = 1L
      )
  }

  pair_key <- function(x) {

    if (
      nrow(x) == 0L ||
      !"peptide_id" %in% names(x) ||
      !"kinase_family" %in% names(x)
    ) {
      return(character())
    }

    paste(
      as.character(x$peptide_id),
      as.character(x$kinase_family),
      sep = "\r"
    )
  }

  require_df(chip_sites, "chip_sites")
  require_df(master_evidence, "master_evidence")
  require_df(family_evidence, "family_evidence")
  require_df(family_mapping, "family_mapping")

  optional_mappings <- list(
    experimental_mapping = experimental_mapping,
    kinase_library_mapping = kinase_library_mapping,
    gps6_mapping = gps6_mapping,
    predictive_concordant_mapping = predictive_concordant_mapping,
    concordant_mapping = concordant_mapping
  )

  supplied_optional_mappings <- optional_mappings[
    !vapply(
      optional_mappings,
      is.null,
      logical(1)
    )
  ]

  for (nm in names(supplied_optional_mappings)) {
    require_df(
      supplied_optional_mappings[[nm]],
      nm
    )
  }

  require_columns(
    chip_sites,
    c(
      "peptide_id"
    ),
    "chip_sites"
  )

  require_columns(
    master_evidence,
    c(
      "peptide_id",
      "kinase_family",
      "evidence_source"
    ),
    "master_evidence"
  )

  require_columns(
    family_evidence,
    c(
      "peptide_id",
      "kinase_family",
      "signor_supported",
      "kinase_library_supported",
      "gps6_supported",
      "predictive_concordant",
      "concordant"
    ),
    "family_evidence"
  )

  require_columns(
    family_mapping,
    c(
      "peptide_id",
      "kinase_family"
    ),
    "family_mapping"
  )

  for (nm in names(supplied_optional_mappings)) {
    require_columns(
      supplied_optional_mappings[[nm]],
      c(
        "peptide_id",
        "kinase_family"
      ),
      nm
    )
  }

  # --------------------------------------------------------------------------
  # Basic dimensions and coverage
  # --------------------------------------------------------------------------

  chip_peptides <- unique(
    as.character(
      chip_sites$peptide_id[
        valid_chr(
          chip_sites$peptide_id
        )
      ]
    )
  )

  family_peptides <- unique(
    as.character(
      family_evidence$peptide_id[
        valid_chr(
          family_evidence$peptide_id
        )
      ]
    )
  )

  mapped_peptides <- unique(
    as.character(
      family_mapping$peptide_id[
        valid_chr(
          family_mapping$peptide_id
        )
      ]
    )
  )

  mapping_objects <- c(
    list(
      chip_sites = chip_sites,
      master_evidence = master_evidence,
      family_evidence = family_evidence,
      family_mapping = family_mapping
    ),
    supplied_optional_mappings
  )

  dimensions <- dplyr::bind_rows(
    lapply(
      names(mapping_objects),
      function(nm) {

        obj <- mapping_objects[[nm]]

        tibble::tibble(
          object = nm,

          rows = nrow(obj),

          unique_peptides = if ("peptide_id" %in% names(obj)) {
            dplyr::n_distinct(
              obj$peptide_id[
                valid_chr(obj$peptide_id)
              ]
            )
          } else {
            NA_integer_
          },

          unique_peptide_family_pairs = if (
            nm == "chip_sites"
          ) {
            NA_integer_
          } else {
            count_unique_pairs(obj)
          }
        )
      }
    )
  )


  # --------------------------------------------------------------------------
  # Duplicate peptide × family rows
  # --------------------------------------------------------------------------

  duplicate_objects <- c(
    list(
      family_evidence = family_evidence,
      family_mapping = family_mapping
    ),
    supplied_optional_mappings
  )

  duplicate_summary <- dplyr::bind_rows(
    lapply(
      names(duplicate_objects),
      function(nm) {
        duplicate_pairs(
          duplicate_objects[[nm]],
          nm
        )
      }
    )
  )


  # --------------------------------------------------------------------------
  # Master linkage
  # --------------------------------------------------------------------------

  master_linked <- valid_chr(
    master_evidence$peptide_id
  )

  master_linkage <- master_evidence |>
    dplyr::mutate(
      linked_to_chip =
        master_linked
    ) |>
    dplyr::group_by(
      .data$evidence_source
    ) |>
    dplyr::summarise(
      total_rows =
        dplyr::n(),

      linked_rows =
        sum(
          .data$linked_to_chip,
          na.rm = TRUE
        ),

      unlinked_rows =
        sum(
          !.data$linked_to_chip,
          na.rm = TRUE
        ),

      linked_percent =
        100 *
          .data$linked_rows /
          .data$total_rows,

      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$total_rows
      )
    )

  # --------------------------------------------------------------------------
  # Family assignment missingness
  # --------------------------------------------------------------------------

  family_missingness <- master_evidence |>
    dplyr::mutate(
      missing_family =
        !valid_chr(
          .data$kinase_family
        )
    ) |>
    dplyr::group_by(
      .data$evidence_source
    ) |>
    dplyr::summarise(
      total_rows =
        dplyr::n(),

      missing_family_rows =
        sum(
          .data$missing_family,
          na.rm = TRUE
        ),

      assigned_family_rows =
        .data$total_rows -
          .data$missing_family_rows,

      assigned_percent =
        100 *
          .data$assigned_family_rows /
          .data$total_rows,

      .groups = "drop"
    )

  # --------------------------------------------------------------------------
  # Evidence-source consistency
  # --------------------------------------------------------------------------

  signor_pairs <- family_evidence |>
    dplyr::filter(
      .data$signor_supported %in% TRUE
    )

  kl_pairs <- family_evidence |>
    dplyr::filter(
      .data$kinase_library_supported %in% TRUE
    )

  gps6_pairs <- family_evidence |>
    dplyr::filter(
      .data$gps6_supported %in% TRUE
    )

  predictive_pairs <- family_evidence |>
    dplyr::filter(
      .data$predictive_concordant %in% TRUE
    )

  concordant_pairs <- family_evidence |>
    dplyr::filter(
      .data$concordant %in% TRUE
    )

  expected_keys_all <- list(
    experimental_mapping =
      unique(
        pair_key(signor_pairs)
      ),

    kinase_library_mapping =
      unique(
        pair_key(kl_pairs)
      ),

    gps6_mapping =
      unique(
        pair_key(gps6_pairs)
      ),

    predictive_concordant_mapping =
      unique(
        pair_key(predictive_pairs)
      ),

    concordant_mapping =
      unique(
        pair_key(concordant_pairs)
      )
  )

  if (length(supplied_optional_mappings) > 0L) {

    source_consistency <- dplyr::bind_rows(
      lapply(
        names(supplied_optional_mappings),
        function(nm) {

          actual <- unique(
            pair_key(
              supplied_optional_mappings[[nm]]
            )
          )

          expected <- expected_keys_all[[nm]]

          tibble::tibble(
            object =
              nm,

            expected_pairs =
              length(expected),

            actual_pairs =
              length(actual),

            missing_expected_pairs =
              length(
                setdiff(
                  expected,
                  actual
                )
              ),

            unexpected_pairs =
              length(
                setdiff(
                  actual,
                  expected
                )
              )
          )
        }
      )
    )

  } else {

    source_consistency <- tibble::tibble(
      object = character(),
      expected_pairs = integer(),
      actual_pairs = integer(),
      missing_expected_pairs = integer(),
      unexpected_pairs = integer()
    )
  }


  # --------------------------------------------------------------------------
  # Logical concordance invariants
  # --------------------------------------------------------------------------

  predictive_flag_invalid <- family_evidence |>
    dplyr::filter(
      .data$predictive_concordant %in% TRUE &
        !(
          .data$kinase_library_supported %in% TRUE &
            .data$gps6_supported %in% TRUE
        )
    )

  concordant_flag_invalid <- family_evidence |>
    dplyr::filter(
      .data$concordant %in% TRUE &
        (
          (
            .data$signor_supported %in% TRUE
          ) +
            (
              .data$kinase_library_supported %in% TRUE
            ) +
            (
              .data$gps6_supported %in% TRUE
            )
        ) < 2L
    )

  flag_checks <- tibble::tibble(
    check = c(
      "predictive_concordant requires KL + GPS6",
      "concordant requires >=2 evidence sources"
    ),
    failing_rows = c(
      nrow(
        predictive_flag_invalid
      ),
      nrow(
        concordant_flag_invalid
      )
    )
  )

  # --------------------------------------------------------------------------
  # Mapping coverage relative to chip
  # --------------------------------------------------------------------------

  mapping_coverage <- tibble::tibble(
    metric = c(
      "Total chip peptides",
      "Peptides represented in family evidence",
      "Peptides represented in primary family mapping",
      "Chip peptides absent from family evidence",
      "Chip peptides absent from primary family mapping"
    ),
    value = c(
      length(chip_peptides),
      length(family_peptides),
      length(mapped_peptides),
      length(
        setdiff(
          chip_peptides,
          family_peptides
        )
      ),
      length(
        setdiff(
          chip_peptides,
          mapped_peptides
        )
      )
    )
  )

  unmapped_chip_peptides <- tibble::tibble(
    peptide_id =
      setdiff(
        chip_peptides,
        mapped_peptides
      )
  )

  # --------------------------------------------------------------------------
  # Overall QC status
  # --------------------------------------------------------------------------

  n_duplicate_mapping_rows <- nrow(
    duplicate_summary |>
      dplyr::filter(
        .data$object != "family_evidence"
      )
  )

  n_source_consistency_failures <- if (nrow(source_consistency) > 0L) {
    sum(
      source_consistency$missing_expected_pairs +
        source_consistency$unexpected_pairs
    )
  } else {
    0L
  }

  n_flag_failures <- sum(
    flag_checks$failing_rows
  )

  # Family evidence can legitimately contain multiple rows per peptide/family
  # if lower-level provenance has not yet been fully collapsed. The exported
  # mapping variants, however, should each be unique at peptide × family level.
  qc_pass <- (
    n_duplicate_mapping_rows == 0L &&
      n_source_consistency_failures == 0L &&
      n_flag_failures == 0L
  )

  qc_status <- tibble::tibble(
    qc_pass =
      qc_pass,

    duplicate_mapping_pair_rows =
      n_duplicate_mapping_rows,

    source_consistency_failures =
      n_source_consistency_failures,

    logical_flag_failures =
      n_flag_failures,

    unlinked_master_rows =
      sum(
        !master_linked
      ),

    unmapped_chip_peptides =
      nrow(
        unmapped_chip_peptides
      ),

    optional_mapping_tables_checked =
      length(
        supplied_optional_mappings
      )
  )

  if (!qc_pass) {
    warning(
      "Mapping QC completed with one or more structural consistency failures. ",
      "Inspect `qc_summary$qc_status`, `qc_summary$duplicates`, ",
      "`qc_summary$source_consistency`, and `qc_summary$flag_checks`.",
      call. = FALSE
    )
  }

  list(
    qc_status =
      qc_status,

    dimensions =
      dimensions,

    mapping_coverage =
      mapping_coverage,

    master_linkage =
      master_linkage,

    family_missingness =
      family_missingness,

    source_consistency =
      source_consistency,

    flag_checks =
      flag_checks,

    duplicates =
      duplicate_summary,

    unmapped_chip_peptides =
      unmapped_chip_peptides
  )
}

