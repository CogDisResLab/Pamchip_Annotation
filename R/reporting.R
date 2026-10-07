# R/reporting.R

#' Write Mapping Outputs
#'
#' Export the primary KRSA-compatible peptide -> kinase-family mapping together
#' with provenance-rich mapping tables, sensitivity mappings, evidence tables,
#' QC summaries, and coverage summaries.
#'
#' The primary deliverable is written in two forms:
#'
#'   1. krsa_family_mapping.tsv
#'      Minimal KRSA-ready two-column table:
#'          Peptide    Kinase
#'
#'   2. krsa_family_mapping_full.tsv
#'      Full primary mapping with all available evidence/provenance columns.
#'
#' @param chip_sites Validated chip-site annotations.
#' @param master_evidence Harmonized source-native evidence table after linkage.
#' @param family_evidence Family-level collapsed evidence table.
#' @param family_mapping Primary KRSA family mapping.
#' @param experimental_mapping PhosphoSIGNOR-only mapping.
#' @param kinase_library_mapping Kinase Library-only mapping.
#' @param gps6_mapping GPS6-only mapping.
#' @param predictive_concordant_mapping GPS6 + Kinase Library concordant mapping.
#' @param concordant_mapping General multi-source concordant mapping.
#' @param qc_summary Named QC list from perform_mapping_qc().
#' @param coverage_summary Named coverage list from calculate_mapping_coverage().
#' @param manifest Pipeline manifest.
#'
#' @return Character vector of all files written. Suitable for a
#'   targets target with format = "file".
write_mapping_outputs <- function(
  chip_sites,
  master_evidence,
  family_evidence,
  family_mapping,
  experimental_mapping,
  kinase_library_mapping,
  gps6_mapping,
  predictive_concordant_mapping,
  concordant_mapping,
  qc_summary,
  coverage_summary,
  manifest
) {

  # ============================================================================
  # Validation
  # ============================================================================

  require_df <- function(x, name) {
    if (!is.data.frame(x)) {
      stop(
        "`", name, "` must be a data.frame or tibble.",
        call. = FALSE
      )
    }
  }

  require_mapping_columns <- function(x, name) {

    require_df(
      x,
      name
    )

    required <- c(
      "peptide_id",
      "kinase_family"
    )

    missing <- setdiff(
      required,
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

  require_df(
    chip_sites,
    "chip_sites"
  )

  require_df(
    master_evidence,
    "master_evidence"
  )

  require_df(
    family_evidence,
    "family_evidence"
  )

  require_mapping_columns(
    family_mapping,
    "family_mapping"
  )

  require_mapping_columns(
    experimental_mapping,
    "experimental_mapping"
  )

  require_mapping_columns(
    kinase_library_mapping,
    "kinase_library_mapping"
  )

  require_mapping_columns(
    gps6_mapping,
    "gps6_mapping"
  )

  require_mapping_columns(
    predictive_concordant_mapping,
    "predictive_concordant_mapping"
  )

  require_mapping_columns(
    concordant_mapping,
    "concordant_mapping"
  )

  if (!is.list(qc_summary)) {
    stop(
      "`qc_summary` must be the named list returned by perform_mapping_qc().",
      call. = FALSE
    )
  }

  if (!is.list(coverage_summary)) {
    stop(
      "`coverage_summary` must be the named list returned by calculate_mapping_coverage().",
      call. = FALSE
    )
  }

  if (is.null(manifest)) {
    manifest <- list()
  }

  # ============================================================================
  # Output directories
  # ============================================================================

  output_dir <-
    manifest$outputs$directory %||%
    manifest$output$directory %||%
    manifest$outputs$dir %||%
    manifest$output$dir %||%
    "results"

  mapping_dir <- file.path(
    output_dir,
    "mappings"
  )

  evidence_dir <- file.path(
    output_dir,
    "evidence"
  )

  qc_dir <- file.path(
    output_dir,
    "qc"
  )

  dir.create(
    mapping_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    evidence_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    qc_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ============================================================================
  # Helpers
  # ============================================================================

  clean_mapping <- function(x) {

    x |>
      dplyr::filter(
        !is.na(.data$peptide_id),
        nzchar(
          trimws(
            as.character(
              .data$peptide_id
            )
          )
        ),
        !is.na(.data$kinase_family),
        nzchar(
          trimws(
            as.character(
              .data$kinase_family
            )
          )
        )
      ) |>
      dplyr::mutate(
        peptide_id =
          as.character(
            .data$peptide_id
          ),

        kinase_family =
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(
                .data$kinase_family
              )
            )
          )
      ) |>
      dplyr::distinct(
        .data$peptide_id,
        .data$kinase_family,
        .keep_all = TRUE
      ) |>
      dplyr::arrange(
        .data$peptide_id,
        .data$kinase_family
      )
  }

  as_krsa_mapping <- function(x) {

    clean_mapping(x) |>
      dplyr::transmute(
        Peptide =
          .data$peptide_id,

        Kinase =
          .data$kinase_family
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(
        .data$Peptide,
        .data$Kinase
      )
  }

  write_tsv_safe <- function(x, path) {

    readr::write_tsv(
      x,
      path,
      na = ""
    )

    if (
      !file.exists(path) ||
      file.info(path)$size <= 0L
    ) {
      stop(
        "Failed to create output file: ",
        path,
        call. = FALSE
      )
    }

    normalizePath(
      path,
      winslash = "/",
      mustWork = TRUE
    )
  }

  write_csv_safe <- function(x, path) {

    readr::write_csv(
      x,
      path,
      na = ""
    )

    if (
      !file.exists(path) ||
      file.info(path)$size <= 0L
    ) {
      stop(
        "Failed to create output file: ",
        path,
        call. = FALSE
      )
    }

    normalizePath(
      path,
      winslash = "/",
      mustWork = TRUE
    )
  }

  write_list_tables <- function(x, directory, prefix) {

    paths <- character()

    for (nm in names(x)) {

      obj <- x[[nm]]

      if (is.data.frame(obj)) {

        path <- file.path(
          directory,
          paste0(
            prefix,
            "_",
            nm,
            ".csv"
          )
        )

        paths <- c(
          paths,
          write_csv_safe(
            obj,
            path
          )
        )
      }
    }

    paths
  }

  # ============================================================================
  # Normalize mapping variants
  # ============================================================================

  family_mapping_clean <-
    clean_mapping(
      family_mapping
    )

  experimental_mapping_clean <-
    clean_mapping(
      experimental_mapping
    )

  kinase_library_mapping_clean <-
    clean_mapping(
      kinase_library_mapping
    )

  gps6_mapping_clean <-
    clean_mapping(
      gps6_mapping
    )

  predictive_concordant_mapping_clean <-
    clean_mapping(
      predictive_concordant_mapping
    )

  concordant_mapping_clean <-
    clean_mapping(
      concordant_mapping
    )

  # ============================================================================
  # Primary deliverables
  # ============================================================================

  primary_krsa_path <- write_tsv_safe(
    as_krsa_mapping(
      family_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_family_mapping.tsv"
    )
  )

  primary_full_path <- write_tsv_safe(
    family_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_family_mapping_full.tsv"
    )
  )

  # Also write a CSV version for convenient manual inspection.
  primary_csv_path <- write_csv_safe(
    family_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_family_mapping_full.csv"
    )
  )

  # ============================================================================
  # Sensitivity / source-specific mappings
  # ============================================================================

  experimental_path <- write_tsv_safe(
    as_krsa_mapping(
      experimental_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_mapping_phosphosignor_only.tsv"
    )
  )

  kinase_library_path <- write_tsv_safe(
    as_krsa_mapping(
      kinase_library_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_mapping_kinase_library_only.tsv"
    )
  )

  gps6_path <- write_tsv_safe(
    as_krsa_mapping(
      gps6_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_mapping_gps6_only.tsv"
    )
  )

  predictive_concordant_path <- write_tsv_safe(
    as_krsa_mapping(
      predictive_concordant_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_mapping_predictive_concordant.tsv"
    )
  )

  concordant_path <- write_tsv_safe(
    as_krsa_mapping(
      concordant_mapping_clean
    ),
    file.path(
      mapping_dir,
      "krsa_mapping_multisource_concordant.tsv"
    )
  )

  # Full provenance-bearing variants are useful for auditing why a given
  # peptide-family pair appears in a sensitivity mapping.
  experimental_full_path <- write_tsv_safe(
    experimental_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_mapping_phosphosignor_only_full.tsv"
    )
  )

  kinase_library_full_path <- write_tsv_safe(
    kinase_library_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_mapping_kinase_library_only_full.tsv"
    )
  )

  gps6_full_path <- write_tsv_safe(
    gps6_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_mapping_gps6_only_full.tsv"
    )
  )

  predictive_concordant_full_path <- write_tsv_safe(
    predictive_concordant_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_mapping_predictive_concordant_full.tsv"
    )
  )

  concordant_full_path <- write_tsv_safe(
    concordant_mapping_clean,
    file.path(
      mapping_dir,
      "krsa_mapping_multisource_concordant_full.tsv"
    )
  )

  # ============================================================================
  # Evidence/provenance tables
  # ============================================================================

  chip_sites_path <- write_tsv_safe(
    chip_sites,
    file.path(
      evidence_dir,
      "validated_chip_sites.tsv"
    )
  )

  master_evidence_path <- write_tsv_safe(
    master_evidence,
    file.path(
      evidence_dir,
      "master_kinase_evidence.tsv"
    )
  )

  family_evidence_path <- write_tsv_safe(
    family_evidence,
    file.path(
      evidence_dir,
      "family_evidence.tsv"
    )
  )

  # ============================================================================
  # QC and coverage tables
  # ============================================================================

  qc_paths <- write_list_tables(
    qc_summary,
    qc_dir,
    "qc"
  )

  coverage_paths <- write_list_tables(
    coverage_summary,
    qc_dir,
    "coverage"
  )

  # ============================================================================
  # Export summary / machine-readable manifest
  # ============================================================================

  export_summary <- tibble::tibble(
    mapping = c(
      "primary",
      "experimental_only",
      "kinase_library_only",
      "gps6_only",
      "predictive_concordant",
      "multisource_concordant"
    ),

    peptide_family_pairs = c(
      nrow(
        as_krsa_mapping(
          family_mapping_clean
        )
      ),
      nrow(
        as_krsa_mapping(
          experimental_mapping_clean
        )
      ),
      nrow(
        as_krsa_mapping(
          kinase_library_mapping_clean
        )
      ),
      nrow(
        as_krsa_mapping(
          gps6_mapping_clean
        )
      ),
      nrow(
        as_krsa_mapping(
          predictive_concordant_mapping_clean
        )
      ),
      nrow(
        as_krsa_mapping(
          concordant_mapping_clean
        )
      )
    ),

    unique_peptides = c(
      dplyr::n_distinct(
        family_mapping_clean$peptide_id
      ),
      dplyr::n_distinct(
        experimental_mapping_clean$peptide_id
      ),
      dplyr::n_distinct(
        kinase_library_mapping_clean$peptide_id
      ),
      dplyr::n_distinct(
        gps6_mapping_clean$peptide_id
      ),
      dplyr::n_distinct(
        predictive_concordant_mapping_clean$peptide_id
      ),
      dplyr::n_distinct(
        concordant_mapping_clean$peptide_id
      )
    ),

    unique_families = c(
      dplyr::n_distinct(
        family_mapping_clean$kinase_family
      ),
      dplyr::n_distinct(
        experimental_mapping_clean$kinase_family
      ),
      dplyr::n_distinct(
        kinase_library_mapping_clean$kinase_family
      ),
      dplyr::n_distinct(
        gps6_mapping_clean$kinase_family
      ),
      dplyr::n_distinct(
        predictive_concordant_mapping_clean$kinase_family
      ),
      dplyr::n_distinct(
        concordant_mapping_clean$kinase_family
      )
    )
  )

  export_summary_path <- write_csv_safe(
    export_summary,
    file.path(
      output_dir,
      "mapping_export_summary.csv"
    )
  )

  all_written <- c(
    primary_krsa_path,
    primary_full_path,
    primary_csv_path,

    experimental_path,
    kinase_library_path,
    gps6_path,
    predictive_concordant_path,
    concordant_path,

    experimental_full_path,
    kinase_library_full_path,
    gps6_full_path,
    predictive_concordant_full_path,
    concordant_full_path,

    chip_sites_path,
    master_evidence_path,
    family_evidence_path,

    qc_paths,
    coverage_paths,
    export_summary_path
  )

  all_written <- unique(
    normalizePath(
      all_written,
      winslash = "/",
      mustWork = TRUE
    )
  )

  message(
    "Mapping outputs written to: ",
    normalizePath(
      output_dir,
      winslash = "/",
      mustWork = TRUE
    )
  )

  message(
    "Primary KRSA mapping: ",
    primary_krsa_path
  )

  all_written
}
