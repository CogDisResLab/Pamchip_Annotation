# R/mapping.R
#
# Evidence integration and KRSA-compatible mapping construction for the
# PamChip annotation pipeline.
#
# Design principle
# ----------------
# PhosphoSIGNOR, GPS6, and The Kinase Library are preserved as separate
# evidence streams because they represent different evidence classes and use
# different native score scales. No weighted cross-source score is calculated.
#
# Pipeline contract
# -----------------
# This file provides the mapping functions expected by _targets.R:
#
#   build_master_evidence_table(
#     chip_sites,
#     harmonized_evidence,
#     manifest
#   )
#
#   collapse_to_kinase_family(
#     evidence,
#     manifest
#   )
#
#   build_krsa_mapping(
#     family_evidence,
#     mapping_type,
#     manifest
#   )
#
# The resulting architecture is:
#
#   source-specific harmonized evidence
#       -> master evidence table
#       -> peptide x kinase-family evidence table
#       -> primary / sensitivity / concordance mappings


# ==============================================================================
# Small internal helpers
# ==============================================================================

.mapping_null_coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    y
  } else {
    x
  }
}


.clean_mapping_character <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}


.clean_mapping_uniprot <- function(x) {
  x <- .clean_mapping_character(x)
  x <- toupper(x)
  x <- sub("-\\d+$", "", x)
  x
}


.clean_mapping_phosphosite <- function(x) {
  x <- .clean_mapping_character(x)
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  x
}


.any_true <- function(x) {
  x <- as.logical(x)
  any(x %in% TRUE, na.rm = TRUE)
}


.max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    NA_real_
  } else {
    max(x)
  }
}


.collapse_character_values <- function(x) {
  x <- .clean_mapping_character(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0L) {
    NA_character_
  } else {
    paste(x, collapse = ";")
  }
}


# ==============================================================================
# Canonical schemas
# ==============================================================================

.empty_master_evidence_table <- function() {

  tibble::tibble(
    peptide_id = character(),
    substrate_uniprot = character(),
    phosphosite = character(),
    peptide_site_position = integer(),
    phosphoacceptor = character(),
    kinase_gene = character(),
    kinase_uniprot = character(),
    kinase_family = character(),
    kinase_superfamily = character(),
    kinase_group = character(),
    source_native_kinase = character(),
    evidence_source = character(),
    evidence_class = character(),
    source_score = double(),
    source_cutoff = double(),
    source_supported = logical(),
    signor_supported = logical(),
    signor_score = double(),
    signor_pmid = character(),
    kinase_library_percentile = double(),
    kinase_library_supported = logical(),
    gps6_node = character(),
    gps6_node_root = character(),
    gps6_node_parent = character(),
    gps6_node_leaf = character(),
    gps6_node_depth = integer(),
    gps6_peptide_window = character(),
    gps6_score = double(),
    gps6_cutoff = double(),
    gps6_supported = logical()
  )
}


.empty_family_evidence_table <- function() {

  tibble::tibble(
    peptide_id = character(),
    substrate_uniprot = character(),
    phosphosite = character(),
    peptide_site_position = integer(),
    phosphoacceptor = character(),
    kinase_family = character(),
    kinase_superfamily = character(),
    kinase_group = character(),
    kinase_genes = character(),
    kinase_uniprots = character(),
    evidence_sources = character(),
    evidence_classes = character(),
    source_record_count = integer(),
    evidence_source_count = integer(),
    predictive_source_count = integer(),
    experimental_supported = logical(),
    signor_supported = logical(),
    kinase_library_supported = logical(),
    gps6_supported = logical(),
    predictive_concordant = logical(),
    concordant = logical(),
    max_signor_score = double(),
    signor_pmids = character(),
    max_kinase_library_percentile = double(),
    max_gps6_score = double(),
    gps6_nodes = character()
  )
}


.ensure_master_columns <- function(data) {

  template <- .empty_master_evidence_table()

  for (column in names(template)) {

    if (!column %in% names(data)) {

      prototype <- template[[column]]

      if (is.character(prototype)) {
        data[[column]] <- rep(NA_character_, nrow(data))
      } else if (is.integer(prototype)) {
        data[[column]] <- rep(NA_integer_, nrow(data))
      } else if (is.logical(prototype)) {
        data[[column]] <- rep(NA, nrow(data))
      } else {
        data[[column]] <- rep(NA_real_, nrow(data))
      }
    }
  }

  data
}


# ==============================================================================
# Chip-site crosswalk for evidence rows lacking peptide IDs
# ==============================================================================

.build_chip_crosswalk <- function(chip_sites) {

  if (is.null(chip_sites) || !is.data.frame(chip_sites) || nrow(chip_sites) == 0L) {
    return(tibble::tibble())
  }

  required <- c(
    "peptide_id",
    "substrate_uniprot",
    "phosphosite"
  )

  if (!all(required %in% names(chip_sites))) {
    return(tibble::tibble())
  }

  crosswalk <- tibble::tibble(
    peptide_id_chip = as.character(chip_sites$peptide_id),
    substrate_uniprot_key = .clean_mapping_uniprot(chip_sites$substrate_uniprot),
    phosphosite_key = .clean_mapping_phosphosite(chip_sites$phosphosite),
    peptide_site_position_chip =
      if ("phosphosite_position" %in% names(chip_sites)) {
        suppressWarnings(as.integer(chip_sites$phosphosite_position))
      } else {
        rep(NA_integer_, nrow(chip_sites))
      },
    phosphoacceptor_chip =
      if ("residue" %in% names(chip_sites)) {
        as.character(chip_sites$residue)
      } else if ("phosphosite_residue" %in% names(chip_sites)) {
        as.character(chip_sites$phosphosite_residue)
      } else {
        rep(NA_character_, nrow(chip_sites))
      },
    res_position_chip =
      if ("res_position" %in% names(chip_sites)) {
        suppressWarnings(as.integer(chip_sites$res_position))
      } else {
        rep(NA_integer_, nrow(chip_sites))
      }
  )

  crosswalk <- dplyr::filter(
    crosswalk,
    !is.na(.data$peptide_id_chip),
    !is.na(.data$substrate_uniprot_key)
  )

  dplyr::distinct(crosswalk)
}


.backfill_chip_linkage <- function(master, chip_sites) {

  crosswalk <- .build_chip_crosswalk(chip_sites)

  if (nrow(crosswalk) == 0L) {
    return(master)
  }

  master$.mapping_row_id <- seq_len(nrow(master))

  unresolved <- is.na(master$peptide_id) |
    !nzchar(trimws(master$peptide_id))

  if (!any(unresolved)) {
    master$.mapping_row_id <- NULL
    return(master)
  }

  # --------------------------------------------------------------------------
  # Pass 1: exact UniProt + phosphosite label match.
  # Only uniquely resolvable chip matches are used.
  # --------------------------------------------------------------------------

  exact_crosswalk <- dplyr::filter(
    crosswalk,
    !is.na(.data$substrate_uniprot_key),
    !is.na(.data$phosphosite_key)
  )

  exact_crosswalk <- dplyr::group_by(
    exact_crosswalk,
    .data$substrate_uniprot_key,
    .data$phosphosite_key
  )

  exact_crosswalk <- dplyr::filter(
    exact_crosswalk,
    dplyr::n_distinct(.data$peptide_id_chip) == 1L
  )

  exact_crosswalk <- dplyr::ungroup(
    exact_crosswalk
  )

  exact_crosswalk <- dplyr::distinct(
    exact_crosswalk,
    .data$substrate_uniprot_key,
    .data$phosphosite_key,
    .keep_all = TRUE
  )

  if (nrow(exact_crosswalk) > 0L) {

    lookup_input <- tibble::tibble(
      .mapping_row_id = master$.mapping_row_id[unresolved],
      substrate_uniprot_key = .clean_mapping_uniprot(master$substrate_uniprot[unresolved]),
      phosphosite_key = .clean_mapping_phosphosite(master$phosphosite[unresolved])
    )

    matched <- dplyr::left_join(
      lookup_input,
      exact_crosswalk,
      by = c("substrate_uniprot_key", "phosphosite_key")
    )

    matched <- dplyr::filter(
      matched,
      !is.na(.data$peptide_id_chip)
    )

    if (nrow(matched) > 0L) {

      idx <- match(
        matched$.mapping_row_id,
        master$.mapping_row_id
      )

      master$peptide_id[idx] <- matched$peptide_id_chip

      fill_position <- is.na(master$peptide_site_position[idx])
      master$peptide_site_position[idx[fill_position]] <-
        matched$peptide_site_position_chip[fill_position]

      fill_residue <- is.na(master$phosphoacceptor[idx])
      master$phosphoacceptor[idx[fill_residue]] <-
        matched$phosphoacceptor_chip[fill_residue]
    }
  }

  # --------------------------------------------------------------------------
  # Pass 2: UniProt + protein residue coordinate, when both tables expose it.
  # This is useful for curated resources that carry numeric positions but not
  # a preformatted phosphosite label.
  # --------------------------------------------------------------------------

  unresolved <- is.na(master$peptide_id) |
    !nzchar(trimws(master$peptide_id))

  if (
    any(unresolved) &&
    "res_position" %in% names(master) &&
    any(!is.na(crosswalk$res_position_chip))
  ) {

    position_crosswalk <- dplyr::filter(
      crosswalk,
      !is.na(.data$substrate_uniprot_key),
      !is.na(.data$res_position_chip)
    )

    position_crosswalk <- dplyr::group_by(
      position_crosswalk,
      .data$substrate_uniprot_key,
      .data$res_position_chip
    )

    position_crosswalk <- dplyr::filter(
      position_crosswalk,
      dplyr::n_distinct(.data$peptide_id_chip) == 1L
    )

    position_crosswalk <- dplyr::ungroup(
      position_crosswalk
    )

    position_crosswalk <- dplyr::distinct(
      position_crosswalk,
      .data$substrate_uniprot_key,
      .data$res_position_chip,
      .keep_all = TRUE
    )

    if (nrow(position_crosswalk) > 0L) {

      lookup_input <- tibble::tibble(
        .mapping_row_id = master$.mapping_row_id[unresolved],
        substrate_uniprot_key = .clean_mapping_uniprot(master$substrate_uniprot[unresolved]),
        res_position_chip = suppressWarnings(as.integer(master$res_position[unresolved]))
      )

      matched <- dplyr::left_join(
        lookup_input,
        position_crosswalk,
        by = c("substrate_uniprot_key", "res_position_chip")
      )

      matched <- dplyr::filter(
        matched,
        !is.na(.data$peptide_id_chip)
      )

      if (nrow(matched) > 0L) {

        idx <- match(
          matched$.mapping_row_id,
          master$.mapping_row_id
        )

        master$peptide_id[idx] <- matched$peptide_id_chip

        fill_position <- is.na(master$peptide_site_position[idx])
        master$peptide_site_position[idx[fill_position]] <-
          matched$peptide_site_position_chip[fill_position]

        fill_residue <- is.na(master$phosphoacceptor[idx])
        master$phosphoacceptor[idx[fill_residue]] <-
          matched$phosphoacceptor_chip[fill_residue]
      }
    }
  }

  master$.mapping_row_id <- NULL
  master
}


# ==============================================================================
# Build master evidence table
# ==============================================================================

#' Build Master Kinase Evidence Table
#'
#' Creates the canonical long-form source-specific evidence table used for
#' downstream kinase-family collapse, coverage/QC reporting, and sensitivity
#' mappings.
#'
#' @param chip_sites Validated chip-site table.
#' @param harmonized_evidence Output of harmonize_kinase_evidence().
#' @param manifest Parsed pipeline manifest. Accepted for pipeline/provenance
#'   compatibility; source thresholds should already have been applied upstream.
#'
#' @return Canonical long-form evidence tibble.
#'
#' @export
build_master_evidence_table <- function(
  chip_sites,
  harmonized_evidence,
  manifest = NULL
) {

  if (!is.data.frame(harmonized_evidence)) {
    stop(
      "`harmonized_evidence` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  if (nrow(harmonized_evidence) == 0L) {
    warning(
      "`harmonized_evidence` contains zero rows. Returning an empty master evidence table.",
      call. = FALSE
    )
    return(.empty_master_evidence_table())
  }

  required_columns <- c(
    "peptide_id",
    "evidence_source",
    "evidence_class"
  )

  missing_required <- setdiff(
    required_columns,
    names(harmonized_evidence)
  )

  if (length(missing_required) > 0L) {
    stop(
      "Harmonized evidence is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  master <- .ensure_master_columns(
    harmonized_evidence
  )

  master <- dplyr::mutate(
    master,
    peptide_id = .clean_mapping_character(.data$peptide_id),
    substrate_uniprot = .clean_mapping_uniprot(.data$substrate_uniprot),
    phosphosite = .clean_mapping_phosphosite(.data$phosphosite),
    peptide_site_position = suppressWarnings(as.integer(.data$peptide_site_position)),
    phosphoacceptor = toupper(.clean_mapping_character(.data$phosphoacceptor)),
    kinase_gene = toupper(.clean_mapping_character(.data$kinase_gene)),
    kinase_uniprot = .clean_mapping_uniprot(.data$kinase_uniprot),
    kinase_family = toupper(.clean_mapping_character(.data$kinase_family)),
    kinase_superfamily = toupper(.clean_mapping_character(.data$kinase_superfamily)),
    kinase_group = toupper(.clean_mapping_character(.data$kinase_group)),
    source_native_kinase = .clean_mapping_character(.data$source_native_kinase),
    evidence_source = .clean_mapping_character(.data$evidence_source),
    evidence_class = .clean_mapping_character(.data$evidence_class),
    source_score = suppressWarnings(as.numeric(.data$source_score)),
    source_cutoff = suppressWarnings(as.numeric(.data$source_cutoff)),
    source_supported = as.logical(.data$source_supported),
    signor_supported = as.logical(.data$signor_supported),
    signor_score = suppressWarnings(as.numeric(.data$signor_score)),
    signor_pmid = .clean_mapping_character(.data$signor_pmid),
    kinase_library_percentile = suppressWarnings(as.numeric(.data$kinase_library_percentile)),
    kinase_library_supported = as.logical(.data$kinase_library_supported),
    gps6_node = .clean_mapping_character(.data$gps6_node),
    gps6_node_root = .clean_mapping_character(.data$gps6_node_root),
    gps6_node_parent = .clean_mapping_character(.data$gps6_node_parent),
    gps6_node_leaf = .clean_mapping_character(.data$gps6_node_leaf),
    gps6_node_depth = suppressWarnings(as.integer(.data$gps6_node_depth)),
    gps6_peptide_window = .clean_mapping_character(.data$gps6_peptide_window),
    gps6_score = suppressWarnings(as.numeric(.data$gps6_score)),
    gps6_cutoff = suppressWarnings(as.numeric(.data$gps6_cutoff)),
    gps6_supported = as.logical(.data$gps6_supported)
  )

  # Recover chip linkage for curated rows that entered harmonization without a
  # peptide ID. Only unambiguous chip matches are filled.
  master <- .backfill_chip_linkage(
    master = master,
    chip_sites = chip_sites
  )

  # Keep supported rows. NA is retained because curated/legacy evidence may not
  # expose an explicit support flag even though its presence itself is evidence.
  master <- dplyr::filter(
    master,
    is.na(.data$source_supported) |
      .data$source_supported
  )

  master <- dplyr::distinct(master)

  master <- dplyr::arrange(
    master,
    .data$peptide_id,
    .data$peptide_site_position,
    .data$kinase_family,
    .data$kinase_gene,
    .data$evidence_source,
    .data$gps6_node
  )

  canonical_columns <- names(
    .empty_master_evidence_table()
  )

  extra_columns <- setdiff(
    names(master),
    canonical_columns
  )

  master <- dplyr::select(
    master,
    dplyr::all_of(canonical_columns),
    dplyr::all_of(extra_columns)
  )

  missing_peptide <- is.na(master$peptide_id) |
    !nzchar(trimws(master$peptide_id))

  if (any(missing_peptide)) {
    warning(
      sum(missing_peptide),
      " master-evidence row(s) remain unlinked to a chip peptide after unambiguous chip-site matching. ",
      "They are retained for provenance but excluded from peptide-level family mappings.",
      call. = FALSE
    )
  }

  unexpected_source <- setdiff(
    unique(stats::na.omit(master$evidence_source)),
    c("PhosphoSIGNOR", "Kinase Library", "GPS6")
  )

  if (length(unexpected_source) > 0L) {
    warning(
      "Unexpected evidence source label(s): ",
      paste(unexpected_source, collapse = ", "),
      call. = FALSE
    )
  }

  master
}


# ==============================================================================
# Collapse master evidence to peptide x kinase-family resolution
# ==============================================================================

#' Collapse Evidence to Kinase-Family Resolution
#'
#' Produces one row per chip peptide and fine-grained kinase family (e.g. AKT, MAPK, RAF). Evidence-source identity
#' is summarized with explicit support flags rather than combined numerical
#' scores.
#'
#' @param evidence Master evidence table.
#' @param manifest Parsed pipeline manifest.
#'
#' @return Family-level evidence tibble.
#'
#' @export
collapse_to_kinase_family <- function(
  evidence,
  manifest = NULL
) {

  if (!is.data.frame(evidence)) {
    stop(
      "`evidence` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  if (nrow(evidence) == 0L) {
    warning(
      "`evidence` contains zero rows. Returning an empty family evidence table.",
      call. = FALSE
    )
    return(.empty_family_evidence_table())
  }

  required <- c(
    "peptide_id",
    "kinase_family",
    "evidence_source",
    "evidence_class"
  )

  missing_required <- setdiff(required, names(evidence))

  if (length(missing_required) > 0L) {
    stop(
      "Master evidence is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  family <- dplyr::filter(
    evidence,
    !is.na(.data$peptide_id),
    nzchar(trimws(.data$peptide_id)),
    !is.na(.data$kinase_family),
    nzchar(trimws(.data$kinase_family))
  )

  if (nrow(family) == 0L) {
    warning(
      "No chip-linked evidence rows with kinase-family assignments were available.",
      call. = FALSE
    )
    return(.empty_family_evidence_table())
  }

  family <- dplyr::mutate(
    family,
    peptide_id = as.character(.data$peptide_id),
    substrate_uniprot = as.character(.data$substrate_uniprot),
    phosphosite = as.character(.data$phosphosite),
    peptide_site_position = suppressWarnings(as.integer(.data$peptide_site_position)),
    phosphoacceptor = as.character(.data$phosphoacceptor),
    kinase_family = toupper(as.character(.data$kinase_family)),
    kinase_superfamily = toupper(as.character(.data$kinase_superfamily)),
    kinase_group = toupper(as.character(.data$kinase_group)),
    evidence_source = as.character(.data$evidence_source),
    evidence_class = as.character(.data$evidence_class)
  )

  # Collapse at the biological unit used for family-level mapping:
  #
  #   peptide x kinase_family
  #
  # Do NOT include kinase_group in the grouping key. GPS6 and Kinase Library
  # can assign different group-level metadata to the same family, and grouping
  # on kinase_group would split otherwise concordant evidence into separate rows.
  family <- dplyr::group_by(
    family,
    .data$peptide_id,
    .data$substrate_uniprot,
    .data$phosphosite,
    .data$peptide_site_position,
    .data$phosphoacceptor,
    .data$kinase_family
  )

  family <- dplyr::summarise(
    family,

    # Broad superfamily/group labels are retained as metadata rather than used
    # as a grouping key. The primary KRSA mapping remains at the fine family
    # level (e.g. AKT rather than AGC).
    kinase_superfamily =
      .collapse_character_values(.data$kinase_superfamily),

    kinase_group =
      .collapse_character_values(.data$kinase_group),

    kinase_genes =
      .collapse_character_values(.data$kinase_gene),

    kinase_uniprots =
      .collapse_character_values(.data$kinase_uniprot),

    evidence_sources =
      .collapse_character_values(.data$evidence_source),

    evidence_classes =
      .collapse_character_values(.data$evidence_class),

    source_record_count =
      as.integer(dplyr::n()),

    signor_supported =
      .any_true(
        (.data$evidence_source == "PhosphoSIGNOR") &
          (is.na(.data$source_supported) | .data$source_supported)
      ),

    kinase_library_supported =
      .any_true(
        (.data$evidence_source == "Kinase Library") &
          (is.na(.data$source_supported) | .data$source_supported)
      ),

    gps6_supported =
      .any_true(
        (.data$evidence_source == "GPS6") &
          (is.na(.data$source_supported) | .data$source_supported)
      ),

    max_signor_score =
      .max_or_na(.data$signor_score),

    signor_pmids =
      .collapse_character_values(.data$signor_pmid),

    max_kinase_library_percentile =
      .max_or_na(.data$kinase_library_percentile),

    max_gps6_score =
      .max_or_na(.data$gps6_score),

    gps6_nodes =
      .collapse_character_values(.data$gps6_node),

    .groups = "drop"
  )

  family <- dplyr::mutate(
    family,

    experimental_supported =
      .data$signor_supported,

    predictive_source_count =
      as.integer(.data$kinase_library_supported) +
      as.integer(.data$gps6_supported),

    evidence_source_count =
      as.integer(.data$signor_supported) +
      as.integer(.data$kinase_library_supported) +
      as.integer(.data$gps6_supported),

    predictive_concordant =
      .data$kinase_library_supported &
      .data$gps6_supported,

    concordant =
      .data$evidence_source_count >= 2L
  )

  family <- dplyr::select(
    family,
    dplyr::all_of(names(.empty_family_evidence_table()))
  )

  dplyr::arrange(
    family,
    .data$peptide_id,
    .data$kinase_family
  )
}


# ==============================================================================
# KRSA-compatible mapping generation
# ==============================================================================

.get_primary_mapping_rule <- function(manifest) {

  if (is.null(manifest)) {
    return("any_source")
  }

  candidates <- list(
    manifest$mapping$primary_rule,
    manifest$mapping$inclusion_rule,
    manifest$mappings$primary_rule,
    manifest$mappings$inclusion_rule
  )

  for (candidate in candidates) {
    if (!is.null(candidate) && length(candidate) > 0L && !is.na(candidate[[1L]])) {
      return(tolower(as.character(candidate[[1L]])))
    }
  }

  "any_source"
}


.primary_mapping_mask <- function(family_evidence, manifest) {

  rule <- .get_primary_mapping_rule(manifest)

  if (rule %in% c("any", "any_source", "union")) {
    return(family_evidence$evidence_source_count >= 1L)
  }

  if (rule %in% c("concordant", "multi_source", "two_sources", ">=2_sources")) {
    return(family_evidence$evidence_source_count >= 2L)
  }

  if (rule %in% c("predictive_concordant", "predictive_agreement")) {
    return(family_evidence$predictive_concordant)
  }

  if (rule %in% c("experimental", "experimental_only", "signor")) {
    return(family_evidence$signor_supported)
  }

  if (rule %in% c("experimental_or_concordant", "signor_or_concordant")) {
    return(
      family_evidence$signor_supported |
        family_evidence$concordant
    )
  }

  warning(
    "Unknown primary mapping rule '",
    rule,
    "'. Falling back to any_source.",
    call. = FALSE
  )

  family_evidence$evidence_source_count >= 1L
}


#' Build a KRSA-Compatible Peptide-to-Kinase-Family Mapping
#'
#' @param family_evidence Output of collapse_to_kinase_family().
#' @param mapping_type One of primary, experimental_only, kinase_library_only,
#'   gps6_only, predictive_concordant, or concordant.
#' @param manifest Parsed pipeline manifest.
#'
#' @return Tibble containing retained peptide-family relationships plus
#'   evidence provenance columns and simple Peptide/Kinase aliases for KRSA; `Kinase` uses the fine-grained kinase_family
#'   interoperability.
#'
#' @export
build_krsa_mapping <- function(
  family_evidence,
  mapping_type = "primary",
  manifest = NULL
) {

  if (!is.data.frame(family_evidence)) {
    stop(
      "`family_evidence` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  mapping_type <- tolower(as.character(mapping_type[[1L]]))

  valid_types <- c(
    "primary",
    "experimental_only",
    "kinase_library_only",
    "gps6_only",
    "predictive_concordant",
    "concordant"
  )

  if (!mapping_type %in% valid_types) {
    stop(
      "Unsupported mapping_type '",
      mapping_type,
      "'. Valid values are: ",
      paste(valid_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (nrow(family_evidence) == 0L) {
    return(
      tibble::tibble(
        peptide_id = character(),
        kinase_family = character(),
        Peptide = character(),
        Kinase = character(),
        mapping_type = character()
      )
    )
  }

  required <- c(
    "peptide_id",
    "kinase_family",
    "signor_supported",
    "kinase_library_supported",
    "gps6_supported",
    "predictive_concordant",
    "concordant",
    "evidence_source_count"
  )

  missing_required <- setdiff(required, names(family_evidence))

  if (length(missing_required) > 0L) {
    stop(
      "Family evidence is missing required column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  keep <- switch(
    mapping_type,

    primary =
      .primary_mapping_mask(
        family_evidence,
        manifest
      ),

    experimental_only =
      family_evidence$signor_supported,

    kinase_library_only =
      family_evidence$kinase_library_supported,

    gps6_only =
      family_evidence$gps6_supported,

    predictive_concordant =
      family_evidence$predictive_concordant,

    concordant =
      family_evidence$concordant
  )

  keep[is.na(keep)] <- FALSE

  mapping <- family_evidence[keep, , drop = FALSE]

  mapping <- dplyr::mutate(
    mapping,
    Peptide = as.character(.data$peptide_id),
    Kinase = as.character(.data$kinase_family),
    mapping_type = mapping_type
  )

  mapping <- dplyr::distinct(
    mapping,
    .data$peptide_id,
    .data$kinase_family,
    .keep_all = TRUE
  )

  mapping <- dplyr::arrange(
    mapping,
    .data$peptide_id,
    .data$kinase_family
  )

  mapping <- dplyr::select(
    mapping,
    "peptide_id",
    "kinase_family",
    "Peptide",
    "Kinase",
    "mapping_type",
    dplyr::everything()
  )

  mapping
}


# ==============================================================================
# Deprecated weighted-ensemble interface
# ==============================================================================

#' Deprecated Weighted Ensemble Mapping Interface
#'
#' The original implementation numerically combined PhosphoSIGNOR, GPS6, and
#' Kinase Library scores. That method is intentionally disabled because the
#' source-native scores are not directly comparable.
#'
#' @export
build_ensemble_mapping <- function(
  chip_sites,
  harmonized_sources,
  weights
) {

  stop(
    "`build_ensemble_mapping()` is deprecated and intentionally disabled. ",
    "Use `build_master_evidence_table()`, `collapse_to_kinase_family()`, and ",
    "`build_krsa_mapping()` instead.",
    call. = FALSE
  )
}
