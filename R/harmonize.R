# R/harmonize.R

#' Harmonize Kinase Evidence Across PhosphoSIGNOR, Kinase Library, and GPS 6.0
#'
#' Standardizes source-specific evidence into a common schema while preserving
#' the native evidence fields required for provenance and downstream family-level
#' collapse.
#'
#' Important design principle:
#' GPS6 model nodes are hierarchical and should NOT automatically be interpreted
#' as individual HGNC kinase symbols. For example:
#'
#'   AGC/Akt
#'   AGC/Akt/AKT1
#'
#' represent different levels of the GPS6 kinase hierarchy.
#'
#' Therefore, this function preserves the full GPS6 node and derives:
#'
#'   kinase_group   = top-level GPS6 group (e.g. AGC)
#'   kinase_family  = second hierarchy level when available (e.g. Akt)
#'   kinase_gene    = deepest node only when the GPS6 node is sufficiently
#'                    specific to represent an individual kinase
#'
#' Final kinase-family standardization/crosswalking can then occur downstream.
#'
#' @param signor_data Processed PhosphoSIGNOR evidence.
#' @param kinase_library_data Processed Kinase Library evidence.
#' @param gps6_data Processed GPS 6.0 evidence from process_gps6_mappings().
#'
#' @return A tibble containing harmonized evidence from all three sources.
#'
#' @export
harmonize_kinase_evidence <- function(
  signor_data,
  kinase_library_data,
  gps6_data
) {

  # ============================================================================
  # Internal helpers
  # ============================================================================

  # Return the first matching column among several possible source schemas.
  #
  # This keeps the harmonizer tolerant of minor naming differences between
  # source-processing functions without silently guessing when a required field
  # is absent.
  get_column <- function(
    data,
    candidates,
    default = NA
  ) {

    hit <- candidates[candidates %in% names(data)]

    if (length(hit) > 0L) {
      return(data[[hit[[1L]]]])
    }

    rep(
      default,
      nrow(data)
    )
  }


  # --------------------------------------------------------------------------
  # Normalize symbols without destroying biologically meaningful identity.
  #
  # Unlike the previous implementation, "/" is NOT blindly converted to "_"
  # before GPS6 hierarchy parsing.
  # --------------------------------------------------------------------------

  clean_symbol <- function(x) {

    x <- as.character(x)

    x <- stringr::str_trim(x)

    x[x == ""] <- NA_character_

    stringr::str_to_upper(x)
  }


  # --------------------------------------------------------------------------
  # Normalize UniProt accessions
  # --------------------------------------------------------------------------

  clean_uniprot <- function(x) {

    x <- as.character(x)

    x <- stringr::str_trim(x)

    x[x == ""] <- NA_character_

    x
  }


  # --------------------------------------------------------------------------
  # Extract an integer residue/site position where possible
  # --------------------------------------------------------------------------

  clean_position <- function(x) {

    x <- as.character(x)

    extracted <- stringr::str_extract(
      x,
      "[0-9]+"
    )

    suppressWarnings(
      as.integer(extracted)
    )
  }


  # --------------------------------------------------------------------------
  # Normalize logical/support fields
  # --------------------------------------------------------------------------

  as_support_flag <- function(x) {

    if (is.logical(x)) {
      return(x)
    }

    if (is.numeric(x)) {
      return(
        dplyr::case_when(
          is.na(x) ~ NA,
          x != 0 ~ TRUE,
          TRUE ~ FALSE
        )
      )
    }

    x <- stringr::str_to_lower(
      stringr::str_trim(
        as.character(x)
      )
    )

    dplyr::case_when(
      x %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y",
        "supported"
      ) ~ TRUE,

      x %in% c(
        "false",
        "f",
        "0",
        "no",
        "n",
        "unsupported"
      ) ~ FALSE,

      TRUE ~ NA
    )
  }


  # --------------------------------------------------------------------------
  # Extract phosphoacceptor residue from a label such as S473 / T308 / Y123
  # --------------------------------------------------------------------------

  residue_from_site <- function(x) {

    x <- stringr::str_to_upper(
      as.character(x)
    )

    residue <- stringr::str_extract(
      x,
      "^[STY]"
    )

    residue
  }


  # --------------------------------------------------------------------------
  # Clean GPS6 hierarchy strings
  # --------------------------------------------------------------------------

  clean_gps6_node <- function(x) {

    x <- as.character(x)

    x <- stringr::str_trim(x)

    x <- stringr::str_replace_all(
      x,
      "^/+|/+$",
      ""
    )

    x[x == ""] <- NA_character_

    x
  }


  # --------------------------------------------------------------------------
  # Extract GPS6 hierarchy level
  # --------------------------------------------------------------------------

  gps6_level <- function(
    node,
    level
  ) {

    node <- clean_gps6_node(node)

    vapply(
      strsplit(
        ifelse(
          is.na(node),
          "",
          node
        ),
        "/",
        fixed = TRUE
      ),
      function(parts) {

        parts <- parts[
          nzchar(parts)
        ]

        if (length(parts) < level) {
          return(NA_character_)
        }

        parts[[level]]
      },
      character(1)
    )
  }


  # --------------------------------------------------------------------------
  # GPS6 deepest node
  # --------------------------------------------------------------------------

  gps6_leaf <- function(node) {

    node <- clean_gps6_node(node)

    vapply(
      strsplit(
        ifelse(
          is.na(node),
          "",
          node
        ),
        "/",
        fixed = TRUE
      ),
      function(parts) {

        parts <- parts[
          nzchar(parts)
        ]

        if (length(parts) == 0L) {
          return(NA_character_)
        }

        parts[[length(parts)]]
      },
      character(1)
    )
  }


  # --------------------------------------------------------------------------
  # GPS6 hierarchy depth
  # --------------------------------------------------------------------------

  gps6_depth <- function(node) {

    node <- clean_gps6_node(node)

    vapply(
      strsplit(
        ifelse(
          is.na(node),
          "",
          node
        ),
        "/",
        fixed = TRUE
      ),
      function(parts) {

        parts <- parts[
          nzchar(parts)
        ]

        if (length(parts) == 0L) {
          return(NA_integer_)
        }

        as.integer(
          length(parts)
        )
      },
      integer(1)
    )
  }


  # ============================================================================
  # 1. PhosphoSIGNOR
  # ============================================================================

  signor_peptide_id <- get_column(
    signor_data,
    c(
      "peptide_id",
      "id",
      "ID"
    )
  )

  signor_uniprot <- get_column(
    signor_data,
    c(
      "substrate_uniprot",
      "uniprot_id",
      "uniprot",
      "substrate_accession"
    )
  )

  signor_position_raw <- get_column(
    signor_data,
    c(
      "peptide_site_position",
      "res_position",
      "residue_position",
      "position"
    )
  )

  signor_site <- get_column(
    signor_data,
    c(
      "phosphosite",
      "site",
      "site_label"
    )
  )

  # If no explicit site label exists, preserve the position rather than
  # inventing a phosphoacceptor residue.
  signor_site <- dplyr::if_else(
    is.na(signor_site) |
      !nzchar(as.character(signor_site)),
    as.character(signor_position_raw),
    as.character(signor_site)
  )

  signor_kinase <- get_column(
    signor_data,
    c(
      "kinase_gene",
      "kinase_symbol",
      "kinase",
      "enzyme"
    )
  )

  signor_kinase_uniprot <- get_column(
    signor_data,
    c(
      "kinase_uniprot",
      "kinase_accession"
    )
  )

  signor_family <- get_column(
    signor_data,
    c(
      "kinase_family",
      "family"
    )
  )

  signor_score <- get_column(
    signor_data,
    c(
      "score",
      "signor_score"
    ),
    default = NA_real_
  )

  signor_pmid <- get_column(
    signor_data,
    c(
      "signor_pmid",
      "pmid",
      "PMID"
    )
  )

  signor_supported <- get_column(
    signor_data,
    c(
      "signor_supported",
      "supported"
    ),
    default = TRUE
  )

  signor_clean <- tibble::tibble(

    peptide_id =
      as.character(signor_peptide_id),

    substrate_uniprot =
      clean_uniprot(signor_uniprot),

    phosphosite =
      as.character(signor_site),

    peptide_site_position =
      clean_position(signor_position_raw),

    phosphoacceptor =
      residue_from_site(signor_site),

    kinase_gene =
      clean_symbol(signor_kinase),

    kinase_uniprot =
      clean_uniprot(signor_kinase_uniprot),

    kinase_family =
      as.character(signor_family),

    kinase_group =
      NA_character_,

    source_native_kinase =
      as.character(signor_kinase),

    evidence_source =
      "PhosphoSIGNOR",

    evidence_class =
      "experimental",

    source_score =
      suppressWarnings(
        as.numeric(signor_score)
      ),

    source_cutoff =
      NA_real_,

    source_supported =
      as_support_flag(signor_supported),

    # ------------------------------------------------------------------------
    # Source-specific fields
    # ------------------------------------------------------------------------

    signor_supported =
      as_support_flag(signor_supported),

    signor_score =
      suppressWarnings(
        as.numeric(signor_score)
      ),

    signor_pmid =
      as.character(signor_pmid),

    kinase_library_percentile =
      NA_real_,

    kinase_library_supported =
      NA,

    gps6_node =
      NA_character_,

    gps6_node_root =
      NA_character_,

    gps6_node_parent =
      NA_character_,

    gps6_node_leaf =
      NA_character_,

    gps6_node_depth =
      NA_integer_,

    gps6_peptide_window =
      NA_character_,

    gps6_score =
      NA_real_,

    gps6_cutoff =
      NA_real_,

    gps6_supported =
      NA
  )


  # ============================================================================
  # 2. Kinase Library
  # ============================================================================

  kl_peptide_id <- get_column(
    kinase_library_data,
    c(
      "peptide_id",
      "id",
      "ID"
    )
  )

  kl_uniprot <- get_column(
    kinase_library_data,
    c(
      "substrate_uniprot",
      "uniprot_id",
      "uniprot"
    )
  )

  kl_position_raw <- get_column(
    kinase_library_data,
    c(
      "peptide_site_position",
      "res_position",
      "residue_position",
      "position"
    )
  )

  kl_site <- get_column(
    kinase_library_data,
    c(
      "phosphosite",
      "site",
      "site_label"
    )
  )

  kl_site <- dplyr::if_else(
    is.na(kl_site) |
      !nzchar(as.character(kl_site)),
    as.character(kl_position_raw),
    as.character(kl_site)
  )

  kl_kinase <- get_column(
    kinase_library_data,
    c(
      "kinase_gene",
      "kinase_symbol",
      "kinase"
    )
  )

  kl_kinase_uniprot <- get_column(
    kinase_library_data,
    c(
      "kinase_uniprot",
      "kinase_accession"
    )
  )

  kl_family <- get_column(
    kinase_library_data,
    c(
      "kinase_family",
      "family"
    )
  )

  kl_percentile <- get_column(
    kinase_library_data,
    c(
      "kinase_library_percentile",
      "percentile",
      "score"
    ),
    default = NA_real_
  )

  kl_supported <- get_column(
    kinase_library_data,
    c(
      "kinase_library_supported",
      "supported"
    ),
    default = NA
  )

  # If a support flag was not explicitly provided, infer it only if the
  # processed Kinase Library table already represents a filtered set.
  #
  # Do NOT impose another cutoff here because the cutoff belongs upstream in
  # run_kinase_library_scoring().
  if (all(is.na(kl_supported))) {
    kl_supported <- rep(
      TRUE,
      nrow(kinase_library_data)
    )
  }

  kl_clean <- tibble::tibble(

    peptide_id =
      as.character(kl_peptide_id),

    substrate_uniprot =
      clean_uniprot(kl_uniprot),

    phosphosite =
      as.character(kl_site),

    peptide_site_position =
      clean_position(kl_position_raw),

    phosphoacceptor =
      residue_from_site(kl_site),

    kinase_gene =
      clean_symbol(kl_kinase),

    kinase_uniprot =
      clean_uniprot(kl_kinase_uniprot),

    kinase_family =
      as.character(kl_family),

    kinase_group =
      NA_character_,

    source_native_kinase =
      as.character(kl_kinase),

    evidence_source =
      "Kinase Library",

    evidence_class =
      "predictive",

    source_score =
      suppressWarnings(
        as.numeric(kl_percentile)
      ),

    source_cutoff =
      NA_real_,

    source_supported =
      as_support_flag(kl_supported),

    # ------------------------------------------------------------------------
    # Source-specific fields
    # ------------------------------------------------------------------------

    signor_supported =
      NA,

    signor_score =
      NA_real_,

    signor_pmid =
      NA_character_,

    kinase_library_percentile =
      suppressWarnings(
        as.numeric(kl_percentile)
      ),

    kinase_library_supported =
      as_support_flag(kl_supported),

    gps6_node =
      NA_character_,

    gps6_node_root =
      NA_character_,

    gps6_node_parent =
      NA_character_,

    gps6_node_leaf =
      NA_character_,

    gps6_node_depth =
      NA_integer_,

    gps6_peptide_window =
      NA_character_,

    gps6_score =
      NA_real_,

    gps6_cutoff =
      NA_real_,

    gps6_supported =
      NA
  )


  # ============================================================================
  # 3. GPS 6.0
  # ============================================================================

  gps6_node <- clean_gps6_node(
    get_column(
      gps6_data,
      c(
        "gps6_node",
        "Kinase",
        "kinase"
      )
    )
  )

  gps6_root <- get_column(
    gps6_data,
    c(
      "gps6_node_root"
    )
  )

  gps6_parent <- get_column(
    gps6_data,
    c(
      "gps6_node_parent"
    )
  )

  gps6_leaf_value <- get_column(
    gps6_data,
    c(
      "gps6_node_leaf"
    )
  )

  gps6_depth_value <- get_column(
    gps6_data,
    c(
      "gps6_node_depth"
    ),
    default = NA_integer_
  )

  # Derive hierarchy metadata if process_gps6_mappings() did not already
  # provide it.
  missing_root <- is.na(gps6_root) |
    !nzchar(as.character(gps6_root))

  gps6_root[missing_root] <-
    gps6_level(
      gps6_node[missing_root],
      1L
    )

  missing_leaf <- is.na(gps6_leaf_value) |
    !nzchar(as.character(gps6_leaf_value))

  gps6_leaf_value[missing_leaf] <-
    gps6_leaf(
      gps6_node[missing_leaf]
    )

  missing_depth <- is.na(gps6_depth_value)

  gps6_depth_value[missing_depth] <-
    gps6_depth(
      gps6_node[missing_depth]
    )

  # --------------------------------------------------------------------------
  # Derive GPS6 family.
  #
  # For nodes such as:
  #
  #   AGC/Akt
  #   AGC/Akt/AKT1
  #
  # level 2 ("Akt") is the useful family/subfamily representation.
  #
  # Root-only nodes such as "AGC" are retained as group-level evidence and do
  # not receive a fabricated kinase_family.
  # --------------------------------------------------------------------------

  gps6_family <- gps6_level(
    gps6_node,
    2L
  )

  # --------------------------------------------------------------------------
  # Derive a candidate individual kinase gene only from sufficiently deep
  # hierarchy nodes.
  #
  # This deliberately avoids interpreting:
  #
  #   AGC/Akt
  #
  # as if "Akt" were a specific HGNC kinase.
  #
  # Individual-node names are still only SOURCE-NATIVE candidates here. A
  # downstream crosswalk can standardize aliases if needed.
  # --------------------------------------------------------------------------

  gps6_gene <- dplyr::if_else(
    gps6_depth_value >= 3L,
    clean_symbol(gps6_leaf_value),
    NA_character_
  )

  gps6_score_value <- get_column(
    gps6_data,
    c(
      "gps6_score",
      "Score",
      "score"
    ),
    default = NA_real_
  )

  gps6_cutoff_value <- get_column(
    gps6_data,
    c(
      "gps6_cutoff",
      "Cutoff",
      "cutoff"
    ),
    default = NA_real_
  )

  gps6_supported_value <- get_column(
    gps6_data,
    c(
      "gps6_supported",
      "Pass",
      "supported"
    ),
    default = TRUE
  )

  gps6_clean <- tibble::tibble(

    peptide_id =
      as.character(
        get_column(
          gps6_data,
          c(
            "peptide_id",
            "ID",
            "id"
          )
        )
      ),

    substrate_uniprot =
      clean_uniprot(
        get_column(
          gps6_data,
          c(
            "substrate_uniprot",
            "uniprot_id",
            "uniprot"
          )
        )
      ),

    phosphosite =
      as.character(
        get_column(
          gps6_data,
          c(
            "phosphosite",
            "site"
          )
        )
      ),

    peptide_site_position =
      clean_position(
        get_column(
          gps6_data,
          c(
            "peptide_site_position",
            "Position",
            "position"
          )
        )
      ),

    phosphoacceptor =
      stringr::str_to_upper(
        as.character(
          get_column(
            gps6_data,
            c(
              "phosphoacceptor",
              "Code",
              "residue"
            )
          )
        )
      ),

    kinase_gene =
      gps6_gene,

    # GPS6 itself does not provide a canonical UniProt kinase accession.
    kinase_uniprot =
      NA_character_,

    kinase_family =
      as.character(gps6_family),

    kinase_group =
      as.character(gps6_root),

    source_native_kinase =
      as.character(gps6_leaf_value),

    evidence_source =
      "GPS6",

    evidence_class =
      "predictive",

    source_score =
      suppressWarnings(
        as.numeric(gps6_score_value)
      ),

    source_cutoff =
      suppressWarnings(
        as.numeric(gps6_cutoff_value)
      ),

    source_supported =
      as_support_flag(gps6_supported_value),

    # ------------------------------------------------------------------------
    # Source-specific fields
    # ------------------------------------------------------------------------

    signor_supported =
      NA,

    signor_score =
      NA_real_,

    signor_pmid =
      NA_character_,

    kinase_library_percentile =
      NA_real_,

    kinase_library_supported =
      NA,

    gps6_node =
      as.character(gps6_node),

    gps6_node_root =
      as.character(gps6_root),

    gps6_node_parent =
      as.character(gps6_parent),

    gps6_node_leaf =
      as.character(gps6_leaf_value),

    gps6_node_depth =
      as.integer(gps6_depth_value),

    gps6_peptide_window =
      as.character(
        get_column(
          gps6_data,
          c(
            "gps6_peptide_window",
            "Peptide"
          )
        )
      ),

    gps6_score =
      suppressWarnings(
        as.numeric(gps6_score_value)
      ),

    gps6_cutoff =
      suppressWarnings(
        as.numeric(gps6_cutoff_value)
      ),

    gps6_supported =
      as_support_flag(gps6_supported_value)
  )


  # ============================================================================
  # 4. Combine Sources
  # ============================================================================

  harmonized <- dplyr::bind_rows(
    signor_clean,
    kl_clean,
    gps6_clean
  )


  # ============================================================================
  # 5. Final normalization
  # ============================================================================

  harmonized <- harmonized %>%

    dplyr::mutate(

      peptide_id =
        as.character(.data$peptide_id),

      substrate_uniprot =
        clean_uniprot(.data$substrate_uniprot),

      kinase_gene =
        clean_symbol(.data$kinase_gene),

      kinase_uniprot =
        clean_uniprot(.data$kinase_uniprot),

      kinase_family =
        dplyr::na_if(
          stringr::str_trim(
            as.character(.data$kinase_family)
          ),
          ""
        ),

      kinase_group =
        dplyr::na_if(
          stringr::str_trim(
            as.character(.data$kinase_group)
          ),
          ""
        ),

      evidence_source =
        as.character(.data$evidence_source),

      evidence_class =
        as.character(.data$evidence_class)
    ) %>%

    # Remove only exact duplicate evidence rows.
    #
    # We intentionally DO NOT collapse evidence from different GPS6 hierarchy
    # levels or different sources here.
    dplyr::distinct() %>%

    dplyr::arrange(
      .data$peptide_id,
      .data$peptide_site_position,
      .data$kinase_family,
      .data$kinase_gene,
      .data$evidence_source,
      .data$gps6_node
    )


  # ============================================================================
  # 6. Basic integrity checks
  # ============================================================================

  if (anyNA(harmonized$peptide_id)) {
    warning(
      "Harmonized evidence contains rows with missing peptide_id.",
      call. = FALSE
    )
  }

  unknown_source <- setdiff(
    unique(harmonized$evidence_source),
    c(
      "PhosphoSIGNOR",
      "Kinase Library",
      "GPS6"
    )
  )

  if (length(unknown_source) > 0L) {
    warning(
      "Unexpected evidence source label(s) after harmonization: ",
      paste(
        unknown_source,
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  harmonized
}


# ==============================================================================
# Backward-compatible wrapper
# ==============================================================================

#' @rdname harmonize_kinase_evidence
#' @export
harmonize_kinase_names <- function(
  signor_data,
  gps6_data,
  kl_data
) {

  harmonize_kinase_evidence(
    signor_data = signor_data,
    kinase_library_data = kl_data,
    gps6_data = gps6_data
  )
}
