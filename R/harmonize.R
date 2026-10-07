# R/harmonize.R

# ==============================================================================
# GPS6 hierarchy crosswalk
# ==============================================================================

#' Build a Canonical Kinase Hierarchy from GPS6 Model Assets
#'
#' Uses the GPS6 model directory structure as the fine-grained kinase hierarchy.
#' For example:
#'
#'   AGC/Akt/AKT1
#'
#' becomes:
#'
#'   kinase_gene        = AKT1
#'   kinase_family      = AKT
#'   kinase_superfamily = AGC
#'
#' Only runnable GPS6 model directories containing comb_LRM.pkl are considered.
#' Nodes with fewer than three hierarchy levels are excluded from the gene-level
#' crosswalk because they represent groups/families rather than individual
#' kinases.
#'
#' @param gps6_model_files Character vector of tracked GPS6 model/resource files.
#' @param pre_dir GPS6 pre directory containing ST/ and optionally Y/.
#'
#' @return A tibble mapping source-native GPS6 kinase leaves to canonical
#'   gene/family/superfamily labels.
#'
#' @export
build_gps6_kinase_hierarchy <- function(
  gps6_model_files,
  pre_dir = "data/pre"
) {

  files <- normalizePath(
    gps6_model_files,
    winslash = "/",
    mustWork = FALSE
  )

  pre_dir <- normalizePath(
    pre_dir,
    winslash = "/",
    mustWork = FALSE
  )

  ensemble_files <- files[
    basename(files) == "comb_LRM.pkl"
  ]

  if (length(ensemble_files) == 0L) {
    stop(
      "No GPS6 comb_LRM.pkl model files were available to build the kinase hierarchy.",
      call. = FALSE
    )
  }

  model_dirs <- dirname(ensemble_files)

  parse_model_dir <- function(model_dir) {

    relative <- sub(
      paste0(
        "^",
        gsub(
          "([][{}()+*^$|\\\\?.])",
          "\\\\\\1",
          pre_dir
        ),
        "/?"
      ),
      "",
      model_dir
    )

    parts <- strsplit(
      relative,
      "/",
      fixed = TRUE
    )[[1L]]

    parts <- parts[
      nzchar(parts)
    ]

    if (length(parts) < 2L) {
      return(NULL)
    }

    kinase_type <- parts[[1L]]
    hierarchy <- parts[-1L]

    if (length(hierarchy) < 3L) {
      return(NULL)
    }

    tibble::tibble(
      kinase_type =
        as.character(kinase_type),

      gps6_hierarchy =
        paste(
          hierarchy,
          collapse = "/"
        ),

      gps6_source_gene =
        toupper(
          hierarchy[[length(hierarchy)]]
        ),

      kinase_gene =
        toupper(
          hierarchy[[length(hierarchy)]]
        ),

      kinase_family =
        toupper(
          hierarchy[[2L]]
        ),

      kinase_superfamily =
        toupper(
          hierarchy[[1L]]
        )
    )
  }

  hierarchy <- dplyr::bind_rows(
    lapply(
      model_dirs,
      parse_model_dir
    )
  )

  if (nrow(hierarchy) == 0L) {
    stop(
      "GPS6 model files were found, but no individual-kinase hierarchy nodes could be derived.",
      call. = FALSE
    )
  }

  hierarchy <- dplyr::distinct(
    hierarchy
  )

  # ============================================================================
  # Kinase Library metadata helpers
  # ============================================================================

  py_path <- Sys.getenv(
    "RETICULATE_PYTHON",
    unset = "/root/.virtualenvs/pamchip-env/bin/python"
  )

  if (file.exists(py_path)) {

    reticulate::use_python(
      py_path,
      required = TRUE
    )

  } else {

    reticulate::use_virtualenv(
      "pamchip-env",
      required = TRUE
    )
  }

  kl <- reticulate::import(
    "kinase_library",
    convert = FALSE
  )

  first_non_missing <- function(x) {

    x <- as.character(x)

    x <- x[
      !is.na(x) &
        nzchar(
          stringr::str_trim(x)
        )
    ]

    if (length(x) == 0L) {
      return(NA_character_)
    }

    x[[1L]]
  }

  extract_info_field <- function(
    info,
    candidates
  ) {

    if (is.null(info)) {
      return(NA_character_)
    }

    info_r <- tryCatch(
      reticulate::py_to_r(info),
      error = function(e) NULL
    )

    if (is.null(info_r)) {
      return(NA_character_)
    }

    if (is.data.frame(info_r)) {

      names_upper <- toupper(
        names(info_r)
      )

      for (candidate in candidates) {

        hit <- which(
          names_upper ==
            toupper(candidate)
        )

        if (length(hit) > 0L) {

          value <- first_non_missing(
            info_r[[hit[[1L]]]]
          )

          if (!is.na(value)) {
            return(value)
          }
        }
      }
    }

    if (!is.null(names(info_r))) {

      names_upper <- toupper(
        names(info_r)
      )

      for (candidate in candidates) {

        hit <- which(
          names_upper ==
            toupper(candidate)
        )

        if (length(hit) > 0L) {

          value <- first_non_missing(
            info_r[[hit[[1L]]]]
          )

          if (!is.na(value)) {
            return(value)
          }
        }
      }
    }

    NA_character_
  }

  clean_alias <- function(x) {

    x <- stringr::str_to_upper(
      stringr::str_trim(
        as.character(x)
      )
    )

    x[
      is.na(x) |
        !nzchar(x)
    ] <- NA_character_

    x
  }

  # ============================================================================
  # Canonicalize GPS6 source-native leaves
  # ============================================================================
  #
  # GPS6 often uses kinase-library-style matrix aliases rather than canonical
  # gene symbols, e.g.:
  #
  #   ERK2 -> MAPK1
  #   PKCA -> PRKCA
  #   RSK2 -> RPS6KA3
  #
  # Query the installed Kinase Library metadata using name_type="matrix".
  # Preserve the original GPS6 leaf for provenance.

  gps6_aliases <- unique(
    hierarchy$gps6_source_gene
  )

  alias_rows <- lapply(
    gps6_aliases,
    function(alias) {

      info <- tryCatch(
        kl$get_kinase_info(
          alias,
          name_type = "matrix"
        ),
        error = function(e) NULL
      )

      canonical_gene <- extract_info_field(
        info,
        c(
          "GENE_NAME",
          "GENE",
          "HGNC",
          "HGNC_SYMBOL",
          "SYMBOL",
          "GENE SYMBOL",
          "GENE NAME"
        )
      )

      matrix_name <- extract_info_field(
        info,
        c(
          "MATRIX_NAME",
          "KINASE"
        )
      )

      display_name <- extract_info_field(
        info,
        c(
          "DISPLAY_NAME"
        )
      )

      canonical_uniprot <- extract_info_field(
        info,
        c(
          "UNIPROT_ID",
          "UNIPROT",
          "UNIPROT_ACCESSION",
          "ACCESSION"
        )
      )

      tibble::tibble(
        gps6_source_gene =
          clean_alias(alias),

        canonical_kinase_gene =
          clean_alias(canonical_gene),

        kl_matrix_name =
          clean_alias(matrix_name),

        kl_display_name =
          clean_alias(display_name),

        canonical_kinase_uniprot =
          {
            value <- stringr::str_trim(
              as.character(canonical_uniprot)
            )

            if (
              length(value) == 0L ||
              is.na(value[[1L]]) ||
              !nzchar(value[[1L]])
            ) {
              NA_character_
            } else {
              value[[1L]]
            }
          }
      )
    }
  )

  alias_crosswalk <- dplyr::bind_rows(
    alias_rows
  )

  hierarchy <- dplyr::left_join(
    hierarchy,
    alias_crosswalk,
    by = "gps6_source_gene"
  )

  hierarchy <- dplyr::mutate(
    hierarchy,

    kinase_gene =
      dplyr::coalesce(
        .data$canonical_kinase_gene,
        .data$kinase_gene
      ),

    gps6_source_gene =
      clean_alias(
        .data$gps6_source_gene
      ),

    kl_matrix_name =
      clean_alias(
        .data$kl_matrix_name
      ),

    kl_display_name =
      clean_alias(
        .data$kl_display_name
      )
  )

  hierarchy <- dplyr::select(
    hierarchy,
    -dplyr::any_of(
      c(
        "canonical_kinase_gene"
      )
    )
  )

  # ============================================================================
  # Retain only unambiguous canonical gene -> family/superfamily assignments
  # ============================================================================

  hierarchy <- dplyr::group_by(
    hierarchy,
    .data$kinase_gene
  )

  hierarchy <- dplyr::filter(
    hierarchy,
    dplyr::n_distinct(
      .data$kinase_family
    ) == 1L,
    dplyr::n_distinct(
      .data$kinase_superfamily
    ) == 1L
  )

  hierarchy <- dplyr::ungroup(
    hierarchy
  )

  hierarchy <- dplyr::distinct(
    hierarchy,
    .data$kinase_gene,
    .keep_all = TRUE
  )

  dplyr::arrange(
    hierarchy,
    .data$kinase_superfamily,
    .data$kinase_family,
    .data$kinase_gene
  )
}


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
#'   kinase_superfamily = top-level GPS6 group (e.g. AGC)
#'   kinase_family      = second hierarchy level when available (e.g. AKT)
#'   kinase_gene    = deepest node only when the GPS6 node is sufficiently
#'                    specific to represent an individual kinase
#'
#' Final kinase-family standardization/crosswalking can then occur downstream.
#'
#' @param signor_data Processed PhosphoSIGNOR evidence.
#' @param kinase_library_data Processed Kinase Library evidence.
#' @param gps6_data Processed GPS 6.0 evidence from process_gps6_mappings().
#' @param gps6_hierarchy Canonical gene-to-family/superfamily crosswalk returned
#'   by build_gps6_kinase_hierarchy().
#'
#' @return A tibble containing harmonized evidence from all three sources.
#'
#' @export
harmonize_kinase_evidence <- function(
  signor_data,
  kinase_library_data,
  gps6_data,
  gps6_hierarchy
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
  # GPS6 hierarchy helpers
  # --------------------------------------------------------------------------

  clean_gps6_node <- function(x) {

    # Protect against factors, list-columns, and other non-character classes.
    if (is.list(x)) {

      x <- vapply(
        x,
        function(value) {

          if (
            length(value) == 0L ||
            all(is.na(value))
          ) {
            return(NA_character_)
          }

          paste(
            as.character(value),
            collapse = "/"
          )
        },
        character(1)
      )
    }

    x <- as.character(x)

    x <- stringr::str_trim(x)

    x <- stringr::str_replace_all(
      x,
      "^/+|/+$",
      ""
    )

    x[
      is.na(x) |
        !nzchar(x)
    ] <- NA_character_

    x
  }


  gps6_parts <- function(node) {

    node <- clean_gps6_node(node)

    # base::strsplit() requires a plain character vector.
    node_for_split <- as.character(node)

    node_for_split[
      is.na(node_for_split)
    ] <- ""

    parts <- strsplit(
      node_for_split,
      "/",
      fixed = TRUE
    )

    lapply(
      parts,
      function(x) {

        x[
          !is.na(x) &
            nzchar(x)
        ]
      }
    )
  }


  # --------------------------------------------------------------------------
  # Extract GPS6 hierarchy level
  # --------------------------------------------------------------------------

  gps6_level <- function(
    node,
    level
  ) {

    if (
      length(level) != 1L ||
      is.na(level) ||
      level < 1L
    ) {
      stop(
        "`level` must be one positive integer.",
        call. = FALSE
      )
    }

    level <- as.integer(level)

    parts <- gps6_parts(node)

    vapply(
      parts,
      function(x) {

        if (length(x) < level) {
          return(NA_character_)
        }

        as.character(
          x[[level]]
        )
      },
      character(1)
    )
  }


  # --------------------------------------------------------------------------
  # GPS6 deepest node
  # --------------------------------------------------------------------------

  gps6_leaf <- function(node) {

    parts <- gps6_parts(node)

    vapply(
      parts,
      function(x) {

        if (length(x) == 0L) {
          return(NA_character_)
        }

        as.character(
          x[[length(x)]]
        )
      },
      character(1)
    )
  }


  # --------------------------------------------------------------------------
  # GPS6 hierarchy depth
  # --------------------------------------------------------------------------

  gps6_depth <- function(node) {

    parts <- gps6_parts(node)

    vapply(
      parts,
      function(x) {

        if (length(x) == 0L) {
          return(NA_integer_)
        }

        as.integer(
          length(x)
        )
      },
      integer(1)
    )
  }


  # --------------------------------------------------------------------------
  # Canonical GPS6 / Kinase Library alias -> fine family -> superfamily lookup
  # --------------------------------------------------------------------------

  if (!is.data.frame(gps6_hierarchy)) {
    stop(
      "`gps6_hierarchy` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  required_hierarchy_columns <- c(
    "kinase_gene",
    "kinase_family",
    "kinase_superfamily",
    "gps6_source_gene"
  )

  missing_hierarchy_columns <- setdiff(
    required_hierarchy_columns,
    names(gps6_hierarchy)
  )

  if (length(missing_hierarchy_columns) > 0L) {
    stop(
      "`gps6_hierarchy` is missing required column(s): ",
      paste(
        missing_hierarchy_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  normalize_lookup_key <- function(x) {

    x <- stringr::str_to_upper(
      stringr::str_trim(
        as.character(x)
      )
    )

    x[
      is.na(x) |
        !nzchar(x)
    ] <- NA_character_

    x
  }

  make_lookup <- function(
    key_column,
    method
  ) {

    if (!key_column %in% names(gps6_hierarchy)) {
      return(
        tibble::tibble(
          lookup_key = character(),
          hierarchy_kinase_gene = character(),
          hierarchy_kinase_family = character(),
          hierarchy_kinase_superfamily = character(),
          hierarchy_mapping_method = character()
        )
      )
    }

    lookup <- dplyr::transmute(
      gps6_hierarchy,

      lookup_key =
        normalize_lookup_key(
          .data[[key_column]]
        ),

      hierarchy_kinase_gene =
        clean_symbol(
          .data$kinase_gene
        ),

      hierarchy_kinase_family =
        stringr::str_to_upper(
          stringr::str_trim(
            as.character(
              .data$kinase_family
            )
          )
        ),

      hierarchy_kinase_superfamily =
        stringr::str_to_upper(
          stringr::str_trim(
            as.character(
              .data$kinase_superfamily
            )
          )
        ),

      hierarchy_mapping_method =
        as.character(method)
    )

    lookup <- dplyr::filter(
      lookup,
      !is.na(.data$lookup_key),
      nzchar(.data$lookup_key)
    )

    # An alias is usable only when it maps unambiguously to one family and one
    # superfamily. This prevents broad/generic aliases from silently assigning
    # an incorrect fine family.
    lookup <- dplyr::group_by(
      lookup,
      .data$lookup_key
    )

    lookup <- dplyr::filter(
      lookup,
      dplyr::n_distinct(
        .data$hierarchy_kinase_family
      ) == 1L,
      dplyr::n_distinct(
        .data$hierarchy_kinase_superfamily
      ) == 1L
    )

    lookup <- dplyr::ungroup(
      lookup
    )

    dplyr::distinct(
      lookup,
      .data$lookup_key,
      .keep_all = TRUE
    )
  }

  gps6_gene_lookup <- make_lookup(
    "kinase_gene",
    "gps6_exact_gene"
  )

  gps6_source_alias_lookup <- make_lookup(
    "gps6_source_gene",
    "gps6_source_alias"
  )

  gps6_matrix_alias_lookup <- make_lookup(
    "kl_matrix_name",
    "gps6_kl_matrix_alias"
  )

  gps6_display_alias_lookup <- make_lookup(
    "kl_display_name",
    "gps6_kl_display_alias"
  )

  resolve_hierarchy_key <- function(
    keys,
    lookup
  ) {

    keys <- normalize_lookup_key(keys)

    match_index <- match(
      keys,
      lookup$lookup_key
    )

    tibble::tibble(
      resolved_gene =
        lookup$hierarchy_kinase_gene[
          match_index
        ],

      resolved_family =
        lookup$hierarchy_kinase_family[
          match_index
        ],

      resolved_superfamily =
        lookup$hierarchy_kinase_superfamily[
          match_index
        ],

      resolved_method =
        lookup$hierarchy_mapping_method[
          match_index
        ]
    )
  }

  apply_gps6_gene_hierarchy <- function(
    data,
    gene_column = "kinase_gene",
    fallback_column = NULL,
    matrix_alias_column = NULL,
    display_alias_column = NULL
  ) {

    n <- nrow(data)

    resolved_gene <- rep(
      NA_character_,
      n
    )

    resolved_family <- rep(
      NA_character_,
      n
    )

    resolved_superfamily <- rep(
      NA_character_,
      n
    )

    resolved_method <- rep(
      NA_character_,
      n
    )

    fill_from <- function(
      keys,
      lookup
    ) {

      candidate <- resolve_hierarchy_key(
        keys,
        lookup
      )

      missing <- is.na(resolved_family)

      resolved_gene[missing] <<-
        candidate$resolved_gene[missing]

      resolved_family[missing] <<-
        candidate$resolved_family[missing]

      resolved_superfamily[missing] <<-
        candidate$resolved_superfamily[missing]

      resolved_method[missing] <<-
        candidate$resolved_method[missing]
    }

    if (gene_column %in% names(data)) {
      fill_from(
        data[[gene_column]],
        gps6_gene_lookup
      )
    }

    if (
      !is.null(fallback_column) &&
      fallback_column %in% names(data)
    ) {
      fill_from(
        data[[fallback_column]],
        gps6_source_alias_lookup
      )
    }

    if (
      !is.null(matrix_alias_column) &&
      matrix_alias_column %in% names(data)
    ) {
      fill_from(
        data[[matrix_alias_column]],
        gps6_matrix_alias_lookup
      )

      # Matrix aliases are also common GPS6 source aliases, so try that lookup
      # as a secondary interpretation of the same value.
      fill_from(
        data[[matrix_alias_column]],
        gps6_source_alias_lookup
      )
    }

    if (
      !is.null(display_alias_column) &&
      display_alias_column %in% names(data)
    ) {
      fill_from(
        data[[display_alias_column]],
        gps6_display_alias_lookup
      )

      fill_from(
        data[[display_alias_column]],
        gps6_source_alias_lookup
      )
    }

    data$hierarchy_kinase_gene <-
      resolved_gene

    data$hierarchy_kinase_family <-
      resolved_family

    data$hierarchy_kinase_superfamily <-
      resolved_superfamily

    data$family_mapping_method <-
      dplyr::coalesce(
        resolved_method,
        "unresolved"
      )

    data$family_mapping_source <-
      dplyr::case_when(
        data$family_mapping_method == "gps6_exact_gene" ~
          "GPS6 hierarchy",

        data$family_mapping_method %in% c(
          "gps6_source_alias",
          "gps6_kl_matrix_alias",
          "gps6_kl_display_alias"
        ) ~
          "GPS6 hierarchy + Kinase Library alias metadata",

        TRUE ~
          "unresolved"
      )

    data
  }


  # --------------------------------------------------------------------------
  # Kinase Library kinome metadata crosswalk
  # --------------------------------------------------------------------------
  #
  # The installed kinase_library package exposes source-native kinase names
  # (e.g. CK1A, CK1D, BRAF) and provides its own kinome metadata through
  # get_kinase_info() / get_kinase_family().
  #
  # We query that metadata once per unique kinase name and use it to populate
  # canonical gene/family/group fields without hard-coding aliases in R.
  #
  # The helper is deliberately defensive because Kinase Library versions may
  # differ slightly in metadata field names.

  build_kinase_library_crosswalk <- function(
    kinase_names,
    name_type = "matrix"
  ) {

    kinase_names <- unique(
      stats::na.omit(
        as.character(kinase_names)
      )
    )

    kinase_names <- kinase_names[
      nzchar(
        stringr::str_trim(kinase_names)
      )
    ]

    if (length(kinase_names) == 0L) {

      return(
        tibble::tibble(
          source_native_kinase = character(),
          kl_canonical_gene = character(),
          kl_matrix_name = character(),
          kl_display_name = character(),
          kl_family = character(),
          kl_group = character(),
          kl_kinase_type = character()
        )
      )
    }

    py_path <- Sys.getenv(
      "RETICULATE_PYTHON",
      unset = "/root/.virtualenvs/pamchip-env/bin/python"
    )

    if (file.exists(py_path)) {

      reticulate::use_python(
        py_path,
        required = TRUE
      )

    } else {

      reticulate::use_virtualenv(
        "pamchip-env",
        required = TRUE
      )
    }

    kl <- reticulate::import(
      "kinase_library",
      convert = FALSE
    )

    first_non_missing <- function(x) {

      x <- as.character(x)

      x <- x[
        !is.na(x) &
          nzchar(
            stringr::str_trim(x)
          )
      ]

      if (length(x) == 0L) {
        return(NA_character_)
      }

      x[[1L]]
    }

    extract_info_field <- function(
      info,
      candidates
    ) {

      if (is.null(info)) {
        return(NA_character_)
      }

      info_r <- tryCatch(
        reticulate::py_to_r(info),
        error = function(e) NULL
      )

      if (is.null(info_r)) {
        return(NA_character_)
      }

      # pandas Series commonly becomes a named vector or one-row data.frame.
      if (is.data.frame(info_r)) {

        names_upper <- toupper(
          names(info_r)
        )

        for (candidate in candidates) {

          hit <- which(
            names_upper ==
              toupper(candidate)
          )

          if (length(hit) > 0L) {

            return(
              first_non_missing(
                info_r[[hit[[1L]]]]
              )
            )
          }
        }
      }

      if (!is.null(names(info_r))) {

        names_upper <- toupper(
          names(info_r)
        )

        for (candidate in candidates) {

          hit <- which(
            names_upper ==
              toupper(candidate)
          )

          if (length(hit) > 0L) {

            return(
              first_non_missing(
                info_r[[hit[[1L]]]]
              )
            )
          }
        }
      }

      NA_character_
    }

    rows <- lapply(
      kinase_names,
      function(kinase_name) {

        info <- tryCatch(
          kl$get_kinase_info(
            kinase_name,
            name_type = name_type
          ),
          error = function(e) NULL
        )

        family <- tryCatch(
          reticulate::py_to_r(
            kl$get_kinase_family(
              kinase_name,
              name_type = name_type
            )
          ),
          error = function(e) NA_character_
        )

        kin_type <- tryCatch(
          reticulate::py_to_r(
            kl$get_kinase_type(
              kinase_name,
              name_type = name_type
            )
          ),
          error = function(e) NA_character_
        )

        canonical_gene <- extract_info_field(
          info,
          c(
            "GENE",
            "GENE_NAME",
            "HGNC",
            "HGNC_SYMBOL",
            "SYMBOL",
            "GENE SYMBOL",
            "GENE NAME"
          )
        )

        matrix_name <- extract_info_field(
          info,
          c(
            "MATRIX_NAME",
            "KINASE"
          )
        )

        display_name <- extract_info_field(
          info,
          c(
            "DISPLAY_NAME"
          )
        )

        group <- extract_info_field(
          info,
          c(
            "GROUP",
            "KINASE_GROUP",
            "KINASE GROUP"
          )
        )

        # Some versions expose only FAMILY plus kinase type. Keep group missing
        # rather than fabricating a Manning group from ser_thr/tyrosine.
        tibble::tibble(
          source_native_kinase =
            as.character(kinase_name),

          kl_canonical_gene =
            clean_symbol(
              canonical_gene
            ),

          kl_matrix_name =
            clean_symbol(
              matrix_name
            ),

          kl_display_name =
            clean_symbol(
              display_name
            ),

          kl_family =
            {
              value <- as.character(family)

              if (
                length(value) == 0L ||
                is.na(value[[1L]]) ||
                !nzchar(
                  stringr::str_trim(
                    value[[1L]]
                  )
                )
              ) {
                NA_character_
              } else {
                value[[1L]]
              }
            },

          kl_group =
            {
              value <- as.character(group)

              if (
                length(value) == 0L ||
                is.na(value[[1L]]) ||
                !nzchar(
                  stringr::str_trim(
                    value[[1L]]
                  )
                )
              ) {
                NA_character_
              } else {
                value[[1L]]
              }
            },

          kl_kinase_type =
            {
              value <- as.character(kin_type)

              if (
                length(value) == 0L ||
                is.na(value[[1L]]) ||
                !nzchar(
                  stringr::str_trim(
                    value[[1L]]
                  )
                )
              ) {
                NA_character_
              } else {
                value[[1L]]
              }
            }
        )
      }
    )

    dplyr::bind_rows(
      rows
    )
  }


  # ============================================================================
  # 1. PhosphoSIGNOR
  # ============================================================================

  # process_signor_mappings() emits canonical substrate/site fields:
  #
  #   substrate_uniprot
  #   phosphosite
  #   res_position
  #   phosphoacceptor
  #
  # Keep backwards-compatible fallbacks for older cached targets, but prefer
  # those canonical fields whenever present.

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
      "substrate_id",
      "uniprot",
      "substrate_accession"
    )
  )

  signor_position_raw <- get_column(
    signor_data,
    c(
      "res_position",
      "peptide_site_position",
      "residue_position",
      "position",
      "RESIDUE"
    )
  )

  signor_site <- get_column(
    signor_data,
    c(
      "phosphosite",
      "site",
      "site_label",
      "RESIDUE"
    )
  )

  signor_residue <- get_column(
    signor_data,
    c(
      "phosphoacceptor",
      "phosphosite_residue",
      "residue"
    )
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

  signor_site <- as.character(
    signor_site
  )

  missing_site <- is.na(signor_site) |
    !nzchar(
      stringr::str_trim(signor_site)
    )

  signor_site[missing_site] <- ifelse(
    !is.na(signor_residue[missing_site]) &
      !is.na(
        clean_position(
          signor_position_raw[missing_site]
        )
      ),
    paste0(
      stringr::str_to_upper(
        as.character(
          signor_residue[missing_site]
        )
      ),
      clean_position(
        signor_position_raw[missing_site]
      )
    ),
    as.character(
      signor_position_raw[missing_site]
    )
  )

  signor_phosphoacceptor <- stringr::str_to_upper(
    as.character(signor_residue)
  )

  missing_acceptor <- is.na(signor_phosphoacceptor) |
    !signor_phosphoacceptor %in% c(
      "S",
      "T",
      "Y"
    )

  signor_phosphoacceptor[missing_acceptor] <-
    residue_from_site(
      signor_site[missing_acceptor]
    )

  # Support legacy SIGNOR labels such as "Ser768", "Thr308", and "Tyr204".
  legacy_acceptor <- dplyr::case_when(
    stringr::str_detect(
      signor_site,
      stringr::regex("^ser", ignore_case = TRUE)
    ) ~ "S",

    stringr::str_detect(
      signor_site,
      stringr::regex("^thr", ignore_case = TRUE)
    ) ~ "T",

    stringr::str_detect(
      signor_site,
      stringr::regex("^tyr", ignore_case = TRUE)
    ) ~ "Y",

    TRUE ~ NA_character_
  )

  signor_phosphoacceptor <- dplyr::coalesce(
    signor_phosphoacceptor,
    legacy_acceptor
  )

  signor_position <- clean_position(
    signor_position_raw
  )

  signor_clean <- tibble::tibble(

    peptide_id =
      as.character(signor_peptide_id),

    substrate_uniprot =
      clean_uniprot(signor_uniprot),

    phosphosite =
      dplyr::if_else(
        !is.na(signor_phosphoacceptor) &
          !is.na(signor_position),
        paste0(
          signor_phosphoacceptor,
          signor_position
        ),
        as.character(signor_site)
      ),

    peptide_site_position =
      signor_position,

    phosphoacceptor =
      signor_phosphoacceptor,

    kinase_gene =
      clean_symbol(signor_kinase),

    kinase_uniprot =
      clean_uniprot(signor_kinase_uniprot),

    kinase_family =
      as.character(signor_family),

    kinase_group =
      NA_character_,

    kinase_superfamily =
      NA_character_,

    source_native_kinase =
      as.character(signor_kinase),

    family_mapping_method =
      NA_character_,

    family_mapping_source =
      NA_character_,

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

  # Only retain S/T/Y site-specific records. Older cached SIGNOR targets may
  # still contain non-phosphorylation or non-site-specific rows; those records
  # cannot contribute to this STK/PTK substrate-mapping layer.
  signor_clean <- dplyr::filter(
    signor_clean,
    !is.na(.data$substrate_uniprot),
    !is.na(.data$peptide_site_position),
    .data$phosphoacceptor %in% c(
      "S",
      "T",
      "Y"
    ),
    !is.na(.data$kinase_gene),
    nzchar(.data$kinase_gene)
  )

  # --------------------------------------------------------------------------
  # Enrich SIGNOR kinase nomenclature through Kinase Library metadata
  # --------------------------------------------------------------------------

  signor_crosswalk_gene <- tryCatch(
    build_kinase_library_crosswalk(
      kinase_names = signor_clean$source_native_kinase,
      name_type = "gene"
    ),
    error = function(e) {
      tibble::tibble()
    }
  )

  signor_crosswalk_matrix <- tryCatch(
    build_kinase_library_crosswalk(
      kinase_names = signor_clean$source_native_kinase,
      name_type = "matrix"
    ),
    error = function(e) {
      tibble::tibble()
    }
  )

  signor_crosswalk <- dplyr::bind_rows(
    signor_crosswalk_gene,
    signor_crosswalk_matrix
  )

  if (nrow(signor_crosswalk) > 0L) {

    first_value <- function(x) {
      x <- stats::na.omit(
        as.character(x)
      )

      x <- x[
        nzchar(
          stringr::str_trim(x)
        )
      ]

      if (length(x) == 0L) {
        NA_character_
      } else {
        x[[1L]]
      }
    }

    signor_crosswalk <- dplyr::group_by(
      signor_crosswalk,
      .data$source_native_kinase
    )

    signor_crosswalk <- dplyr::summarise(
      signor_crosswalk,

      kl_canonical_gene =
        first_value(
          .data$kl_canonical_gene
        ),

      kl_matrix_name =
        first_value(
          .data$kl_matrix_name
        ),

      kl_display_name =
        first_value(
          .data$kl_display_name
        ),

      kl_family =
        first_value(
          .data$kl_family
        ),

      kl_group =
        first_value(
          .data$kl_group
        ),

      .groups = "drop"
    )

    signor_clean <- dplyr::left_join(
      signor_clean,
      signor_crosswalk,
      by = "source_native_kinase"
    )

    signor_clean <- dplyr::mutate(
      signor_clean,

      # Prefer a canonical KL gene when one exists. This is important for
      # SIGNOR aliases such as ERK1/2-style labels and source-native names.
      kinase_gene =
        dplyr::coalesce(
          .data$kl_canonical_gene,
          .data$kinase_gene
        ),

      kinase_group =
        dplyr::coalesce(
          dplyr::na_if(
            stringr::str_trim(
              as.character(
                .data$kinase_group
              )
            ),
            ""
          ),
          .data$kl_group,
          .data$kl_family
        )
    )

  } else {

    # Ensure the alias columns exist so the common hierarchy resolver below can
    # be called without branching.
    signor_clean$kl_matrix_name <- NA_character_
    signor_clean$kl_display_name <- NA_character_
  }

  # --------------------------------------------------------------------------
  # Map SIGNOR kinase identities onto the shared GPS6 fine-family hierarchy
  # --------------------------------------------------------------------------

  signor_clean <- apply_gps6_gene_hierarchy(
    signor_clean,
    gene_column = "kinase_gene",
    fallback_column = "source_native_kinase",
    matrix_alias_column = "kl_matrix_name",
    display_alias_column = "kl_display_name"
  )

  signor_clean <- dplyr::mutate(
    signor_clean,

    kinase_family =
      .data$hierarchy_kinase_family,

    kinase_superfamily =
      dplyr::coalesce(
        .data$hierarchy_kinase_superfamily,
        dplyr::na_if(
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(
                .data$kinase_group
              )
            )
          ),
          ""
        )
      ),

    kinase_group =
      .data$kinase_superfamily
  )

  # Remove temporary crosswalk columns only AFTER hierarchy resolution.
  signor_clean <- dplyr::select(
    signor_clean,
    -dplyr::any_of(
      c(
        "kl_canonical_gene",
        "kl_matrix_name",
        "kl_display_name",
        "kl_family",
        "kl_group",
        "hierarchy_kinase_gene",
        "hierarchy_kinase_family",
        "hierarchy_kinase_superfamily"
      )
    )
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

  kl_family_source <- get_column(
    kinase_library_data,
    c(
      "kinase_family",
      "family"
    )
  )

  kl_group_source <- get_column(
    kinase_library_data,
    c(
      "kinase_group",
      "group"
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

  # If a support flag was not explicitly provided, infer support only from the
  # fact that run_kinase_library_scoring() already returned a filtered table.
  # Do NOT impose a second percentile cutoff here.
  if (all(is.na(kl_supported))) {

    kl_supported <- rep(
      TRUE,
      nrow(kinase_library_data)
    )
  }

  # Build one metadata lookup per unique Kinase Library source-native kinase.
  kl_crosswalk <- build_kinase_library_crosswalk(
    kinase_names = kl_kinase,
    name_type = "matrix"
  )

  kl_base <- tibble::tibble(

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

    source_native_kinase =
      as.character(kl_kinase),

    source_kinase_gene =
      clean_symbol(kl_kinase),

    source_kinase_uniprot =
      clean_uniprot(kl_kinase_uniprot),

    source_kinase_family =
      as.character(kl_family_source),

    source_kinase_group =
      as.character(kl_group_source),

    kinase_library_percentile =
      suppressWarnings(
        as.numeric(kl_percentile)
      ),

    kinase_library_supported =
      as_support_flag(kl_supported)
  )

  kl_base <- dplyr::left_join(
    kl_base,
    kl_crosswalk,
    by = "source_native_kinase"
  )

  kl_base <- dplyr::mutate(
    kl_base,
    canonical_kinase_gene =
      dplyr::coalesce(
        .data$kl_canonical_gene,
        .data$source_kinase_gene
      )
  )

  kl_base <- apply_gps6_gene_hierarchy(
    kl_base,
    gene_column = "canonical_kinase_gene",
    fallback_column = "source_native_kinase",
    matrix_alias_column = "kl_matrix_name",
    display_alias_column = "kl_display_name"
  )

  kl_clean <- dplyr::transmute(
    kl_base,

    peptide_id =
      .data$peptide_id,

    substrate_uniprot =
      .data$substrate_uniprot,

    phosphosite =
      .data$phosphosite,

    peptide_site_position =
      .data$peptide_site_position,

    phosphoacceptor =
      .data$phosphoacceptor,

    # Prefer Kinase Library's canonical gene name, then map that gene into the
    # finer GPS6 family hierarchy. This converts broad KL labels such as AGC
    # into family-level labels such as AKT when the gene can be resolved.
    kinase_gene =
      .data$canonical_kinase_gene,

    kinase_uniprot =
      .data$source_kinase_uniprot,

    # Fine family comes only from the shared GPS6 hierarchy. The Kinase
    # Library FAMILY field is retained at the broader superfamily level below.
    kinase_family =
      .data$hierarchy_kinase_family,

    kinase_superfamily =
      dplyr::coalesce(
        .data$hierarchy_kinase_superfamily,
        dplyr::na_if(
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(.data$kl_family)
            )
          ),
          ""
        ),
        dplyr::na_if(
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(.data$kl_group)
            )
          ),
          ""
        ),
        dplyr::na_if(
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(.data$source_kinase_group)
            )
          ),
          ""
        )
      ),

    # Backward-compatible alias for broad hierarchy level.
    kinase_group =
      .data$kinase_superfamily,

    source_native_kinase =
      .data$source_native_kinase,

    family_mapping_method =
      .data$family_mapping_method,

    family_mapping_source =
      .data$family_mapping_source,

    evidence_source =
      "Kinase Library",

    evidence_class =
      "predictive",

    source_score =
      .data$kinase_library_percentile,

    source_cutoff =
      NA_real_,

    source_supported =
      .data$kinase_library_supported,

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
      .data$kinase_library_percentile,

    kinase_library_supported =
      .data$kinase_library_supported,

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
      stringr::str_to_upper(
        as.character(gps6_family)
      ),

    kinase_superfamily =
      stringr::str_to_upper(
        as.character(gps6_root)
      ),

    # Backward-compatible alias for broad hierarchy level.
    kinase_group =
      stringr::str_to_upper(
        as.character(gps6_root)
      ),

    source_native_kinase =
      as.character(gps6_leaf_value),

    family_mapping_method =
      "gps6_native_hierarchy",

    family_mapping_source =
      "GPS6 hierarchy",

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

  harmonized <- dplyr::mutate(
    harmonized,

    peptide_id =
      as.character(.data$peptide_id),

    substrate_uniprot =
      clean_uniprot(.data$substrate_uniprot),

    kinase_gene =
      clean_symbol(.data$kinase_gene),

    kinase_uniprot =
      clean_uniprot(.data$kinase_uniprot),

    kinase_superfamily =
      dplyr::na_if(
        stringr::str_to_upper(
          stringr::str_trim(
            as.character(.data$kinase_superfamily)
          )
        ),
        ""
      ),

    # Keep kinase_group as a compatibility alias for kinase_superfamily.
    kinase_group =
      dplyr::coalesce(
        .data$kinase_superfamily,
        dplyr::na_if(
          stringr::str_to_upper(
            stringr::str_trim(
              as.character(.data$kinase_group)
            )
          ),
          ""
        )
      ),

    kinase_family =
      dplyr::na_if(
        stringr::str_to_upper(
          stringr::str_trim(
            as.character(.data$kinase_family)
          )
        ),
        ""
      ),

    evidence_source =
      as.character(.data$evidence_source),

    evidence_class =
      as.character(.data$evidence_class),

    family_mapping_method =
      as.character(.data$family_mapping_method),

    family_mapping_source =
      as.character(.data$family_mapping_source)
  )

  # Remove only exact duplicate evidence rows.
  #
  # We intentionally DO NOT collapse evidence from different GPS6 hierarchy
  # levels or different sources here.
  harmonized <- dplyr::distinct(
    harmonized
  )

  harmonized <- dplyr::arrange(
    harmonized,
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
  kl_data,
  gps6_hierarchy
) {

  harmonize_kinase_evidence(
    signor_data = signor_data,
    kinase_library_data = kl_data,
    gps6_data = gps6_data,
    gps6_hierarchy = gps6_hierarchy
  )
}