# R/kinase_taxonomy.R

#' Default KinHub / OpenKinome kinase taxonomy URL
#'
#' The OpenKinome KinHub kinase list is a machine-readable human kinase
#' reference containing HGNC gene names, UniProt accessions, Manning-style
#' groups, families, and subfamilies.
KINHUB_TAXONOMY_URL <- paste0(
  "https://raw.githubusercontent.com/",
  "openkinome/kinodata/master/data/KinHubKinaseList.csv"
)


#' Download and cache the KinHub human kinase taxonomy
#'
#' @param cache_dir Directory in which to cache the downloaded CSV.
#' @param force_update Logical; if TRUE, download even when a cached file exists.
#' @param url Source URL. Defaults to the OpenKinome KinHub CSV.
#' @param filename Cached filename.
#'
#' @return Path to the cached CSV file.
fetch_kinhub_taxonomy <- function(
  cache_dir = file.path("data", "external", "kinase_taxonomy"),
  force_update = FALSE,
  url = KINHUB_TAXONOMY_URL,
  filename = "KinHubKinaseList.csv"
) {

  dir.create(
    cache_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  out <- file.path(
    cache_dir,
    filename
  )

  if (
    file.exists(out) &&
    !isTRUE(force_update)
  ) {
    message(
      "Loading cached KinHub kinase taxonomy from: ",
      out
    )

    return(out)
  }

  message(
    "Downloading KinHub kinase taxonomy from: ",
    url
  )

  response <- tryCatch(
    {
      request <- httr2::request(url) |>
        httr2::req_user_agent(
          "Pamchip_Annotation kinase taxonomy downloader"
        ) |>
        httr2::req_retry(
          max_tries = 5L
        ) |>
        httr2::req_timeout(
          seconds = 60
        )

      httr2::req_perform(request)
    },
    error = function(e) {

      if (file.exists(out)) {
        warning(
          "KinHub taxonomy download failed; using cached file instead: ",
          conditionMessage(e),
          call. = FALSE
        )

        return(NULL)
      }

      stop(
        "Failed to download KinHub kinase taxonomy and no cached copy exists: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.null(response)) {
    return(out)
  }

  if (httr2::resp_status(response) < 200L ||
      httr2::resp_status(response) >= 300L) {
    stop(
      "KinHub taxonomy download returned HTTP status ",
      httr2::resp_status(response),
      ".",
      call. = FALSE
    )
  }

  writeBin(
    httr2::resp_body_raw(response),
    out
  )

  if (
    !file.exists(out) ||
    file.info(out)$size <= 0L
  ) {
    stop(
      "KinHub taxonomy download did not create a valid file: ",
      out,
      call. = FALSE
    )
  }

  out
}


#' Curated gene-level taxonomy overrides
#'
#' Explicitly fills known gaps in the KinHub/OpenKinome human kinase table.
#' These rows are intentionally small, transparent, and version-controlled.
#' They take precedence over KinHub if the same canonical gene is ever added
#' upstream, so changes are reviewable rather than silently altering mappings.
#'
#' @return Tibble with the same canonical taxonomy columns used downstream.
kinase_taxonomy_overrides <- function() {

  tibble::tribble(
    ~kinase_gene, ~kinase_uniprot, ~kinase_superfamily, ~kinase_family,
    ~kinase_subfamily, ~kinhub_xname, ~manning_name, ~kinhub_kinase_name,
    ~taxonomy_source, ~taxonomy_source_url,

    "MAP3K21", "Q5TCX8", "TKL", "MLK",
    NA_character_, NA_character_, "MLK4", "Mitogen-activated protein kinase kinase kinase 21",
    "Pamchip_Annotation curated override", NA_character_,

    "GRK2", "P25098", "AGC", "GRK",
    NA_character_, NA_character_, "GRK2", "G protein-coupled receptor kinase 2",
    "Pamchip_Annotation curated override", NA_character_,

    "GRK3", "P35626", "AGC", "GRK",
    NA_character_, NA_character_, "GRK3", "G protein-coupled receptor kinase 3",
    "Pamchip_Annotation curated override", NA_character_
  )
}


#' Standardize the KinHub human kinase taxonomy
#'
#' @param kinhub_file Path to KinHubKinaseList.csv.
#'
#' @return Tibble with one canonical row per kinase gene and columns:
#'   kinase_gene, kinase_uniprot, kinase_family, kinase_superfamily,
#'   kinase_subfamily, and source-native identifiers.
build_kinase_taxonomy <- function(
  kinhub_file
) {

  if (
    length(kinhub_file) != 1L ||
    !file.exists(kinhub_file)
  ) {
    stop(
      "`kinhub_file` must be exactly one valid file path.",
      call. = FALSE
    )
  }

  raw <- readr::read_csv(
    kinhub_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  # KinHub currently uses non-breaking spaces (U+00A0) in some CSV headers
  # such as "Manning Name" and "HGNC Name". Normalize all header whitespace
  # before validating or selecting columns so the parser is robust to those
  # Unicode formatting differences.
  normalize_header <- function(x) {

    x <- as.character(x)

    # Convert non-breaking spaces and other repeated whitespace to plain spaces.
    x <- stringr::str_replace_all(
      x,
      "\\u00A0",
      " "
    )

    x <- stringr::str_squish(x)

    x
  }

  names(raw) <- normalize_header(
    names(raw)
  )

  required <- c(
    "Manning Name",
    "HGNC Name",
    "Group",
    "Family",
    "SubFamily",
    "UniprotID"
  )

  missing <- setdiff(
    required,
    names(raw)
  )

  if (length(missing) > 0L) {
    stop(
      "KinHub taxonomy is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  clean_chr <- function(x) {

    x <- stringr::str_trim(
      as.character(x)
    )

    x[
      is.na(x) |
        x == "" |
        toupper(x) %in% c(
          "NA",
          "N/A",
          "NONE",
          "NULL"
        )
    ] <- NA_character_

    x
  }

  clean_upper <- function(x) {
    x <- clean_chr(x)
    stringr::str_to_upper(x)
  }

  taxonomy <- raw |>
    dplyr::transmute(

      kinase_gene =
        clean_upper(
          .data[["HGNC Name"]]
        ),

      kinase_uniprot =
        clean_upper(
          .data[["UniprotID"]]
        ),

      kinase_superfamily =
        clean_upper(
          .data[["Group"]]
        ),

      kinase_family =
        clean_upper(
          .data[["Family"]]
        ),

      kinase_subfamily =
        clean_upper(
          .data[["SubFamily"]]
        ),

      kinhub_xname =
        clean_chr(
          if ("xName" %in% names(raw)) {
            .data[["xName"]]
          } else {
            NA_character_
          }
        ),

      manning_name =
        clean_chr(
          .data[["Manning Name"]]
        ),

      kinhub_kinase_name =
        clean_chr(
          if ("Kinase Name" %in% names(raw)) {
            .data[["Kinase Name"]]
          } else {
            NA_character_
          }
        ),

      taxonomy_source =
        "KinHub/OpenKinome",

      taxonomy_source_url =
        KINHUB_TAXONOMY_URL
    ) |>
    dplyr::filter(
      !is.na(.data$kinase_gene),
      nzchar(.data$kinase_gene)
    ) |>
    dplyr::distinct()

  # --------------------------------------------------------------------------
  # KinHub contains explicit second-domain rows for multi-domain kinases.
  #
  # Examples in the current KinHub table include:
  #
  #   JAK1_b   / Domain2_JAK1
  #   RSK1_b   / Domain2_RSK3
  #   GCN2_b   / Domain2_GCN2
  #
  # These are kinase-domain classifications, not contradictory gene-level
  # annotations. For the canonical gene -> family taxonomy used by this
  # pipeline, retain the primary kinase-domain row and exclude explicit
  # Domain2 / *_b rows from the gene-level collapse.
  #
  # Secondary-domain rows are still available in the cached raw KinHub CSV
  # for provenance and future domain-level analyses.
  # --------------------------------------------------------------------------

  taxonomy <- taxonomy |>
    dplyr::mutate(
      is_secondary_domain =
        stringr::str_detect(
          dplyr::coalesce(
            .data$kinhub_xname,
            ""
          ),
          stringr::regex("_b$", ignore_case = TRUE)
        ) |
        stringr::str_detect(
          dplyr::coalesce(
            .data$manning_name,
            ""
          ),
          stringr::regex("^Domain2_", ignore_case = TRUE)
        )
    )

  primary_taxonomy <- taxonomy |>
    dplyr::filter(
      !.data$is_secondary_domain
    )

  # After explicit secondary-domain rows are removed, canonical human genes
  # should have a single family / superfamily assignment. If conflicts remain,
  # fail loudly because those would represent a genuinely ambiguous taxonomy
  # rather than a known multi-domain representation.
  conflicts <- primary_taxonomy |>
    dplyr::group_by(
      .data$kinase_gene
    ) |>
    dplyr::summarise(
      n_family =
        dplyr::n_distinct(
          .data$kinase_family[
            !is.na(.data$kinase_family)
          ]
        ),

      n_superfamily =
        dplyr::n_distinct(
          .data$kinase_superfamily[
            !is.na(.data$kinase_superfamily)
          ]
        ),

      .groups = "drop"
    ) |>
    dplyr::filter(
      .data$n_family > 1L |
        .data$n_superfamily > 1L
    )

  if (nrow(conflicts) > 0L) {
    stop(
      "KinHub contains conflicting PRIMARY-domain family/group assignments for ",
      nrow(conflicts),
      " kinase gene(s): ",
      paste(
        utils::head(
          conflicts$kinase_gene,
          20L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  taxonomy <- primary_taxonomy |>
    dplyr::group_by(
      .data$kinase_gene
    ) |>
    dplyr::summarise(

      kinase_uniprot =
        first_non_missing_taxonomy(
          .data$kinase_uniprot
        ),

      kinase_superfamily =
        first_non_missing_taxonomy(
          .data$kinase_superfamily
        ),

      kinase_family =
        first_non_missing_taxonomy(
          .data$kinase_family
        ),

      kinase_subfamily =
        first_non_missing_taxonomy(
          .data$kinase_subfamily
        ),

      kinhub_xname =
        first_non_missing_taxonomy(
          .data$kinhub_xname
        ),

      manning_name =
        first_non_missing_taxonomy(
          .data$manning_name
        ),

      kinhub_kinase_name =
        first_non_missing_taxonomy(
          .data$kinhub_kinase_name
        ),

      taxonomy_source =
        first_non_missing_taxonomy(
          .data$taxonomy_source
        ),

      taxonomy_source_url =
        first_non_missing_taxonomy(
          .data$taxonomy_source_url
        ),

      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$kinase_superfamily,
      .data$kinase_family,
      .data$kinase_gene
    )

  # --------------------------------------------------------------------------
  # Apply explicit curated coverage overrides.
  #
  # Override rows take precedence over KinHub for the same canonical gene.
  # This currently fills three known Kinase Library coverage gaps:
  #
  #   MAP3K21 -> MLK -> TKL
  #   GRK2    -> GRK -> AGC
  #   GRK3    -> GRK -> AGC
  #
  # Keeping these exceptions here makes them explicit and auditable rather
  # than burying them in source-specific harmonization code.
  # --------------------------------------------------------------------------

  overrides <- kinase_taxonomy_overrides() |>
    dplyr::mutate(
      dplyr::across(
        c(
          "kinase_gene",
          "kinase_uniprot",
          "kinase_superfamily",
          "kinase_family",
          "kinase_subfamily"
        ),
        ~ {
          x <- as.character(.x)
          x[!is.na(x)] <- stringr::str_to_upper(
            stringr::str_trim(x[!is.na(x)])
          )
          x
        }
      )
    )

  taxonomy <- taxonomy |>
    dplyr::filter(
      !.data$kinase_gene %in% overrides$kinase_gene
    ) |>
    dplyr::bind_rows(
      overrides
    ) |>
    dplyr::arrange(
      .data$kinase_superfamily,
      .data$kinase_family,
      .data$kinase_gene
    )

  # Final gene-level invariant: exactly one row per canonical kinase gene.
  duplicate_genes <- taxonomy |>
    dplyr::count(
      .data$kinase_gene,
      name = "n"
    ) |>
    dplyr::filter(
      .data$n != 1L
    )

  if (nrow(duplicate_genes) > 0L) {
    stop(
      "Final kinase taxonomy does not contain exactly one row per gene for: ",
      paste(
        duplicate_genes$kinase_gene,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  taxonomy
}


# Internal helper used by build_kinase_taxonomy().
first_non_missing_taxonomy <- function(x) {

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


#' Merge the external human kinase taxonomy with the GPS6 alias hierarchy
#'
#' This produces a hierarchy table that is directly compatible with the
#' existing harmonize_kinase_evidence() machinery.
#'
#' External KinHub taxonomy is authoritative for canonical human
#' gene -> family -> superfamily classification. GPS6 rows are retained because
#' they contain source-native aliases required to resolve GPS6 model names and
#' legacy kinase aliases.
#'
#' @param kinase_taxonomy Standardized external taxonomy from
#'   build_kinase_taxonomy().
#' @param gps6_hierarchy Existing GPS6 hierarchy from
#'   build_gps6_kinase_hierarchy().
#'
#' @return Combined hierarchy tibble.
combine_kinase_taxonomies <- function(
  kinase_taxonomy,
  gps6_hierarchy
) {

  required_external <- c(
    "kinase_gene",
    "kinase_family",
    "kinase_superfamily"
  )

  missing_external <- setdiff(
    required_external,
    names(kinase_taxonomy)
  )

  if (length(missing_external) > 0L) {
    stop(
      "`kinase_taxonomy` is missing required column(s): ",
      paste(missing_external, collapse = ", "),
      call. = FALSE
    )
  }

  required_gps6 <- c(
    "kinase_gene",
    "kinase_family",
    "kinase_superfamily",
    "gps6_source_gene"
  )

  missing_gps6 <- setdiff(
    required_gps6,
    names(gps6_hierarchy)
  )

  if (length(missing_gps6) > 0L) {
    stop(
      "`gps6_hierarchy` is missing required column(s): ",
      paste(missing_gps6, collapse = ", "),
      call. = FALSE
    )
  }

  external_rows <- kinase_taxonomy |>
    dplyr::transmute(

      kinase_type =
        NA_character_,

      gps6_hierarchy =
        NA_character_,

      # Setting the source gene to the canonical gene lets the existing
      # source-alias lookup resolve a canonical external taxonomy row without
      # requiring changes to harmonize.R.
      gps6_source_gene =
        as.character(.data$kinase_gene),

      kinase_gene =
        as.character(.data$kinase_gene),

      kinase_family =
        as.character(.data$kinase_family),

      kinase_superfamily =
        as.character(.data$kinase_superfamily),

      canonical_kinase_uniprot =
        as.character(.data$kinase_uniprot),

      kl_matrix_name =
        NA_character_,

      kl_display_name =
        NA_character_,

      hierarchy_source =
        as.character(.data$taxonomy_source)
    )

  gps6_rows <- gps6_hierarchy

  for (column in c(
    "canonical_kinase_uniprot",
    "kl_matrix_name",
    "kl_display_name"
  )) {
    if (!column %in% names(gps6_rows)) {
      gps6_rows[[column]] <- NA_character_
    }
  }

  gps6_rows <- gps6_rows |>
    dplyr::mutate(
      hierarchy_source =
        "GPS6"
    ) |>
    dplyr::select(
      dplyr::all_of(
        names(external_rows)
      )
    )

  combined <- dplyr::bind_rows(
    external_rows,
    gps6_rows
  )

  # --------------------------------------------------------------------------
  # Enforce external taxonomy for canonical genes.
  #
  # GPS6 remains authoritative only for source-native aliases. For a canonical
  # human gene that exists in KinHub, family/group values are overwritten with
  # the external taxonomy assignment.
  # --------------------------------------------------------------------------

  canonical_lookup <- kinase_taxonomy |>
    dplyr::select(
      "kinase_gene",
      external_family = "kinase_family",
      external_superfamily = "kinase_superfamily",
      external_uniprot = "kinase_uniprot"
    )

  combined <- combined |>
    dplyr::left_join(
      canonical_lookup,
      by = "kinase_gene"
    ) |>
    dplyr::mutate(

      kinase_family =
        dplyr::coalesce(
          .data$external_family,
          .data$kinase_family
        ),

      kinase_superfamily =
        dplyr::coalesce(
          .data$external_superfamily,
          .data$kinase_superfamily
        ),

      canonical_kinase_uniprot =
        dplyr::coalesce(
          .data$external_uniprot,
          .data$canonical_kinase_uniprot
        )
    ) |>
    dplyr::select(
      -"external_family",
      -"external_superfamily",
      -"external_uniprot"
    )

  # Preserve multiple GPS6 aliases for a canonical gene, but eliminate exact
  # duplicate hierarchy rows.
  combined |>
    dplyr::distinct() |>
    dplyr::arrange(
      .data$kinase_superfamily,
      .data$kinase_family,
      .data$kinase_gene,
      .data$gps6_source_gene
    )
}
