#' Fetch SIGNOR Kinase-Substrate Interactions
#'
#' Downloads the complete SIGNOR dataset directly via bulk release or API fallback.
#'
#' @param cache_dir Directory to cache raw downloaded files
#' @param force_update Logical, if TRUE re-downloads data even if cached file exists
#' @return A tibble of raw SIGNOR interactions
fetch_signor_data <- function(cache_dir = "data/external/signor", force_update = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, "signor_raw.rds")

  # Standard 29-column schema definitions for SIGNOR downloads
  signor_standard_cols <- c(
    "ENTITYA", "TYPEA", "IDA", "DATABASEA",
    "ENTITYB", "TYPEB", "IDB", "DATABASEB",
    "EFFECT", "MECHANISM", "RESIDUE", "SEQUENCE",
    "TAX_ID", "CELL_DATA", "TISSUE_DATA", "MODULATOR",
    "METHOD", "DIRECT", "NOTES", "ANNOTATOR",
    "PMID", "DIRECT_BOOL", "SENTENCE", "SCORE",
    "SIGNOR_ID", "CONFIDENCE", "EXTRA_1", "EXTRA_2", "EXTRA_3"
  )

  signor_data <- NULL

  # Load from cache if present and valid
  if (file.exists(cache_file) && !force_update) {
    message("Loading cached SIGNOR data from: ", cache_file)
    signor_data <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  }

  cache_invalid <- is.null(signor_data) || 
    nrow(signor_data) == 0 || 
    any(startsWith(colnames(signor_data), "...")) ||
    !("ENTITYA" %in% colnames(signor_data) || "IDA" %in% colnames(signor_data))

  if (cache_invalid) {
    message("Fetching fresh SIGNOR interaction data from source...")

    primary_url <- "https://signor.uniroma2.it/getData.php?organism=9606"
    fallback_url <- "https://signor.uniroma2.it/downloads/SIGNOR_all_data_signed.tsv"

    signor_data <- tryCatch({
      req <- httr2::request(primary_url) %>%
        httr2::req_user_agent("Mozilla/5.0 (PamChipAnnotationPipeline/1.0; R/httr2)") %>%
        httr2::req_timeout(30)
      
      resp <- httr2::req_perform(req)
      raw_text <- httr2::resp_body_string(resp)

      if (nchar(trimws(raw_text)) > 0 && !startsWith(trimws(raw_text), "{")) {
        # Wrapped in I() to adhere to readr >= 2.2.0 literal string requirements
        readr::read_tsv(
          I(raw_text),
          col_names = TRUE,
          show_col_types = FALSE
        )
      } else {
        jsonlite::fromJSON(raw_text) %>% tibble::as_tibble()
      }
    }, error = function(e) {
      warning("Primary SIGNOR fetch failed: ", e$message, ". Attempting fallback...")
      tryCatch({
        req_fb <- httr2::request(fallback_url) %>%
          httr2::req_user_agent("Mozilla/5.0 (PamChipAnnotationPipeline/1.0)") %>%
          httr2::req_timeout(30)
        resp_fb <- httr2::req_perform(req_fb)
        readr::read_tsv(
          I(httr2::resp_body_string(resp_fb)),
          col_names = TRUE,
          show_col_types = FALSE
        )
      }, error = function(e_fb) {
        stop("All SIGNOR download attempts failed: ", e_fb$message)
      })
    })
  }

  # Guard against 0-row results
  if (nrow(signor_data) == 0) {
    stop("Downloaded SIGNOR dataset contains 0 rows.")
  }

  # Repair headers if auto-generated or shifted into row 1
  if (any(startsWith(colnames(signor_data), "...")) || !"ENTITYA" %in% colnames(signor_data)) {
    first_row_vals <- unlist(signor_data[1, ], use.names = FALSE)
    if (any(c("protein", "UNIPROT", "ENTITYA") %in% first_row_vals)) {
      signor_data <- signor_data[-1, ]
    }
    n_cols <- min(ncol(signor_data), length(signor_standard_cols))
    colnames(signor_data)[1:n_cols] <- signor_standard_cols[1:n_cols]
  }

  saveRDS(signor_data, cache_file)
  return(signor_data)
}


#' Process and Filter SIGNOR Mappings for PamChip Target Sites
#'
#' @param signor_raw Raw data frame fetched from PhosphoSIGNOR
#' @param target_sites Processed PamChip sites target tibble
#' @return Cleaned tibble of curated kinase-substrate interactions
process_signor_mappings <- function(signor_raw, target_sites) {

  empty_result <- function() {
    tibble::tibble(
      peptide_id = character(),
      id = character(),
      substrate_uniprot = character(),
      uniprot_id = character(),
      phosphosite = character(),
      res_position = integer(),
      phosphoacceptor = character(),
      kinase_symbol = character(),
      signor_pmid = character(),
      signor_id = character(),
      mechanism = character(),
      signor_sequence = character(),
      score = numeric(),
      source = character()
    )
  }

  if (
    is.null(signor_raw) ||
    length(colnames(signor_raw)) == 0L ||
    nrow(signor_raw) == 0L
  ) {
    warning(
      "signor_raw is empty or contains no columns. Returning empty mapping table.",
      call. = FALSE
    )
    return(empty_result())
  }

  if (
    is.null(target_sites) ||
    nrow(target_sites) == 0L
  ) {
    warning(
      "target_sites is empty. Returning empty SIGNOR mapping table.",
      call. = FALSE
    )
    return(empty_result())
  }

  cols <- colnames(signor_raw)

  first_existing_name <- function(candidates) {
    hit <- candidates[candidates %in% cols]
    if (length(hit) == 0L) {
      return(NA_character_)
    }
    hit[[1L]]
  }

  substrate_col <- first_existing_name(
    c(
      "substrate_id",
      "IDB",
      "ENTITYB_ID",
      "SUBSTRATE_ID",
      "substrate_uniprot",
      "uniprot_id",
      "UNIPROT_ID"
    )
  )

  kinase_col <- first_existing_name(
    c(
      "kinase_symbol",
      "ENTITYA",
      "kinase_id",
      "KINASE",
      "IDA"
    )
  )

  residue_col <- first_existing_name(
    c(
      "RESIDUE",
      "residue",
      "phosphosite",
      "site"
    )
  )

  mechanism_col <- first_existing_name(
    c(
      "MECHANISM",
      "mechanism"
    )
  )

  sequence_col <- first_existing_name(
    c(
      "SEQUENCE",
      "sequence",
      "site_sequence"
    )
  )

  pmid_col <- first_existing_name(
    c(
      "PMID",
      "pmid",
      "signor_pmid"
    )
  )

  signor_id_col <- first_existing_name(
    c(
      "SIGNOR_ID",
      "signor_id"
    )
  )

  if (is.na(substrate_col)) {
    stop(
      "Could not locate a SIGNOR substrate UniProt column. Available columns: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(kinase_col)) {
    stop(
      "Could not locate a SIGNOR upstream-entity column. Available columns: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(residue_col)) {
    stop(
      "Could not locate a SIGNOR residue column. Available columns: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.na(mechanism_col)) {
    stop(
      "Could not locate a SIGNOR mechanism column. Available columns: ",
      paste(cols, collapse = ", "),
      call. = FALSE
    )
  }

  clean_uniprot <- function(x) {
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

  normalize_residue_name <- function(x) {
    x <- stringr::str_trim(as.character(x))

    residue <- dplyr::case_when(
      stringr::str_detect(
        x,
        stringr::regex("^ser", ignore_case = TRUE)
      ) ~ "S",

      stringr::str_detect(
        x,
        stringr::regex("^thr", ignore_case = TRUE)
      ) ~ "T",

      stringr::str_detect(
        x,
        stringr::regex("^tyr", ignore_case = TRUE)
      ) ~ "Y",

      stringr::str_detect(
        stringr::str_to_upper(x),
        "^[STY][0-9]+"
      ) ~ stringr::str_sub(
        stringr::str_to_upper(x),
        1L,
        1L
      ),

      TRUE ~ NA_character_
    )

    residue
  }

  extract_position <- function(x) {
    suppressWarnings(
      as.integer(
        stringr::str_extract(
          as.character(x),
          "[0-9]+"
        )
      )
    )
  }

  raw_substrate <- signor_raw[[substrate_col]]
  raw_kinase <- signor_raw[[kinase_col]]
  raw_residue <- signor_raw[[residue_col]]
  raw_mechanism <- signor_raw[[mechanism_col]]

  raw_sequence <- if (!is.na(sequence_col)) {
    signor_raw[[sequence_col]]
  } else {
    rep(NA_character_, nrow(signor_raw))
  }

  raw_pmid <- if (!is.na(pmid_col)) {
    signor_raw[[pmid_col]]
  } else {
    rep(NA_character_, nrow(signor_raw))
  }

  raw_signor_id <- if (!is.na(signor_id_col)) {
    signor_raw[[signor_id_col]]
  } else {
    rep(NA_character_, nrow(signor_raw))
  }

  signor_df <- tibble::tibble(
    peptide_id = NA_character_,
    id = NA_character_,

    substrate_uniprot =
      clean_uniprot(raw_substrate),

    uniprot_id =
      clean_uniprot(raw_substrate),

    residue_raw =
      as.character(raw_residue),

    phosphoacceptor =
      normalize_residue_name(raw_residue),

    res_position =
      extract_position(raw_residue),

    phosphosite =
      dplyr::if_else(
        !is.na(normalize_residue_name(raw_residue)) &
          !is.na(extract_position(raw_residue)),
        paste0(
          normalize_residue_name(raw_residue),
          extract_position(raw_residue)
        ),
        NA_character_
      ),

    kinase_symbol =
      stringr::str_trim(
        as.character(raw_kinase)
      ),

    mechanism =
      stringr::str_to_lower(
        stringr::str_trim(
          as.character(raw_mechanism)
        )
      ),

    signor_sequence =
      as.character(raw_sequence),

    signor_pmid =
      as.character(raw_pmid),

    signor_id =
      as.character(raw_signor_id),

    source =
      "SIGNOR",

    score =
      1.0
  )

  # --------------------------------------------------------------------------
  # Restrict SIGNOR to the evidence class this pipeline actually models:
  # phosphorylation of serine/threonine/tyrosine residues.
  #
  # Do not use DIRECT / DIRECT_BOOL here. The currently observed SIGNOR export
  # contains publication-like values in that field, so it is not safe as a
  # boolean filter without separately validating the upstream download schema.
  # --------------------------------------------------------------------------

  signor_df <- dplyr::filter(
    signor_df,
    .data$mechanism == "phosphorylation",
    !is.na(.data$substrate_uniprot),
    !is.na(.data$res_position),
    .data$phosphoacceptor %in% c("S", "T", "Y"),
    !is.na(.data$kinase_symbol),
    nzchar(.data$kinase_symbol)
  )

  # --------------------------------------------------------------------------
  # Restrict to proteins represented on the current PamChip.
  #
  # Exact phosphosite-to-peptide linking is intentionally deferred to the
  # master-evidence linker. This preserves all curated phosphorylation evidence
  # on array proteins while keeping the site join explicit and auditable.
  # --------------------------------------------------------------------------

  target_uniprots <- target_sites %>%
    dplyr::transmute(
      substrate_uniprot =
        clean_uniprot(
          dplyr::coalesce(
            as.character(.data$substrate_uniprot),
            as.character(.data$uniprot_id)
          )
        )
    ) %>%
    dplyr::filter(
      !is.na(.data$substrate_uniprot)
    ) %>%
    dplyr::distinct() %>%
    dplyr::pull(.data$substrate_uniprot)

  signor_df <- dplyr::filter(
    signor_df,
    .data$substrate_uniprot %in% target_uniprots
  )

  # Keep exact curated records. Multiple kinases can legitimately regulate the
  # same substrate/site, and the same site can have multiple publications.
  signor_df <- dplyr::distinct(
    signor_df
  )

  if (nrow(signor_df) == 0L) {
    warning(
      "No SIGNOR S/T/Y phosphorylation records matched PamChip substrate proteins.",
      call. = FALSE
    )
  }

  signor_df
}

