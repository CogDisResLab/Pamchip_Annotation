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
  
  if (is.null(signor_raw) || length(colnames(signor_raw)) == 0 || nrow(signor_raw) == 0) {
    warning("signor_raw is empty or contains no columns. Returning empty mapping table.")
    return(
      tibble::tibble(
        id = character(),
        uniprot_id = character(),
        kinase_symbol = character(),
        score = numeric(),
        source = character()
      )
    )
  }

  cols <- colnames(signor_raw)

  # Flexible column detection supporting SIGNOR standard headers
  sub_col <- dplyr::case_when(
    "IDB"           %in% cols ~ "IDB",
    "ENTITYB_ID"    %in% cols ~ "ENTITYB_ID",
    "substrate_id"  %in% cols ~ "substrate_id",
    "SUBSTRATE_ID"  %in% cols ~ "SUBSTRATE_ID",
    "uniprot_id"    %in% cols ~ "uniprot_id",
    "UNIPROT_ID"    %in% cols ~ "UNIPROT_ID",
    "ENTITYB"       %in% cols ~ "ENTITYB",
    TRUE                      ~ NA_character_
  )

  kin_col <- dplyr::case_when(
    "ENTITYA"       %in% cols ~ "ENTITYA",
    "kinase_symbol" %in% cols ~ "kinase_symbol",
    "IDA"           %in% cols ~ "IDA",
    "kinase_id"     %in% cols ~ "kinase_id",
    "KINASE"        %in% cols ~ "KINASE",
    TRUE                      ~ NA_character_
  )

  if (is.na(sub_col)) {
    warning("Could not locate a substrate ID column in signor_raw. Available columns: ", paste(cols, collapse = ", "))
    return(
      tibble::tibble(
        id = character(),
        uniprot_id = character(),
        kinase_symbol = character(),
        score = numeric(),
        source = character()
      )
    )
  }

  signor_df <- signor_raw %>%
    dplyr::rename(substrate_id = dplyr::all_of(sub_col))

  if (!is.na(kin_col) && kin_col != "kinase_symbol") {
    signor_df <- signor_df %>% dplyr::rename(kinase_symbol = dplyr::all_of(kin_col))
  }

  target_uniprots <- target_sites %>%
    dplyr::pull(uniprot_id) %>%
    unique() %>%
    stats::na.omit()

  # Safely handle substrate matching without vector length mismatch inside filter()
  if ("IDB" %in% colnames(signor_df)) {
    processed_signor <- signor_df %>%
      dplyr::filter(substrate_id %in% target_uniprots | IDB %in% target_uniprots)
  } else {
    processed_signor <- signor_df %>%
      dplyr::filter(substrate_id %in% target_uniprots)
  }

  processed_signor <- processed_signor %>%
    dplyr::mutate(
      source = "SIGNOR",
      score = 1.0
    )

  return(processed_signor)
}
