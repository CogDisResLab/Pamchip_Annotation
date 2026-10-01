#' Parse Raw PamChip CSV Layout Files and Extract Phosphopeptides
#'
#' @param chip_files Character vector of file paths in data/raw/
#' @param manifest List object containing config settings from manifest.yml
#' @return A tidy tibble of standardized PamChip target sites
parse_pamchip_annotations <- function(chip_files, manifest) {
  if (length(chip_files) == 0) {
    stop("No PamChip layout files found in data/raw/")
  }
  
  parsed_data <- purrr::map_dfr(chip_files, function(file_path) {
    ext <- tools::file_ext(file_path)
    if (ext == "csv") {
      df <- readr::read_csv(file_path, show_col_types = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      df <- readr::read_tsv(file_path, show_col_types = FALSE)
    } else {
      stop("Unsupported file extension: ", ext)
    }
    
    # 1. Filter out alignment grid control spots (#REF)
    df_clean <- df %>%
      dplyr::filter(!is.na(ID), ID != "#REF") %>%
      dplyr::transmute(
        chip_type        = manifest$chip$type %||% "STK",
        id               = as.character(ID),
        uniprot_id       = as.character(UniprotAccession),
        description      = as.character(Description),
        raw_sequence     = as.character(Sequence),
        # Remove phospho-annotations like (pS), (pT), (pY) to get clean sequence
        clean_sequence   = stringr::str_replace_all(raw_sequence, "\\(p[STY]\\)", ""),
        ser_positions    = as.character(Ser),
        thr_positions    = as.character(Thr)
      )
    
    return(df_clean)
  })
  
  # 2. Extract and expand site positions from bracketed strings [e.g., "[303, 304]"]
  expanded_sites <- parsed_data %>%
    tidyr::pivot_longer(
      cols = c(ser_positions, thr_positions),
      names_to = "residue_type",
      values_to = "pos_string"
    ) %>%
    dplyr::mutate(
      residue = dplyr::if_else(residue_type == "ser_positions", "S", "T"),
      pos_clean = stringr::str_extract_all(pos_string, "\\d+")
    ) %>%
    tidyr::unnest(pos_clean) %>%
    dplyr::transmute(
      chip_type,
      id,
      uniprot_id,
      description,
      raw_sequence,
      sequence_fragment = clean_sequence,
      residue,
      res_position = as.numeric(pos_clean)
    ) %>%
    dplyr::distinct(id, uniprot_id, res_position, .keep_all = TRUE)
  
  return(expanded_sites)
}
