#' Harmonize Kinase Nomenclature and UniProt Identifiers Across Sources
#'
#' @param signor_data Processed SIGNOR data
#' @param gps6_data Processed GPS 6.0 data
#' @param kl_data Processed Kinase Library data
#' @return A unified tibble with harmonized kinase names and standardized IDs
harmonize_kinase_names <- function(signor_data, gps6_data, kl_data) {
  
  # Helper to clean kinase symbol strings
  clean_symbol <- function(symbols) {
    symbols %>%
      stringr::str_to_upper() %>%
      stringr::str_replace_all("[/-]", "_") %>%
      stringr::str_trim()
  }
  
  signor_clean <- signor_data %>%
    dplyr::mutate(kinase_clean = clean_symbol(kinase_symbol)) %>%
    dplyr::select(id, substrate_uniprot, res_position, kinase_clean, source, score)
  
  gps6_clean <- gps6_data %>%
    dplyr::mutate(kinase_clean = clean_symbol(kinase_symbol)) %>%
    dplyr::select(id, substrate_uniprot = uniprot_id, res_position, kinase_clean, source, score)
  
  kl_clean <- kl_data %>%
    dplyr::mutate(kinase_clean = clean_symbol(kinase_symbol)) %>%
    dplyr::select(id, substrate_uniprot = uniprot_id, res_position, kinase_clean, source, score = percentile)
  
  harmonized <- dplyr::bind_rows(signor_clean, gps6_clean, kl_clean)
  
  return(harmonized)
}
