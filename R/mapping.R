#' Build Weighted Ensemble Kinase-Substrate Mapping
#'
#' @param chip_sites Parsed PamChip sites tibble
#' @param harmonized_sources Harmonized predictions/curated entries from harmonize_kinase_names()
#' @param weights List of weights for each evidence source (from manifest.yml)
#' @return Final integrated ensemble score matrix
build_ensemble_mapping <- function(chip_sites, harmonized_sources, weights) {
  
  w_signor <- weights$signor %||% 0.50
  w_gps6   <- weights$gps6   %||% 0.25
  w_kl     <- weights$kinase_library %||% 0.25
  
  # Pivot sources to wide format per peptide-kinase pair
  wide_scores <- harmonized_sources %>%
    dplyr::group_by(id, substrate_uniprot, res_position, kinase_clean, source) %>%
    dplyr::summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = source,
      values_from = max_score,
      values_fill = 0
    )
  
  # Ensure all source columns exist
  if (!"SIGNOR" %in% colnames(wide_scores)) wide_scores$SIGNOR <- 0
  if (!"GPS6" %in% colnames(wide_scores)) wide_scores$GPS6 <- 0
  
  # Calculate integrated score
  ensemble <- wide_scores %>%
    dplyr::mutate(
      ensemble_score = (SIGNOR * w_signor) + (`GPS6.0` * w_gps6) + (KinaseLibrary * w_kl),
      evidence_count = (SIGNOR > 0) + (`GPS6.0` > 0) + (KinaseLibrary > 0)
    ) %>%
    dplyr::arrange(dplyr::desc(ensemble_score)) %>%
    dplyr::inner_join(chip_sites, by = c("id", "substrate_uniprot" = "uniprot_id", "res_position"))
  
  return(ensemble)
}
