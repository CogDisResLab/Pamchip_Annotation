#' Perform Quality Control and Audit Checks on Mappings
#'
#' @param chip_sites Parsed PamChip input sites
#' @param final_map Final integrated ensemble mapping
#' @return A list summarizing QC metrics and target coverage
perform_mapping_qc <- function(chip_sites, final_map) {
  
  total_peptides <- nrow(chip_sites)
  mapped_peptides <- dplyr::n_distinct(final_map$id)
  coverage_pct <- (mapped_peptides / total_peptides) * 100
  
  evidence_breakdown <- final_map %>%
    dplyr::count(evidence_count, name = "pair_count")
  
  source_overlaps <- tibble::tibble(
    total_pairs = nrow(final_map),
    signor_supported = sum(final_map$SIGNOR > 0),     gps6_supported = sum(final_map$`GPS6.0` > 0),
    kinase_library_supported = sum(final_map$KinaseLibrary > 0)
  )
  
  qc_results <- list(
    total_input_peptides = total_peptides,
    mapped_peptides = mapped_peptides,
    coverage_percentage = round(coverage_pct, 2),
    evidence_breakdown = evidence_breakdown,
    source_overlaps = source_overlaps
  )
  
  message(sprintf("QC Complete: %.2f%% PamChip peptide coverage (%d/%d)", 
                  coverage_pct, mapped_peptides, total_peptides))
  
  return(qc_results)
}
