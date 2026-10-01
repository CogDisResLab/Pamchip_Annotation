#' Generate Summary Tables for Quarto Report
#'
#' @param final_map Final ensemble mapping tibble
#' @param qc_results QC list from perform_mapping_qc()
#' @return A list of summary tables ready for reporting display
generate_report_tables <- function(final_map, qc_results) {
  
  top_kinases <- final_map %>%
    dplyr::group_by(kinase_clean) %>%
    dplyr::summarise(
      target_peptides = dplyr::n_distinct(id),
      mean_ensemble_score = mean(ensemble_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::slice_max(order_by = target_peptides, n = 20)
  
  coverage_summary <- tibble::tibble(
    Metric = c("Total Input Peptides", "Mapped Peptides", "Coverage (%)"),
    Value = c(
      qc_results$total_input_peptides,
      qc_results$mapped_peptides,
      paste0(qc_results$coverage_percentage, "%")
    )
  )
  
  return(list(
    top_kinases = top_kinases,
    coverage_summary = coverage_summary
  ))
}
