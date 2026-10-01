#' Run Kinase Library Scoring via Python Reticulate Bindings
#'
#' @param sites Parsed PamChip sites tibble
#' @param percentile_cutoff Numeric, minimum percentile score threshold (0-1)
#' @param cache_dir Directory to store intermediate scoring results
#' @return A processed tibble of Kinase Library score predictions
run_kinase_library_scoring <- function(sites, percentile_cutoff = 0.90, cache_dir = "data/external/kinase_library") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------------------------
  # 1. Standardize Sequence and ID Columns from parsed_chip_sites
  # ------------------------------------------------------------------------------
  if (!"sequence_15mer" %in% colnames(sites)) {
    if ("sequence_fragment" %in% colnames(sites)) {
      sites <- sites %>% dplyr::rename(sequence_15mer = sequence_fragment)
    } else if ("raw_sequence" %in% colnames(sites)) {
      sites <- sites %>% dplyr::rename(sequence_15mer = raw_sequence)
    } else {
      stop("Could not locate a valid peptide sequence column in parsed_chip_sites.")
    }
  }

  if (!"id" %in% colnames(sites)) {
    if ("ID" %in% colnames(sites)) {
      sites <- sites %>% dplyr::rename(id = ID)
    } else if ("site_id" %in% colnames(sites)) {
      sites <- sites %>% dplyr::rename(id = site_id)
    }
  }

  # ------------------------------------------------------------------------------
  # 2. Configure Reticulate & Python Virtual Environment
  # ------------------------------------------------------------------------------
  py_path <- Sys.getenv("RETICULATE_PYTHON", unset = "/root/.virtualenvs/pamchip-env/bin/python")

  if (file.exists(py_path)) {
    reticulate::use_python(py_path, required = TRUE)
  } else {
    # Fallback to the named virtualenv if path isn't explicit
    reticulate::use_virtualenv("pamchip-env", required = TRUE)
  }
  
  kl <- reticulate::import("kinase_library")

  # ------------------------------------------------------------------------------
  # 3. Extract Valid Non-Empty Sequences for Scoring
  # ------------------------------------------------------------------------------
  sequences <- sites$sequence_15mer
  ids <- sites$id

  results_list <- list()
  
  message("Scoring ", length(sequences), " sequences with The Kinase Library...")

  for (i in seq_along(sequences)) {
    seq <- sequences[i]
    chip_id <- ids[i]
    
    # Skip NA or empty sequence strings
    if (is.na(seq) || nchar(trimws(seq)) == 0) {
      next
    }
    
    tryCatch({
      # Call Kinase Library scoring function on 15-mer sequence
      scores_df <- kl$score_sequence(seq)
      
      # Convert pandas DataFrame or Series to R tibble
      r_scores <- reticulate::py_to_r(scores_df) %>%
        tibble::rownames_to_column(var = "kinase_symbol") %>%
        dplyr::rename(score = 2) %>%
        dplyr::mutate(
          id = chip_id,
          sequence_15mer = seq,
          percentile = dplyr::percent_rank(score),
          source = "KinaseLibrary"
        ) %>%
        dplyr::filter(percentile >= percentile_cutoff)
      
      results_list[[i]] <- r_scores
    }, error = function(e) {
      warning("Failed to score sequence ", seq, " for ID ", chip_id, ": ", e$message)
    })
  }
  
  # ------------------------------------------------------------------------------
  # 4. Aggregate & Join Predictions back with Input Sites
  # ------------------------------------------------------------------------------
  all_scores <- dplyr::bind_rows(results_list)

  if (nrow(all_scores) == 0) {
    warning("No sequence predictions met the criteria or succeeded during scoring.")
    return(tibble::tibble())
  }

  final_scores <- dplyr::inner_join(all_scores, sites, by = c("id", "sequence_15mer"))
  
  return(final_scores)
}
