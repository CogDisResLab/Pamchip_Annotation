#' Run Kinase Library Scoring via Python Reticulate Bindings
#'
#' Scores PamChip peptide phosphosites with the current kinase_library Python API.
#'
#' The current Kinase Library package exposes a Substrate class rather than a
#' top-level score_sequence() function. A substrate is represented as a peptide
#' sequence in which the phosphorylated residue is lower-case, e.g.:
#'
#'   RFIGRRQsLIEDARK
#'
#' for an S phosphosite at position 8.
#'
#' @param sites Parsed/validated PamChip sites tibble.
#' @param percentile_cutoff Numeric minimum Kinase Library percentile threshold
#'   on the 0-100 scale. Default = 90.
#' @param cache_dir Directory to store intermediate scoring results.
#'
#' @return A tibble containing Kinase Library predictions joined to chip
#'   phosphosite metadata.
#'
#' @export
run_kinase_library_scoring <- function(
  sites,
  percentile_cutoff = 90,
  cache_dir = "data/external/kinase_library"
) {

  dir.create(
    cache_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ============================================================================
  # 1. Validate input
  # ============================================================================

  if (!is.data.frame(sites)) {
    stop(
      "`sites` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  if (nrow(sites) == 0L) {
    warning(
      "`sites` contains zero rows. Returning empty Kinase Library table.",
      call. = FALSE
    )

    return(
      tibble::tibble(
        peptide_id = character(),
        substrate_uniprot = character(),
        phosphosite = character(),
        peptide_site_position = integer(),
        phosphoacceptor = character(),
        kinase_symbol = character(),
        kinase_library_score = double(),
        kinase_library_score_rank = double(),
        kinase_library_percentile = double(),
        kinase_library_percentile_rank = double(),
        kinase_library_supported = logical(),
        kinase_library_sequence = character(),
        evidence_source = character(),
        evidence_class = character()
      )
    )
  }

  if (
    !is.numeric(percentile_cutoff) ||
    length(percentile_cutoff) != 1L ||
    is.na(percentile_cutoff) ||
    percentile_cutoff < 0 ||
    percentile_cutoff > 100
  ) {
    stop(
      "`percentile_cutoff` must be one numeric value between 0 and 100.",
      call. = FALSE
    )
  }

  first_existing <- function(
    data,
    candidates,
    required = TRUE,
    label = NULL
  ) {

    hit <- candidates[
      candidates %in% names(data)
    ]

    if (length(hit) > 0L) {
      return(hit[[1L]])
    }

    if (required) {

      if (is.null(label)) {
        label <- paste(
          candidates,
          collapse = " / "
        )
      }

      stop(
        "Could not locate required ",
        label,
        " column. Available columns: ",
        paste(
          names(data),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    NULL
  }

  peptide_col <- first_existing(
    sites,
    c(
      "peptide_id",
      "id",
      "ID"
    ),
    label = "peptide identifier"
  )

  sequence_col <- first_existing(
    sites,
    c(
      "peptide_sequence",
      "sequence_fragment",
      "sequence_15mer",
      "raw_sequence"
    ),
    label = "peptide sequence"
  )

  position_col <- first_existing(
    sites,
    c(
      "phosphosite_position",
      "peptide_site_position"
    ),
    label = "peptide-local phosphosite position"
  )

  residue_col <- first_existing(
    sites,
    c(
      "residue",
      "phosphosite_residue",
      "phosphoacceptor"
    ),
    label = "phosphoacceptor residue"
  )

  uniprot_col <- first_existing(
    sites,
    c(
      "substrate_uniprot",
      "uniprot_id",
      "UniprotAccession"
    ),
    required = FALSE
  )

  phosphosite_col <- first_existing(
    sites,
    c(
      "phosphosite",
      "site"
    ),
    required = FALSE
  )

  # ============================================================================
  # 2. Configure Reticulate
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

  # ============================================================================
  # 3. Helper: encode phosphosite for Kinase Library
  # ============================================================================

  make_kinase_library_sequence <- function(
    sequence,
    phosphosite_position,
    phosphoacceptor
  ) {

    if (
      is.na(sequence) ||
      !nzchar(
        trimws(
          as.character(sequence)
        )
      )
    ) {
      return(NA_character_)
    }

    sequence <- toupper(
      trimws(
        as.character(sequence)
      )
    )

    phosphosite_position <- suppressWarnings(
      as.integer(
        phosphosite_position
      )
    )

    phosphoacceptor <- toupper(
      as.character(
        phosphoacceptor
      )
    )

    if (
      is.na(phosphosite_position) ||
      phosphosite_position < 1L ||
      phosphosite_position > nchar(sequence)
    ) {
      return(NA_character_)
    }

    if (
      is.na(phosphoacceptor) ||
      !phosphoacceptor %in% c(
        "S",
        "T",
        "Y"
      )
    ) {
      return(NA_character_)
    }

    observed <- substr(
      sequence,
      phosphosite_position,
      phosphosite_position
    )

    if (observed != phosphoacceptor) {

      stop(
        "Phosphosite residue mismatch while preparing Kinase Library input: ",
        "sequence='",
        sequence,
        "', position=",
        phosphosite_position,
        ", expected=",
        phosphoacceptor,
        ", observed=",
        observed,
        ".",
        call. = FALSE
      )
    }

    prefix <- if (phosphosite_position > 1L) {
      substr(
        sequence,
        1L,
        phosphosite_position - 1L
      )
    } else {
      ""
    }

    suffix <- if (phosphosite_position < nchar(sequence)) {
      substr(
        sequence,
        phosphosite_position + 1L,
        nchar(sequence)
      )
    } else {
      ""
    }

    paste0(
      prefix,
      tolower(
        phosphoacceptor
      ),
      suffix
    )
  }

  # ============================================================================
  # 4. Prepare one row per unique chip phosphosite
  # ============================================================================

  scoring_sites <- tibble::tibble(

    peptide_id =
      as.character(
        sites[[peptide_col]]
      ),

    sequence =
      as.character(
        sites[[sequence_col]]
      ),

    peptide_site_position =
      suppressWarnings(
        as.integer(
          sites[[position_col]]
        )
      ),

    phosphoacceptor =
      toupper(
        as.character(
          sites[[residue_col]]
        )
      ),

    substrate_uniprot =
      if (!is.null(uniprot_col)) {
        as.character(
          sites[[uniprot_col]]
        )
      } else {
        rep(
          NA_character_,
          nrow(sites)
        )
      },

    phosphosite =
      if (!is.null(phosphosite_col)) {
        as.character(
          sites[[phosphosite_col]]
        )
      } else {
        rep(
          NA_character_,
          nrow(sites)
        )
      }
  )

  scoring_sites <- dplyr::mutate(
    scoring_sites,

    kinase_library_sequence =
      mapply(
        make_kinase_library_sequence,
        sequence,
        peptide_site_position,
        phosphoacceptor,
        USE.NAMES = FALSE
      )
  )

  scoring_sites <- dplyr::filter(
    scoring_sites,
    !is.na(.data$kinase_library_sequence)
  )

  scoring_sites <- dplyr::distinct(
    scoring_sites,
    .data$peptide_id,
    .data$substrate_uniprot,
    .data$phosphosite,
    .data$peptide_site_position,
    .data$phosphoacceptor,
    .data$kinase_library_sequence,
    .keep_all = TRUE
  )

  if (nrow(scoring_sites) == 0L) {

    warning(
      "No chip phosphosites could be converted to Kinase Library substrate ",
      "sequences.",
      call. = FALSE
    )

    return(
      tibble::tibble(
        peptide_id = character(),
        substrate_uniprot = character(),
        phosphosite = character(),
        peptide_site_position = integer(),
        phosphoacceptor = character(),
        kinase_symbol = character(),
        kinase_library_score = double(),
        kinase_library_score_rank = double(),
        kinase_library_percentile = double(),
        kinase_library_percentile_rank = double(),
        kinase_library_supported = logical(),
        kinase_library_sequence = character(),
        evidence_source = character(),
        evidence_class = character()
      )
    )
  }

  message(
    "Scoring ",
    nrow(scoring_sites),
    " unique chip phosphosites with The Kinase Library..."
  )

  # ============================================================================
  # 5. Score each substrate with current Kinase Library API
  # ============================================================================

  results_list <- vector(
    "list",
    nrow(scoring_sites)
  )

  failures <- list()

  for (i in seq_len(nrow(scoring_sites))) {

    row <- scoring_sites[
      i,
      ,
      drop = FALSE
    ]

    substrate_sequence <- row$kinase_library_sequence[[1L]]

    if (
      i == 1L ||
      i %% 25L == 0L ||
      i == nrow(scoring_sites)
    ) {
      message(
        "Kinase Library substrate ",
        i,
        "/",
        nrow(scoring_sites),
        ": ",
        row$peptide_id[[1L]]
      )
    }

    result <- tryCatch(
      {

        substrate <- kl$Substrate(
  substrate_sequence,
  phos_pos = as.integer(
    row$peptide_site_position[[1L]]
  )
)

        prediction_py <- substrate$predict(
          sort_by = "percentile"
        )

        prediction <- reticulate::py_to_r(
          prediction_py
        )

        if (!is.data.frame(prediction)) {
          prediction <- as.data.frame(
            prediction
          )
        }

        if (nrow(prediction) == 0L) {
          return(NULL)
        }

        prediction <- tibble::rownames_to_column(
          prediction,
          var = "kinase_symbol"
        )

        # Normalize Python column names.
        normalized_names <- names(prediction)

        normalized_names <- gsub(
          "[^A-Za-z0-9]+",
          "_",
          normalized_names
        )

        normalized_names <- gsub(
          "^_|_$",
          "",
          normalized_names
        )

        normalized_names <- tolower(
          normalized_names
        )

        names(prediction) <- normalized_names

        find_prediction_col <- function(candidates) {

          hit <- candidates[
            candidates %in% names(prediction)
          ]

          if (length(hit) == 0L) {
            return(NULL)
          }

          hit[[1L]]
        }

        score_col <- find_prediction_col(
          c(
            "score",
            "log2_score"
          )
        )

        score_rank_col <- find_prediction_col(
          c(
            "score_rank",
            "rank_score"
          )
        )

        percentile_col <- find_prediction_col(
          c(
            "percentile"
          )
        )

        percentile_rank_col <- find_prediction_col(
          c(
            "percentile_rank",
            "rank_percentile"
          )
        )

        if (is.null(percentile_col)) {
          stop(
            "Kinase Library predict() output does not contain a percentile ",
            "column. Available columns: ",
            paste(
              names(prediction),
              collapse = ", "
            ),
            call. = FALSE
          )
        }

        result <- tibble::tibble(

          peptide_id =
            row$peptide_id[[1L]],

          substrate_uniprot =
            row$substrate_uniprot[[1L]],

          phosphosite =
            row$phosphosite[[1L]],

          peptide_site_position =
            row$peptide_site_position[[1L]],

          phosphoacceptor =
            row$phosphoacceptor[[1L]],

          kinase_symbol =
            as.character(
              prediction$kinase_symbol
            ),

          kinase_library_score =
            if (!is.null(score_col)) {
              suppressWarnings(
                as.numeric(
                  prediction[[score_col]]
                )
              )
            } else {
              rep(
                NA_real_,
                nrow(prediction)
              )
            },

          kinase_library_score_rank =
            if (!is.null(score_rank_col)) {
              suppressWarnings(
                as.numeric(
                  prediction[[score_rank_col]]
                )
              )
            } else {
              rep(
                NA_real_,
                nrow(prediction)
              )
            },

          kinase_library_percentile =
            suppressWarnings(
              as.numeric(
                prediction[[percentile_col]]
              )
            ),

          kinase_library_percentile_rank =
            if (!is.null(percentile_rank_col)) {
              suppressWarnings(
                as.numeric(
                  prediction[[percentile_rank_col]]
                )
              )
            } else {
              rep(
                NA_real_,
                nrow(prediction)
              )
            },

          kinase_library_sequence =
            substrate_sequence
        )

        result <- dplyr::mutate(
          result,

          kinase_library_supported =
            !is.na(.data$kinase_library_percentile) &
              .data$kinase_library_percentile >= percentile_cutoff,

          evidence_source =
            "Kinase Library",

          evidence_class =
            "predictive"
        )

        result <- dplyr::filter(
          result,
          .data$kinase_library_supported
        )

        result
      },

      error = function(e) {

        failures[[length(failures) + 1L]] <<- tibble::tibble(
          peptide_id = row$peptide_id[[1L]],
          kinase_library_sequence = substrate_sequence,
          error = conditionMessage(e)
        )

        NULL
      }
    )

    results_list[[i]] <- result
  }

  # ============================================================================
  # 6. Aggregate results
  # ============================================================================

  all_scores <- dplyr::bind_rows(
    results_list
  )

  if (length(failures) > 0L) {

    failure_table <- dplyr::bind_rows(
      failures
    )

    failure_file <- file.path(
      cache_dir,
      "kinase_library_failures.csv"
    )

    readr::write_csv(
      failure_table,
      failure_file
    )

    warning(
      nrow(failure_table),
      " Kinase Library substrate(s) failed scoring. Details written to ",
      failure_file,
      ".",
      call. = FALSE
    )
  }

  if (nrow(all_scores) == 0L) {

    warning(
      "Kinase Library completed but no predictions met percentile cutoff ",
      percentile_cutoff,
      ".",
      call. = FALSE
    )

    return(
      tibble::tibble(
        peptide_id = character(),
        substrate_uniprot = character(),
        phosphosite = character(),
        peptide_site_position = integer(),
        phosphoacceptor = character(),
        kinase_symbol = character(),
        kinase_library_score = double(),
        kinase_library_score_rank = double(),
        kinase_library_percentile = double(),
        kinase_library_percentile_rank = double(),
        kinase_library_supported = logical(),
        kinase_library_sequence = character(),
        evidence_source = character(),
        evidence_class = character()
      )
    )
  }

  all_scores <- dplyr::distinct(
    all_scores
  )

  all_scores <- dplyr::arrange(
    all_scores,
    .data$peptide_id,
    .data$peptide_site_position,
    dplyr::desc(
      .data$kinase_library_percentile
    ),
    .data$kinase_symbol
  )

  # ============================================================================
  # 7. Cache processed results
  # ============================================================================

  output_file <- file.path(
    cache_dir,
    "kinase_library_predictions.csv"
  )

  readr::write_csv(
    all_scores,
    output_file
  )

  message(
    "Kinase Library scoring complete: ",
    format(
      nrow(all_scores),
      big.mark = ","
    ),
    " predictions retained at percentile >= ",
    percentile_cutoff,
    "."
  )

  all_scores
}
