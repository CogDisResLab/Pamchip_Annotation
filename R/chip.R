# R/chip.R
#
# Chip-layout parsing and validation for the PamChip annotation pipeline.
#
# This module converts a vendor-supplied chip layout into a canonical schema
# used by:
#
#   - PhosphoSIGNOR matching
#   - Kinase Library scoring
#   - GPS 6.0 FASTA generation
#   - evidence harmonization
#   - KRSA family-level mapping
#
# Important distinction
# ---------------------
# PamGene layout files may contain BOTH:
#
#   1. protein-level phosphosite coordinates
#      e.g. S473 in the source protein
#
#   2. peptide-local phosphoacceptor positions
#      e.g. position 8 within a 15-mer peptide
#
# These are not interchangeable.
#
# The canonical output therefore preserves:
#
#   res_position
#       Protein / UniProt residue coordinate when supplied by the chip layout.
#
#   phosphosite_position
#       Position of the phosphoacceptor within the clean peptide sequence.
#
# GPS6 uses phosphosite_position to restrict predictions back to the intended
# array phosphoacceptor, while SIGNOR and other experimental resources can use
# res_position for protein-site matching.


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Null-coalescing operator
#'
#' Defined locally in case another sourced R file has not already defined it.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    y
  } else {
    x
  }
}


#' Read one chip-layout file
#'
#' @param file_path Path to CSV, TSV, or TXT layout file.
#' @return A tibble containing the raw layout.
.read_chip_layout_file <- function(file_path) {

  if (!file.exists(file_path)) {
    stop(
      "Chip-layout file does not exist: ",
      file_path,
      call. = FALSE
    )
  }

  ext <- tolower(
    tools::file_ext(file_path)
  )

  if (ext == "csv") {

    return(
      readr::read_csv(
        file_path,
        show_col_types = FALSE
      )
    )
  }

  if (ext %in% c("tsv", "txt")) {

    return(
      readr::read_tsv(
        file_path,
        show_col_types = FALSE
      )
    )
  }

  stop(
    "Unsupported chip-layout file extension: .",
    ext,
    ". Supported formats are CSV, TSV, and TXT.",
    call. = FALSE
  )
}


#' Parse integer positions from vendor fields
#'
#' Handles values such as:
#'
#'   473
#'   "[473]"
#'   "[473, 474]"
#'   "473,474"
#'   "S473"
#'
#' @param x Input value.
#' @return Sorted unique integer vector.
.parse_position_field <- function(x) {

  if (
    length(x) == 0L ||
    is.na(x) ||
    !nzchar(trimws(as.character(x)))
  ) {
    return(integer())
  }

  x <- as.character(x)

  matches <- stringr::str_extract_all(
    x,
    "\\d+"
  )[[1L]]

  if (length(matches) == 0L) {
    return(integer())
  }

  values <- suppressWarnings(
    as.integer(matches)
  )

  values <- values[
    !is.na(values) &
      values > 0L
  ]

  sort(
    unique(values)
  )
}


#' Remove PamGene phosphosite markup from a peptide sequence
#'
#' Examples:
#'
#'   RRR(pS)TAA -> RRRSTAA
#'   AAA(pT)BBB -> AAATBBB
#'
#' @param x Raw sequence.
#' @return Clean uppercase amino-acid sequence.
.clean_chip_sequence <- function(x) {

  x <- as.character(x)

  x <- stringr::str_to_upper(x)

  x <- stringr::str_replace_all(
    x,
    "\\s+",
    ""
  )

  x <- stringr::str_replace_all(
    x,
    "\\(P([STY])\\)",
    "\\1"
  )

  x
}


#' Locate explicitly annotated phosphoacceptors in a raw sequence
#'
#' Converts positions in a sequence containing (pS)/(pT)/(pY) markup into
#' positions in the cleaned peptide sequence.
#'
#' For example:
#'
#'   raw:   ABC(pS)DE(pT)FG
#'   clean: ABCSDETFG
#'
#' returns:
#'
#'   residue local_position
#'   S       4
#'   T       7
#'
#' @param raw_sequence Raw annotated peptide sequence.
#' @return Tibble with residue and local_position.
.extract_annotated_local_sites <- function(raw_sequence) {

  if (
    length(raw_sequence) == 0L ||
    is.na(raw_sequence) ||
    !nzchar(raw_sequence)
  ) {
    return(
      tibble::tibble(
        residue = character(),
        local_position = integer()
      )
    )
  }

  raw_sequence <- stringr::str_to_upper(
    stringr::str_replace_all(
      as.character(raw_sequence),
      "\\s+",
      ""
    )
  )

  # Find all explicit PamGene-style phosphosite tokens.
  locations <- stringr::str_locate_all(
    raw_sequence,
    "\\(P[STY]\\)"
  )[[1L]]

  if (nrow(locations) == 0L) {
    return(
      tibble::tibble(
        residue = character(),
        local_position = integer()
      )
    )
  }

  tokens <- stringr::str_sub(
    raw_sequence,
    locations[, "start"],
    locations[, "end"]
  )

  residues <- stringr::str_extract(
    tokens,
    "[STY]"
  )

  # Position in clean sequence =
  # number of real amino acids occurring before the token + 1.
  #
  # Each previous "(pX)" token represents exactly one amino acid after cleaning.
  local_positions <- vapply(
    seq_len(nrow(locations)),
    function(i) {

      prefix <- if (locations[i, "start"] > 1L) {
        stringr::str_sub(
          raw_sequence,
          1L,
          locations[i, "start"] - 1L
        )
      } else {
        ""
      }

      clean_prefix <- .clean_chip_sequence(
        prefix
      )

      nchar(clean_prefix) + 1L
    },
    integer(1)
  )

  tibble::tibble(
    residue = residues,
    local_position = local_positions
  )
}


#' Locate residue occurrences in an unannotated peptide
#'
#' This is a fallback only. If a peptide has exactly one occurrence of a
#' residue type represented by the vendor Ser/Thr field, that position can be
#' assigned unambiguously.
#'
#' @param sequence Clean peptide sequence.
#' @param residue One of S, T, or Y.
#' @return Integer positions within the peptide.
.find_residue_positions <- function(
  sequence,
  residue
) {

  if (
    is.na(sequence) ||
    !nzchar(sequence)
  ) {
    return(integer())
  }

  chars <- strsplit(
    sequence,
    "",
    fixed = TRUE
  )[[1L]]

  which(
    chars == residue
  )
}


#' Match protein positions to peptide-local annotated positions
#'
#' PamGene Ser/Thr fields contain protein-level positions. The raw peptide
#' sequence may separately mark which S/T residue is phosphorylated using
#' (pS)/(pT).
#'
#' This function combines those two representations.
#'
#' @param raw_sequence Vendor sequence.
#' @param clean_sequence Sequence after removal of phosphosite markup.
#' @param ser_field Vendor Ser field.
#' @param thr_field Vendor Thr field.
#' @return Tibble containing one row per intended phosphosite.
.build_site_table <- function(
  raw_sequence,
  clean_sequence,
  ser_field,
  thr_field
) {

  protein_ser <- .parse_position_field(
    ser_field
  )

  protein_thr <- .parse_position_field(
    thr_field
  )

  annotated <- .extract_annotated_local_sites(
    raw_sequence
  )

  annotated_ser <- annotated %>%
    dplyr::filter(.data$residue == "S") %>%
    dplyr::pull(.data$local_position)

  annotated_thr <- annotated %>%
    dplyr::filter(.data$residue == "T") %>%
    dplyr::pull(.data$local_position)

  # --------------------------------------------------------------------------
  # Fallback for layouts where Sequence lacks explicit (pS)/(pT) markup.
  #
  # Only infer a local position when the assignment is unambiguous.
  # --------------------------------------------------------------------------

  if (
    length(protein_ser) > 0L &&
    length(annotated_ser) == 0L
  ) {

    candidate_ser <- .find_residue_positions(
      clean_sequence,
      "S"
    )

    if (
      length(protein_ser) == 1L &&
      length(candidate_ser) == 1L
    ) {
      annotated_ser <- candidate_ser
    }
  }

  if (
    length(protein_thr) > 0L &&
    length(annotated_thr) == 0L
  ) {

    candidate_thr <- .find_residue_positions(
      clean_sequence,
      "T"
    )

    if (
      length(protein_thr) == 1L &&
      length(candidate_thr) == 1L
    ) {
      annotated_thr <- candidate_thr
    }
  }

  rows <- list()


  # --------------------------------------------------------------------------
  # Serine
  # --------------------------------------------------------------------------

  if (
    length(protein_ser) > 0L ||
    length(annotated_ser) > 0L
  ) {

    n <- max(
      length(protein_ser),
      length(annotated_ser)
    )

    rows[["S"]] <- tibble::tibble(

      residue =
        rep(
          "S",
          n
        ),

      res_position =
        c(
          protein_ser,
          rep(
            NA_integer_,
            n - length(protein_ser)
          )
        ),

      phosphosite_position =
        c(
          annotated_ser,
          rep(
            NA_integer_,
            n - length(annotated_ser)
          )
        )
    )
  }


  # --------------------------------------------------------------------------
  # Threonine
  # --------------------------------------------------------------------------

  if (
    length(protein_thr) > 0L ||
    length(annotated_thr) > 0L
  ) {

    n <- max(
      length(protein_thr),
      length(annotated_thr)
    )

    rows[["T"]] <- tibble::tibble(

      residue =
        rep(
          "T",
          n
        ),

      res_position =
        c(
          protein_thr,
          rep(
            NA_integer_,
            n - length(protein_thr)
          )
        ),

      phosphosite_position =
        c(
          annotated_thr,
          rep(
            NA_integer_,
            n - length(annotated_thr)
          )
        )
    )
  }


  if (length(rows) == 0L) {

    return(
      tibble::tibble(
        residue = character(),
        res_position = integer(),
        phosphosite_position = integer()
      )
    )
  }

  dplyr::bind_rows(
    rows
  )
}


# ==============================================================================
# Public API: parse_chip_annotations()
# ==============================================================================

#' Parse and Standardize a Chip Layout
#'
#' Reads one chip-layout file and converts it to the canonical schema expected
#' by the downstream annotation pipeline.
#'
#' @param chip_file Path to chip-layout file.
#' @param manifest Parsed manifest.yml configuration.
#'
#' @return A tibble with one row per intended phosphoacceptor.
#'
#' @export
parse_chip_annotations <- function(
  chip_file,
  manifest
) {

  if (
    length(chip_file) != 1L ||
    is.na(chip_file) ||
    !nzchar(chip_file)
  ) {
    stop(
      "`chip_file` must be exactly one valid file path.",
      call. = FALSE
    )
  }

  raw <- .read_chip_layout_file(
    chip_file
  )


  # ==========================================================================
  # Validate vendor columns
  # ==========================================================================

  required_columns <- c(
    "ID",
    "Sequence"
  )

  missing_required <- setdiff(
    required_columns,
    names(raw)
  )

  if (length(missing_required) > 0L) {
    stop(
      "Chip layout is missing required column(s): ",
      paste(
        missing_required,
        collapse = ", "
      ),
      ". Available columns: ",
      paste(
        names(raw),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # Optional vendor columns.
  #
  # Missing columns are created so that the parser remains usable with
  # different chip-layout versions.
  optional_columns <- c(
    "Row",
    "Col",
    "Ser",
    "Thr",
    "UniprotAccession",
    "Description",
    "Xoff",
    "Yoff"
  )

  for (column in optional_columns) {

    if (!column %in% names(raw)) {
      raw[[column]] <- NA
    }
  }


  # ==========================================================================
  # Remove non-peptide/control spots
  # ==========================================================================

  raw <- raw %>%

    dplyr::filter(
      !is.na(.data$ID),
      nzchar(
        stringr::str_trim(
          as.character(.data$ID)
        )
      ),
      .data$ID != "#REF"
    )


  if (nrow(raw) == 0L) {
    stop(
      "No peptide rows remained after removal of control/empty spots.",
      call. = FALSE
    )
  }


  # ==========================================================================
  # Parse each chip row
  # ==========================================================================

  parsed_rows <- vector(
    "list",
    nrow(raw)
  )


  for (i in seq_len(nrow(raw))) {

    peptide_id <- as.character(
      raw$ID[[i]]
    )

    raw_sequence <- as.character(
      raw$Sequence[[i]]
    )

    clean_sequence <- .clean_chip_sequence(
      raw_sequence
    )

    site_table <- .build_site_table(
      raw_sequence = raw_sequence,
      clean_sequence = clean_sequence,
      ser_field = raw$Ser[[i]],
      thr_field = raw$Thr[[i]]
    )


    # ------------------------------------------------------------------------
    # A chip row without a usable S/T annotation is retained once.
    #
    # Validation will flag it downstream rather than silently deleting it.
    # ------------------------------------------------------------------------

    if (nrow(site_table) == 0L) {

      site_table <- tibble::tibble(
        residue = NA_character_,
        res_position = NA_integer_,
        phosphosite_position = NA_integer_
      )
    }


    parsed_rows[[i]] <- site_table %>%

      dplyr::mutate(

        chip_type =
          manifest$chip$type %||%
          "STK",

        chip_name =
          manifest$chip$name %||%
          basename(chip_file),

        chip_file =
          basename(chip_file),

        row =
          suppressWarnings(
            as.integer(raw$Row[[i]])
          ),

        col =
          suppressWarnings(
            as.integer(raw$Col[[i]])
          ),

        peptide_id =
          peptide_id,

        # Legacy alias retained because some existing source-specific code may
        # still use `id`.
        id =
          peptide_id,

        substrate_uniprot =
          as.character(
            raw$UniprotAccession[[i]]
          ),

        # Legacy alias retained for existing functions.
        uniprot_id =
          as.character(
            raw$UniprotAccession[[i]]
          ),

        description =
          as.character(
            raw$Description[[i]]
          ),

        raw_sequence =
          raw_sequence,

        peptide_sequence =
          clean_sequence,

        # Legacy alias retained for existing Kinase Library / other functions.
        sequence_fragment =
          clean_sequence,

        phosphosite_residue =
          .data$residue,

        phosphosite =
          dplyr::if_else(
            !is.na(.data$residue) &
              !is.na(.data$res_position),
            paste0(
              .data$residue,
              .data$res_position
            ),
            NA_character_
          ),

        ser_positions_raw =
          as.character(
            raw$Ser[[i]]
          ),

        thr_positions_raw =
          as.character(
            raw$Thr[[i]]
          ),

        x_offset =
          suppressWarnings(
            as.numeric(raw$Xoff[[i]])
          ),

        y_offset =
          suppressWarnings(
            as.numeric(raw$Yoff[[i]])
          )
      ) %>%

      dplyr::select(
        .data$chip_type,
        .data$chip_name,
        .data$chip_file,
        .data$row,
        .data$col,
        .data$peptide_id,
        .data$id,
        .data$substrate_uniprot,
        .data$uniprot_id,
        .data$description,
        .data$raw_sequence,
        .data$peptide_sequence,
        .data$sequence_fragment,
        .data$phosphosite_residue,
        .data$residue,
        .data$res_position,
        .data$phosphosite,
        .data$phosphosite_position,
        .data$ser_positions_raw,
        .data$thr_positions_raw,
        .data$x_offset,
        .data$y_offset
      )
  }


  parsed <- dplyr::bind_rows(
    parsed_rows
  )


  # ==========================================================================
  # Remove exact duplicates only
  # ==========================================================================

  parsed <- parsed %>%

    dplyr::distinct(

      .data$peptide_id,
      .data$substrate_uniprot,
      .data$residue,
      .data$res_position,
      .data$phosphosite_position,

      .keep_all = TRUE
    ) %>%

    dplyr::arrange(
      .data$row,
      .data$col,
      .data$peptide_id,
      .data$res_position
    )


  parsed
}


# ==============================================================================
# Public API: validate_chip_annotations()
# ==============================================================================

#' Validate Standardized Chip Annotation
#'
#' Performs structural and biological sanity checks before downstream resources
#' are queried.
#'
#' Fatal problems stop the pipeline. Potentially incomplete annotations generate
#' warnings so that they are visible without silently removing chip content.
#'
#' @param chip_sites Output of parse_chip_annotations().
#' @param manifest Parsed manifest.yml configuration.
#'
#' @return Validated chip-site tibble.
#'
#' @export
validate_chip_annotations <- function(
  chip_sites,
  manifest
) {

  if (!is.data.frame(chip_sites)) {
    stop(
      "`chip_sites` must be a data.frame or tibble.",
      call. = FALSE
    )
  }


  if (nrow(chip_sites) == 0L) {
    stop(
      "Chip annotation contains zero rows.",
      call. = FALSE
    )
  }


  # ==========================================================================
  # Required canonical columns
  # ==========================================================================

  required_columns <- c(
    "peptide_id",
    "peptide_sequence",
    "substrate_uniprot",
    "residue",
    "res_position",
    "phosphosite_position"
  )

  missing_columns <- setdiff(
    required_columns,
    names(chip_sites)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Parsed chip annotation is missing required canonical column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Peptide identifiers
  # ==========================================================================

  bad_id <- is.na(chip_sites$peptide_id) |
    !nzchar(
      stringr::str_trim(
        as.character(
          chip_sites$peptide_id
        )
      )
    )

  if (any(bad_id)) {
    stop(
      sum(bad_id),
      " parsed chip row(s) have missing peptide identifiers.",
      call. = FALSE
    )
  }


  # ==========================================================================
  # Sequences
  # ==========================================================================

  chip_sites <- chip_sites %>%

    dplyr::mutate(
      peptide_sequence =
        .clean_chip_sequence(
          .data$peptide_sequence
        ),

      sequence_fragment =
        .data$peptide_sequence
    )


  bad_sequence <- is.na(
    chip_sites$peptide_sequence
  ) |
    !nzchar(
      chip_sites$peptide_sequence
    )


  if (any(bad_sequence)) {
    stop(
      sum(bad_sequence),
      " parsed chip row(s) have missing peptide sequences.",
      call. = FALSE
    )
  }


  valid_aa <- grepl(
    "^[ACDEFGHIKLMNPQRSTVWY*]+$",
    chip_sites$peptide_sequence
  )


  if (!all(valid_aa)) {

    bad_ids <- unique(
      chip_sites$peptide_id[
        !valid_aa
      ]
    )

    stop(
      "Invalid amino-acid characters were found in peptide sequence(s): ",
      paste(
        utils::head(
          bad_ids,
          20L
        ),
        collapse = ", "
      ),
      if (length(bad_ids) > 20L) {
        " ..."
      } else {
        ""
      },
      call. = FALSE
    )
  }


  # ==========================================================================
  # Duplicate peptide IDs with conflicting sequences
  # ==========================================================================

  conflicts <- chip_sites %>%

    dplyr::distinct(
      .data$peptide_id,
      .data$peptide_sequence
    ) %>%

    dplyr::count(
      .data$peptide_id,
      name = "n_sequences"
    ) %>%

    dplyr::filter(
      .data$n_sequences > 1L
    )


  if (nrow(conflicts) > 0L) {

    stop(
      "The same peptide ID maps to multiple peptide sequences: ",
      paste(
        conflicts$peptide_id,
        collapse = ", "
      ),
      ". GPS6 FASTA headers must uniquely identify one sequence.",
      call. = FALSE
    )
  }


  # ==========================================================================
  # Phosphoacceptor residues
  # ==========================================================================

  invalid_residue <- !is.na(
    chip_sites$residue
  ) &
    !chip_sites$residue %in% c(
      "S",
      "T",
      "Y"
    )


  if (any(invalid_residue)) {

    stop(
      "Invalid phosphoacceptor residue(s) detected: ",
      paste(
        unique(
          chip_sites$residue[
            invalid_residue
          ]
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Protein-level residue coordinates
  # ==========================================================================

  invalid_protein_position <- !is.na(
    chip_sites$res_position
  ) &
    chip_sites$res_position <= 0


  if (any(invalid_protein_position)) {

    stop(
      "Protein residue positions must be positive integers.",
      call. = FALSE
    )
  }


  # ==========================================================================
  # Peptide-local residue coordinates
  # ==========================================================================

  invalid_local_position <- !is.na(
    chip_sites$phosphosite_position
  ) &
    (
      chip_sites$phosphosite_position <= 0L |
        chip_sites$phosphosite_position >
          nchar(
            chip_sites$peptide_sequence
          )
    )


  if (any(invalid_local_position)) {

    bad <- chip_sites[
      invalid_local_position,
      c(
        "peptide_id",
        "peptide_sequence",
        "phosphosite_position"
      )
    ]

    stop(
      "One or more peptide-local phosphosite positions fall outside the ",
      "peptide sequence. Example: ",
      paste(
        utils::capture.output(
          print(
            utils::head(
              bad,
              5L
            )
          )
        ),
        collapse = " "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Confirm local position corresponds to expected amino acid
  # ==========================================================================

  checkable <- !is.na(
    chip_sites$phosphosite_position
  ) &
    !is.na(
      chip_sites$residue
    )


  observed_residue <- rep(
    NA_character_,
    nrow(chip_sites)
  )


  observed_residue[checkable] <- mapply(

    FUN = function(sequence, position) {

      substr(
        sequence,
        position,
        position
      )
    },

    sequence =
      chip_sites$peptide_sequence[
        checkable
      ],

    position =
      chip_sites$phosphosite_position[
        checkable
      ],

    USE.NAMES = FALSE
  )


  residue_mismatch <- checkable &
    observed_residue != chip_sites$residue


  if (any(residue_mismatch)) {

    bad <- chip_sites[
      residue_mismatch,
      c(
        "peptide_id",
        "peptide_sequence",
        "residue",
        "phosphosite_position"
      )
    ]


    stop(
      "Peptide-local phosphosite annotation does not match the amino acid ",
      "present in the peptide sequence. Example: ",
      paste(
        utils::capture.output(
          print(
            utils::head(
              bad,
              5L
            )
          )
        ),
        collapse = " "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Missing phosphosite-localization warnings
  # ==========================================================================

  no_local_site <- is.na(
    chip_sites$phosphosite_position
  )


  if (any(no_local_site)) {

    affected <- unique(
      chip_sites$peptide_id[
        no_local_site
      ]
    )


    warning(
      length(affected),
      " peptide ID(s) lack an unambiguous peptide-local phosphosite ",
      "position. GPS6 predictions for those peptides cannot be restricted ",
      "to the intended chip phosphoacceptor and should not contribute to the ",
      "final GPS6 mapping unless resolved. Example(s): ",
      paste(
        utils::head(
          affected,
          10L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Missing protein-coordinate warnings
  # ==========================================================================

  no_protein_site <- is.na(
    chip_sites$res_position
  )


  if (any(no_protein_site)) {

    affected <- unique(
      chip_sites$peptide_id[
        no_protein_site
      ]
    )


    warning(
      length(affected),
      " peptide ID(s) lack a protein-level residue position. ",
      "Experimental kinase-substrate matching may be unavailable for these ",
      "entries. Example(s): ",
      paste(
        utils::head(
          affected,
          10L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Missing UniProt warnings
  # ==========================================================================

  no_uniprot <- is.na(
    chip_sites$substrate_uniprot
  ) |
    !nzchar(
      stringr::str_trim(
        as.character(
          chip_sites$substrate_uniprot
        )
      )
    )


  if (any(no_uniprot)) {

    affected <- unique(
      chip_sites$peptide_id[
        no_uniprot
      ]
    )


    warning(
      length(affected),
      " peptide ID(s) lack a substrate UniProt accession. ",
      "PhosphoSIGNOR matching may therefore be unavailable for these entries. ",
      "Example(s): ",
      paste(
        utils::head(
          affected,
          10L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # ==========================================================================
  # Final deterministic ordering
  # ==========================================================================

  chip_sites <- chip_sites %>%

    dplyr::arrange(
      .data$row,
      .data$col,
      .data$peptide_id,
      .data$res_position,
      .data$phosphosite_position
    )


  chip_sites
}


# ==============================================================================
# Backward-compatible wrapper
# ==============================================================================

#' Parse Raw PamChip Layout Files
#'
#' Backward-compatible wrapper around parse_chip_annotations().
#'
#' The previous pipeline accepted a vector of chip-layout files and combined
#' them. The current targets pipeline tracks one manifest-selected layout file.
#' This wrapper is retained so older scripts do not immediately break.
#'
#' @param chip_files Character vector of chip-layout files.
#' @param manifest Parsed manifest configuration.
#'
#' @return Combined standardized chip annotation.
#'
#' @export
parse_pamchip_annotations <- function(
  chip_files,
  manifest
) {

  if (length(chip_files) == 0L) {

    stop(
      "No chip-layout files were supplied.",
      call. = FALSE
    )
  }


  purrr::map_dfr(
    chip_files,
    function(chip_file) {

      parse_chip_annotations(
        chip_file = chip_file,
        manifest = manifest
      )
    }
  )
}