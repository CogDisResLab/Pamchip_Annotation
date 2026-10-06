# R/gps6.R
#
# Helpers for integrating the local GPS 6.0 runner with the targets pipeline.
#
# Expected targets usage:
#
#   chip_fasta_file <- export_chip_fasta(
#     sites = validated_chip_sites,
#     output_path = "data/external/gps6/chip_peptides.fasta"
#   )
#
#   gps6_processed <- process_gps6_mappings(
#     gps6_file = gps6_raw_output,
#     target_sites = validated_chip_sites,
#     manifest = manifest
#   )
#
# Design principles:
#   - Export the full validated chip as ONE FASTA file.
#   - Preserve the chip peptide ID exactly as the FASTA header.
#   - Collapse exact duplicate peptide ID/sequence rows before FASTA export.
#   - Parse the local GPS6 TSV without collapsing the GPS6 hierarchy.
#   - Restrict GPS6 predictions back to the phosphoacceptor positions represented
#     by the chip whenever those positions can be recovered from the chip layout.
#   - Preserve GPS6 node, score, cutoff, pass/support status, and source identity
#     for downstream harmonization and family-level collapse.


# ==============================================================================
# Internal utilities
# ==============================================================================

.gps6_first_existing_col <- function(data, candidates, required = FALSE, label = NULL) {
  hit <- candidates[candidates %in% names(data)]

  if (length(hit) > 0L) {
    return(hit[[1L]])
  }

  if (isTRUE(required)) {
    if (is.null(label)) {
      label <- paste(candidates, collapse = ", ")
    }

    stop(
      "Could not identify the required ",
      label,
      " column. Tried: ",
      paste(candidates, collapse = ", "),
      ". Available columns: ",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  NULL
}


.gps6_clean_sequence <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  x
}


.gps6_parse_integer_positions <- function(x) {
  # Accept numeric values or strings such as:
  #   "8"
  #   "8,10"
  #   "8; 10"
  #   "S8"
  #   "S8/T10"
  #
  # Returns a sorted unique integer vector. Zero/negative values are ignored.
  if (length(x) == 0L || all(is.na(x))) {
    return(integer())
  }

  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]

  if (length(x) == 0L) {
    return(integer())
  }

  pieces <- unlist(
    regmatches(
      x,
      gregexpr("[0-9]+", x, perl = TRUE)
    ),
    use.names = FALSE
  )

  if (length(pieces) == 0L) {
    return(integer())
  }

  out <- suppressWarnings(as.integer(pieces))
  out <- out[!is.na(out) & out > 0L]

  sort(unique(out))
}


.gps6_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x)) {
    return(!is.na(x) & x != 0)
  }

  y <- tolower(trimws(as.character(x)))

  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE

  out
}


.gps6_node_depth <- function(x) {
  x <- gsub("^/+|/+$", "", as.character(x))
  ifelse(
    is.na(x) | x == "",
    NA_integer_,
    lengths(strsplit(x, "/", fixed = TRUE))
  )
}


.gps6_node_leaf <- function(x) {
  x <- gsub("^/+|/+$", "", as.character(x))

  vapply(
    strsplit(x, "/", fixed = TRUE),
    function(parts) {
      if (length(parts) == 0L) {
        return(NA_character_)
      }
      parts[[length(parts)]]
    },
    character(1)
  )
}


.gps6_node_parent <- function(x) {
  x <- gsub("^/+|/+$", "", as.character(x))

  vapply(
    strsplit(x, "/", fixed = TRUE),
    function(parts) {
      if (length(parts) <= 1L) {
        return(NA_character_)
      }
      paste(parts[-length(parts)], collapse = "/")
    },
    character(1)
  )
}


.gps6_node_root <- function(x) {
  x <- gsub("^/+|/+$", "", as.character(x))

  vapply(
    strsplit(x, "/", fixed = TRUE),
    function(parts) {
      if (length(parts) == 0L) {
        return(NA_character_)
      }
      parts[[1L]]
    },
    character(1)
  )
}


.gps6_expected_chip_sites <- function(target_sites) {
  # Build one row per intended phosphoacceptor represented by the chip.
  #
  # Preferred input:
  #   peptide_id + phosphosite_position (+ phosphosite_residue)
  #
  # Fallback for PamGene-style layouts:
  #   ID + Ser + Thr
  #
  # The Ser/Thr columns are interpreted as peptide-local positions if they
  # contain positive integer values.

  peptide_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "peptide_id",
      "PeptideID",
      "peptideID",
      "peptide",
      "ID",
      "id"
    ),
    required = TRUE,
    label = "peptide identifier"
  )

  sequence_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "peptide_sequence",
      "sequence",
      "Sequence",
      "seq"
    ),
    required = FALSE
  )

  uniprot_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "substrate_uniprot",
      "uniprot",
      "uniprot_accession",
      "UniprotAccession",
      "UniProtAccession",
      "UniProt"
    ),
    required = FALSE
  )

  phosphosite_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "phosphosite",
      "site",
      "site_label",
      "phosphosite_label"
    ),
    required = FALSE
  )

  position_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "phosphosite_position",
      "site_position",
      "peptide_site_position",
      "position"
    ),
    required = FALSE
  )

  residue_col <- .gps6_first_existing_col(
    target_sites,
    c(
      "phosphosite_residue",
      "site_residue",
      "residue"
    ),
    required = FALSE
  )

  ser_col <- .gps6_first_existing_col(
    target_sites,
    c("Ser", "ser", "SER"),
    required = FALSE
  )

  thr_col <- .gps6_first_existing_col(
    target_sites,
    c("Thr", "thr", "THR"),
    required = FALSE
  )

  out <- vector("list", nrow(target_sites))
  out_i <- 0L

  for (i in seq_len(nrow(target_sites))) {
    peptide_id <- as.character(target_sites[[peptide_col]][[i]])

    if (is.na(peptide_id) || !nzchar(peptide_id)) {
      next
    }

    sequence <- if (!is.null(sequence_col)) {
      .gps6_clean_sequence(target_sites[[sequence_col]][[i]])
    } else {
      NA_character_
    }

    substrate_uniprot <- if (!is.null(uniprot_col)) {
      as.character(target_sites[[uniprot_col]][[i]])
    } else {
      NA_character_
    }

    site_label <- if (!is.null(phosphosite_col)) {
      as.character(target_sites[[phosphosite_col]][[i]])
    } else {
      NA_character_
    }

    # --------------------------------------------------------------------------
    # Preferred canonical position column
    # --------------------------------------------------------------------------

    if (!is.null(position_col)) {
      positions <- .gps6_parse_integer_positions(
        target_sites[[position_col]][[i]]
      )

      if (length(positions) > 0L) {
        for (pos in positions) {
          residue <- if (!is.null(residue_col)) {
            toupper(as.character(target_sites[[residue_col]][[i]]))
          } else if (!is.na(sequence) && nchar(sequence) >= pos) {
            substr(sequence, pos, pos)
          } else {
            NA_character_
          }

          out_i <- out_i + 1L
          out[[out_i]] <- data.frame(
            peptide_id = peptide_id,
            peptide_site_position = as.integer(pos),
            phosphoacceptor = residue,
            phosphosite = site_label,
            substrate_uniprot = substrate_uniprot,
            stringsAsFactors = FALSE
          )
        }

        next
      }
    }

    # --------------------------------------------------------------------------
    # PamGene layout fallback: Ser / Thr columns
    # --------------------------------------------------------------------------

    ser_positions <- if (!is.null(ser_col)) {
      .gps6_parse_integer_positions(target_sites[[ser_col]][[i]])
    } else {
      integer()
    }

    thr_positions <- if (!is.null(thr_col)) {
      .gps6_parse_integer_positions(target_sites[[thr_col]][[i]])
    } else {
      integer()
    }

    if (length(ser_positions) > 0L) {
      for (pos in ser_positions) {
        out_i <- out_i + 1L
        out[[out_i]] <- data.frame(
          peptide_id = peptide_id,
          peptide_site_position = as.integer(pos),
          phosphoacceptor = "S",
          phosphosite = site_label,
          substrate_uniprot = substrate_uniprot,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(thr_positions) > 0L) {
      for (pos in thr_positions) {
        out_i <- out_i + 1L
        out[[out_i]] <- data.frame(
          peptide_id = peptide_id,
          peptide_site_position = as.integer(pos),
          phosphoacceptor = "T",
          phosphosite = site_label,
          substrate_uniprot = substrate_uniprot,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- out[seq_len(out_i)]

  if (length(out) == 0L) {
    return(
      data.frame(
        peptide_id = character(),
        peptide_site_position = integer(),
        phosphoacceptor = character(),
        phosphosite = character(),
        substrate_uniprot = character(),
        stringsAsFactors = FALSE
      )
    )
  }

  out <- dplyr::bind_rows(out)

  # If no protein-level phosphosite label is available, preserve an explicit
  # peptide-local site label rather than pretending the peptide position is a
  # UniProt residue coordinate.
  missing_site <- is.na(out$phosphosite) | !nzchar(out$phosphosite)

  out$phosphosite[missing_site] <- paste0(
    out$phosphoacceptor[missing_site],
    out$peptide_site_position[missing_site],
    "_peptide"
  )

  dplyr::distinct(
    out,
    peptide_id,
    peptide_site_position,
    phosphoacceptor,
    .keep_all = TRUE
  )
}


# ==============================================================================
# Public API: export_chip_fasta()
# ==============================================================================

export_chip_fasta <- function(
  sites,
  output_path
) {
  if (!is.data.frame(sites)) {
    stop("`sites` must be a data.frame or tibble.", call. = FALSE)
  }

  if (nrow(sites) == 0L) {
    stop("Cannot export GPS6 FASTA: `sites` contains zero rows.", call. = FALSE)
  }

  peptide_col <- .gps6_first_existing_col(
    sites,
    c(
      "peptide_id",
      "PeptideID",
      "peptideID",
      "peptide",
      "ID",
      "id"
    ),
    required = TRUE,
    label = "peptide identifier"
  )

  sequence_col <- .gps6_first_existing_col(
    sites,
    c(
      "peptide_sequence",
      "sequence",
      "Sequence",
      "seq"
    ),
    required = TRUE,
    label = "peptide sequence"
  )

  fasta_tbl <- data.frame(
    peptide_id = as.character(sites[[peptide_col]]),
    sequence = .gps6_clean_sequence(sites[[sequence_col]]),
    stringsAsFactors = FALSE
  )

  if (anyNA(fasta_tbl$peptide_id) || any(!nzchar(fasta_tbl$peptide_id))) {
    stop(
      "Cannot export GPS6 FASTA: one or more peptide identifiers are missing.",
      call. = FALSE
    )
  }

  if (anyNA(fasta_tbl$sequence) || any(!nzchar(fasta_tbl$sequence))) {
    stop(
      "Cannot export GPS6 FASTA: one or more peptide sequences are missing.",
      call. = FALSE
    )
  }

  valid_sequence <- grepl("^[A-Z*]+$", fasta_tbl$sequence)

  if (!all(valid_sequence)) {
    bad <- unique(fasta_tbl$peptide_id[!valid_sequence])

    stop(
      "Cannot export GPS6 FASTA: invalid characters were found in sequence(s) ",
      "for peptide ID(s): ",
      paste(utils::head(bad, 20L), collapse = ", "),
      if (length(bad) > 20L) " ..." else "",
      call. = FALSE
    )
  }

  # GPS6's FASTA parser stores sequences in a dictionary keyed by FASTA ID.
  # Therefore duplicate IDs with conflicting sequences would silently overwrite
  # one another and must be rejected here.
  conflicting <- fasta_tbl |>
    dplyr::distinct(peptide_id, sequence) |>
    dplyr::count(peptide_id, name = "n_sequences") |>
    dplyr::filter(n_sequences > 1L)

  if (nrow(conflicting) > 0L) {
    stop(
      "Cannot export GPS6 FASTA because the same peptide ID maps to multiple ",
      "different sequences: ",
      paste(conflicting$peptide_id, collapse = ", "),
      call. = FALSE
    )
  }

  fasta_tbl <- fasta_tbl |>
    dplyr::distinct(peptide_id, sequence)

  output_path <- normalizePath(
    output_path,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    dirname(output_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  con <- file(
    output_path,
    open = "wt"
  )

  on.exit(close(con), add = TRUE)

  for (i in seq_len(nrow(fasta_tbl))) {
    writeLines(
      c(
        paste0(">", fasta_tbl$peptide_id[[i]]),
        fasta_tbl$sequence[[i]]
      ),
      con = con
    )
  }

  output_path
}


# ==============================================================================
# Public API: process_gps6_mappings()
# ==============================================================================

process_gps6_mappings <- function(
  gps6_file,
  target_sites,
  manifest = NULL
) {
  if (!file.exists(gps6_file)) {
    stop(
      "GPS6 output file does not exist: ",
      gps6_file,
      call. = FALSE
    )
  }

  if (!is.data.frame(target_sites)) {
    stop(
      "`target_sites` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  gps <- readr::read_tsv(
    gps6_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  required <- c(
    "ID",
    "Position",
    "Code",
    "Kinase",
    "Peptide",
    "Score",
    "Cutoff"
  )

  missing <- setdiff(
    required,
    names(gps)
  )

  if (length(missing) > 0L) {
    stop(
      "GPS6 output is missing required column(s): ",
      paste(missing, collapse = ", "),
      ". Available columns: ",
      paste(names(gps), collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(gps) == 0L) {
    return(
      tibble::tibble(
        peptide_id = character(),
        peptide_site_position = integer(),
        phosphoacceptor = character(),
        phosphosite = character(),
        substrate_uniprot = character(),
        gps6_node = character(),
        gps6_node_root = character(),
        gps6_node_parent = character(),
        gps6_node_leaf = character(),
        gps6_node_depth = integer(),
        gps6_peptide_window = character(),
        gps6_score = double(),
        gps6_cutoff = double(),
        gps6_supported = logical(),
        evidence_source = character(),
        evidence_class = character()
      )
    )
  }

  gps <- gps |>
    dplyr::transmute(
      peptide_id = as.character(.data$ID),
      peptide_site_position = as.integer(.data$Position),
      phosphoacceptor = toupper(as.character(.data$Code)),
      gps6_node = gsub(
        "^/+|/+$",
        "",
        as.character(.data$Kinase)
      ),
      gps6_peptide_window = as.character(.data$Peptide),
      gps6_score = as.numeric(.data$Score),
      gps6_cutoff = suppressWarnings(as.numeric(.data$Cutoff)),
      gps6_supported = if ("Pass" %in% names(gps)) {
        .gps6_as_logical(gps$Pass)
      } else {
        # For h/m/l runs, gps6_local.py writes only retained predictions.
        # If Pass is absent, infer support from the numeric cutoff when possible.
        ifelse(
          is.na(suppressWarnings(as.numeric(gps$Cutoff))),
          TRUE,
          as.numeric(gps$Score) >= as.numeric(gps$Cutoff)
        )
      }
    )

  if (anyNA(gps$peptide_site_position)) {
    stop(
      "GPS6 output contains non-integer or missing Position values.",
      call. = FALSE
    )
  }

  gps$gps6_node_root <- .gps6_node_root(gps$gps6_node)
  gps$gps6_node_parent <- .gps6_node_parent(gps$gps6_node)
  gps$gps6_node_leaf <- .gps6_node_leaf(gps$gps6_node)
  gps$gps6_node_depth <- .gps6_node_depth(gps$gps6_node)

  # ---------------------------------------------------------------------------
  # Restrict GPS6 predictions to phosphoacceptor positions represented by chip
  # ---------------------------------------------------------------------------

  expected <- .gps6_expected_chip_sites(target_sites)

  if (nrow(expected) > 0L) {
    input_ids <- unique(
      as.character(
        target_sites[[
          .gps6_first_existing_col(
            target_sites,
            c(
              "peptide_id",
              "PeptideID",
              "peptideID",
              "peptide",
              "ID",
              "id"
            ),
            required = TRUE,
            label = "peptide identifier"
          )
        ]]
      )
    )

    annotated_ids <- unique(expected$peptide_id)
    missing_annotation_ids <- setdiff(input_ids, annotated_ids)

    if (length(missing_annotation_ids) > 0L) {
      warning(
        length(missing_annotation_ids),
        " chip peptide ID(s) had no recoverable phosphosite position and ",
        "therefore cannot contribute GPS6 mappings. Example(s): ",
        paste(
          utils::head(missing_annotation_ids, 10L),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    gps <- gps |>
      dplyr::inner_join(
        expected,
        by = c(
          "peptide_id",
          "peptide_site_position",
          "phosphoacceptor"
        )
      )
  } else {
    warning(
      "No chip phosphosite positions could be recovered from `target_sites`. ",
      "GPS6 predictions will therefore be retained for every S/T/Y site found ",
      "within each exported peptide sequence. This is less specific than ",
      "filtering to the intended chip phosphoacceptor.",
      call. = FALSE
    )

    gps$phosphosite <- paste0(
      gps$phosphoacceptor,
      gps$peptide_site_position,
      "_peptide"
    )

    uniprot_col <- .gps6_first_existing_col(
      target_sites,
      c(
        "substrate_uniprot",
        "uniprot",
        "uniprot_accession",
        "UniprotAccession",
        "UniProtAccession",
        "UniProt"
      ),
      required = FALSE
    )

    if (!is.null(uniprot_col)) {
      peptide_col <- .gps6_first_existing_col(
        target_sites,
        c(
          "peptide_id",
          "PeptideID",
          "peptideID",
          "peptide",
          "ID",
          "id"
        ),
        required = TRUE,
        label = "peptide identifier"
      )

      uniprot_lookup <- target_sites |>
        dplyr::transmute(
          peptide_id = as.character(.data[[peptide_col]]),
          substrate_uniprot = as.character(.data[[uniprot_col]])
        ) |>
        dplyr::distinct()

      gps <- gps |>
        dplyr::left_join(
          uniprot_lookup,
          by = "peptide_id"
        )
    } else {
      gps$substrate_uniprot <- NA_character_
    }
  }

  # ---------------------------------------------------------------------------
  # Add source metadata
  # ---------------------------------------------------------------------------

  gps <- gps |>
    dplyr::mutate(
      evidence_source = "GPS6",
      evidence_class = "predictive"
    ) |>
    dplyr::select(
      .data$peptide_id,
      .data$substrate_uniprot,
      .data$phosphosite,
      .data$peptide_site_position,
      .data$phosphoacceptor,
      .data$gps6_node,
      .data$gps6_node_root,
      .data$gps6_node_parent,
      .data$gps6_node_leaf,
      .data$gps6_node_depth,
      .data$gps6_peptide_window,
      .data$gps6_score,
      .data$gps6_cutoff,
      .data$gps6_supported,
      .data$evidence_source,
      .data$evidence_class
    ) |>
    dplyr::distinct()

  # Preserve deterministic ordering for reproducible diffs/output.
  gps <- gps |>
    dplyr::arrange(
      .data$peptide_id,
      .data$peptide_site_position,
      .data$gps6_node
    )

  gps
}
