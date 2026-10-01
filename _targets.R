# _targets.R
library(targets)
library(tarchetypes)

# ------------------------------------------------------------------------------
# 1. Global Target Options & Reticulate Environment
# ------------------------------------------------------------------------------
# Explicitly configure reticulate to use the shared-library-enabled virtualenv
Sys.setenv(RETICULATE_PYTHON = "/root/.virtualenvs/pamchip-env/bin/python")

tar_option_set(
  packages = c(
    "yaml",
    "dplyr",
    "purrr",
    "readr",
    "stringr",
    "tidyr",
    "httr2",
    "jsonlite",
    "reticulate"
  ),
  format = "rds"
)

# ------------------------------------------------------------------------------
# 2. Source R Functions from R/ Directory
# ------------------------------------------------------------------------------
# Sources chip.R, gps6.R, harmonize.R, kinase_library.R, mapping.R, qc.R, reporting.R, signor.R
lapply(list.files("R", full.names = TRUE, pattern = "\\.R$"), source)

# ------------------------------------------------------------------------------
# 3. Pipeline Targets Plan
# ------------------------------------------------------------------------------
list(
  # ============================================================================
  # Step A: Configuration & PamChip Input File Tracking
  # ============================================================================
  tar_target(
    name = manifest_file,
    command = "config/manifest.yml",
    format = "file"
  ),
  
  tar_target(
    name = manifest,
    command = yaml::read_yaml(manifest_file)
  ),
  
  # Tracks PamGene array layout files in data/raw/
  tar_target(
    name = raw_chip_files,
    command = list.files(
      path = file.path("data", "raw"),
      pattern = "\\.(csv|tsv|txt|xlsx)$",
      full.names = TRUE
    ),
    format = "file"
  ),

  # ============================================================================
  # Step B: PamChip Peptide Annotation Parsing
  # ============================================================================
  # Extracts PamGene peptide IDs, standardizes 15-mer sequence windows, and maps UniProt IDs
  tar_target(
    name = parsed_chip_sites,
    command = parse_pamchip_annotations(raw_chip_files, manifest)
  ),

  # ============================================================================
  # Step C: Database Query - PhosphoSIGNOR REST API
  # ============================================================================
  tar_target(
    name = signor_raw,
    command = fetch_signor_data(
      cache_dir = file.path("data", "external", "signor"),
      force_update = manifest$sources$signor$force_update %||% FALSE
    )
  ),
  
  tar_target(
    name = signor_processed,
    command = process_signor_mappings(
      signor_raw = signor_raw, 
      target_sites = parsed_chip_sites
    )
  ),

  # ============================================================================
  # Step D: Predictive Tool 1 - GPS 6.0 Standalone Executable
  # ============================================================================
  # Formats PamChip 15-mer peptide windows into FASTA
  tar_target(
    name = chip_fasta_file,
    command = export_chip_fasta(
      sites = parsed_chip_sites,
      output_path = file.path("data", "external", "gps6", "chip_peptides.fasta")
    ),
    format = "file"
  ),
  
  # Runs the /opt/gps6 binary installed in Docker
  tar_target(
    name = gps6_raw_output,
    command = run_gps6_cli(
      fasta_path = chip_fasta_file,
      out_dir = file.path("data", "external", "gps6"),
      threshold = manifest$sources$gps6$threshold %||% "medium"
    ),
    format = "file"
  ),
  
  tar_target(
    name = gps6_processed,
    command = parse_gps6_output(
      gps6_file = gps6_raw_output, 
      target_sites = parsed_chip_sites
    )
  ),

  # ============================================================================
  # Step E: Predictive Tool 2 - The Kinase Library (Python Module via Reticulate)
  # ============================================================================
  # Scores 15-mer sequence windows using Python kinase-library package
  tar_target(
    name = kinase_library_processed,
    command = run_kinase_library_scoring(
      sites = parsed_chip_sites,
      percentile_cutoff = manifest$sources$kinase_library$cutoff %||% 0.90,
      cache_dir = file.path("data", "external", "kinase_library")
    )
  ),

  # ============================================================================
  # Step F: Entity Harmonization & Weighted Ensemble Mapping
  # ============================================================================
  # Standardizes Kinase UniProt IDs and HGNC symbols across SIGNOR, GPS 6.0, and Kinase Library
  tar_target(
    name = harmonized_mappings,
    command = harmonize_kinase_names(
      signor_data = signor_processed,
      gps6_data = gps6_processed,
      kl_data = kinase_library_processed
    )
  ),

  # Combines predictions and curated interactions into a final score matrix
  tar_target(
    name = final_kinase_substrate_map,
    command = build_ensemble_mapping(
      chip_sites = parsed_chip_sites,
      harmonized_sources = harmonized_mappings,
      weights = manifest$mapping_weights
    )
  ),

  # ============================================================================
  # Step G: Quality Control Checks
  # ============================================================================
  tar_target(
    name = qc_summary,
    command = perform_mapping_qc(
      chip_sites = parsed_chip_sites,
      final_map = final_kinase_substrate_map
    )
  ),

  # ============================================================================
  # Step H: Quarto HTML Report Generation
  # ============================================================================
  tar_quarto(
    name = mapping_report,
    path = "report/mapping_report.qmd",
    extra_files = c("config/manifest.yml"),
    quiet = FALSE
  )
)
