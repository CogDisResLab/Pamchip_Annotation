# _targets.R
library(targets)
library(tarchetypes)

# ------------------------------------------------------------------------------
# 1. Global Target Options & Python Environment
# ------------------------------------------------------------------------------

# Reticulate is used by The Kinase Library helper functions.
# GPS6 itself is executed as a standalone Python process with system2().
Sys.setenv(
  RETICULATE_PYTHON = "/root/.virtualenvs/pamchip-env/bin/python"
)

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

# Expected helper files include:
#   chip.R
#   gps6.R
#   harmonize.R
#   kinase_library.R
#   kinase_taxonomy.R
#   mapping.R
#   qc.R
#   reporting.R
#   signor.R
#
# gps6.R should provide:
#   export_chip_fasta()
#   process_gps6_mappings()
#
# harmonize.R should accept GPS6 as a third evidence source in
# harmonize_kinase_evidence().
lapply(
  list.files(
    path = "R",
    full.names = TRUE,
    pattern = "\\.R$"
  ),
  source
)

# ------------------------------------------------------------------------------
# 3. Pipeline Targets Plan
# ------------------------------------------------------------------------------

list(

  # ============================================================================
  # Step A: Configuration & Chip Input Tracking
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

  # Track the chip-layout file specified in the manifest.
  #
  # The pipeline should support any chip layout that can be standardized by
  # parse_chip_annotations() / validate_chip_annotations().
  tar_target(
    name = raw_chip_file,
    command = manifest$chip$file,
    format = "file"
  ),

  # ============================================================================
  # Step B: Chip Annotation Parsing & Standardization
  # ============================================================================

  tar_target(
    name = parsed_chip_sites,
    command = parse_chip_annotations(
      chip_file = raw_chip_file,
      manifest = manifest
    )
  ),

  tar_target(
    name = validated_chip_sites,
    command = validate_chip_annotations(
      chip_sites = parsed_chip_sites,
      manifest = manifest
    )
  ),

  # ============================================================================
  # Step C: Experimental Evidence - PhosphoSIGNOR
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
      target_sites = validated_chip_sites
    )
  ),

  # ============================================================================
  # Step D: Predictive Evidence - The Kinase Library
  # ============================================================================

  tar_target(
  name = kinase_library_processed,
  command = run_kinase_library_scoring(
    sites = validated_chip_sites,
    percentile_cutoff =
      manifest$sources$kinase_library$percentile_cutoff %||% 90,
    cache_dir = file.path(
      "data",
      "external",
      "kinase_library"
    )
  )
),

  # ============================================================================
  # Step E: External Human Kinase Taxonomy
  # ============================================================================

  # KinHub/OpenKinome provides a machine-readable human kinase catalog with
  # canonical HGNC gene names, UniProt accessions, Manning-style groups,
  # families, and subfamilies. This is the primary gene -> family taxonomy.
  tar_target(
    name = kinhub_taxonomy_file,
    command = fetch_kinhub_taxonomy(
      cache_dir = file.path(
        "data",
        "external",
        "kinase_taxonomy"
      ),
      force_update =
        manifest$sources$kinase_taxonomy$force_update %||%
          FALSE,
      url =
        manifest$sources$kinase_taxonomy$url %||%
          KINHUB_TAXONOMY_URL
    ),
    format = "file"
  ),

  tar_target(
    name = kinase_taxonomy,
    command = build_kinase_taxonomy(
      kinhub_file = kinhub_taxonomy_file
    )
  ),

  # ============================================================================
  # Step F: Predictive Evidence - GPS 6.0
  # ============================================================================

  # Track the local GPS6 runner itself so changes to the Python implementation
  # invalidate downstream GPS6 targets.
  tar_target(
    name = gps6_script,
    command =
      manifest$sources$gps6$script %||%
        "bin/GPS6.0/cmd/gps6_local.py",
    format = "file"
  ),

  # Track all GPS6 feature resources and model files as explicit dependencies.
  tar_target(
    name = gps6_model_files,
    command = {
      pre_dir <-
        manifest$sources$gps6$pre_dir %||%
        "data/pre"

      files <- list.files(
        path = pre_dir,
        recursive = TRUE,
        full.names = TRUE,
        all.files = FALSE
      )

      files <- files[file.info(files)$isdir %in% FALSE]

      if (length(files) == 0L) {
        stop(
          "No GPS6 model/resource files were found under: ",
          pre_dir
        )
      }

      files
    },
    format = "file"
  ),

  # Derive an explicit gene -> fine family -> broad superfamily crosswalk
  # from the tracked GPS6 model hierarchy. Example:
  #
  #   AKT1 -> AKT -> AGC
  #
  # This lets Kinase Library and PhosphoSIGNOR evidence be harmonized to the
  # same fine family level used by GPS6 while preserving the broader group.
  tar_target(
    name = gps6_hierarchy,
    command = build_gps6_kinase_hierarchy(
      gps6_model_files = gps6_model_files,
      pre_dir =
        manifest$sources$gps6$pre_dir %||%
          "data/pre"
    )
  ),


  # Merge the independent KinHub taxonomy with GPS6's source-native hierarchy.
  #
  # KinHub is authoritative for canonical gene -> fine family -> superfamily
  # classification. GPS6-specific aliases remain in the table so the existing
  # harmonizer can resolve model names such as ERK2, PKCA, etc.
  tar_target(
    name = kinase_hierarchy,
    command = combine_kinase_taxonomies(
      kinase_taxonomy = kinase_taxonomy,
      gps6_hierarchy = gps6_hierarchy
    )
  ),

  # Export the entire validated chip to one FASTA file.
  #
  # The former <=9-sequence batching was only needed for the unstable web server.
  tar_target(
    name = chip_fasta_file,
    command = export_chip_fasta(
      sites = validated_chip_sites,
      output_path = file.path(
        "data",
        "external",
        "gps6",
        "chip_peptides.fasta"
      )
    ),
    format = "file"
  ),

  # Run local GPS6 across every complete runnable model node when model_scope
  # is ALL. The Python runner retains the full GPS6 hierarchy as provenance.
  tar_target(
    name = gps6_raw_output,
    command = {
      gps6_script
      gps6_model_files
      chip_fasta_file

      python_exe <-
        manifest$sources$gps6$python %||%
        Sys.getenv(
          "RETICULATE_PYTHON",
          unset = "python"
        )

      threshold <-
        manifest$sources$gps6$threshold %||%
        "h"

      model_scope <-
        manifest$sources$gps6$model_scope %||%
        "ALL"

      pre_dir <-
        manifest$sources$gps6$pre_dir %||%
        "data/pre"

      out <- file.path(
        "data",
        "external",
        "gps6",
        "gps6_predictions.tsv"
      )

      dir.create(
        dirname(out),
        recursive = TRUE,
        showWarnings = FALSE
      )

      status <- system2(
        command = python_exe,
        args = c(
          "-u",
          gps6_script,
          threshold,
          model_scope,
          chip_fasta_file,
          out,
          "--pre-dir",
          pre_dir
        ),
        stdout = "",
        stderr = ""
      )

      if (!identical(status, 0L)) {
        stop(
          "GPS6 local prediction failed with exit status ",
          status,
          "."
        )
      }

      if (!file.exists(out)) {
        stop(
          "GPS6 completed without creating the expected output file: ",
          out
        )
      }

      out
    },
    format = "file"
  ),

  # Parse raw GPS6 output into a canonical evidence table while preserving:
  # peptide ID, phosphosite, full GPS6 node, score, cutoff, and support flag.
  tar_target(
    name = gps6_processed,
    command = process_gps6_mappings(
      gps6_file = gps6_raw_output,
      target_sites = validated_chip_sites,
      manifest = manifest
    )
  ),

  # ============================================================================
  # Step G: Harmonization
  # ============================================================================

  # `kinase_hierarchy` combines:
  #
  #   1. KinHub/OpenKinome canonical human gene -> family -> superfamily
  #   2. GPS6 source-native aliases needed to interpret GPS6 model names
  #
  # The current harmonizer already accepts this hierarchy schema through its
  # `gps6_hierarchy` argument, so no additional harmonizer API is required.
  tar_target(
    name = harmonized_evidence,
    command = harmonize_kinase_evidence(
      signor_data = signor_processed,
      kinase_library_data = kinase_library_processed,
      gps6_data = gps6_processed,
      gps6_hierarchy = kinase_hierarchy
    )
  ),

  # ============================================================================
  # Step H: Master Evidence Table
  # ============================================================================

  # Preserve all source-native evidence before family collapse, including:
  #
  #   peptide_id
  #   substrate_uniprot
  #   phosphosite
  #   kinase_gene
  #   kinase_uniprot
  #   kinase_family
  #   kinase_superfamily
  #
  #   kinase_library_percentile
  #   kinase_library_supported
  #
  #   gps6_node
  #   gps6_score
  #   gps6_cutoff
  #   gps6_supported
  #
  #   signor_supported
  #   signor_pmid
  #
  #   evidence_class
  #   evidence_sources
  tar_target(
    name = master_evidence,
    command = build_master_evidence_table(
      chip_sites = validated_chip_sites,
      harmonized_evidence = harmonized_evidence,
      manifest = manifest
    )
  ),

  # ============================================================================
  # Step I: Family-Level Collapse
  # ============================================================================

  tar_target(
    name = family_evidence,
    command = collapse_to_kinase_family(
      evidence = master_evidence,
      manifest = manifest
    )
  ),

  # ============================================================================
  # Step J: KRSA-Compatible Mapping Generation
  # ============================================================================

  # Primary family-level mapping. Inclusion logic is defined in manifest.yml.
  tar_target(
    name = family_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "primary",
      manifest = manifest
    )
  ),

  # PhosphoSIGNOR-only sensitivity mapping.
  tar_target(
    name = experimental_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "experimental_only",
      manifest = manifest
    )
  ),

  # Kinase-Library-only sensitivity mapping.
  tar_target(
    name = kinase_library_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "kinase_library_only",
      manifest = manifest
    )
  ),

  # GPS6-only sensitivity mapping.
  tar_target(
    name = gps6_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "gps6_only",
      manifest = manifest
    )
  ),

  # Requires support from BOTH predictive sources: Kinase Library and GPS6.
  tar_target(
    name = predictive_concordant_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "predictive_concordant",
      manifest = manifest
    )
  ),

  # General multi-source concordance mapping. The exact rule should be defined
  # in manifest.yml / build_krsa_mapping(), e.g. >=2 distinct evidence sources.
  tar_target(
    name = concordant_mapping,
    command = build_krsa_mapping(
      family_evidence = family_evidence,
      mapping_type = "concordant",
      manifest = manifest
    )
  ),

  # ============================================================================
  # Step K: Quality Control & Coverage
  # ============================================================================

    tar_target(
    name = qc_summary,
    command = perform_mapping_qc(
      chip_sites = validated_chip_sites,
      master_evidence = master_evidence,
      family_evidence = family_evidence,
      family_mapping = family_mapping,
      experimental_mapping = experimental_mapping,
      kinase_library_mapping = kinase_library_mapping,
      gps6_mapping = gps6_mapping,
      predictive_concordant_mapping = predictive_concordant_mapping,
      concordant_mapping = concordant_mapping
    )
  ),

  tar_target(
    name = coverage_summary,
    command = calculate_mapping_coverage(
      chip_sites = validated_chip_sites,
      master_evidence = master_evidence,
      family_evidence = family_evidence
    )
  ),

  # ============================================================================
  # Step L: Output Assets
  # ============================================================================

  tar_target(
    name = output_files,
    command = write_mapping_outputs(
      chip_sites = validated_chip_sites,
      master_evidence = master_evidence,
      family_evidence = family_evidence,
      family_mapping = family_mapping,
      experimental_mapping = experimental_mapping,
      kinase_library_mapping = kinase_library_mapping,
      gps6_mapping = gps6_mapping,
      predictive_concordant_mapping = predictive_concordant_mapping,
      concordant_mapping = concordant_mapping,
      qc_summary = qc_summary,
      coverage_summary = coverage_summary,
      manifest = manifest
    ),
    format = "file"
  ),

  # ============================================================================
  # Step M: Quarto HTML Report
  # ============================================================================

  tar_quarto(
    name = mapping_report,
    path = "report/mapping_report.qmd",
    extra_files = c("config/manifest.yml"),
    quiet = FALSE
  )
)
