#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# Split chip-layout peptide sequences into FASTA files containing <= 9 peptides
#
# Usage:
#   bash scripts/make_fasta_batches.sh data/raw/87202-Array-Layout.csv
#
# Optional output directory:
#   bash scripts/make_fasta_batches.sh \
#       data/raw/87202-Array-Layout.csv \
#       data/fasta_batches
# ==============================================================================

INPUT_CSV="${1:-}"
OUTPUT_DIR="${2:-data/fasta_batches}"

# Change these if your CSV uses different column names.
PEPTIDE_ID_COLUMN="ID"
SEQUENCE_COLUMN="Sequence"

MAX_PER_FASTA=9

# ------------------------------------------------------------------------------
# Validate input
# ------------------------------------------------------------------------------

if [[ -z "${INPUT_CSV}" ]]; then
    echo "ERROR: No input CSV supplied."
    echo
    echo "Usage:"
    echo "  bash $0 <chip_layout.csv> [output_directory]"
    exit 1
fi

if [[ ! -f "${INPUT_CSV}" ]]; then
    echo "ERROR: Input file does not exist:"
    echo "  ${INPUT_CSV}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Remove FASTA files from previous runs so stale batches are not retained.
rm -f "${OUTPUT_DIR}"/*.fasta

# ------------------------------------------------------------------------------
# Create FASTA batches
# ------------------------------------------------------------------------------

python3 - \
    "${INPUT_CSV}" \
    "${OUTPUT_DIR}" \
    "${PEPTIDE_ID_COLUMN}" \
    "${SEQUENCE_COLUMN}" \
    "${MAX_PER_FASTA}" <<'PYTHON'

import csv
import os
import re
import sys

input_csv = sys.argv[1]
output_dir = sys.argv[2]
id_column = sys.argv[3]
sequence_column = sys.argv[4]
max_per_file = int(sys.argv[5])


def clean_sequence(sequence):
    """
    Standardize peptide sequence for FASTA output.

    - removes whitespace
    - converts to uppercase
    - removes common phosphorylation markers such as '*'
    """
    sequence = sequence.strip().upper()
    sequence = re.sub(r"\s+", "", sequence)
    sequence = sequence.replace("*", "")
    return sequence


def clean_identifier(identifier):
    """
    Make the peptide identifier safe for a FASTA header.
    """
    identifier = str(identifier).strip()
    identifier = re.sub(r"\s+", "_", identifier)
    return identifier


with open(input_csv, newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)

    if reader.fieldnames is None:
        raise RuntimeError("CSV contains no header row.")

    missing = [
        col
        for col in (id_column, sequence_column)
        if col not in reader.fieldnames
    ]

    if missing:
        print("ERROR: Required CSV column(s) not found:", file=sys.stderr)

        for col in missing:
            print(f"  - {col}", file=sys.stderr)

        print("\nAvailable columns:", file=sys.stderr)

        for col in reader.fieldnames:
            print(f"  - {col}", file=sys.stderr)

        sys.exit(1)

    peptides = []

    for row_number, row in enumerate(reader, start=2):
        peptide_id = clean_identifier(row.get(id_column, ""))
        sequence = clean_sequence(row.get(sequence_column, ""))

        if not peptide_id:
            print(
                f"WARNING: skipping row {row_number}: missing peptide ID",
                file=sys.stderr
            )
            continue

        if not sequence:
            print(
                f"WARNING: skipping {peptide_id}: missing sequence",
                file=sys.stderr
            )
            continue

        if not re.fullmatch(r"[A-Z]+", sequence):
            print(
                f"WARNING: skipping {peptide_id}: "
                f"sequence contains unexpected characters: {sequence}",
                file=sys.stderr
            )
            continue

        peptides.append((peptide_id, sequence))


if not peptides:
    raise RuntimeError("No valid peptide sequences were found.")


# --------------------------------------------------------------------------
# Write batches
# --------------------------------------------------------------------------

batch_number = 0

for start in range(0, len(peptides), max_per_file):
    batch_number += 1
    batch = peptides[start:start + max_per_file]

    output_file = os.path.join(
        output_dir,
        f"peptides_batch_{batch_number:02d}.fasta"
    )

    with open(output_file, "w", encoding="utf-8") as out:
        for peptide_id, sequence in batch:
            out.write(f">{peptide_id}\n")
            out.write(f"{sequence}\n")

    print(
        f"Wrote {output_file}: "
        f"{len(batch)} peptide{'s' if len(batch) != 1 else ''}"
    )


print()
print(f"Total peptides: {len(peptides)}")
print(f"Total FASTA files: {batch_number}")
print(f"Maximum peptides per FASTA: {max_per_file}")

PYTHON

echo
echo "FASTA generation complete."
echo "Output directory: ${OUTPUT_DIR}"
