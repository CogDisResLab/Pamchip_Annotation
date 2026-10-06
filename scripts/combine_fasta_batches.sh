#!/bin/bash

# Define the output file name
OUTPUT_FILE="combined.fasta"

# Remove the output file if it already exists to prevent duplicate entries on multiple runs
rm -f "$OUTPUT_FILE"

# Concatenate all FASTA files in the current folder into the output file
cat *.fasta > "$OUTPUT_FILE"

# Count the number of headers (lines starting with '>') in the combined file
HEADER_COUNT=$(grep -c "^>" "$OUTPUT_FILE")

echo "Successfully combined files into: $OUTPUT_FILE"
echo "Total number of peptides/sequences: $HEADER_COUNT"