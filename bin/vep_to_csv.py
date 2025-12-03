#!/usr/bin/env python3
"""
vep_to_csv.py

Convert a VEP-annotated VCF (with CSQ field) into a simple flat CSV table
for downstream inspection (e.g. Excel, R, Python).

For each variant, we output:
- Basic VCF fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER
- Selected CSQ / VEP fields for the FIRST transcript annotation

This is for training / portfolio use only.
"""

import argparse
import csv
import gzip
from pathlib import Path
from typing import List, Dict, Optional


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def open_maybe_gzip(path: Path):
    """
    Open a file that may be gzipped (.gz) or plain text.
    Returns a text-mode file handle.
    """
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def parse_csq_header(lines: List[str]) -> List[str]:
    """
    Parse the CSQ header line to extract the field order.

    Looks for a line starting with '##INFO=<ID=CSQ' and containing 'Format: '.
    The format string after 'Format: ' is a pipe-delimited list of field names.

    Returns:
        A list of CSQ field names in order, or an empty list if not found.
    """
    for line in lines:
        if line.startswith("##INFO=<ID=CSQ"):
            # Example:
            # ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...">
            if "Format: " not in line:
                continue
            fmt = line.split("Format: ")[1].rstrip(">\n").strip().strip('"')
            fields = fmt.split("|")
            return fields
    return []


def extract_info_field(info_str: str, key: str) -> Optional[str]:
    """
    From an INFO string like 'AC=1;AF=0.5;CSQ=..', extract the value for a key.

    Returns:
        The value string (without 'key=') or None if not found.
    """
    for entry in info_str.split(";"):
        if entry.startswith(key + "="):
            return entry.split("=", 1)[1]
    return None


# ---------------------------------------------------------------------------
# Main conversion logic
# ---------------------------------------------------------------------------

def vep_vcf_to_csv(vep_vcf: Path, out_csv: Path) -> None:
    """
    Read a VEP-annotated VCF and write a CSV file with selected fields.
    """

    # Collect header lines first to parse CSQ format
    header_lines: List[str] = []
    with open_maybe_gzip(vep_vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line)
            else:
                # Stop once we hit the column header line (#CHROM ...)
                break

    csq_fields = parse_csq_header(header_lines)

    # Choose a subset of VEP fields to export.
    # These must exist in csq_fields; if they don't, they will just be blank.
    selected_csq_fields = [
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature",
        "BIOTYPE",
        "EXON",
        "HGVSc",
        "HGVSp",
        "Existing_variation",
    ]

    # Build mapping from CSQ field name -> index
    csq_index: Dict[str, int] = {name: i for i, name in enumerate(csq_fields)}

    # Prepare CSV header
    base_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    header = base_columns + selected_csq_fields

    with open_maybe_gzip(vep_vcf) as fh_in, out_csv.open("w", newline="") as fh_out:
        writer = csv.writer(fh_out)
        writer.writerow(header)

        for line in fh_in:
            # Skip header lines
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue  # malformed line

            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]

            # Basic VCF fields
            base_values = [chrom, pos, vid, ref, alt, qual, filt]

            # Default CSQ-based values are empty
            csq_values_out = [""] * len(selected_csq_fields)

            csq_value = extract_info_field(info, "CSQ")
            if csq_value:
                # Take only the FIRST transcript annotation
                first_anno = csq_value.split(",")[0]
                csq_parts = first_anno.split("|")

                for j, field_name in enumerate(selected_csq_fields):
                    idx = csq_index.get(field_name)
                    if idx is not None and idx < len(csq_parts):
                        csq_values_out[j] = csq_parts[idx]

            writer.writerow(base_values + csq_values_out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert a VEP-annotated VCF (with CSQ) to CSV."
    )
    parser.add_argument(
        "--vep_vcf",
        required=True,
        type=Path,
        help="Path to VEP-annotated VCF (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Path to output CSV file",
    )
    args = parser.parse_args()

    vep_vcf_to_csv(args.vep_vcf, args.out)


if __name__ == "__main__":
    main()
