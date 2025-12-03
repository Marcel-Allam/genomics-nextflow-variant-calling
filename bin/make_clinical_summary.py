#!/usr/bin/env python3
"""
make_clinical_summary.py

Generate a simplified clinical-style Markdown report using:

- VEP-annotated VCF (for high-impact variant summary)
- bcftools stats on the RAW VCF (for variant QC)
- samtools flagstat / stats (for alignment QC)
- samtools coverage (for coverage QC)

This is for training / portfolio purposes only and is NOT validated
for clinical decision-making.
"""

import argparse
import gzip
import re
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple


# ---------------------------
# Parsing helper functions
# ---------------------------

def parse_bcftools_stats(path: Path) -> Dict[str, Optional[float]]:
    """
    Parse bcftools stats output and extract:
    - total_variants
    - snps
    - indels
    - ts_tv_ratio
    """
    stats = {
        "total_variants": None,
        "snps": None,
        "indels": None,
        "ts_tv_ratio": None,
    }

    if not path.is_file():
        return stats

    with path.open() as fh:
        for line in fh:
            if line.startswith("SN"):
                # "SN\t0\tnumber of records:\t1234"
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    label = parts[2].strip().lower()
                    value = parts[3].strip()
                    if "number of records" in label:
                        stats["total_variants"] = float(value)
                    elif "number of snps" in label:
                        stats["snps"] = float(value)
                    elif "number of indels" in label:
                        stats["indels"] = float(value)

            elif line.startswith("TSTV"):
                # "TSTV\tall\tTS:TV\tts_count\t tv_count\t ratio"
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    try:
                        stats["ts_tv_ratio"] = float(parts[-1])
                    except ValueError:
                        pass

    return stats


def parse_flagstat(path: Path) -> Dict[str, Optional[float]]:
    """
    Parse samtools flagstat output to extract:
    - total_reads
    - mapped_percent
    - properly_paired_percent
    """
    stats = {
        "total_reads": None,
        "mapped_percent": None,
        "properly_paired_percent": None,
    }

    if not path.is_file():
        return stats

    with path.open() as fh:
        for line in fh:
            # Example lines:
            # "1234 + 0 in total (QC-passed reads + QC-failed reads)"
            # "1200 + 0 mapped (97.24% : N/A)"
            # "1180 + 0 properly paired (95.61% : N/A)"
            line = line.strip()

            if "in total" in line and "QC-passed" in line:
                # Extract total reads at start of line
                try:
                    total = int(line.split()[0])
                    stats["total_reads"] = float(total)
                except ValueError:
                    pass

            elif "mapped (" in line:
                # Extract percent inside parentheses
                m = re.search(r"\(([\d\.]+)%", line)
                if m:
                    stats["mapped_percent"] = float(m.group(1))

            elif "properly paired (" in line:
                m = re.search(r"\(([\d\.]+)%", line)
                if m:
                    stats["properly_paired_percent"] = float(m.group(1))

    return stats


def parse_samtools_stats_for_duplicates(path: Path) -> Dict[str, Optional[float]]:
    """
    Parse samtools stats output to get approximate duplicate counts if available.
    This is optional; if not found, values remain None.
    """
    stats = {
        "reads_duplicated": None,
    }

    if not path.is_file():
        return stats

    with path.open() as fh:
        for line in fh:
            if line.startswith("SN") and "reads duplicated:" in line:
                # SN 0 reads duplicated: 123
                parts = line.strip().split(":")
                try:
                    dup = int(parts[-1].strip())
                    stats["reads_duplicated"] = float(dup)
                except ValueError:
                    pass

    return stats


def parse_coverage(path: Path) -> Dict[str, Optional[float]]:
    """
    Parse samtools coverage output to estimate:
    - mean_coverage (avg of 'meandepth' column)
    - mean_breadth (avg of 'coverage' column, fraction 0-1)

    This is a simplification but sufficient for demonstration.
    """
    stats = {
        "mean_coverage": None,
        "mean_breadth": None,
    }

    if not path.is_file():
        return stats

    total_depth = 0.0
    total_cov = 0.0
    n = 0

    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            try:
                cov_fraction = float(parts[5])  # "coverage" column
                mean_depth = float(parts[6])    # "meandepth" column
            except ValueError:
                continue

            total_cov += cov_fraction
            total_depth += mean_depth
            n += 1

    if n > 0:
        stats["mean_coverage"] = total_depth / n
        stats["mean_breadth"] = total_cov / n

    return stats


def parse_vep_impact(vep_vcf: Path, max_records: int = 5000) -> Dict[str, Optional[int]]:
    """
    Parse a VEP-annotated VCF to count HIGH and MODERATE impact variants.

    To keep it simple and fast:
    - We parse the CSQ header to find the index of the IMPACT field.
    - Then scan up to `max_records` variant lines and count IMPACT=HIGH/MODERATE.

    For larger cohorts, this could be expanded or optimised.
    """
    stats = {
        "high_impact_variants": 0,
        "moderate_impact_variants": 0,
        "records_scanned": 0,
    }

    if not vep_vcf.is_file():
        return stats
    
def get_tool_versions() -> Dict[str, str]:
    """
    Try to capture versions of key tools for audit-style reporting.
    If any command fails or is unavailable, 'NA' is used.
    """
    versions: Dict[str, str] = {}

    def run(cmd, key):
        try:
            out = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
            first_line = out.strip().splitlines()[0]
            versions[key] = first_line.strip()
        except Exception:
            versions[key] = "NA"

    run(["nextflow", "-version"], "nextflow")
    run(["samtools", "--version"], "samtools")
    run(["bcftools", "--version"], "bcftools")
    run(["vep", "--version"], "vep")

    return versions



    # Determine if file is gzipped
    open_func = gzip.open if vep_vcf.suffix == ".gz" else open

    impact_index = None

    with open_func(vep_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                # Extract the FORMAT description
                # e.g. ... Format: Allele|Consequence|IMPACT|...
                m = re.search(r"Format: ([^\">]+)", line)
                if m:
                    fields = m.group(1).strip().split("|")
                    for i, f in enumerate(fields):
                        if f.strip().upper() == "IMPACT":
                            impact_index = i
                            break
            elif line.startswith("#"):
                # Other header lines
                continue
            else:
                # Variant line
                if impact_index is None:
                    # No CSQ header parsed; cannot count impact reliably
                    continue

                if stats["records_scanned"] >= max_records:
                    break

                stats["records_scanned"] += 1

                cols = line.strip().split("\t")
                if len(cols) < 8:
                    continue
                info = cols[7]

                csq_match = re.search(r"CSQ=([^;]+)", info)
                if not csq_match:
                    continue

                csq_entries = csq_match.group(1).split(",")
                for entry in csq_entries:
                    fields = entry.split("|")
                    if len(fields) <= impact_index:
                        continue
                    impact = fields[impact_index].upper()
                    if impact == "HIGH":
                        stats["high_impact_variants"] += 1
                    elif impact == "MODERATE":
                        stats["moderate_impact_variants"] += 1

    return stats


# ---------------------------
# QC decision logic
# ---------------------------

def qc_decision(
    align_qc: Dict[str, Optional[float]],
    cov_qc: Dict[str, Optional[float]],
    var_qc: Dict[str, Optional[float]],
) -> Tuple[str, str]:
    """
    Very simplified QC decision rules.

    Returns:
      (overall_status, explanation_string)
    """
    reasons = []
    status = "PASS"

    mapped = align_qc.get("mapped_percent")
    proper = align_qc.get("properly_paired_percent")
    mean_cov = cov_qc.get("mean_coverage")
    breadth = cov_qc.get("mean_breadth")
    tstv = var_qc.get("ts_tv_ratio")

    # Thresholds are illustrative, not validated.
    if mapped is not None and mapped < 90.0:
        status = "FAIL"
        reasons.append(f"Mapping rate below threshold: {mapped:.1f}% < 90%")

    if proper is not None and proper < 80.0:
        status = "FAIL"
        reasons.append(
            f"Properly paired reads below threshold: {proper:.1f}% < 80%"
        )

    if mean_cov is not None and mean_cov < 30.0:
        status = "FAIL"
        reasons.append(
            f"Mean coverage below threshold: {mean_cov:.1f}x < 30x"
        )

    if breadth is not None and breadth < 0.90:
        status = "FAIL"
        reasons.append(
            f"Breadth of coverage below threshold: {breadth*100:.1f}% < 90%"
        )

    if tstv is not None and tstv < 1.5:
        # For exome-like data, Ts/Tv ratio much lower than this can indicate issues.
        if status == "PASS":
            status = "WARN"
        reasons.append(
            f"Ts/Tv ratio is lower than expected: {tstv:.2f} < 1.5"
        )

    if not reasons:
        reasons.append("All key QC metrics met the example thresholds.")

    explanation = "\n".join(f"- {r}" for r in reasons)
    return status, explanation


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a simplified clinical-style Markdown summary."
    )
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--vep_vcf", required=True, help="VEP-annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("--vcf_stats", required=True, help="bcftools stats output for RAW VCF")
    parser.add_argument("--flagstat", required=True, help="samtools flagstat output")
    parser.add_argument("--align_stats", required=True, help="samtools stats output")
    parser.add_argument("--idxstats", required=True, help="samtools idxstats output")
    parser.add_argument("--coverage", required=True, help="samtools coverage output")
    parser.add_argument("--out", required=True, help="Output Markdown file")

    args = parser.parse_args()

    sample_id = args.sample

    vcf_stats_path = Path(args.vcf_stats)
    flagstat_path = Path(args.flagstat)
    align_stats_path = Path(args.align_stats)
    idxstats_path = Path(args.idxstats)
    coverage_path = Path(args.coverage)
    vep_vcf_path = Path(args.vep_vcf)

    # Parse QC files
    variant_qc = parse_bcftools_stats(vcf_stats_path)
    align_flagstat_qc = parse_flagstat(flagstat_path)
    align_dup_qc = parse_samtools_stats_for_duplicates(align_stats_path)
    coverage_qc = parse_coverage(coverage_path)

    # VEP impact parsing (defensive: ensure we always have a dict)
    vep_qc = parse_vep_impact(vep_vcf_path) or {
        "high_impact_variants": 0,
        "moderate_impact_variants": 0,
        "records_scanned": 0,
    }
    
    # Merge alignment QC dictionaries
    align_qc = {**align_flagstat_qc, **align_dup_qc}

    # Make QC decision
    overall_qc_status, qc_explanation = qc_decision(
        align_qc=align_qc,
        cov_qc=coverage_qc,
        var_qc=variant_qc,
    )

    # Build Markdown report
    out_path = Path(args.out)
    with out_path.open("w") as out:
        out.write(f"# Clinical-style Variant Summary – {sample_id}\n\n")
        out.write("> **For training / portfolio use only – NOT for clinical decision-making**\n\n")


        # SECTION 0: Summary Overview 
        out.write("## 0. Summary Overview\n\n")
        out.write(f"- **Overall QC Status:** {overall_qc_status}\n")
        out.write("- **Total reads:** " +
                  (f"{align_qc.get('total_reads'):.0f}" if align_qc.get('total_reads') is not None else "NA") + "\n")
        out.write("- **Mapped reads (%):** " +
                  (f"{align_qc.get('mapped_percent'):.2f}%" if align_qc.get('mapped_percent') is not None else "NA") + "\n")
        out.write("- **Properly paired (%):** " +
                  (f"{align_qc.get('properly_paired_percent'):.2f}%" if align_qc.get('properly_paired_percent') is not None else "NA") + "\n")
        out.write("- **Mean coverage:** " +
                  (f"{coverage_qc.get('mean_coverage'):.2f}x" if coverage_qc.get('mean_coverage') is not None else "NA") + "\n")
        out.write("- **Breadth ≥1x (%):** " +
                  (f"{coverage_qc.get('mean_breadth')*100:.2f}%" if coverage_qc.get('mean_breadth') is not None else "NA") + "\n")
        out.write("- **Total variants:** " +
                  (f"{variant_qc.get('total_variants'):.0f}" if variant_qc.get('total_variants') is not None else "NA") + "\n")
        out.write("- **Ts/Tv ratio:** " +
                  (f"{variant_qc.get('ts_tv_ratio'):.2f}" if variant_qc.get('ts_tv_ratio') is not None else "NA") + "\n")
        out.write("- **High-impact variants:** " +
                  f"{vep_qc.get('high_impact_variants')}" + "\n")
        out.write("\n---\n\n")

        # SECTION 1: Sample and pipeline overview
        out.write("## 1. Sample and Pipeline Overview\n\n")
        out.write(f"- **Sample ID:** `{sample_id}`\n")
        out.write("- **Reference genome:** GRCh38 (as configured in the pipeline)\n")
        out.write("- **Variant caller:** bcftools mpileup + call\n")
        out.write("- **Annotation:** Ensembl VEP (offline cache, GRCh38)\n\n")

        # SECTION 2: Overall QC status
        out.write("## 2. Overall QC Status\n\n")
        out.write(f"- **QC Status:** **{overall_qc_status}**\n\n")
        out.write("**QC Rationale:**\n\n")
        out.write(qc_explanation + "\n\n")

        # SECTION 3: Alignment QC
        out.write("## 3. Alignment QC (samtools)\n\n")
        out.write("| Metric | Value |\n")
        out.write("|--------|-------|\n")
        out.write(f"| Total reads | {align_qc.get('total_reads') if align_qc.get('total_reads') is not None else 'NA'} |\n")
        out.write(f"| Mapped reads (%) | {align_qc.get('mapped_percent'):.2f}% |\n" if align_qc.get("mapped_percent") is not None else "| Mapped reads (%) | NA |\n")
        out.write(f"| Properly paired reads (%) | {align_qc.get('properly_paired_percent'):.2f}% |\n" if align_qc.get("properly_paired_percent") is not None else "| Properly paired reads (%) | NA |\n")
        out.write(f"| Reads duplicated | {align_qc.get('reads_duplicated') if align_qc.get('reads_duplicated') is not None else 'NA'} |\n")
        out.write("\n")

        # SECTION 4: Coverage QC
        out.write("## 4. Coverage QC (samtools coverage)\n\n")
        out.write("| Metric | Value |\n")
        out.write("|--------|-------|\n")
        if coverage_qc.get("mean_coverage") is not None:
            out.write(f"| Mean coverage (x) | {coverage_qc['mean_coverage']:.2f} |\n")
        else:
            out.write("| Mean coverage (x) | NA |\n")

        if coverage_qc.get("mean_breadth") is not None:
            out.write(f"| Mean breadth of coverage (%) | {coverage_qc['mean_breadth']*100:.2f}% |\n")
        else:
            out.write("| Mean breadth of coverage (%) | NA |\n")

        out.write("\n")
        out.write("_Note: Coverage metrics are simplified for demonstration and do not replace formal panel/region-level coverage analysis._\n\n")

        # SECTION 5: Variant QC (RAW VCF – bcftools stats)\n
        out.write("## 5. Variant QC (bcftools stats on RAW VCF)\n\n")
        out.write("| Metric | Value |\n")
        out.write("|--------|-------|\n")
        out.write(f"| Total variants | {variant_qc.get('total_variants') if variant_qc.get('total_variants') is not None else 'NA'} |\n")
        out.write(f"| SNPs | {variant_qc.get('snps') if variant_qc.get('snps') is not None else 'NA'} |\n")
        out.write(f"| INDELs | {variant_qc.get('indels') if variant_qc.get('indels') is not None else 'NA'} |\n")
        if variant_qc.get("ts_tv_ratio") is not None:
            out.write(f"| Ts/Tv ratio | {variant_qc['ts_tv_ratio']:.2f} |\n")
        else:
            out.write("| Ts/Tv ratio | NA |\n")
        out.write("\n")

        # SECTION 6: VEP Annotation Summary
        out.write("## 6. VEP Annotation Summary (approximate)\n\n")
        out.write("| Metric | Value |\n")
        out.write("|--------|-------|\n")
        out.write(f"| High-impact variants (IMPACT=HIGH, first {vep_qc['records_scanned']} records scanned) | {vep_qc['high_impact_variants']} |\n")
        out.write(f"| Moderate-impact variants (IMPACT=MODERATE, first {vep_qc['records_scanned']} records scanned) | {vep_qc['moderate_impact_variants']} |\n")
        out.write("\n")
        out.write("_Note: Impact counts are based on scanning a subset of variants for demonstration purposes._\n\n")

        # SECTION 7: Interpretation Notes / Disclaimer
        out.write("## 7. Interpretation Notes and Disclaimer\n\n")
        out.write("- This report is automatically generated by a training pipeline.\n")
        out.write("- It is intended to demonstrate bioinformatics workflow design only.\n")
        out.write("- It has not been validated for clinical use and must **not** be used\n")
        out.write("  to guide patient management or clinical decision-making.\n")

    print(f"Clinical-style summary written to: {out_path}")
    

if __name__ == '__main__':
    main()
