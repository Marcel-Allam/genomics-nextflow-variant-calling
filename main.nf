// main.nf
//
// Genomics – Variant Calling Pipeline (Nextflow DSL2)
//
// Overview:
//   For each sample (paired-end FASTQ), this pipeline:
//     1) Performs raw read QC with FastQC
//     2) Optionally trims reads using fastp (quality-based trimming)
//     3) Aligns reads to a reference genome with bwa-mem2
//     4) Sorts and indexes BAM files using samtools
//     5) Generates alignment and coverage QC metrics
//     6) Calls variants using bcftools mpileup + call
//     7) Runs variant-level QC with bcftools stats
//     8) Annotates variants using VEP (offline cache, GRCh38)
//     9) Exports annotated variants to CSV
//    10) Builds a clinical-style Markdown/HTML report per sample
//    11) Aggregates QC metrics using MultiQC
//
// Notes:
//   - This is a portfolio pipeline, not a clinical pipeline.
//   - It demonstrates a realistic structure similar to pipelines used in
//     cancer / WES workflows.
//
// Requirements:
//   - Input FASTQ(s) under: data/fastq/*_{R1,R2}.fastq.gz
//   - Reference FASTA at:   data/reference/genome.fa
//   - Conda env in:         envs/genomics_env.yml
//   - VEP offline cache installed in: ~/.vep (e.g. homo_sapiens/GRCh38)

nextflow.enable.dsl = 2

// -----------------------------
// PARAMETERS (with fallbacks)
// -----------------------------

// Paired-end reads glob
params.reads     = params.reads     ?: "data/fastq/*_{R1,R2}.fastq.gz"

// Reference FASTA (must be pre-indexed for bwa-mem2 + samtools + bcftools)
params.reference = params.reference ?: "data/reference/genome.fa"

// Output directory
params.outdir    = params.outdir    ?: "results"

// Trimming behaviour
params.trim_reads = params.trim_reads instanceof Boolean ? params.trim_reads : false
params.trim_qual  = params.trim_qual ?: 20  // 20 = standard, 30 = strict

// Threads per process
params.cpus       = params.cpus ?: 4

// Ensure main output directory exists
new File(params.outdir).mkdirs()
// Ensure main output directory exists
new File(params.outdir).mkdirs()

// -----------------------------
// CHANNELS
// -----------------------------
//
// 1) Paired-end FASTQ files
//
// This creates tuples of:
//   ( sample_id, [R1.fastq.gz, R2.fastq.gz] )
// based on the pattern defined in params.reads.

Channel
    .fromFilePairs( params.reads, flat: false )
    .set { read_pairs_ch }

// -----------------------------
// PROCESSES
// -----------------------------

// -----------------------------------------------------------------------------
// Process: FASTQC
// Purpose:
//   Perform raw read quality control using FastQC. This step provides an initial
//   assessment of sequencing quality, adapter contamination, GC distribution,
//   and read length profiles. Raw QC is always run prior to trimming or alignment.
//
// Notes:
//   - Output is published to results/fastqc for direct inspection.
//   - No trimming is performed here; modern aligners handle soft-clipping.
// -----------------------------------------------------------------------------
process FASTQC {

    // Label task with sample ID for clear logging and traceability.
    tag "${sample_id}"

    // Export FastQC reports (HTML + zipped data) to a structured output dir.
    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

    input:
    // Paired FASTQ files grouped as: (sample_id, [R1.fastq.gz, R2.fastq.gz])
    tuple val(sample_id), path(reads)

    output:
    // FastQC emits multiple files; wildcard captures both HTML and .zip outputs.
    path "*_fastqc.*"

    """
    # Run FastQC on both paired-end read files. FastQC automatically identifies
    # per-base quality issues, adapter content, and sequence duplication levels.
    fastqc \
        --outdir . \
        ${reads.join(' ')}
    """
}

// -----------------------------------------------------------------------------
// Process: TRIM_READS
// Purpose:
//   Conditionally trim paired-end FASTQ files using fastp when
//   `--trim_reads true` is specified. Trimming is *not* applied by default,
//   consistent with modern variant-calling workflows where aligner soft-clipping
//   typically handles low-quality tails and adapter contamination.
//
// What this step performs:
//   - Automatic adapter detection for paired-end reads
//   - Removal of low-quality bases at 5' and 3' ends
//   - Optional strict quality enforcement via `--trim_qual`
//       Q20 ≈ 1% error rate (standard)
//       Q30 ≈ 0.1% error rate (strict)
//
// Notes:
//   - Trimming is intentionally user-controlled, as over-trimming can reduce
//     effective read length and coverage, impacting downstream variant calling.
//   - Output files are renamed consistently to <sample>_R1/R2.trimmed.fastq.gz.
// -----------------------------------------------------------------------------
process TRIM_READS {

    // Tag workflow execution logs with the sample ID.
    tag "${sample_id}"

    // Publish trimmed FASTQs to a dedicated results directory.
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    input:
    // Receive paired FASTQs packaged as: (sample_id, [R1, R2])
    tuple val(sample_id), path(reads)

    output:
    // Emit trimmed R1 and R2 FASTQs under consistent names.
    tuple val(sample_id), path("${sample_id}_R{1,2}.trimmed.fastq.gz")

    """
    # Perform adapter trimming and base-quality filtering with fastp.
    # The threshold is user-defined via --trim_qual (default Q20).
    # Higher thresholds (e.g., Q30) enforce stricter trimming at the cost
    # of shorter reads.
    fastp \
      --in1 ${reads[0]} \
      --in2 ${reads[1]} \
      --out1 ${sample_id}_R1.trimmed.fastq.gz \
      --out2 ${sample_id}_R2.trimmed.fastq.gz \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_front_window_size 4 \
      --cut_front_mean_quality ${params.trim_qual} \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality ${params.trim_qual} \
      --qualified_quality_phred ${params.trim_qual} \
      --thread ${params.cpus}
    """
}

// -----------------------------------------------------------------------------
// Process: ALIGN_AND_SORT
// Purpose:
//   Align paired-end reads to the reference genome using bwa-mem2, convert the
//   output to BAM, and produce a coordinate-sorted, indexed BAM suitable for
//   downstream QC and variant calling.
//
// Rationale:
//   - bwa-mem2 provides high-performance alignment optimised for WGS/WES data.
//   - All GRCh38 reference indexes (bwa-mem2, samtools, bcftools, picard)
//     are expected to be pre-built outside the workflow to avoid expensive
//     indexing steps inside Nextflow.
//   - Sorting and indexing at this stage ensures compatibility with tools such
//     as samtools, bcftools, and VEP, which require random-access BAM input.
//
// Notes:
//   - Trimming does not affect this step: both raw and trimmed reads are valid
//     inputs, depending on user configuration.
//   - Downstream QC modules (flagstat, stats, idxstats, coverage) consume the
//     sorted BAM emitted here.
// -----------------------------------------------------------------------------
process ALIGN_AND_SORT {

    // Label run logs and trace files by sample for clarity.
    tag "${sample_id}"

    // Publish sorted BAM + index to a dedicated results directory.
    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true

    input:
    // Receive paired reads (either raw or trimmed): (sample_id, [R1, R2])
    tuple val(sample_id), path(reads)

    output:
    // Emit a coordinate-sorted BAM (<sample>.sorted.bam)
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    """
    # Perform alignment with bwa-mem2 using the pre-indexed GRCh38/chr22 reference.
    # Output is piped directly into samtools to avoid intermediate file creation.
    bwa-mem2 mem -t ${params.cpus} ${params.reference} ${reads[0]} ${reads[1]} \
        | samtools view -bS - \
        | samtools sort -o ${sample_id}.sorted.bam -

    # Index the BAM for downstream tools.
    samtools index ${sample_id}.sorted.bam
    """
}

// -----------------------------------------------------------------------------
// Process: ALIGN_QC
// Purpose:
//   Generate alignment-level QC metrics from the coordinate-sorted BAM produced
//   during read alignment. These metrics provide insight into mapping quality,
//   read pairing, duplication, and potential contamination—all essential for
//   validating sequencing integrity in both research and clinical workflows.
//
// QC Outputs:
//   - <sample>.flagstat.txt       : High-level mapping summary (mapped %, paired %, duplicates).
//   - <sample>.samtools_stats.txt : Comprehensive statistics (read lengths, error rates, insert sizes).
//   - <sample>.idxstats.txt       : Per-contig read counts, useful for detecting
//                                   off-target alignment or reference contamination.
//
// Rationale:
//   - These QC metrics are routinely required in diagnostic labs and large-scale
//     sequencing projects for auditability and troubleshooting.
//   - Keeping this step modular allows MultiQC to aggregate results automatically.
// -----------------------------------------------------------------------------
process ALIGN_QC {

    // Tag logs/trace entries with the sample identifier.
    tag "${sample_id}"

    // Publish all QC artefacts to a dedicated results directory.
    publishDir "${params.outdir}/alignment_qc", mode: 'copy', overwrite: true

    input:
    // Receive sorted BAM from ALIGN_AND_SORT.
    tuple val(sample_id), path(sorted_bam)

    output:
    // Emit the three QC output files as a structured tuple.
    tuple val(sample_id),
          path("${sample_id}.flagstat.txt"),
          path("${sample_id}.samtools_stats.txt"),
          path("${sample_id}.idxstats.txt")

    script:
    """
    # Generate a high-level alignment summary.
    samtools flagstat ${sorted_bam} > ${sample_id}.flagstat.txt

    # Produce detailed alignment metrics including error profiles and insert sizes.
    samtools stats ${sorted_bam} > ${sample_id}.samtools_stats.txt

    # Calculate per-contig read distribution to assess uniformity and contamination.
    samtools idxstats ${sorted_bam} > ${sample_id}.idxstats.txt
    """
}

// -----------------------------------------------------------------------------
// Process: COVERAGE_QC
// Purpose:
//   Compute coverage-based quality metrics from the coordinate-sorted BAM.
//   Coverage is one of the most critical QC dimensions in variant calling,
//   influencing sensitivity, callable regions, and downstream interpretation.
//
// Outputs:
//   - <sample>.coverage.txt : Summary depth and breadth metrics, including
//                             mean coverage and % of bases above key thresholds.
//   - <sample>.depth.txt    : Per-base depth across the reference. Useful for
//                             panel/WES uniformity checks but can be large for WGS.
//
// Rationale:
//   - Coverage uniformity and depth are primary drivers of variant detection power.
//   - samtools coverage provides fast, robust high-level metrics.
//   - samtools depth enables fine-grained inspection (e.g., dropout regions,
//     capture performance, or CNV signal).
// -----------------------------------------------------------------------------
process COVERAGE_QC {

    // Attach sample ID to execution logs for easy debugging and traceability.
    tag "${sample_id}"

    // Publish outputs to a dedicated QC directory for structured reporting.
    publishDir "${params.outdir}/coverage_qc", mode: 'copy', overwrite: true

    input:
    // Receive sorted BAM from ALIGN_AND_SORT.
    tuple val(sample_id), path(sorted_bam)

    // Emit sample ID and two core coverage QC artefacts:
    // - <sample>.coverage.txt : High-level summary metrics (mean depth, breadth)
    // - <sample>.depth.txt    : Base-level depth across the reference
    output:
    tuple val(sample_id),
          path("${sample_id}.coverage.txt"),
          path("${sample_id}.depth.txt")

    script:
    """
    # Generate high-level coverage summary: mean depth, breadth at thresholds,
    # and distribution across the reference.
    samtools coverage ${sorted_bam} > ${sample_id}.coverage.txt

    # Generate per-base depth. This provides fine-resolution coverage profiles
    # but produces large files for whole-genome datasets.
    samtools depth ${sorted_bam} > ${sample_id}.depth.txt
    """
}

// -----------------------------------------------------------------------------
// Process: CALL_VARIANTS
// Purpose:
//   Generate a raw, unfiltered VCF using bcftools mpileup + call. This provides
//   an initial set of SNP and INDEL candidates based on the aligned reads.
//   Variant annotation and downstream filtering steps consume this output.
//
// Rationale:
//   - bcftools mpileup computes genotype likelihoods from the BAM file.
//   - bcftools call performs multiallelic variant calling (-mv), producing a
//     compact bgzipped VCF suitable for indexing and further analysis.
//   - Using the reference directly via `params.reference` avoids unnecessary
//     channel duplication.
//
// Notes:
//   - This is a research-grade caller configuration appropriate for test datasets.
//   - Clinical-grade variant calling generally involves recalibration and
//     specialised models (e.g., GATK HaplotypeCaller, Strelka2, DeepVariant).
// -----------------------------------------------------------------------------
process CALL_VARIANTS {

    // Log and trace this process under the sample ID.
    tag "${sample_id}"

    // Publish raw VCF outputs to the standard VCF results directory.
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true

    input:
    // Receive coordinate-sorted BAM from ALIGN_AND_SORT.
    tuple val(sample_id), path(sorted_bam)

    output:
    // Emit a bgzipped raw VCF (<sample>.raw.vcf.gz) for annotation and QC.
    tuple val(sample_id), path("${sample_id}.raw.vcf.gz")

    """
    # Compute genotype likelihoods via mpileup and call SNPs/INDELs using bcftools.
    # -f supplies the reference FASTA (must be indexed with samtools).
    bcftools mpileup -f ${params.reference} ${sorted_bam} \
        | bcftools call -mv -Oz -o ${sample_id}.raw.vcf.gz

    # Index the VCF for random-access queries (required by bcftools and VEP).
    bcftools index ${sample_id}.raw.vcf.gz
    """
}

// -----------------------------------------------------------------------------
// Process: ANNOTATE_VEP
// Purpose:
//   Functionally annotate called variants using Ensembl VEP in offline mode.
//   This step enriches raw VCF records with gene-level, transcript-level, and
//   consequence-level information required for interpretation and reporting.
//
// Rationale:
//   - Offline mode with a local VEP cache (~/.vep) avoids repeated network calls
//     and ensures reproducibility across runs.
//   - Using VCF as both input and output format keeps the annotation compatible
//     with bcftools, downstream filters, and reporting layers.
//   - Explicit format specification (`--format vcf`) and streaming via STDOUT
//     avoids common VEP I/O issues and temporary large files.
//
// Assumptions:
//   - A valid VEP cache for `homo_sapiens` / `GRCh38` is installed under ~/.vep.
//   - The input VCF is bgzipped and indexed (produced by CALL_VARIANTS).
// -----------------------------------------------------------------------------
process ANNOTATE_VEP {

    // Attach sample ID to logs and trace output.
    tag "${sample_id}"

    // Publish annotated VCF back into the main VCF results directory.
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true

    // Receive the raw VCF from CALL_VARIANTS.
    input:
    tuple val(sample_id), path(vcf)

    // Emit the VEP-annotated variant call file:
    //   <sample>.annotated.vcf.gz
    //       A bgzipped VCF enriched with gene, transcript, consequence,
    //       population frequency, and functional impact annotations.
    //       Indexed for compatibility with bcftools and downstream reporting.
    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz")

    """
     # Expand the compressed VCF to a plain-text file for VEP consumption.
    zcat ${vcf} > input.vcf

    # Run Ensembl VEP in offline mode against the local cache.
    # --everything enables a comprehensive annotation set (gene, consequence,
    #   protein change, population frequencies, and more), which is useful for
    #   exploratory and training pipelines.
    vep \
        --format vcf \
        --input_file input.vcf \
        --output_file STDOUT \
        --vcf \
        --cache \
        --offline \
        --species homo_sapiens \
        --assembly GRCh38 \
        --everything \
        --force_overwrite \
        --quiet \
        | bgzip -c > ${sample_id}.annotated.vcf.gz

    # Index the annotated VCF for efficient querying and downstream filtering.
    bcftools index ${sample_id}.annotated.vcf.gz
    """
}

// -----------------------------------------------------------------------------
// Process: EXPORT_VARIANTS_CSV
// Purpose:
//   Convert the VEP-annotated VCF into a flat, tabular CSV containing key
//   variant-level annotations. This representation is easier to review in
//   downstream analytical tools such as Python, R, or Excel, and is often used
//   in exploratory variant interpretation workflows.
//
// Rationale:
//   - VCF is a compact, structured format optimised for genomic tooling but is
//     difficult to inspect manually.
//   - Exporting to CSV enables rapid filtering, plotting, and integration with
//     reporting pipelines.
//   - A custom helper script (`vep_to_csv.py`) extracts relevant INFO fields,
//     standardises column names, and writes a tidy tabular file.
//
// Input:
//   - (sample_id, annotated_vcf_gz) emitted by ANNOTATE_VEP.
//
// Output:
//   - <sample>.variants.csv placed under results/variants_csv/
//     containing a clean, analysis-ready variant table.
// -----------------------------------------------------------------------------
process EXPORT_VARIANTS_CSV {

    // Tag logs and run metadata with the sample identifier.
    tag "${sample_id}"

    // Store all CSV exports in a structured output directory.
    publishDir "${params.outdir}/variants_csv", mode: 'copy', overwrite: true

    // Receive the annotated VCF from ANNOTATE_VEP.
    input:
    tuple val(sample_id), path(annotated_vcf_gz)

    // Emit a per-sample CSV summarising VEP annotations:
        //   <sample>.variants.csv
        //       A flattened, analysis-ready table containing variant position,
        //       gene annotations, predicted consequences, protein effects,
        //       population frequencies, and any custom VEP-derived fields.
    output:
    tuple val(sample_id), path("${sample_id}.variants.csv")

    script:
    """
    # Convert the VEP-annotated VCF to a tidy CSV using a helper script.
    # This ensures compatibility with downstream R / Python analysis pipelines.
    ${projectDir}/bin/vep_to_csv.py \\
        --vep_vcf ${annotated_vcf_gz} \\
        --out ${sample_id}.variants.csv
    """
}

// -----------------------------------------------------------------------------
// Process: VARIANT_QC
// Purpose:
//   Generate variant-level QC metrics from the raw, unannotated VCF using
//   bcftools stats. This step quantifies variant composition and distribution,
//   providing an essential quality assessment before annotation or filtering.
//
// Rationale:
//   - Variant QC should be applied to the *raw* VCF to capture properties of the
//     caller output without the influence of annotation layers.
//   - bcftools stats produces a standardised metrics file summarising SNP/INDEL
//     counts, transition/transversion ratio (Ts/Tv), heterozygosity, and
//     per-chromosome variant distribution—key indicators of dataset quality.
//   - These metrics are often used as acceptance criteria in research and
//     diagnostic environments.
//
// Inputs:
//   - sample_id  : unique sample identifier
//   - vcf_gz     : compressed raw VCF produced by CALL_VARIANTS
//
// Outputs:
//   - <sample>.vcfstats.txt : variant-level QC metrics used in reporting,
//                             MultiQC aggregation, and downstream interpretation.
// -----------------------------------------------------------------------------
process VARIANT_QC {

    // Track execution logs under the corresponding sample name.
    tag "${sample_id}"

    // Store QC outputs in a dedicated directory for downstream summarisation.
    publishDir "${params.outdir}/variant_qc", mode: 'copy', overwrite: true

    input:
    // Receive raw VCF from CALL_VARIANTS.
    tuple val(sample_id), path(vcf_gz)

    output:
    // Emit a structured QC summary:
        //   <sample>.vcfstats.txt
        //       Contains SNP/INDEL counts, Ts/Tv ratio, genotype tallies,
        //       per-chromosome variant distribution, and additional caller metrics.
        //       Designed for both manual review and automated MultiQC collection.
    tuple val(sample_id),
          path("${sample_id}.vcfstats.txt")

    script:
    """
    # Generate comprehensive variant QC metrics.
    # This file is included automatically in MultiQC reports.
    bcftools stats ${vcf_gz} > ${sample_id}.vcfstats.txt

    # Optional: enable visualisation of stat profiles via plot-vcfstats.
    # Disabled by default to avoid LaTeX/TeX dependencies in lightweight setups.
    #
    # mkdir -p ${sample_id}.vcfstats_plots
    # plot-vcfstats \\
    #   -p ${sample_id}.vcfstats_plots \\
    #   ${sample_id}.vcfstats.txt
    """
}

// -----------------------------------------------------------------------------
// Process: CLINICAL_SUMMARY
// Purpose:
//   Generate a consolidated QC and variant interpretation report in Markdown
//   and standalone HTML format. This mirrors a simplified clinical-style variant
//   summary integrating multiple QC layers and annotated variant data.
//
// Report includes:
//   - Alignment QC        : mapping quality, error profiles, read distribution
//   - Coverage QC         : depth and breadth metrics, per-base depth features
//   - Variant QC          : SNP/INDEL composition, Ts/Tv ratio, distribution
//   - VEP Annotation      : functional impact summaries, high-impact variants
//
// Rationale:
//   - Clinical and translational genomics workflows routinely compile integrated
//     QC + variant summaries for human review (e.g., MDT meetings, analyst sign-off).
//   - Generating a structured Markdown file enables flexible downstream rendering.
//   - An HTML report provides shareable, user-friendly visualisation.
//
// Notes:
//   - This workflow is designed for training and portfolio demonstration only,
//     not for diagnostic use. It omits required clinical validation and SOPs.
// -----------------------------------------------------------------------------
process CLINICAL_SUMMARY {

    // Attach sample ID to logs and execution metadata.
    tag "${sample_id}"

    // Publish Markdown + HTML report into a dedicated reporting directory.
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

    input:
    // Annotated VCF from ANNOTATE_VEP
    tuple val(sample_id), path(annotated_vcf_gz)

    // Variant QC from VARIANT_QC (bcftools stats on RAW VCF)
    tuple val(sample_id_qc), path(vcfstats_txt)

    // Alignment QC metrics from ALIGN_QC
    tuple val(sample_id_align),
          path(flagstat_txt),
          path(align_stats_txt),
          path(idxstats_txt)

    // Coverage QC metrics from COVERAGE_QC
    tuple val(sample_id_cov),
          path(coverage_txt),
          path(depth_txt)

    output:
     // Emit Markdown + HTML report:
        //   <sample>_clinical_summary.md
        //   <sample>_clinical_summary.html
        // These reports aggregate QC layers and annotated variant summaries
        // into a human-readable overview resembling clinical reporting formats.
    path("${sample_id}_clinical_summary.*")

    when:
    // Ensure all input sample identifiers match before generating report.
    sample_id == sample_id_qc &&
    sample_id == sample_id_align &&
    sample_id == sample_id_cov

       script:
    """
    # Generate a unified Markdown report summarising QC metrics and
    # annotated variant findings for downstream interpretation.
    ${projectDir}/bin/make_clinical_summary.py \
        --sample ${sample_id} \
        --vep_vcf ${annotated_vcf_gz} \
        --vcf_stats ${vcfstats_txt} \
        --flagstat ${flagstat_txt} \
        --align_stats ${align_stats_txt} \
        --idxstats ${idxstats_txt} \
        --coverage ${coverage_txt} \
        --out ${sample_id}_clinical_summary.md

    # Convert Markdown into a fully self-contained HTML document using pandoc
    # for portable review (e.g., for analysts, supervisors, or demonstration).
    pandoc ${sample_id}_clinical_summary.md \
        -f markdown -t html -s \
        -o ${sample_id}_clinical_summary.html
    """
}

// -----------------------------------------------------------------------------
// Process: MULTIQC
// Purpose:
//   Generate a unified MultiQC report by aggregating QC metrics produced across
//   the workflow (FastQC, alignment QC, coverage QC, variant QC, etc.).
//
// Rationale:
//   - MultiQC automatically discovers supported output formats within the
//     results directory and compiles them into an integrated HTML report.
//   - A dependency on the CLINICAL_SUMMARY output ensures this step executes
//     only after all upstream QC modules have completed.
//   - Consolidated QC reporting is standard practice in clinical, research,
//     and production sequencing pipelines.
//
// Inputs:
//   - At least one path from CLINICAL_SUMMARY is required solely as a
//     synchronization point. The contents themselves are not used by MultiQC.
//
// Output:
//   - multiqc_report.html : An aggregated QC dashboard stored under
//     results/multiqc/, suitable for human review and archiving.
// -----------------------------------------------------------------------------
process MULTIQC {

    // Use a fixed tag since this process aggregates multiple samples.
    tag "multiqc"

    // Place the final HTML report into a dedicated directory.
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    // Upstream dependency only; file contents not consumed by script.
    path(summary_reports)

    output:
    // Emit the MultiQC HTML dashboard:
    //   multiqc_report.html
    // Provides a consolidated visual summary of QC across all samples.
    path("multiqc_report.html")

    script:
    """
    # Run MultiQC over the full results directory rather than per-process
    # work directories. This ensures all QC artefacts are detected and included.
    multiqc ${projectDir}/${params.outdir} --outdir .
    """
}

// -----------------------------------------------------------------------------
// WORKFLOW (Nextflow DSL2)
// Purpose:
//   Define the orchestration of all pipeline processes. Channels from earlier
//   steps feed directly into downstream modules, forming a reproducible,
//   modular sequencing → alignment → QC → variant calling → annotation pipeline.
//
// Notes:
//   - Processes behave like pure functions: channels in → channels out.
//   - QC is executed in parallel wherever possible.
//   - Ordering is enforced implicitly through channel dependencies.
// -----------------------------------------------------------------------------

workflow {

    // -------------------------------------------------------------------------
    // 1) Raw FASTQ quality control (always executed)
    //    Provides initial metrics for assessing read quality prior to trimming
    //    or alignment. Raw QC is essential for determining whether trimming
    //    should be enabled.
    // -------------------------------------------------------------------------
    fastqc_raw_ch = FASTQC(read_pairs_ch)

    
    // -------------------------------------------------------------------------
    // 2) Conditional trimming (fastp)
    //    Trimming is applied *only when requested* via --trim_reads true.
    //    Modern aligners soft-clip adapters/low-quality tails, so trimming is
    //    optional and dataset-dependent.
    // -------------------------------------------------------------------------
    reads_for_alignment_ch = params.trim_reads ? TRIM_READS(read_pairs_ch)
                                               : read_pairs_ch

    
    // -------------------------------------------------------------------------
    // 3) Alignment and primary QC
    // -------------------------------------------------------------------------
    // 3a) Align reads to reference genome → sorted BAM.
    sorted_bam_ch = ALIGN_AND_SORT(reads_for_alignment_ch)

    // 3b) Generate alignment QC (flagstat, stats, idxstats).
    align_qc_ch = ALIGN_QC(sorted_bam_ch)

    // 3c) Coverage QC (samtools coverage / depth)
    coverage_qc_ch = COVERAGE_QC(sorted_bam_ch)


    // -------------------------------------------------------------------------
    // 4) Variant calling (bcftools)
    //    Produces a raw, unfiltered VCF suitable for QC and downstream annotation.
    // -------------------------------------------------------------------------
    raw_vcf_ch = CALL_VARIANTS(sorted_bam_ch)


    // -------------------------------------------------------------------------
    // 5) Variant QC and annotation
    // -------------------------------------------------------------------------
    
    // 5a) Variant-level QC on raw VCF (best practice before annotation).
    variant_qc_ch = VARIANT_QC(raw_vcf_ch)

    // 5b) Annotate variants with VEP (functional consequences, gene effects).
    annotated_vcf_ch = ANNOTATE_VEP(raw_vcf_ch)


    // -------------------------------------------------------------------------
    // 6) Export annotated variants to CSV (flat table)
    //    Useful for downstream R/Python analysis or human inspection.
    // -------------------------------------------------------------------------
    variants_csv_ch = EXPORT_VARIANTS_CSV(annotated_vcf_ch)


    // -------------------------------------------------------------------------
    // 7) Clinical-style QC + variant summary report
    //    Integrates:
    //       - alignment QC
    //       - coverage QC
    //       - variant QC
    //       - VEP annotation summaries
    //    Produces *.md and *.html reports.
    // -------------------------------------------------------------------------
     clinical_summary_ch = CLINICAL_SUMMARY(annotated_vcf_ch, variant_qc_ch, align_qc_ch, coverage_qc_ch)


    // -------------------------------------------------------------------------
    // 8) MultiQC aggregation
    //    Consolidates all QC outputs across samples into a single HTML report.
    // -------------------------------------------------------------------------
    multiqc_ch = MULTIQC(clinical_summary_ch)

 
    // -------------------------------------------------------------------------
    // Final Outputs Summary:
    //   - results/fastqc/         : FASTQC HTML and zip reports
    //   - results/trimmed/        : trimmed FASTQs (if trimming enabled)
    //   - results/bam/            : sorted BAM + BAM index
    //   - results/alignment_qc/   : flagstat, stats, idxstats
    //   - results/coverage_qc/    : depth + coverage metrics
    //   - results/vcf/            : raw + annotated VCFs
    //   - results/variant_qc/     : bcftools variant QC summaries
    //   - results/variants_csv/   : tabular variant summaries
    //   - results/reports/        : Markdown + HTML clinical-style reports
    //   - results/multiqc/        : aggregated MultiQC HTML dashboard
    // -------------------------------------------------------------------------
}
