# Data Flow Specification — HCV Quasispecies Pipeline

> Notation: `$SID` = sample_id, `$GT` = genotype label (e.g., `1a`, `2b`; per-branch outputs always use a single subtype).

## Step 5.1 — Sample input parsing
- **Module:** `INPUT_CHECK` (`modules/local/input_check.nf`)
- **Input:** `params.input` (CSV samplesheet)
- **Output:** Channel of `[meta, fastq_path]` tuples; `meta = [id: ..., metadata: ...]`
- **Tool:** `bin/check_samplesheet.py`
- **Pass:** Valid `sample_id` (regex `^[A-Za-z0-9._-]+$`), unique IDs, one FASTQ per sample, FASTQ readable.
- **Fail:** Any violation → workflow errors with descriptive message.
- **Resources:** 1 CPU, 1 GB, 5 min.

## Step 5.2 — Raw QC
- **Module:** `NANOPLOT_RAW`, `NANOQ_RAW`
- **Output:** `qc/raw/${SID}/NanoPlot-report.html`, `qc/raw/${SID}/nanoq.json`
- **Commands:**
  ```
  NanoPlot --fastq ${reads} -o nanoplot_raw --tsv_stats --no_static
  nanoq -i ${reads} --json -o ${SID}.nanoq.json --report
  ```
- **Resources:** 2 CPU, 4 GB, 15 min.
- **Fail:** Empty FASTQ → `EMPTY_INPUT` flag, skip sample.

## Step 5.3 — Read filtering
- **Module:** `CHOPPER`
- **Output:** `${SID}.filtered.fastq.gz`, `chopper.log`
- **Command:**
  ```
  zcat ${reads} | chopper -q ${params.min_qual} \
      --minlength ${params.min_length} --maxlength ${params.max_length} \
      --threads ${task.cpus} 2> chopper.log \
    | pigz -p ${task.cpus} > ${SID}.filtered.fastq.gz
  ```
- **Resources:** 4 CPU, 4 GB, 20 min.
- **Fail:** Output empty → `ALL_READS_FILTERED` flag, skip sample.

## Step 5.4 — Host depletion
- **Module:** `HOST_DEPLETE_MINIMAP2`
- **Output:** `${SID}.hostdep.fastq.gz`, `${SID}.host_stats.json`
- **Command:**
  ```
  minimap2 -ax map-ont -t ${task.cpus} ${host_index} ${reads} \
    | samtools view -@ ${task.cpus} -b -f 4 - \
    | samtools fastq -@ ${task.cpus} - \
    | pigz -p ${task.cpus} > ${SID}.hostdep.fastq.gz
  ```
- **Resources:** 16 CPU, 32 GB, 90 min.
- **Fail:** `%host ≥ 99.5%` → `NO_VIRAL_READS_LIKELY` warning; pipeline continues flagged.

## Step 5.5 — Post-host QC
- **Module:** `NANOQ_POSTHOST`
- **Output:** `qc/posthost/${SID}/nanoq.json`
- **Resources:** 1 CPU, 1 GB, 5 min.

## Step 5.6 — HCV panel indexing
- **Module:** `MINIMAP2_INDEX_PANEL`
- **Output:** `panel.mmi`, `panel.fasta`, `panel.fasta.fai`
- **Command:**
  ```
  minimap2 -x map-ont -t ${task.cpus} -d panel.mmi ${fasta}
  samtools faidx ${fasta}
  ```
- **Resources:** 4 CPU, 8 GB, 10 min.
- **Note:** Cached for whole run via Nextflow content hashing.

## Step 5.7 — Round 1 competitive mapping
- **Module:** `MINIMAP2_ROUND1`
- **Output:** `${SID}.round1.bam`, `.bai`
- **Command:**
  ```
  minimap2 -ax map-ont -t ${task.cpus} --secondary=no -N 5 -Y --MD --eqx \
    -R "@RG\tID:${SID}\tSM:${SID}\tPL:ONT" \
    ${panel.mmi} ${reads} \
    | samtools sort -@ ${task.cpus} -O bam -o ${SID}.round1.bam -
  samtools index -@ ${task.cpus} ${SID}.round1.bam
  ```
- **Resources:** 16 CPU, 32 GB, 60 min.
- **Fail (gate):** <100 primary mapped reads → emit `NO_HCV_DETECTED`, terminate this sample.

## Step 5.8 — Genotype classification + mixed-infection detection
- **Module:** `GENOTYPE_CLASSIFY` (`bin/classify_genotype.py`)
- **Output:**
  - `${SID}.genotype_assignments.tsv` (read_id, ref_id, subtype, genotype, AS, mapq)
  - `${SID}.genotype_summary.json` (see schema in `docs/configuration.md`)
- **Script logic:**
  1. Walk primary alignments only (skip secondary/supplementary).
  2. Extract best-hit reference (highest AS); parse genotype/subtype from FASTA header with regex `^([0-9]+[a-z]?[a-z]?)_`.
  3. Compute fractions; flag mixed if non-primary genotype ≥ `min_secondary_fraction`.
  4. Compute ambiguous pool: reads whose top two hits differ by ΔAS < `params.ambiguous_delta_as`.
- **Resources:** 2 CPU, 4 GB, 15 min.

## Step 5.9 — Per-genotype read partitioning (mixed only)
- **Module:** `PARTITION_READS` (`bin/partition_reads.py`)
- **Output:** `${SID}.${GT}.reads.fastq.gz` per genotype; `${SID}.ambiguous.fastq.gz` (QC only)
- **Note:** Single-genotype samples skip this step; use `${SID}.hostdep.fastq.gz` directly.
- **Resources:** 2 CPU, 4 GB, 20 min.

## Step 5.10 — Consensus build
- **Module:** `BUILD_CONSENSUS`
- **Sub-steps:**
  1. Identify dominant panel reference per genotype (most reads).
  2. Re-map per-genotype reads to dominant reference:
     ```
     minimap2 -ax map-ont -t ${task.cpus} -Y --MD --eqx ${dom_ref} ${reads} \
       | samtools sort -@ ${task.cpus} -O bam -o ${GT}.dom.bam -
     samtools index ${GT}.dom.bam
     ```
  3. High-confidence variant call (`bcftools mpileup | bcftools call --ploidy 1`). Keep DP ≥ 10 and (AD[1]/DP) ≥ 0.5.
  4. Generate mask BED (positions <10× coverage via `mosdepth --quantize 0:1:10:`).
  5. Apply consensus, rename header to `${SID}_${GT}_consensus`, re-index:
     ```
     bcftools consensus -f ${dom_ref} -m mask.bed -o ${GT}.consensus.fasta ${GT}.consensus.vcf.gz
     samtools faidx ${GT}.consensus.fasta
     minimap2 -x map-ont -d ${GT}.consensus.mmi ${GT}.consensus.fasta
     ```
- **Output:** `consensus/${GT}/consensus.fasta`, `.fai`, `.mmi`, `consensus.vcf.gz`, `mask.bed`
- **Resources:** 8 CPU, 16 GB, 30 min.
- **Pass:** Single contig, length [4000, 11000] bp, %N < 30%. Fail: `LOW_COVERAGE_CONSENSUS` flag.

## Step 5.11 — Round 2 mapping
- **Module:** `MINIMAP2_ROUND2`
- **Command:** Same as Round 1 but without `-N 5 --secondary=no` (single reference).
- **Output:** `mapping/${GT}/round2.bam`, `.bai`
- **Resources:** 16 CPU, 16 GB, 30 min.
- **Pass:** >90% primary-mapped (else `LOW_MAPPING_RATE` flag).

## Step 5.12 — BAM post-processing for LoFreq
- **Module:** `LOFREQ_PREPROCESS`
- **Output:** `round2.iq.bam`, `.bai`
- **Commands:**
  ```
  lofreq indelqual --dindel -f ${consensus} -o round2.iq.bam round2.bam
  samtools index round2.iq.bam
  lofreq alnqual -b round2.iq.bam ${consensus} > round2.iq.alnq.bam || cp round2.iq.bam round2.iq.alnq.bam
  samtools index round2.iq.alnq.bam
  ```
- **Resources:** 4 CPU, 8 GB, 30 min.

## Step 5.13 — Coverage QC
- **Module:** `MOSDEPTH_R2`
- **Output:** `qc/coverage/${GT}/mosdepth.summary.txt`, per-window BED
- **Command:** `mosdepth -t ${task.cpus} -n --fast-mode --by 100 ${GT}.r2 round2.bam`
- **Resources:** 4 CPU, 4 GB, 10 min.
- **Pass:** Mean coverage ≥ `params.min_mean_coverage` (default 100). If <100, set `LOW_COVERAGE`, skip DEVIDER.

## Step 5.14 — rasusa subsample for LoFreq
- **Module:** `RASUSA_LOFREQ`
- **Output:** `lofreq_input.fastq.gz`, then remapped to produce `lofreq.iq.bam`
- **Command:**
  ```
  rasusa reads --coverage ${params.lofreq_max_depth} --genome-size 9646 \
    --seed 42 -o lofreq_input.fastq.gz ${reads}
  ```
- **Resources:** 4 CPU, 8 GB, 20 min total (subsample + remap + indelqual).

## Step 5.15 — LoFreq variant calling
- **Module:** `LOFREQ_CALL`
- **Output:** `variants/${GT}/lofreq.vcf.gz`, `.tbi`
- **Command:**
  ```
  lofreq call-parallel --pp-threads ${task.cpus} \
    --call-indels --min-mq 20 --min-bq 7 --min-cov 20 --sig 0.01 \
    -f ${consensus} -o lofreq.vcf lofreq.iq.bam
  bgzip lofreq.vcf && tabix -p vcf lofreq.vcf.gz
  ```
- **Resources:** 16 CPU, 16 GB, 90 min.

## Step 5.16 — Variant filtering
- **Module:** `VARIANT_FILTER`
- **Output:** `lofreq.filtered.vcf.gz`, `variants.tsv`
- **Command:**
  ```
  bcftools view -i 'INFO/AF >= ${params.min_report_af} & INFO/DP >= ${params.min_variant_depth}' \
    -Oz -o lofreq.filtered.vcf.gz lofreq.vcf.gz
  tabix -p vcf lofreq.filtered.vcf.gz
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/DP\t%INFO/SB\n' \
    lofreq.filtered.vcf.gz > variants.tsv
  ```
- **Resources:** 1 CPU, 2 GB, 5 min.

## Step 5.17 — Optional Clair3 corroboration
- **Module:** `CLAIR3_CORROBORATE` (opt-in via `--run_clair3`)
- **Output:** `variants/${GT}/clair3/merge_output.vcf.gz`, `clair3_concordance.tsv`
- **Model:** `r1041_e82_400bps_hac_v520`
- **Resources:** 16 CPU, 32 GB, 60 min.

## Step 5.18 — rasusa subsample for DEVIDER
- **Module:** `RASUSA_DEVIDER`
- **Output:** `devider_input.fastq.gz` → remapped → `devider.bam`
- **Command:**
  ```
  rasusa reads --coverage ${params.devider_max_depth} --genome-size 9646 \
    --seed 43 -o devider_input.fastq.gz ${reads}
  minimap2 -ax map-ont -t ${task.cpus} -Y --MD --eqx ${consensus_mmi} devider_input.fastq.gz \
    | samtools sort -@ ${task.cpus} -O bam -o devider.bam -
  samtools index devider.bam
  ```
- **Resources:** 8 CPU, 16 GB, 30 min.

## Step 5.19 — DEVIDER haplotype reconstruction
- **Module:** `DEVIDER_RUN`
- **Input:** `devider.bam`, `${GT}.consensus.fasta`, `lofreq.filtered.vcf.gz` (DEVIDER phases SNPs from this VCF; it does not call variants itself)
- **Output:** `haplotypes/${GT}/devider/` — parse by directory listing, not hard-coded filenames
- **Command (DEVIDER v0.0.1 — verify flags against `devider --help` at implementation time):**
  ```
  devider \
    -b devider.bam \
    -r ${consensus} \
    -v lofreq.filtered.vcf.gz \
    -o devider_out \
    -O \
    -t ${task.cpus} \
    --preset nanopore-r10 \
    --min-cov ${params.devider_min_cov} \
    --min-abund ${params.devider_min_abund} \
    --output-reads \
    --allele-output
  ```
  Critical notes:
  - Use `--preset nanopore-r10` (**not** `nanopore-r9` — wrong chemistry — **not** `ont` — does not exist).
  - `-O` overwrites output dir (required for Nextflow work-dir reruns).
  - `--output-reads` emits the haplotype-tagged BAM needed by Step 5.20.
  - `--allele-output` writes nucleotide alleles (not 0/1 codes); required by stitching/reporting scripts.
  - v0.0.1 has **no `--merge-windows` flag** — cross-region stitching is Step 5.20.
- **Resources:** 16 CPU, 64 GB, 4 hours (worst case).
- **Fail (graceful):** Non-zero exit → emit empty dir + `devider.failed`; do not fail the sample.

## Step 5.20 — Haplotype stitching
- **Module:** `STITCH_HAPLOTYPES` (`bin/stitch_haplotypes.py`)
- **Output:** `haplotypes/${GT}/merged_haplotypes.fasta`, `stitching_report.json`
- **Logic:** Walk DEVIDER haplotype regions in order. For each adjacent pair, find reads in the haplotype-tagged BAM that span the junction. Link if ≥ `params.stitch_min_reads` reads support a specific haplotype-X → haplotype-Y concatenation (Hamming distance over SNV positions). Emit separately where no spanning reads exist.
- **Resources:** 4 CPU, 8 GB, 30 min.

## Step 5.21 — Per-sample report
- **Module:** `SAMPLE_REPORT` (`bin/render_sample_report.py` + Jinja2)
- **Output:** `reports/${SID}_summary.html`, `reports/${SID}_summary.json`
- **Resources:** 1 CPU, 2 GB, 5 min.

## Step 5.22 — MultiQC aggregation
- **Module:** `MULTIQC`
- **Output:** `reports/multiqc_report.html`
- **Command:** `multiqc -f -o . ${all_qc_dirs}`
- **Resources:** 2 CPU, 4 GB, 10 min.
