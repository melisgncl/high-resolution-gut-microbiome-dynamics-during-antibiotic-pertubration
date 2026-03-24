# Paenibacillus Spectinomycin Resistance — Execution Roadmap

**Author:** Senior bioinformatics review, Melis Gencel project
**Date:** 2026-03-22
**Status:** APPROVED MASTER REFERENCE — no code written; plan only
**Precursor document:** `analysis/paenibacillus_resistance_plan/resistance_analysis_plan.md`

---

## Preamble: What You Are Walking Into

Three isolates are usable: M5 (103x, 1 contig), M6 (101x, 1 contig), M8 (84x, 1 contig). M7 is 9x, 7 fragmented contigs with each contig hitting a different species by mash — exclude from all resistance mechanism conclusions until re-sequenced.

The central biological question has three competing explanations:
1. **Intrinsic resistance** (H4): Paenibacillus macerans was always resistant; the PCR positive is contamination.
2. **De novo mutation** (H1/H2): rpsE or 16S helix 34 mutated in vivo under spectinomycin selection.
3. **HGT from E. coli** (H3): the aadA/speR gene transferred from the inoculum strain.

These are ordered by how quickly they can be falsified. The roadmap is designed to eliminate hypotheses as fast as possible, not to test all hypotheses with equal effort.

**The critical path bottleneck is I6.** Almost every decisive test (rpsE vs I6, 16S vs I6, size anomaly, species confirmation) requires the P. macerans I6 reference genome (GCF_000172175.2, ATCC 7068). Phase 0 and Phase 1 do not require it. Everything else does. Download it before anything else.

---

## Data Map (exact paths)

```
data/sequence_data/
├── TPHP3P_fastq/
│   ├── TPHP3P_1_M5__+_.fastq.gz          # 1.17 Gb, 191,702 reads, rN50=35 kb
│   ├── TPHP3P_2_M6__+_.fastq.gz          # 1.03 Gb, 145,608 reads
│   ├── TPHP3P_3_M7__+_.fastq.gz          # 69.9 Mb, 8,693 reads  ← UNRELIABLE
│   └── TPHP3P_4_M8__+_.fastq.gz          # 637 Mb, 87,841 reads
└── TPHP3P_results/
    ├── TPHP3P_1_M5__+_/
    │   ├── annotation/TPHP3P_1_M5__+_.fna     # assembly (single contig, 7.47 Mb)
    │   ├── annotation/TPHP3P_1_M5__+_.faa     # all proteins
    │   ├── annotation/TPHP3P_1_M5__+_.gff3    # annotation features
    │   ├── annotation/TPHP3P_1_M5__+_.tsv     # annotation table
    │   ├── annotation/TPHP3P_1_M5__+_.hypotheticals.faa
    │   ├── annotation/TPHP3P_1_M5__+_.hypotheticals.tsv
    │   └── annotation/TPHP3P_1_M5__+_.gbff
    ├── TPHP3P_2_M6__+_/   (same structure)
    ├── TPHP3P_3_M7__+_/   (same structure — 7 contigs, contaminated)
    └── TPHP3P_4_M8__+_/   (same structure)

Reference to download:
    data/references/
    ├── GCF_000172175.2_ASM17217v2_genomic.fna   # P. macerans I6 (ATCC 7068)
    ├── GCF_000172175.2_ASM17217v2_protein.faa   # I6 proteins
    └── GCF_000172175.2_ASM17217v2_genomic.gff   # I6 annotation

Analysis outputs to:
    results/genomics/
    ├── alignments/          # BAM files from minimap2
    ├── variants/            # Clair3 VCF outputs
    ├── whole_genome_align/  # nucmer/dnadiff outputs
    ├── amr_search/          # AMR database results
    ├── comparative/         # pan-genome, synteny
    ├── hgt_signals/         # ICE, GI, CRISPR outputs
    └── figures/             # plots
```

---

## Phase 0: Pre-flight Checks (No New Tools — Do This Session 1)

**Scientific goal:** Eliminate hypotheses using data already in hand. These require only shell commands on existing files, no alignment, no downloading. Results in this phase reorder the priority of all downstream work.

**Rationale for doing this before anything else:** Three of the most important hypotheses (H10, H4 partial, H3 partial) can be partially addressed within 30 minutes using the existing 16S R analysis results and the Bakta GFF3 files. Do not start downloading genomes or running alignments until this is done.

---

### 0.1 Check Nanopore chemistry from FASTQ headers

**Why it must be first:** All Nanopore-specific variant calling (Clair3) requires the correct basecaller model. Using the wrong model inflates false SNP rates. This is a dataset property you read once and use everywhere.

**Input:** `data/sequence_data/TPHP3P_fastq/TPHP3P_1_M5__+_.fastq.gz`

**Operation:** Read the first read header. It will contain fields like:
- `model_version_id=dna_r9.4.1_e8.1_...` (R9, Guppy)
- `model_version_id=dna_r10.4.1_e8.2_...` (R10, Dorado)
- Or just `basecall_model_version_id=...`

**Extract:** `zcat <FASTQ> | head -1` — the header line of the first read.

**Decision table:**

| Header contains | Chemistry | Clair3 model |
|-----------------|-----------|--------------|
| `r9.4.1` + `sup` | R9.4.1 SUP | `r941_prom_sup_g5014` |
| `r9.4.1` + `hac` | R9.4.1 HAC | `r941_prom_hac_g5014` |
| `r10.4.1` + `sup` + `400bps` | R10.4.1 SUP | `r1041_e82_400bps_sup_v430` |
| No model in header (older) | R9.4.1 default | `r941_prom_high_g360` |

**Save:** Record the chemistry and Clair3 model in `results/genomics/variants/nanopore_chemistry.txt`. Every downstream Clair3 call references this.

**Failure mode:** FASTQ has no model field (very old basecaller run). In that case, use `medaka consensus` instead of Clair3 — medaka has its own internal model selection logic.

---

### 0.2 Test H10/H4 with existing 16S data

**Why it must be early:** H10 (ecological release) is the hardest hypothesis to falsify with genomics alone — it requires the population-level view, which we already have. H4 (intrinsic resistance) predicts Paenibacillus expands in control mice too. This is the fastest test in the entire project.

**Input:** The existing 16S family-level tables already generated:
- `results/tables/16S/family/c_m1_family.csv`
- `results/tables/16S/family/c_m2_family.csv` (and c_m3, c_m4 if generated)
- The colonized mice family tables (m1–m8)

**Operation:** In R (using existing tidyverse/ggplot2 infrastructure), extract Paenibacillaceae relative abundance over time for all 12 mice. Plot:
1. Paenibacillaceae RA vs time, colonized mice (m1–m8)
2. Paenibacillaceae RA vs time, control mice (c_m1–c_m4)

**Decision criteria:**

| Pattern | Interpretation |
|---------|---------------|
| Paenibacillaceae expands in colonized mice but NOT in controls | E. coli colonization or its interaction with spectinomycin is driving expansion → H10 weakened, resistance hypothesis strengthened |
| Paenibacillaceae expands in BOTH colonized and control mice | Spectinomycin alone drives the expansion, regardless of E. coli → H10 is viable (ecological release from spectinomycin-sensitive taxa), H4 (intrinsic resistance) gains support |
| Paenibacillaceae does NOT expand in any group | The isolates came from rare, not abundant Paenibacillus → the expansion narrative needs reexamination |
| Expansion timing precedes E. coli peak | E. coli is not suppressing Paenibacillus → ecological interaction is unlikely cause |

**Save:** `results/genomics/figures/paenibacillus_RA_all_mice.pdf`

**This result determines how strongly to pursue H3/H9 (HGT) vs H4 (intrinsic).**

---

### 0.3 Cross-isolate rpsE comparison using annotation tables

**Why now:** The Bakta .tsv files for M5, M6, M8 all exist. The rpsE protein sequences can be extracted from the .faa files without any alignment. A 5-minute protein sequence extraction and visual alignment of three ~160 aa sequences tells you immediately whether rpsE varies between mice — which collapses or expands the H1 hypothesis before running a single alignment.

**Input:**
- `data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.faa`
- `data/sequence_data/TPHP3P_results/TPHP3P_2_M6__+_/annotation/TPHP3P_2_M6__+_.faa`
- `data/sequence_data/TPHP3P_results/TPHP3P_4_M8__+_/annotation/TPHP3P_4_M8__+_.faa`

**Operation:** Extract rpsE protein sequences. From the M5 annotation, rpsE is `NKFIDM_05734`. Check whether M6 and M8 annotation tables have the same locus tag pattern and extract their rpsE equivalents. Align the three sequences (MAFFT or muscle, pairwise is fine for three 160-aa sequences).

**Decision criteria:**

| Pattern | Implication |
|---------|-------------|
| All three rpsE identical | Either H1 is ancestral (in I6 or common ancestor) or H1 is false; compare to I6 in Phase 2 to distinguish |
| All three rpsE differ from each other | Three independent mutations → strongest possible evidence for convergent evolution under selection — prioritize H1 immediately |
| Two identical, one different | Partial convergence or one mutant lineage |

**Save:** `results/genomics/comparative/rpsE_cross_isolate_alignment.fa`

---

### 0.4 In silico PCR: resolve the PCR contradiction

**Why now:** The PCR positive result is the foundation of the HGT claim and the central unresolved contradiction (databases found nothing). Before running any wet genomics analysis of H3, determine what the PCR primers actually detect.

**Prerequisite:** You need the PCR primer sequences from the original study or methods. If they are not documented, this step cannot proceed computationally — add it to the wet lab list (Phase 7) and note it as a blocker.

**If primers are available:**

**Input:** Primer sequences + assembly FNA files for M5/M6/M8

**Operation:** Run `isPCR` (Jim Kent's in silico PCR tool, part of UCSC toolkit) or BLAST primer sequences against the assemblies. Command pattern: `isPCR assembly.fna primers.txt stdout -tileSize=7 -minPerfect=7 -minGood=12 -maxSize=3000`

**Decision criteria:**

| Result | Interpretation |
|--------|---------------|
| Primers amplify a product in M5/M6/M8 | A sequence matching the primer pair exists in Paenibacillus — even if databases missed it. Extract and BLAST this product to identify it. This is the sequence to investigate for H3. |
| Primers fail to amplify in M5/M6/M8 | The PCR-positive result was not detecting Paenibacillus DNA. Either E. coli contamination in the culture (most likely) or primer cross-reactivity with something not present in the assembly. |

**If primers not available:** Skip this step, flag that the PCR cannot be validated in silico, and add "obtain primer sequences from original methods" to Phase 7 wet lab list. This is a critical gap.

---

## Phase 1: Assembly Integrity Validation

**Scientific goal:** Confirm that the assemblies are accurate enough to support single-nucleotide-level conclusions. No SNP call, no comparative genomics, no resistance mutation claim is valid if the assembly has errors in the region of interest.

**Must precede:** Phases 2–5 entirely. If an assembly fails here, all conclusions from that isolate are restricted to structural variants only (which are more robust to assembly error).

**Inputs:** FASTQ reads (all 4) + assembled FNA (all 4)

---

### 1.1 NanoPlot read QC

**Tool:** NanoPlot v1.41+ (not NanoStat — NanoPlot produces better visualizations and is more actively maintained)

**Parameters:**
```
NanoPlot --fastq <FASTQ.gz> \
         --outdir results/genomics/QC/<sample> \
         --threads 8 \
         --plots hex dot \
         --N50 \
         --title "<sample>"
```

**One run per sample** (M5, M6, M7, M8).

**Metrics and thresholds:**

| Metric | Pass | Caution | Fail |
|--------|------|---------|------|
| Median Q-score | ≥ Q12 | Q10–Q12 | < Q10 |
| Read N50 | ≥ 10 kb | 5–10 kb | < 5 kb |
| % reads > 10 kb | ≥ 40% | 20–40% | < 20% |
| Unexpected bimodal length distribution | Not present | — | Present → barcoding artifact |

**For M7 specifically:** Expected to fail median Q threshold and show very short reads consistent with degraded DNA or inadequate yield. This confirms the exclude-from-analysis decision.

**Failure mode specific to this dataset:** M8 has 85x raw but only 84x assembled — nearly identical, suggesting minimal quality filtering was applied during assembly. Check whether the basecalling was done in Guppy SUP (super accurate) or HAC (high accuracy) — this is detectable from the FASTQ headers (Phase 0.1) and determines the per-base error rate. If HAC: base accuracy ~97-98% (Q17–19 median). If SUP: ~99% (Q20+).

**Save:** `results/genomics/QC/<sample>/NanoPlot-report.html` (one per isolate)

---

### 1.2 Read-to-assembly mapping and coverage validation

**This is the most important step in Phase 1.** A high-quality assembly with poor coverage at specific loci is useless for calling mutations at those loci.

**Tool:** minimap2 v2.26+

**Parameters:**
```
minimap2 -ax map-ont \
         --secondary=no \
         -t 8 \
         <assembly.fna> \
         <reads.fastq.gz> \
| samtools sort -o <sample>.sorted.bam -@ 4 \
&& samtools index <sample>.sorted.bam
```

**Justification of `-ax map-ont`:** This preset sets `-k15 -w10 -A2 -B4 -O4,24 -E2,1 -g5000 -r2000`. It is tuned for high-error Nanopore reads. Do NOT use `-ax sr` (short reads) or `-ax asm5` (assembly-to-assembly) for read mapping.

**Why `--secondary=no` here:** For coverage analysis of unique regions, secondary alignments at repetitive regions inflate apparent coverage and are misleading. Use separate run with `--secondary=yes -N5` for rRNA phasing (Phase 3.3).

**Run for M5, M6, M8 only. Do not invest compute on M7.**

**Coverage analysis:** After mapping:
```
samtools depth -d 0 -a <sample>.sorted.bam > <sample>_depth.txt
```
`-d 0`: disables the 8000x default cap — critical because rRNA operons will appear at ~500x coverage (5 copies × 100x) if collapsed in assembly.
`-a`: report all positions including zero coverage.

**Loci to check manually** (report coverage mean ± SD for each window):

| Locus | M5 coordinates | Expected depth | Flag if |
|-------|---------------|----------------|---------|
| rpsE | 6,418,100–6,418,588 | ~100x | < 50x (assembly gap) or > 200x (duplication) |
| rsmI | 51,713–52,606 | ~100x | < 50x or > 200x |
| Each 16S rRNA gene | Extract 5 coordinates from GFF3 | ~500x (5 copies collapsed) | < 400x (some copies missing) or uniform ~100x (copies resolved separately) |
| uS5 acetyltransferase | 5,830,084–5,830,644 | ~100x | < 50x |
| All 151 transposase-annotated regions | Batch | ~100x each | Any 2x+ elevation signals tandem duplication |

**Decision criteria:**

| Coverage at locus | Action |
|-------------------|--------|
| 80–150x for rpsE/rsmI | Proceed — variant calls reliable |
| < 50x at rpsE | Flag: assembly gap or mismapping; validate rpsE by Sanger (Phase 7) |
| > 200x at rpsE | Flag: possible tandem duplication → pursue H6 immediately |
| ~500x at all five 16S loci | Confirms all 5 copies are collapsed into one assembly region (expected) |
| ~100x at each of 5 distinct 16S loci | Confirms assembly has resolved copies separately — Nanopore long read advantage. Phase 3.3 rRNA phasing is less urgent. |
| 0x anywhere in the assembly | Misassembly. Contact the assembly pipeline team. |

**Save:**
- `results/genomics/alignments/M5.sorted.bam` (+ .bai)
- `results/genomics/alignments/M5_depth.txt`
- `results/genomics/alignments/M5_locus_coverage_summary.tsv` (manual or scripted table of the loci above)

---

### 1.3 Verify rRNA operon count and structure

**Tool:** Parse M5 GFF3 directly — no new tool required.

**Input:** `data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.gff3`

**Operation:** Count all features annotated as `rRNA` in the GFF3. Confirm: exactly 25 rRNA genes (5 × 16S + 5 × 23S + 5 × 5S). Extract the coordinates of all 5 16S rRNA genes. These coordinates are used in every downstream rRNA analysis.

**Expected structure of a ribosomal operon:**
```
16S rRNA (1,542 bp) → internal transcribed spacer → 23S rRNA (2,900 bp) → 5S rRNA (115 bp)
```
The five operons should be spread around the chromosome. Note their coordinates.

**Failure mode:** If fewer than 25 rRNA features are annotated, the assembly may have collapsed some operons. This does not invalidate the analysis but means rRNA phasing will be essential (you cannot assume all 5 copies are represented).

**Save:** `results/genomics/comparative/M5_rRNA_coordinates.bed` — BED format with all 25 rRNA features; the 5 16S entries are the primary target.

---

## Phase 2: Reference Genome Establishment

**Scientific goal:** Establish the correct reference genome and determine its quality. Every comparative variant call is relative to the reference; using the wrong reference or a draft with gaps makes all SNP calls uninterpretable.

**Must precede:** Phases 3 and 4. Phase 0 and 1 can proceed in parallel with this phase.

---

### 2.1 Download and assess I6 (P. macerans ATCC 7068, GCF_000172175.2)

**Operation:** Download from NCBI FTP:
```
ncbi-datasets genome download accession GCF_000172175.2 \
    --include genome,protein,gff3 \
    --filename data/references/GCF_000172175.2.zip
```
Or via Entrez utilities.

**Immediately check the assembly quality:**
- Count sequences in the FASTA: `grep -c "^>" GCF_000172175.2*.fna`
- If 1 sequence: closed chromosome — ideal reference
- If > 1 sequence: note the number of contigs. If > 100 contigs, I6 is a low-quality draft. **In that case:** use M5 as the reference for M6 and M8 comparisons (M5 is single circular, 103x, high confidence — the best available reference for within-species comparisons).

**Decision table:**

| I6 assembly quality | Consequence | Action |
|--------------------|-------------|--------|
| 1 contig (closed) | Gold standard reference | Use I6 as reference for all comparisons |
| 2–10 contigs | Partial draft | Use I6 but flag gaps; supplement with M5-vs-I6 for gap regions |
| > 10 contigs | Draft assembly | Use M5 as the reference for M6/M8; use I6 only for species-level comparisons, not SNP calling |
| I6 not available / deprecated | Major gap | Use the closest available closed P. macerans genome from NCBI; document the substitution |

**Note:** As of 2025, GCF_000172175.2 has assembly level "Scaffold" — it is likely a draft with multiple contigs. Anticipate using M5 as the primary reference for cross-isolate SNP calling.

**Save:** `data/references/GCF_000172175.2_genomic.fna`, `data/references/GCF_000172175.2_protein.faa`, `data/references/GCF_000172175.2_genomic.gff`

---

### 2.2 Confirm species identity (FastANI)

**Tool:** FastANI v1.34+

**Why FastANI over mash:** Mash uses random k-mer sketches against a database — the low 79% identity for M5 against P. woosongensis is almost certainly a database gap (P. macerans I6 is not in the mash chromosome reference database). FastANI does a direct whole-genome alignment and is the gold standard for ANI. The claimed 99% FastANI identity to I6 needs explicit confirmation.

**Parameters:**
```
fastANI -q data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna \
        -r data/references/GCF_000172175.2_genomic.fna \
        -o results/genomics/comparative/M5_vs_I6_fastANI.txt \
        --fragLen 3000 \
        --minFraction 0.2
```

`--fragLen 3000`: default is 3000. For Nanopore assemblies that may have occasional indels from homopolymer errors, keeping default is correct — short fragments would be overly sensitive to local errors.
`--minFraction 0.2`: lowered from default 0.2 (already default). If I6 is a draft with many contigs, coverage may be lower — 0.2 means 20% of the query must be covered for a result.

**Run for M5, M6, M8 vs I6. Also run M5 vs M6 vs M8 pairwise.**

**Decision criteria:**

| ANI | Interpretation |
|-----|---------------|
| ≥ 99% M5 vs I6 | Same species, I6 is the correct reference |
| 95–99% | Same species boundary but significant divergence; verify using MASH-SCREEN across all NCBI P. macerans assemblies |
| < 95% | Different species — I6 is the wrong reference. Must do BLAST sweep of all Paenibacillaceae genomes. Do not proceed with I6 as reference. |
| M5 vs M6 and M5 vs M8 > 99.9% | Nearly clonal isolates — confirms shared ancestry, makes cross-isolate SNP comparison meaningful |

**Save:** `results/genomics/comparative/fastANI_results.txt`

---

### 2.3 Whole-genome alignment — quantify the 1.77 Mb anomaly

**Tool:** MUMmer4 (nucmer + dnadiff + mummerplot)

**Why nucmer over Mauve:** Mauve is GUI-based, slower, and less scriptable. nucmer + dnadiff gives the same information in a pipeline-friendly format. Syri on top of nucmer output handles structural rearrangements.

**Parameters:**
```
nucmer --maxgap=500 \
       --mincluster=100 \
       --prefix=results/genomics/whole_genome_align/M5_vs_I6 \
       data/references/GCF_000172175.2_genomic.fna \
       data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna
```

`--maxgap=500`: allows up to 500 bp gaps in a single alignment block. Increase to `--maxgap=2000` if expecting large IS element insertions that break synteny locally.
`--mincluster=100`: minimum cluster size in bp. Default is fine.

Then:
```
dnadiff --prefix=results/genomics/whole_genome_align/M5_vs_I6 \
        data/references/GCF_000172175.2_genomic.fna \
        data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna
```

dnadiff output (`.report` file) gives: % identity, insertion lengths, unaligned regions. The unaligned regions in M5 relative to I6 are the extra 1.77 Mb.

**Additionally, run Syri** to classify structural rearrangements:
```
syri -c M5_vs_I6.delta \
     -r data/references/GCF_000172175.2_genomic.fna \
     -q <M5_assembly.fna> \
     --prefix M5_vs_I6_syri
```

**What to extract from the dnadiff .report:**
1. Total aligned bases (M5 side): should be ~5.7 Mb if I6 covers 5.7 Mb
2. Total unaligned bases in M5: should be ~1.77 Mb
3. Largest unaligned block in M5: if > 500 kb, this is a candidate large mobile element
4. Number of insertions > 1 kb in M5 vs I6: a few large insertions vs many small insertions are very different scenarios

**Save:** `results/genomics/whole_genome_align/M5_vs_I6.report` (dnadiff summary), `results/genomics/whole_genome_align/M5_vs_I6_syri.out` (structural variants), and a derived file `results/genomics/whole_genome_align/M5_novel_regions.bed` — the coordinates in M5 that have no I6 synteny (these are the priority regions for all downstream HGT analysis).

**Critical downstream dependency:** The novel regions BED file is the input to Phase 5 (HGT signal analysis). Every HGT and ICE tool runs on these coordinates.

**Failure mode:** If I6 is a draft with many contigs, the synteny plot will be fragmented and the "unaligned" regions will include I6 assembly gaps, not just genuine novel sequence. In this case, supplement nucmer with a read-based approach: map M5 reads to I6 and check which I6 regions have zero coverage (true gaps in I6 vs regions where M5 is simply divergent).

---

## Phase 3: Targeted Locus Investigation — The Critical Path

**Scientific goal:** Answer the three decisive questions (H1, H2, H3, H4) as directly as possible. This phase tests the hypotheses that can be falsified or confirmed with surgical precision.

**Must come before Phase 4** because Phase 4 (pan-genome, whole-genome comparative) is only necessary if this phase returns ambiguous results. If rpsE is mutated in all three mice and the mutation perfectly maps to the spectinomycin-binding interface, Phase 4 remains informative but not decisive for the primary conclusion.

**All of Phase 3 runs on M5, M6, M8. M7 is excluded.**

---

### 3.1 rpsE comparison: M5/M6/M8 vs I6 (H1 critical test)

**The single most decisive computation in the entire project.**

**Input:**
- rpsE protein sequence from M5 (.faa, locus NKFIDM_05734)
- rpsE protein sequences from M6 and M8 (extract from their .faa files by searching for rpsE / ribosomal protein S5)
- rpsE protein from I6 (extract from GCF_000172175.2_protein.faa — search for "ribosomal protein S5")
- Published spectinomycin-sensitive rpsE sequences for comparison: E. coli K12 (NP_416366.1), B. subtilis 168 (NP_388428.1), H. influenzae Rd (NP_438861.1), S. pneumoniae R6 (NP_357882.1)

**Tool:** MAFFT v7.5+ (protein mode)
```
mafft --auto \
      --reorder \
      results/genomics/comparative/rpsE_all.faa \
      > results/genomics/comparative/rpsE_aligned.faa
```

**Do not use Clustal for this** — MAFFT is more accurate for divergent bacterial S5 sequences and handles gaps better.

**What to look for — spectinomycin contact residues in S5:**
Spectinomycin contacts S5 primarily through helix α1 and the loop preceding β1. In E. coli numbering:
- K45 (contacts C1063 of 16S)
- R54 (contacts A1064)
- Q15, G16, G17 (N-terminal loop)

Map these positions onto the aligned P. macerans sequence. Are they conserved? Are they altered in our isolates relative to I6?

**Decision criteria:**

| rpsE result | Implication | Next step |
|-------------|-------------|-----------|
| M5/M6/M8 rpsE = I6 rpsE (identical) | H1 false for acquired mutation — but could still be intrinsic if both I6 and our isolates differ from E. coli at contact residues | Test H4: compare to E. coli/sensitive organisms |
| M5/M6/M8 all differ from I6 at same contact-residue position | H1 supported — ancestral mutation pre-dating colonization, but still conferring resistance | Compare to I6; if I6 also has the variant → natural P. macerans polymorphism at this position |
| M5, M6, M8 each have a DIFFERENT contact-residue mutation | Convergent evolution — strongest possible evidence for positive selection in vivo. Immediately flag as major finding. | Still compare to I6 to ensure none match I6 |
| Mutations present but NOT at contact residues | H1 weakened but not excluded — structural mapping needed | AlphaFold2 or PDB mapping in Phase 3.1.2 |
| M5/M6/M8 have rpsE identical to each other AND to I6 BUT all differ from E. coli at contact positions | H4 supported (intrinsic) — natural P. macerans divergence from drug target | Confirm MIC of I6 needed (wet lab) |

**Save:** `results/genomics/comparative/rpsE_aligned.faa`, `results/genomics/comparative/rpsE_contact_residues.txt` (manually annotated table of positions and amino acids in each isolate vs I6 vs E. coli)

---

### 3.2 16S helix 34 comparison: M5/M6/M8 vs I6 (H2 critical test)

**Context:** The spectinomycin binding site on 16S rRNA involves the helix 34 loop, specifically positions C1063, A1064, A1192, A1193 in E. coli numbering. P. macerans 16S numbering will be offset from E. coli — alignment is required to identify the equivalent positions.

**Input:**
- Extract 16S rRNA sequences from M5 GFF3 coordinates (from Phase 1.3 coordinates)
- Extract I6 16S rRNA from its GFF3 (if annotated) or use barrnap on the I6 assembly
- Reference E. coli K12 16S (NR_145189.1) — download from NCBI

**Tool:** MAFFT for alignment
```
mafft --auto results/genomics/comparative/16S_helix34_all.fna \
      > results/genomics/comparative/16S_helix34_aligned.fna
```

**After alignment:** Extract the helix 34 window (approximately E. coli positions 1063–1200) using the alignment coordinates. This is a ~140-nt region.

**Key question within 16S:** Since there are 5 copies of 16S in the M5 assembly, and the assembly likely collapsed them (see Phase 1.2), the sequence you extract from the FNA represents a consensus of all 5 copies. Variants present in only 1/5 copies will appear as minority variants at ~20% frequency and may be below the threshold for assembly into the consensus.

**To detect minority variants within the 5 16S copies:** Do NOT rely on the assembly consensus. Go back to the read-level mapping (Phase 1.2 BAM files). Examine the per-base frequency at helix 34 positions in the BAM pileup:
```
samtools mpileup -r <contig_1:16S_start-16S_end> \
                 -f <assembly.fna> \
                 -q 20 -Q 20 \
                 <sample>.sorted.bam
```
`-q 20 -Q 20`: minimum mapping quality 20 and base quality 20. This filters most Nanopore systematic errors.

**Interpret the pileup:**
- If a variant is at 100% frequency → all 5 copies carry it (or the assembly already shows it)
- If at ~80% frequency → 4/5 copies
- If at ~20% frequency → 1/5 copies (likely a recent mutation that has not spread to all operons)
- If at 5-10% with no Q20 reads → sequencing error, not a real variant

**Critical Nanopore failure mode for 16S:** The rRNA secondary structure causes polymerase stalling. Expect elevated error rates at stems and loops of helix 34. A variant at 10-15% frequency with moderate quality scores could be a genuine 1/5-copy mutation OR a Nanopore homopolymer error. The only way to distinguish: (a) check if it is present in all three mice (real), or (b) Sanger sequencing of individual 16S copies (Phase 7).

**Save:** `results/genomics/comparative/16S_helix34_aligned.fna`, `results/genomics/comparative/16S_helix34_pileup_M5.txt` (samtools mpileup output at helix 34 region)

---

### 3.3 rRNA copy phasing with long reads (Nanopore advantage)

**When to do this:** Only if Phase 3.2 reveals a variant at < 100% frequency at helix 34. If all five copies appear identical (variant at 100% or 0%), phasing is not informative.

**Goal:** Determine whether the helix 34 variant is in all 5 copies, some copies, or only one.

**Tool:** minimap2 with secondary alignments enabled

**Parameters:**
```
minimap2 -ax map-ont \
         --secondary=yes \
         -N 5 \
         -p 0.8 \
         --cs \
         -t 8 \
         <assembly.fna> \
         <reads.fastq.gz> \
| samtools view -F 256 -b - \
| samtools sort -o M5_rRNA_phasing.sorted.bam -@ 4
samtools index M5_rRNA_phasing.sorted.bam
```

`--secondary=yes -N5`: Allow up to 5 secondary alignments per read (needed because 16S copies are nearly identical and many reads will map equally well to all 5).
`-p 0.8`: Report secondary alignments with score ≥ 80% of primary score.
`--cs`: Output the cs tag (CIGAR-like detailed alignment). This is essential for identifying which copy of the 16S gene a read is most likely from — based on its flanking sequence context.

**Phasing strategy:** Extract reads that cover the helix 34 region. For each such read, examine the flanking sequence outside the 16S gene. Because the 5 operons are at different chromosomal positions with unique flanking sequences, a read long enough to extend past the 16S gene into flanking sequence can be unambiguously assigned to one operon copy. Cluster these reads by flanking sequence identity → determine base call at the helix 34 position per copy.

**Practical check:** With rN50 = 35 kb (M5) and 16S gene = 1.5 kb, most reads spanning a 16S gene will extend 10+ kb on each side — more than enough to identify unique flanking context for each operon.

**Save:** `results/genomics/alignments/M5_rRNA_phasing.sorted.bam`, `results/genomics/comparative/M5_16S_copy_helix34_variants.tsv` (one row per operon copy, base call at each helix 34 position)

---

### 3.4 aadA/speR search — resolve the PCR contradiction (H3 critical test)

**This is the most important negative result to establish.**

#### 3.4.1 Read-level mapping to aadA (most sensitive)

**Do NOT rely solely on the assembly** for this test. If aadA is present at low copy number or is divergent enough to misassemble, it may be absent from the assembly but present in the reads.

**Input:**
- E. coli K12 aadA sequence (obtain from the E. coli K12 inoculum genome, or use GenBank M97202.1 as the reference aadA type sequence)
- All three FASTQ files (M5, M6, M8)

**Tool:** minimap2
```
minimap2 -ax map-ont \
         --secondary=no \
         -t 8 \
         <aadA_reference.fna> \
         <reads.fastq.gz> \
| samtools view -F 4 -b - \
| samtools sort -o <sample>_vs_aadA.sorted.bam
```

`samtools view -F 4`: keep only mapped reads (discard unmapped).

**Check the count of mapped reads:**
```
samtools flagstat <sample>_vs_aadA.sorted.bam
```

**Decision criteria:**

| Reads mapping to aadA | Interpretation |
|-----------------------|---------------|
| 0 mapped reads | aadA sequence is not present in the isolate at any detectable level. H3 is false for aadA specifically. |
| 1–10 reads mapped | Ambiguous — check the alignment quality (MAPQ) and percent identity. If MAPQ < 20, these are likely false alignments to a divergent region. |
| > 10 reads mapped at > 90% identity | aadA-like sequence is present in the reads. Extract these reads and assemble or BLAST to determine the context. |
| > 10 reads mapped at 70–90% identity | Distant homolog present. Could be a native P. macerans ANT family gene with some similarity. Check whether these reads also map elsewhere in the P. macerans assembly. |

#### 3.4.2 Assembly-level BLAST search for aadA and broader ANT family

**Tools:** blastn (nucleotide) and tblastn (protein query vs nucleotide database)

**Commands pattern:**
- `blastn -query aadA.fna -subject M5_assembly.fna -task blastn-short -perc_identity 70 -qcov_hsp_perc 50`
- `tblastn -query ANT9_protein.faa -db M5_assembly.fna -evalue 1e-5` (using the ANT(9) protein from P. aeruginosa or S. aureus as query)

**Also run abricate at reduced thresholds:**
```
abricate --db card --minid 50 --mincov 50 <assembly.fna>
abricate --db resfinder --minid 50 --mincov 50 <assembly.fna>
abricate --db ncbi --minid 50 --mincov 50 <assembly.fna>
```

Default thresholds are 80% identity / 80% coverage — this will miss divergent homologs. 50/50 is aggressive; manually review all hits above this threshold.

#### 3.4.3 HMMER search for nucleotidyltransferase fold (H8 — novel gene)

**Tool:** HMMER3 + Pfam database

**Relevant Pfam models:**
- `PF05347`: Aminoglycoside 3'-phosphotransferase (APH family — not the primary target but related)
- `PF13612`: ANT/nucleotidyltransferase fold — covers all ANT-family enzymes including spectinomycin adenylyltransferases
- `PF02110`: AAC(3) acetyltransferase
- `PF13508`: AAC(6') acetyltransferase

**Command:**
```
hmmsearch --cut_nc \
          --tblout results/genomics/amr_search/M5_hmmer_amr.tblout \
          <Pfam-A.hmm> \
          data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.faa
```

`--cut_nc`: Use noise cutoff threshold (most sensitive — will include borderline hits). Then filter by E-value < 1e-5 and manually review.

**For any HMMER hit with E-value < 1e-5 matching ANT/AAC family:**
1. Note the locus tag (e.g., NKFIDM_XXXXX)
2. Check in the Bakta annotation table (.tsv) — what was it called?
3. Examine the neighboring genes in the GFF3 — is it flanked by IS elements? (IS element context is HGT hallmark)
4. BLAST the protein against NCBI nr at high sensitivity to find closest relatives and their taxonomic distribution

**Save:** `results/genomics/amr_search/M5_aadA_mapping.bam`, `results/genomics/amr_search/M5_abricate_card_50.txt`, `results/genomics/amr_search/M5_hmmer_amr.tblout`

---

### 3.5 rsmI and methyltransferase analysis (H5)

**Prerequisite for this phase:** I6 must be downloaded (Phase 2.1) and annotated (use barrnap for rRNA; for protein annotation, use the GFF from NCBI or re-annotate with Bakta).

**Operation:**
1. Extract rsmI protein from M5 (NKFIDM_00044) and compare to I6 rsmI (search I6 protein FAA for "rsmI" or "16S rRNA methyltransferase" or "2'-O-methyltransferase")
2. If rsmI found in I6: align M5 vs I6 rsmI protein — note all amino acid differences
3. Compare the promoter regions: extract 300 bp upstream of rsmI in M5 and in I6; pairwise align with `needle` (EMBOSS)
4. Check rsmI copy number in M5: `grep -i "rsmI\|00044\|C1402" annotation.tsv` — one hit expected

**Also search for:**
- RlmAII (G1405 methylation, Pfam PF05175) — aminoglycoside resistance in mycobacteria, may exist in Paenibacillus
- Cfr-like methyltransferases (A2503 methylation, Pfam PF08241)

**HMMER search for methyltransferase domains:**
```
hmmsearch --cut_nc \
          -E 1e-5 \
          --tblout results/genomics/amr_search/M5_methyltransferases.tblout \
          PF05175.hmm PF08241.hmm PF01523.hmm \
          <M5.faa>
```

**Decision:** If rsmI is absent from I6 but present in M5/M6/M8 → acquired enzyme, strong H5 signal. If identical in I6 and M5 → native gene, not the resistance explanation. If promoter differs → investigate transcriptional regulation as a potential mechanism (requires RNA-seq to confirm — Phase 7).

---

### 3.6 Variant calling — whole genome (SNP/indel inventory)

**This is a supporting step, not a decisive one.** The decisive comparisons (rpsE, 16S) are done above. Whole-genome variant calling is needed to: (a) provide a comprehensive list of all M5/M6/M8 vs I6 differences, (b) identify IS element insertions that disrupt genes, (c) find large indels in resistance loci.

**Tool selection rationale:**
- **Do NOT use Snippy for Nanopore data** — Snippy uses BWA-MEM, which is designed for short reads. BWA-MEM miscounts indels in Nanopore long reads due to error model mismatch.
- **Use Clair3** — Nanopore-native neural network variant caller. Trained on Nanopore error profiles. Best precision/recall for SNP calling from long reads as of 2025.
- **Alternatively, if Clair3 model is unavailable for the chemistry:** Use `medaka variant` as a fallback. Less accurate but chemistry-agnostic.

**Prerequisite:** Clair3 model from Phase 0.1.

**Command pattern:**
```
run_clair3.sh \
    --bam_fn=results/genomics/alignments/M5.sorted.bam \
    --ref_fn=data/references/GCF_000172175.2_genomic.fna \
    --threads=8 \
    --platform=ont \
    --model_path=<model_from_phase_0.1> \
    --output=results/genomics/variants/M5_vs_I6/ \
    --include_all_ctgs \
    --no_phasing_for_fa \
    --min_coverage=8
```

`--min_coverage=8`: minimum coverage to call a variant — with our ~100x samples, this is very conservative and will produce very high sensitivity.
`--include_all_ctgs`: call on all contigs (since M5 is a single contig, this is fine).
`--no_phasing_for_fa`: skip haplotype phasing (we have haploid bacterial genomes).

**Note:** Clair3 requires the BAM to be aligned to the **same reference** as the `--ref_fn` argument. This means you need to run minimap2 mapping M5 **reads** against I6 (not against the M5 assembly). This is different from the Phase 1.2 mapping (which mapped M5 reads against the M5 assembly for coverage validation).

**Additional mapping run required for Clair3:**
```
minimap2 -ax map-ont --secondary=no -t 8 \
         data/references/GCF_000172175.2_genomic.fna \
         data/sequence_data/TPHP3P_fastq/TPHP3P_1_M5__+_.fastq.gz \
| samtools sort -o results/genomics/alignments/M5_vs_I6.sorted.bam
samtools index results/genomics/alignments/M5_vs_I6.sorted.bam
```

**After Clair3:** Filter the VCF for high-confidence calls:
```
bcftools filter -i 'QUAL>20 && INFO/DP>20' \
                results/genomics/variants/M5_vs_I6/merge_output.vcf.gz \
                > results/genomics/variants/M5_vs_I6_filtered.vcf
```

**Annotate variants with SnpEff** (optional but highly useful): annotates each SNP with the gene it falls in and its predicted effect (missense, nonsense, synonymous, upstream, intronic).

**Save:** `results/genomics/variants/M5_vs_I6_filtered.vcf`, (same for M6, M8)

---

## Phase 4: Genome-Wide Comparative Analysis

**Scientific goal:** Understand the full landscape of what distinguishes our isolates from I6 — pan-genome, structural variation, mobile element census, size anomaly investigation. This provides context for any resistance finding and identifies whether the 1.77 Mb anomaly is relevant.

**Can proceed in parallel with Phase 3** once the reference is established (Phase 2 complete). Phase 3 is the critical path for the main conclusion; Phase 4 provides supporting evidence and explanation.

---

### 4.1 Identify and annotate novel regions (the 1.77 Mb anomaly)

**Input:** `results/genomics/whole_genome_align/M5_novel_regions.bed` (from Phase 2.3)

**This is the most important output of Phase 4.** The 1.77 Mb of novel sequence relative to I6 is either:
1. A large integrated mobile element (ICE, GI, mega-transposon)
2. Expansion of existing gene families (transposons, rRNA operons)
3. A sequenced plasmid that circularized into the assembly
4. Genuine strain-specific accessory genome (large gene clusters for nitrogen fixation, polysaccharide metabolism, etc.)
5. Combination of the above

**Operation:** For each novel region identified by nucmer/dnadiff:
1. Extract the coordinates and the sequence
2. Determine its gene content from the Bakta annotation (overlap with GFF3 features)
3. Compute the GC content in 1 kb windows: `seqkit sliding --window 1000 --step 200 <region.fna> | seqkit fx2tab --gc`
4. Check if the region is flanked by IS elements (look at Bakta GFF3 for flanking transposase annotations)
5. Check if the region is flanked by tRNA genes (classic GI integration sites) in the GFF3
6. Calculate the cumulative size of the 151 IS elements: `number_of_IS × average_IS_length` — how much of the 1.77 Mb do IS elements alone account for?

**Genomic island prediction:**
**Tool:** IslandViewer4 (web server: http://www.pathogenomics.sfu.ca/islandviewer/) OR Alien Hunter (command line)

Alien Hunter is preferred for command line pipelines:
```
alien_hunter <M5_assembly.fna> \
             results/genomics/hgt_signals/M5_alien_hunter_output
```

Cross-reference Alien Hunter predictions with the novel regions BED file — do the predicted GIs overlap the novel regions?

**Save:**
- `results/genomics/comparative/M5_novel_regions_annotated.tsv` (BED + gene content + GC content per region)
- `results/genomics/hgt_signals/M5_alien_hunter_output`

---

### 4.2 Pan-genome analysis: M5, M6, M8 (+ I6 if available in GFF3 format)

**Tool:** Panaroo v1.3+ (preferred over Roary)

**Why Panaroo over Roary:** Panaroo corrects for annotation errors that inflate the accessory genome — particularly relevant here because Bakta and the I6 NCBI annotation will use different gene-calling thresholds. Panaroo also produces a gene presence/absence matrix that is directly useful for identifying isolate-specific genes.

**Input:** GFF3 files (Bakta output for M5/M6/M8; NCBI GFF for I6 — note: you may need to convert NCBI GFF3 to Prokka-compatible format for Panaroo, or run Prokka on I6 to get a consistent annotation)

**Parameters:**
```
panaroo -i M5.gff3 M6.gff3 M8.gff3 I6.gff3 \
        -o results/genomics/comparative/panaroo_output \
        --clean-mode moderate \
        -t 8 \
        --aligner mafft \
        --remove-invalid-genes
```

`--clean-mode moderate`: balances between strict (misses genuine accessory genes) and sensitive (keeps noise).

**Output of interest:**
- `gene_presence_absence.csv`: 0/1 matrix for all genes across all 4 isolates
- Genes present in M5/M6/M8 but absent from I6: these are acquired genes (regardless of whether by HGT, de novo gene gain, or I6 assembly gap)
- Genes present in only one isolate: lineage-specific (check for resistance gene candidates)

**Save:** `results/genomics/comparative/panaroo_output/` (full directory)

---

### 4.3 IS element census and insertion site analysis

**Input:** M5 Bakta GFF3 (151 transposase annotations)

**Operation:**
1. Extract all transposase annotations from GFF3: filter features where product contains "transposase" or "IS4" or similar
2. Classify IS families: IS4, IS1, IS3, IS5, IS256, etc. — use the Bakta product annotations
3. For each IS element insertion: what gene does it fall in (if any)? What is the flanking gene structure?
4. Compare IS element positions between M5, M6, M8 using the whole-genome alignment from Phase 2.3. IS elements in identical positions in all three isolates = ancestral. IS elements in only one isolate = transposed during or after colonization.

**Why IS insertions matter for resistance:** IS element insertion into a repressor gene (marR, AcrR, MexR) de-represses an efflux pump constitutively. IS insertion upstream of a gene can create a new promoter. These are common mechanisms of resistance emergence without acquiring a new gene.

**Tool for automated IS classification:** ISfinder (web) or ISEScan (command line):
```
isescan.py --seqfile <M5_assembly.fna> \
           --output results/genomics/comparative/M5_isescan/ \
           --nthread 4
```

Cross-reference ISEScan output with the 151 Bakta transposase annotations to confirm family assignments.

**Save:** `results/genomics/comparative/M5_IS_element_census.tsv` (IS family, coordinates, flanking genes)

---

### 4.4 Structural variant calling: Sniffles2

**Why:** Clair3 (Phase 3.6) calls small variants (SNPs, small indels). Large structural variants (inversions, translocations, large insertions) require long-read-specific SV callers.

**Tool:** Sniffles2 v2.x

```
sniffles --input results/genomics/alignments/M5_vs_I6.sorted.bam \
         --vcf results/genomics/variants/M5_vs_I6_SVs.vcf \
         --reference data/references/GCF_000172175.2_genomic.fna \
         --threads 8 \
         --minsvlen 50 \
         --minsupport 5
```

`--minsvlen 50`: minimum SV length in bp. Default is 50.
`--minsupport 5`: minimum reads supporting an SV call. At 100x, 5 reads is very conservative — consider increasing to 10.

**What to look for:**
- Large insertions in M5 relative to I6 (the 1.77 Mb should appear as insertion SVs)
- Inversions near resistance loci (inversions flanked by IS elements are sometimes used as phase-variable resistance switches)
- Deletions in M5 vs I6 (if M5 has deleted a susceptibility gene — unlikely for spectinomycin but possible)

**Save:** `results/genomics/variants/M5_vs_I6_SVs.vcf`

---

## Phase 5: HGT Signal Analysis

**Scientific goal:** Determine the mechanism by which any resistance determinant was transferred (if H3 or novel-gene hypotheses are supported), and assess whether the 1.77 Mb anomaly itself is a transferred element. This phase is contingent on Phase 3 results — only pursue fully if there is evidence for an acquired gene.

**Can be deprioritized if:** Phase 3 clearly shows H1 (rpsE mutation in all three mice, not present in I6, structurally maps to drug-binding interface) AND Phase 3.4 shows no aadA reads. In that case, Phase 5 can be reduced to a brief confirmatory check.

---

### 5.1 ICE/integrative element search

**Tools (in order of preference):**
1. **CONJscan** — detects conjugative systems by HMM profiles of core conjugation proteins (T4SS). Command line, most rigorous:
   ```
   MacSyFinder --sequence-db <M5.faa> \
               --db-type gembase \
               -d CONJscan/ \
               -o results/genomics/hgt_signals/M5_conjscan
   ```
2. **ICEfinder** — web server (http://db-mml.sjtu.edu.cn/ICEfinder/). Detects ICEs, IMEs, AICEs. Submit M5 assembly FNA.
3. **PhageBoost or PHASTER** for prophage detection (web: https://phaster.ca/)

**Focus on the novel regions** (Phase 4.1 BED file). Specifically: does the 1.77 Mb novel sequence contain any T4SS or relaxase gene? A T4SS + relaxase + oriT = conjugative element.

**Save:** `results/genomics/hgt_signals/M5_conjscan/`, `results/genomics/hgt_signals/M5_ICEfinder_result.txt` (manual copy from web)

---

### 5.2 CRISPR spacer analysis

**Input:** M5 Bakta GFF3 (11 CRISPR arrays annotated)

**Operation:**
1. Extract all CRISPR spacer sequences from the GFF3 (features annotated as `repeat_region` with `rpt_type=direct` or CRISPR-specific tags from Bakta)
2. Alternatively, use CRISPRdetect de novo on the M5 assembly to get spacer sequences with confidence scores
3. BLAST all spacers against:
   - E. coli phage genomes (download from NCBI RefSeq viral database: `ncbi-datasets download virus --include genome --host-taxon "Escherichia coli"`)
   - E. coli plasmid sequences (PLSDB database)
   - The full NCBI nucleotide database (nt) for any spacer with no phage/plasmid hit — it may match an ICE

**Interpret:** CRISPR spacers matching E. coli phages → Paenibacillus has historically been exposed to E. coli phages (implies the two species coexisted in an environment where E. coli phages were present → consistent with a shared gut environment). Spacers matching E. coli plasmids → even stronger evidence of prior cross-species DNA transfer.

**Save:** `results/genomics/hgt_signals/M5_CRISPR_spacers.fna`, `results/genomics/hgt_signals/M5_CRISPR_blast_vs_ecoli.txt`

---

### 5.3 Codon usage analysis of novel regions

**Only if** a candidate resistance gene is identified in Phase 3.4 or Phase 5.1.

**Tool:** CodonW or EMBOSS `cusp`

**Operation:** For each gene in the novel regions that has an E. coli homolog (or is a candidate resistance gene):
1. Calculate codon usage table (RSCU values)
2. Compare to: (a) P. macerans genome average RSCU (from all M5 CDSs), (b) E. coli K12 genome RSCU (from reference)
3. A gene with E. coli-like codon usage embedded in a Paenibacillus genome is a strong HGT indicator

**Save:** `results/genomics/hgt_signals/codon_usage_analysis.tsv`

---

### 5.4 Phylogenetic placement of any candidate gene

**Only if** a candidate resistance gene is identified.

**Tool:** MAFFT (alignment) + IQ-TREE2 (phylogeny)

**Procedure:**
1. BLAST the candidate gene protein against NCBI nr → download the top 50 hits with their taxonomy
2. Add the M5 candidate sequence to the alignment
3. `mafft --auto candidate_and_homologs.faa > candidate_aligned.faa`
4. `iqtree2 -s candidate_aligned.faa -m TEST -B 1000 --prefix results/genomics/hgt_signals/candidate_gene_tree`
5. Visualize with iTOL (web) or ggtree (R)

**Decision:** If the M5 candidate branches within an E. coli clade → HGT from E. coli. If it branches within Paenibacillus → native gene, not HGT.

---

## Phase 6: Decision Matrix and Reporting

**Operate this phase after Phases 3 and 4 are complete** (Phase 5 may still be in progress for some sub-analyses).

### 6.1 Decision matrix

Apply this logic sequentially:

**Step 1: Is H4 (intrinsic resistance) true?**
- rpsE M5 = rpsE I6? AND 16S helix 34 M5 = I6? AND Paenibacillus expands in control mice?
- If YES to all three: the isolates may be intrinsically resistant. Halt resistance mechanism analysis. Recommend MIC of I6 (Phase 7.1) before making any further claim.
- If NO to at least one: proceed to Step 2.

**Step 2: Is H3 (aadA/HGT gene) present?**
- Any reads map to aadA at >90% identity? OR HMMER finds ANT/AAC family hit in M5?
- If YES: extract the candidate, check its phylogeny (Phase 5.4), check the flanking context for mobile element signatures
- If NO: H3 is false. Proceed to Step 3.

**Step 3: Is H1 (rpsE mutation) supported?**
- rpsE M5/M6/M8 ≠ I6 at contact residues? AND mutations map to drug-binding interface?
- If YES in all three mice (same mutation): ancestral rpsE resistance mutation (pre-colonization or in the common ancestor)
- If YES with different mutations per mouse: convergent evolution — strongest possible claim for positive selection
- If YES in some, NO in others: partial support

**Step 4: Is H2 (16S helix 34 mutation) supported?**
- Variant at helix 34 positions in the pileup? Confirmed across multiple reads at >20% frequency?
- If YES: determine copy number (Phase 3.3). Report with copy fraction.

**Step 5: Residual ambiguity**
- If Steps 1–4 all negative or ambiguous: fall back to the inventory of findings:
  - Novel 1.77 Mb regions: what do they contain?
  - IS element insertions: any in resistance-relevant loci?
  - rsmI promoter: any differences from I6?
- These secondary findings can support a narrative even without a clear primary mechanism.

### 6.2 Required comparisons for any published claim

No claim about resistance mechanism can be published without:
1. Comparison to I6 at every relevant locus (rules out intrinsic resistance)
2. Cross-isolate comparison (M5/M6/M8) — shared vs. independent mutations
3. Structural mapping of any rpsE mutation to the drug-binding interface (PDB 4WFA or equivalent)
4. Resolution of the PCR contradiction (in silico PCR or repeat experiment)

---

## Phase 7: Wet Lab Priority Queue

Ordered by information value / cost ratio:

| Priority | Experiment | What it resolves |
|----------|-----------|-----------------|
| **CRITICAL** | MIC determination for P. macerans I6 ATCC stock | H4 (intrinsic resistance) — if I6 MIC > clinical breakpoint, the genomics analysis is moot |
| **CRITICAL** | Repeat PCR with controls (M5/M6/M8 + I6 negative + E. coli positive + blank) | Resolves the PCR contradiction — validates or invalidates H3 |
| **CRITICAL** | Sanger sequencing of rpsE and 16S helix 34 from M5/M6/M8 colonies | Validates Nanopore SNP calls independent of assembly errors |
| **HIGH** | Obtain original PCR primer sequences from methods | Enables in silico PCR (Phase 0.4 — prerequisite) |
| **MEDIUM** | Re-sequence M7 with R10.4.1 + ≥30x target | M7 is currently unusable — re-sequencing would give 4/4 isolates instead of 3/4 |
| **MEDIUM** | Express M5 rpsE in spectinomycin-sensitive E. coli; measure MIC | Proves M5 rpsE mutation confers resistance (H1 functional validation) |
| **LOW** | RNA-seq M5 vs I6 in identical conditions | Tests H5 (rsmI overexpression) and H7 (efflux) — expensive, pursue only if promoter mutations found |
| **LOW** | Transfer experiment: Paenibacillus × E. coli cross-culture | Tests H9 (conjugative transfer) — only relevant if ICE found |

---

## Phase 8: Outputs and File Organization

**Convention:** All result files go under `results/genomics/`. All scripts go under `analysis/paenibacillus_resistance/` (scripts to be created in a future session after user approval).

```
results/genomics/
├── QC/
│   ├── M5_nanoplot/   M6_nanoplot/   M8_nanoplot/   M7_nanoplot/
├── alignments/
│   ├── M5_vs_M5asm.sorted.bam   (read-to-own-assembly, for coverage)
│   ├── M5_vs_I6.sorted.bam      (read-to-I6-reference, for Clair3)
│   ├── M5_rRNA_phasing.sorted.bam
│   ├── M5_vs_aadA.sorted.bam    (read-to-aadA, for H3 test)
│   ├── (same for M6, M8)
├── variants/
│   ├── nanopore_chemistry.txt
│   ├── M5_vs_I6/               (Clair3 output directory)
│   ├── M5_vs_I6_filtered.vcf
│   ├── M5_vs_I6_SVs.vcf        (Sniffles2)
│   ├── (same for M6, M8)
├── whole_genome_align/
│   ├── M5_vs_I6.report         (dnadiff summary)
│   ├── M5_vs_I6_syri.out
│   ├── M5_novel_regions.bed    ← CRITICAL INPUT TO PHASE 5
├── comparative/
│   ├── rpsE_aligned.faa
│   ├── rpsE_contact_residues.txt
│   ├── 16S_helix34_aligned.fna
│   ├── 16S_helix34_pileup_M5.txt
│   ├── M5_16S_copy_helix34_variants.tsv
│   ├── M5_rRNA_coordinates.bed
│   ├── M5_IS_element_census.tsv
│   ├── M5_novel_regions_annotated.tsv
│   ├── panaroo_output/
│   ├── fastANI_results.txt
├── amr_search/
│   ├── M5_abricate_card_50.txt
│   ├── M5_hmmer_amr.tblout
│   ├── M5_methyltransferases.tblout
│   ├── M5_aadA_mapping_flagstat.txt
├── hgt_signals/
│   ├── M5_alien_hunter_output
│   ├── M5_conjscan/
│   ├── M5_CRISPR_spacers.fna
│   ├── M5_CRISPR_blast_vs_ecoli.txt
│   ├── codon_usage_analysis.tsv
│   ├── candidate_gene_tree/    (if applicable)
└── figures/
    ├── paenibacillus_RA_all_mice.pdf   (Phase 0.2)
    └── (additional plots as generated)

data/references/
├── GCF_000172175.2_genomic.fna
├── GCF_000172175.2_protein.faa
├── GCF_000172175.2_genomic.gff
└── aadA_reference.fna            (from E. coli K12 inoculum genome or GenBank M97202.1)
```

---

## Execution Order Summary

This is the dependency graph compressed into a sequence:

```
Session 1 (no downloads, no alignments):
  Phase 0.1 → Chemistry detection (read FASTQ header)
  Phase 0.2 → Paenibacillus 16S RA in controls (R analysis, existing data)
  Phase 0.3 → rpsE cross-isolate comparison from .faa files (MAFFT, 3 sequences)
  Phase 0.4 → In silico PCR (if primer sequences available)
  Phase 1.3 → rRNA operon count from GFF3 (grep, no alignment)

  --- This session determines: how likely is H4? Are M5/M6/M8 rpsE identical? ---

Session 2 (download I6, run alignments):
  Phase 2.1 → Download and assess I6 (ncbi-datasets)
  Phase 2.2 → FastANI M5/M6/M8 vs I6 (fast, ~5 min per pair)
  Phase 2.3 → nucmer + dnadiff → novel regions BED
  Phase 1.1 → NanoPlot QC (run in background, takes ~30 min per sample)
  Phase 1.2 → minimap2 read-to-assembly mapping (run in background)

  --- After Session 2: reference confirmed, novel regions identified ---

Session 3 (targeted locus analysis):
  Phase 3.1 → rpsE M5/M6/M8 vs I6 alignment (MAFFT) + contact residue mapping
  Phase 3.2 → 16S helix 34 M5/M6/M8 vs I6 (MAFFT + samtools mpileup)
  Phase 3.4 → minimap2 reads vs aadA + abricate + HMMER
  Phase 3.5 → rsmI comparison M5 vs I6

  --- After Session 3: primary hypotheses resolved ---

Session 4 (contextual analysis — only if needed):
  Phase 3.3 → rRNA phasing (only if helix 34 variant at <100%)
  Phase 3.6 → Clair3 variant calling (whole genome)
  Phase 4.1 → Novel regions annotation and GC analysis
  Phase 4.2 → Panaroo pan-genome
  Phase 4.3 → IS element census (ISEScan)
  Phase 4.4 → Sniffles2 SV calling

Session 5 (HGT mechanism — only if H3 or novel gene evidence):
  Phase 5.1 → CONJscan + ICEfinder
  Phase 5.2 → CRISPR spacer analysis
  Phase 5.3 → Codon usage of novel regions
  Phase 5.4 → Phylogenetic tree of candidate gene

Final: Phase 6 decision matrix → report
```

---

## Critical Failure Modes and Dataset-Specific Risks

| Risk | Likely to occur? | Mitigation |
|------|-----------------|------------|
| I6 is a draft assembly (many contigs) | HIGH — GCF_000172175.2 is labeled "Scaffold" level | Confirmed in Phase 2.1; switch to M5 as reference for SNP calling if so |
| rpsE or 16S mutation is in a homopolymer run → Nanopore miscall | MEDIUM — rRNA genes are AT-rich | Use only reads with Q≥20; flag any variant in homopolymer context; validate with Sanger |
| The PCR primer sequences are not documented | HIGH — often not in methods | Must obtain from original authors; without this, H3 cannot be tested in silico |
| M5/M6/M8 rpsE all identical to I6 AND all lack aadA | MEDIUM | Fall back to: IS element disruptions, rsmI, efflux, structural variants; acknowledge mechanism unknown |
| CRISPR spacers match E. coli phages, strengthening HGT case | POSSIBLE | Good outcome — include as supporting evidence |
| The 1.77 Mb anomaly is entirely accounted for by IS element expansion | POSSIBLE (~151 IS × ~1.5 kb = ~225 kb — accounts for only ~13%) | This leaves ~1.5 Mb unexplained; a large GI or integrated plasmid is likely |
| Clair3 model not available for the detected chemistry | LOW if R9.4.1 or R10.4.1 standard | Use medaka variant as fallback |
| M7 re-sequencing reveals M7 is genuinely a different species | MEDIUM — mash hits to 7 species | Excludes M7 permanently from the resistance analysis; reduces n to 3 |
| Panaroo fails due to inconsistent annotation format (Bakta vs NCBI GFF3) | MEDIUM | Re-annotate I6 with Bakta using the same database version for consistency |

---

*This document supersedes all analysis planning documents. The only decision that should modify it is an explicit finding from one of the decision points above. Any new session should begin by reading the Execution Order Summary and determining which session it corresponds to.*
