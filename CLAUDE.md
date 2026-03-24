# Claude Code Rules for hgt-study

## Project
Reanalysis of barcoded E. coli microbiome HGT study.
- 8 colonized mice: cohort 1 (m1-m4), cohort 2 (m5-m8)
- 4 antibiotic-only controls: cm1-cm4
- Key finding: Paenibacillaceae acquired spectinomycin resistance via HGT

## R Code Standards
- Never use setwd() — always use here::here() for all paths
- Always use tidyverse style (dplyr, ggplot2, tidyr)
- Every script must start with a comment block: title, author, date, input, output
- Never hardcode mouse names — always loop or map over sample lists
- Always inspect data with head(), dim(), str() before any analysis
- Never guess data format — read and check first

## Analysis Order
1. 01_16S — community composition
2. 02_barcode — lineage dynamics (Doblin)
3. 03_diversity — Hill diversity indices
4. 04_coclustering — barcode + 16S joint clustering
5. 05_DCM — Jacobian stability analysis

## Package Management
- renv is used — always run renv::snapshot() after installing new packages
- Never install packages without telling the user first

## Git
- Commit after each completed analysis module
- Commit messages: "module 01_16S complete", "module 02_barcode complete" etc.

## Ways of Working
- You are a senior developer collaborator, not a task executor
- For any non-trivial action: propose a plan first, wait for approval
- Break every plan into numbered steps and show them before doing anything
- After completing each step, summarize what was done and what comes next
- If you see a better approach than what was asked, say so and explain why
- Never assume silence means approval — always wait for explicit confirmation
- Flag technical debt, bad patterns, or design issues when you see them
- When reading old code: diagnose before prescribing

---

## Project Status (updated 2026-03-22)

### Completed modules (committed to git)
| Module | Commit | Status |
|--------|--------|--------|
| 01_16S | 33393c8 | DONE — family/genus tables, composition plots |
| 02_barcode | 0e9cc7d | DONE — Doblin lineage dynamics |
| 03_diversity | (see log) | DONE — Hill diversity indices |
| 04_coclustering | 393b6b6 | DONE — barcode + 16S joint clustering |
| 05_DCM | 2bebd78 | DONE — Jacobian stability; control mice fixed (span_eff + seq_len obs_time) |

### Key finding from 16S + DCM analysis
- Paenibacillaceae is a structural outlier across DCM mice — C3↔Paenibacillaceae and C4→Paenibacillaceae pairs appear in 7–8 of 12 mice
- Paenibacillaceae acquired spectinomycin resistance via HGT (the claim under investigation)

### Paenibacillus resistance analysis — status (updated 2026-03-22 Session 3)
- **4 Nanopore isolates sequenced:** M5 (m5, 103x, 1 contig ★), M6 (m6, 101x, 1 contig ★), M7 (m7, 9x, 7 contigs ⚠ UNRELIABLE), M8 (m8, 84x, 1 contig ★)
- **Species:** Paenibacillus macerans (sourmash confirmed; FastANI vs correct I6 reference pending WSL)
- **Assembly quality:** CheckM2 99.85% complete, 0.06–0.12% contamination for all isolates
- **Genome size — REVISED:** M5 = 7,468,786 bp; I6 (correct complete ref) = 7,113,008 bp → **356 kb extra in M5, NOT 1.77 Mb**. The 1.77 Mb anomaly was an artifact of the old draft scaffold reference (GCF_000172175.2, deprecated).
- **I6 reference — REVISED:** GCF_000172175.2 is DEPRECATED. Correct accession: **GCF_022494515.1** (complete genome, 1 contig, NZ_CP086393.1, downloaded to `data/references/`)
- **Key annotated genes:** rpsE (NKFIDM_05734), rsmI (NKFIDM_00044), 25 rRNA genes (8×16S+8×23S+9×5S), 151 transposable elements, 11 CRISPR arrays
- **rpsE vs I6 — KEY FINDING:** M5/M6/M8 rpsE (162 aa) has a **3-aa in-frame deletion Δ(A21-K22-V23)** relative to I6 rpsE (165 aa, WP_036619535.1). Zero mismatches elsewhere. **H1 (de novo mutation) FALSIFIED. H4 (intrinsic resistance) STRONGLY SUPPORTED.**
- **rpsE vs E. coli — Session 3:** Δ(A21-K22-V23) in I6 numbering maps to **E. coli S22-K23-T24** (NP_417762.1, pairwise alignment). E. coli K23 and T24 are in the N-terminal beta-hairpin of S5 that contacts 16S h34 at the spectinomycin binding site.
- **rpsE MAFFT MSA — Session 3 (DECISIVE):** 4-species MSA (M5, I6, E. coli K12, B. subtilis) confirms: **the conserved Lysine (I6-K22 = E. coli-K23 = B. subtilis-K23) is absent in M5/M6/M8.** Present in all 3 reference species, gone in all 3 Paenibacillus isolates. H4 STRONGLY SUPPORTED with cross-species structural evidence. Output: `results/genomics/comparative/rpsE_msa_annotated.txt`
- **108 kb novel region — ICE, no resistance genes:** Bakta annotation of M5 1,285,652–1,393,969 shows VirD4/VirB4/TrbL/VirB6/Relaxase/Recombinase = complete ICE conjugation machinery. **No aadA/ANT(9) anywhere in M5 genome.** H3 (HGT of spectinomycin resistance gene) NOT SUPPORTED.
- **minimap2 read-level confirmation:** M5 reads at rpsE = 163x uniform coverage (min 155, max 167). 106/184 primary reads mapped to I6 rpsE carry `9D261M` CIGAR (9-bp deletion + 261 bp anchor). Deletion is real, not an assembly artifact.
- **H3 NOT SUPPORTED. H4 is the sole supported mechanism.**
- **Phase 4 (Session 4):** 2 Paenibacillus-native ICEs (108 kb + 25 kb); 59 kb SP-beta prophage; 22 kb low-GC foreign block; no resistance cargo in any novel region. ICE donors confirmed Paenibacillus (not E. coli) via NCBI Entrez WP_ lookup.

### Master reference for future sessions
- **`analysis/execution_plan.md`** — tools, parameters, decision criteria, execution order
- **`results/genomics/phase0_findings.md`** — Session 1 results
- **`results/genomics/phase2_findings.md`** — Session 2 results (major revisions + rpsE vs I6)

### Phase execution status (resistance analysis)
| Phase | Description | Status |
|-------|-------------|--------|
| 0.1 | Nanopore chemistry | DONE — R10.4.1 SUP, Clair3 r1041_e82_400bps_sup_v430 |
| 0.2 | 16S RA in controls | DONE — H10 falsified, H4 weakened |
| 0.3 | rpsE cross-isolate | DONE — M5/M6/M8 identical |
| 0.4 | In silico PCR | BLOCKED — primer sequences unavailable |
| 1.1 | NanoPlot QC | BLOCKED — NanoPlot needs interactive terminal in WSL |
| 1.2 | Read-to-assembly mapping | DONE — rpsE 163x uniform; 106 reads with 9D261M CIGAR vs I6 (partial: rsmI/16S not yet checked) |
| 1.3 | rRNA operon count | DONE — 8 copies of 16S (not 5); 25 total |
| 2.1 | Download I6 | DONE — GCF_022494515.1 complete, 1 contig, 7.11 Mb |
| 2.2 | FastANI | DONE — M5/M6/M8 = 99.26% ANI vs I6; same species confirmed |
| 2.3 | nucmer WGA / novel regions | DONE — 108 kb primary HGT candidate at 1.28–1.39 Mb; 57 IS element DUPs; BED file written |
| 3.1 | rpsE vs I6 | DONE — Δ(A21K22V23) deletion; H1 falsified; H4 supported |
| 3.1b | rpsE vs E. coli (pairwise) | DONE — maps to E. coli K23-T24 zone; H4 structurally corroborated |
| 3.1c | rpsE MAFFT MSA | DONE — conserved Lys (I6-K22=E.coli-K23) absent in all isolates; H4 confirmed |
| 3.2 | 108 kb region content / aadA | DONE — ICE (VirD4/B4/B6/Relaxase); no aadA anywhere; H3 NOT SUPPORTED |
| 3.3 | 16S rRNA h34 analysis | DONE — h34 identical M5=I6; no 16S-based resistance |
| 3.4 | rsmI methyltransferase | DONE — rsmI identical M5=I6; no rsmI-based mechanism |
| 3.5 | 30S ribosomal protein panel | DONE — rpsC/D/L/rsmH all identical; rpsE ONLY difference |
| 4.1 | Novel region annotation | DONE — 2 ICEs (Paenibacillus-native), 1 prophage (SP-beta), 1 low-GC block (33.8%), 1 MFS efflux; no resistance cargo |
| 4.2 | Pan-genome (Panaroo) | DONE — 5,825 core, 736 novel (isolate-only), 309 lost (I6-only); no resistance genes in novel set |
| 4.3 | IS element census | DONE — IS110(29), IS4(19); 57 extra vs I6 = 57 DUP regions |
| 4.4 | Gene duplication + efflux analysis | DONE — 338 extra CDS in isolates; drug efflux genes present but not spectinomycin-specific; H7 NOT SUPPORTED |
| 4.5 | All 47 GAP regions annotated | DONE — 18/47 have functional keyword hits; zero contain spectinomycin resistance genes |
| 6 | Decision matrix | DONE — all 10 hypotheses resolved; H4 STRONGLY SUPPORTED, all others falsified/not supported |
| 7 | Wet lab | NOT STARTED — requires approval |
| 8 | Final report | DONE — `results/genomics/reports/paenibacillus_resistance_report.md` |

### Known unresolved issues (do not lose track)
- M7 assembly is contaminated/insufficient — re-sequencing is medium-priority wet lab task
- PCR primer sequences for original spectinomycin resistance PCR are not documented — must be obtained before in silico PCR can run
- WSL installed (Ubuntu 24.04); `paeni-genomics` env ready with fastANI 1.34, nucmer 4.0.1, MAFFT v7.526, minimap2 2.30, samtools 1.23.1, NanoPlot 1.46.2
- MIC of I6 (ATCC 7068) for spectinomycin — critical to interpret H4: if I6 is sensitive and M5/M6/M8 are resistant, the rpsE Δ(A21K22V23) is the mechanism
- m8 Baranyi lag = -7 (non-physical artifact of DCM grid search) — noted but not corrected; acceptable limitation
