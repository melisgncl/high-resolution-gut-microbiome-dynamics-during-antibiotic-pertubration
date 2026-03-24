# Session 5 Plan — Pan-Genome, Gene Duplication, Final Report
**Date:** 2026-03-23
**Continues from:** Session 4 (commit 9fa728d)

---

## Objectives

1. **Phase 4.2 — Panaroo pan-genome analysis** (M5/M6/M8 vs I6)
2. **Gene inventory + duplication analysis** (ancestral vs acquired, efflux pumps, H7)
3. **Annotate ALL 47 GAP novel regions** from Bakta TSV (extend Phase 4.1)
4. **Write comprehensive final report** with all 10 hypotheses (H1–H10)

---

## Step-by-step execution plan

### Phase 4.2: Panaroo

| Step | Action | Input | Output |
|------|--------|-------|--------|
| 1 | Install Panaroo 1.6.0 in WSL `paeni-genomics` | conda | panaroo binary |
| 2 | Verify I6 GFF has embedded FASTA; if not, append from FNA | `GCF_022494515.1_genomic.gff` + `.fna` | I6 GFF3 ready |
| 3 | Run Panaroo `--clean-mode strict -t 4` on M5/M6/M8/I6 GFF3 | 4 GFF3 files | `panaroo_output/` |
| 4 | Parse `gene_presence_absence.csv` for novel/lost genes | Panaroo output | gene lists |

### Gene inventory + duplication

| Step | Action | Input | Output |
|------|--------|-------|--------|
| 5 | Count total CDS in M5/M6/M8 (Bakta TSV) and I6 (NCBI GFF) | TSV + GFF | gene counts |
| 6 | Find duplicated genes within each genome by product name + >90% seq identity | Bakta FAA + I6 FAA | duplication table |
| 7 | Compare: shared duplications (ancestral) vs isolate-only (acquired) | Step 6 output | classification |
| 8 | Efflux pump deep-dive: count, duplication status, IS proximity (<5 kb) | Bakta TSV + BED | H7 assessment |
| 9 | Cross-reference duplicated genes against 47 GAP novel regions | BED + duplication table | HGT vs tandem |

### Novel region annotation (all 47 GAPs)

| Step | Action | Input | Output |
|------|--------|-------|--------|
| 10 | Extend Phase 4.1 to all 47 GAP regions from Bakta TSV | `M5_novel_regions.bed` + Bakta TSV | full annotation |
| 11 | Keyword screen for efflux/ICE/plasmid/AMR across all regions | Step 10 | flagged regions |

### Final report

| Step | Action | Output |
|------|--------|--------|
| 12 | Write `results/genomics/reports/paenibacillus_resistance_report.md` | Final report |
| 13 | Include all 10 hypotheses (H1–H10) with status + evidence | Hypothesis table |

---

## Data files available

- **M5 GFF3:** `data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.gff3`
- **M6 GFF3:** `data/sequence_data/TPHP3P_results/TPHP3P_2_M6__+_/annotation/TPHP3P_2_M6__+_.gff3`
- **M8 GFF3:** `data/sequence_data/TPHP3P_results/TPHP3P_4_M8__+_/annotation/TPHP3P_4_M8__+_.gff3`
- **I6 GFF:** `data/references/GCF_022494515.1_genomic.gff`
- **I6 FNA:** `data/references/GCF_022494515.1_genomic.fna`
- **I6 FAA:** `data/references/GCF_022494515.1_protein.faa`
- **M5 Bakta TSV:** `data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv`
- **M5 novel regions BED:** `results/genomics/comparative/M5_novel_regions.bed`
- **M5 FAA:** `data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.faa`

## Decisions made

- Install Panaroo in **existing** `paeni-genomics` env (not new env)
- Use **Bakta TSV annotation** for novel region characterization (not NCBI BLAST)
- Include **all 10 hypotheses** (H1–H10) in final report
- Exclude M7 (unreliable 9x/7-contig assembly)
