# Phase 0 Findings — Session 1 (2026-03-22)

## 0.1 Nanopore chemistry — COMPLETE

- Chemistry: **R10.4.1, e8.2 pore, 400 bps, Dorado SUP v4.3.0**
- GPU: NVIDIA A100 80GB PCIe
- Sequencing date: 2025-08-21
- Applies to: all four isolates (same run: 250821LV_BACT)
- **Clair3 model: `r1041_e82_400bps_sup_v430`**
- Implication: R10.4.1 SUP has ~Q20+ median accuracy with dramatically better
  homopolymer resolution than R9.4.1. The risk of homopolymer-driven false variants
  at rpsE and 16S rRNA is substantially reduced vs the planning-phase assumption.

---

## 0.2 Paenibacillaceae 16S dynamics — COMPLETE

**Result: Paenibacillaceae does NOT expand in antibiotic-only control mice.**

| Group | Max RA | Last-timepoint mean RA |
|-------|--------|------------------------|
| Colonized cohort 1 (m1–m4) | 53–60% | 46–55% |
| Colonized cohort 2 (m5–m8) | 62–67% | 56–59% |
| Antibiotic-only controls (c_m1–c_m4) | 0.3–1.4% | 0–0.6% |

- Controls max: 1.38% (c_m1), 1.04% (c_m2), 0.33% (c_m3), 0.97% (c_m4)
- Colonized mice: all end experiment >46% Paenibacillaceae

**Decision (from execution_plan.md):**
- **H10 (ecological release, spectinomycin alone) is FALSIFIED.** 40-50x difference
  between colonized and control mice means spectinomycin alone does not drive the
  expansion. E. coli colonization is required.
- **H4 (intrinsic resistance) is WEAKENED.** If P. macerans were simply intrinsically
  resistant and filling an antibiotic-cleared niche, expansion would occur in controls
  too. It does not. Something about the E. coli × spectinomycin interaction is
  necessary for the massive Paenibacillaceae bloom.
- **The resistance/HGT hypothesis is strengthened** — Paenibacillus expanding only
  in E. coli-colonized mice is consistent with HGT from E. coli or with selection of
  a resistance mutation in a competitive context requiring E. coli presence.

Figures: `results/genomics/figures/paenibacillus_RA_all_mice.pdf/.png`
         `results/genomics/figures/paenibacillus_RA_controls_only.pdf/.png`

---

## 0.3 rpsE cross-isolate comparison — COMPLETE

**Result: rpsE is IDENTICAL across M5, M6, and M8.**

```
>M5_rpsE
MRVDPNTLELTERVVNINRVVKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTIPHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVNATLEGLSRLKRAEDVAKLRGKTVEELLG

>M6_rpsE  [identical]
>M8_rpsE  [identical]
```

- Length: 163 amino acids
- Locus tags: M5=NKFIDM_05734, M6=KMIGAE_05734, M8=ILIHOP_05734
- Coordinates: all three at contig_1:6,418,100–6,418,588, minus strand
- No RefSeq/WP_ accession in Bakta annotation — confirmed divergent from reference databases

**Decision:**
- **Convergent evolution (three independent rpsE mutations) is RULED OUT.**
  All three mice carry identical rpsE — no evidence of parallel adaptation.
- **Cannot yet distinguish** between:
  (a) rpsE carries an ancestral resistance mutation (present in the M5/M6/M8 common
      ancestor before colonization — compare to I6 in Phase 2/3)
  (b) rpsE is identical to I6 and is NOT the resistance mechanism (H1 false)
- **Priority for Phase 2/3:** rpsE vs I6 comparison is the single most decisive
  analysis remaining.

File: `results/genomics/comparative/rpsE_isolates.faa`

---

## 0.4 In silico PCR — BLOCKED

PCR primer sequences for the original spectinomycin resistance PCR are not documented
anywhere in the project. Cannot proceed with in silico PCR.

**Action required (wet lab / data retrieval):**
- Obtain primer sequences from original experimental methods
- Add to `data/references/spectinomycin_PCR_primers.txt` when available
- Then re-run Phase 0.4 (isPCR or BLAST primers against M5/M6/M8 assemblies)

**This is a critical blocker for H3 validation.**

---

## 1.3 rRNA operon count — COMPLETE (correction)

**CORRECTION: Previous planning documents stated 5 rRNA operons. Actual count: 8.**

| Feature | Count |
|---------|-------|
| 16S rRNA | **8** |
| 23S rRNA | 8 |
| 5S rRNA | 9 |
| **Total** | **25** |

**16S rRNA copy positions:**
```
Copy 1: 30,224–31,777    (+strand)   operon cluster 1 (chromosome start)
Copy 2: 119,334–120,887  (+strand)
Copy 3: 176,589–178,142  (+strand)
Copy 4: 1,161,778–1,163,331 (+strand)  isolated in mid-chromosome
Copy 5: 5,664,298–5,665,863 (-strand)  *** likely in novel 1.77 Mb region ***
Copy 6: 6,207,588–6,209,153 (-strand)  *** likely in novel 1.77 Mb region ***
Copy 7: 6,399,261–6,400,814 (-strand)  *** likely in novel 1.77 Mb region ***
Copy 8: 6,475,824–6,477,377 (-strand)  *** likely in novel 1.77 Mb region ***
```

**Key observation:** Copies 5–8 are clustered at positions 5.66–6.48 Mb on the minus
strand. Since I6 is ~5.7 Mb total, these four copies likely fall within or adjacent to
the 1.77 Mb extra genomic content. This will be confirmed by the nucmer alignment in
Phase 2.3. If copies 5–8 are in the novel region, P. macerans I6 likely has only 4 rRNA
operons; our isolate acquired 4 additional copies as part of a large horizontal transfer.

**Consequence for helix 34 analysis:**
- A mutation in 1/8 copies appears at ~12.5% frequency in read pileups (not 20%)
- With R10.4.1 at 100x coverage, 12.5% is well above the noise floor
- But the extra 4 copies in the novel region will need to be phased separately

BED file: `results/genomics/comparative/M5_rRNA_16S_coordinates.bed` (coordinates only,
locus_tag column blank due to awk parsing issue — IDs: NKFIDM_00026/00104/00158/00994/
05005/05512/05708/05794)

---

## Session 1 Status

| Step | Status | Key result |
|------|--------|-----------|
| 0.1 Chemistry | DONE | R10.4.1 SUP, Clair3 model r1041_e82_400bps_sup_v430 |
| 0.2 16S RA dynamics | DONE | Paenibacillaceae expands only in colonized mice; H10 falsified |
| 0.3 rpsE cross-isolate | DONE | Identical in M5/M6/M8; compare to I6 next |
| 0.4 In silico PCR | BLOCKED | Primer sequences unavailable |
| 1.3 rRNA count | DONE | 8 copies (correction from planning-assumed 5); 4 likely in novel region |

## Next session: Phase 2 (Session 2)

1. Download I6 (GCF_000172175.2) — assess contig count (draft vs. closed)
2. FastANI M5/M6/M8 vs I6 — confirm 99% ANI
3. nucmer + dnadiff — quantify the 1.77 Mb novel region; produce `M5_novel_regions.bed`
4. NanoPlot QC (can run in background while nucmer runs)
5. minimap2 read-to-assembly mapping for coverage validation (background)

The rpsE vs I6 comparison (Phase 3.1) is unlocked as soon as I6 is downloaded.
