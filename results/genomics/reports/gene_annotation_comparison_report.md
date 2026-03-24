# Gene Annotation Comparison Report: P. macerans Isolates vs I6 Reference

**Author:** Melis Gencel
**Date:** 2026-03-23
**Source:** Panaroo pan-genome analysis (strict mode) on M5, M6, M8, and I6

---

## 1. Pan-Genome Structure

| Category | Gene families | % of total |
|----------|--------------|------------|
| Core (all 4 genomes) | 5,825 | 84.8% |
| Novel (isolates only) | 736 | 10.7% |
| Lost (I6 only) | 309 | 4.5% |
| **Total** | **6,870** | **100%** |

All 736 novel genes are present in all three isolates (M5, M6, M8) — zero isolate-specific genes exist. This confirms the colonising population is clonal. The 309 lost genes are present in the I6 reference (ATCC 7068) but absent from all three isolates.

---

## 2. Functional Breakdown of Novel Genes (736 gene families)

Genes acquired by the isolates relative to the I6 reference, classified by function:

| Functional category | Count | % | In novel regions | Key observations |
|---------------------|-------|---|-----------------|-----------------|
| Other (named, miscellaneous) | 233 | 31.7% | 156 | SLH domains, ATPases, mobilisation proteins |
| Hypothetical / DUF | 201 | 27.3% | 105 | DUF5704, DUF554, DUF3991, DUF4411 — function unknown |
| Metabolism / enzymes | 88 | 12.0% | 57 | SDR oxidoreductases, sulfotransferases, kinases |
| IS transposases | 44 | 6.0% | 24 | IS110 (dominant), IS4, Mu-type; drives 57 DUP regions |
| Regulatory | 40 | 5.4% | 26 | MerR, MarR, HTH regulators, response regulators |
| Transport | 37 | 5.0% | 30 | Sugar ABC, MFS, branched-chain amino acid permeases |
| Phage | 31 | 4.2% | 15 | SP-beta prophage structural proteins (capsid, tail, portal) |
| CRISPR | 15 | 2.0% | 15 | Complete type I-B and I-C CRISPR-Cas systems |
| Conjugation / ICE | 12 | 1.6% | 12 | VirD4, TrbL/VirB6, Relaxase, MobA — two complete ICEs |
| Integrase / recombinase | 8 | 1.1% | 8 | Site-specific integrases at ICE/prophage boundaries |
| Toxin-antitoxin | 7 | 1.0% | 6 | HicA/HicB TA systems — plasmid/ICE maintenance |
| Cell envelope / motility | 7 | 1.0% | 2 | LysM peptidoglycan-binding, chemotaxis proteins |
| Sporulation | 6 | 0.8% | 5 | GerKA/KB germination proteins, Cse60 |
| Ribosomal / translation | 5 | 0.7% | 3 | Extra copies of L33, S14; tRNA amidotransferase |
| **Drug efflux** | **2** | **0.3%** | **2** | **MFS drug:H+ antiporter-2, MATE efflux — NOT spectinomycin-specific** |

### Key insights from novel genes:

**IS element proliferation is the dominant genomic change.** The 44 novel IS transposases (plus their flanking genes) account for much of the 356 kb size difference. IS110 family is the most expanded, with 29 copies in M5 vs 7 in I6 (22 extra). These elements have inserted throughout the chromosome, creating the 57 DUP regions.

**Two complete ICEs provide conjugation capability.** ICE #1 (108 kb, Brevibacillus agri origin) and ICE #2 (25 kb, Lacrimispora saccharolytica origin) carry full conjugation machinery (VirD4, VirB4, TrbL/VirB6, Relaxase) but zero resistance cargo. The HicA/HicB toxin-antitoxin systems likely serve as ICE maintenance modules.

**SP-beta prophage is a known Paenibacillus phage.** The 59 kb prophage encodes 31 novel phage-related proteins (capsid, tail, portal, terminase) and matches Paenibacillus polymyxa phage SP-beta at 72%/66% identity/coverage.

**Two novel CRISPR-Cas systems (type I-B and I-C)** were acquired, with 15 genes. These provide adaptive immunity against foreign DNA — potentially protecting the ICEs from phage interference.

**Only 2 drug efflux genes** are novel — an MFS drug:H+ antiporter-2 and a MATE efflux family protein. Neither is spectinomycin-specific. This rules out efflux-based resistance as the explanation (H7 NOT SUPPORTED).

---

## 3. Functional Breakdown of Lost Genes (309 gene families)

Genes present in I6 but absent from all three isolates:

| Functional category | Count | % | Key observations |
|---------------------|-------|---|-----------------|
| Hypothetical / DUF | 89 | 28.8% | DUF5643, DUF4179, DUF1788 — unknown function |
| Other (named, misc) | 68 | 22.0% | Ankyrin repeats, condensation domains, VOC family |
| Metabolism / enzymes | 53 | 17.2% | SDR oxidoreductases, glycosyltransferases, methyltransferases |
| Transport | 45 | 14.6% | Carbohydrate ABC permeases, MFS, PTS sugar transporters |
| Regulatory | 21 | 6.8% | HTH regulators, PadR, MarR, MerR families |
| CRISPR | 15 | 4.9% | I6 has its own CRISPR-Cas system that differs from the isolates' |
| IS transposases | 12 | 3.9% | I6-specific transposases (different IS families) |
| Sporulation | 2 | 0.6% | GerPD/PB germination proteins |
| Toxin-antitoxin | 2 | 0.6% | Fic/DOC TA system (different from isolates' HicA/B) |
| Phage | 1 | 0.3% | Holin-like toxin |
| Integrase | 1 | 0.3% | Tyrosine recombinase/integrase |

### Key insights from lost genes:

The lost gene set is dominated by hypotheticals (28.8%) and metabolic functions (17.2%). The loss of 45 transport genes may reflect adaptation to the gut environment where certain substrates are unavailable. Notably, I6 and the isolates have **different CRISPR-Cas systems** (15 genes each, non-overlapping) — this represents a CRISPR system replacement event, not simple gain or loss.

No resistance-relevant genes are among the 309 lost genes.

---

## 4. Enrichment Analysis: What is Over-Represented in Novel vs Core?

Categories enriched (>2x proportion) in the novel gene set compared to the core genome:

| Category | % in novel | % in core | Enrichment |
|----------|-----------|----------|------------|
| Hypothetical / DUF | 27.3% | 11.6% | **2.4x** |
| IS transposase | 6.0% | 1.3% | **4.6x** |
| Phage | 4.2% | 0.9% | **4.7x** |
| CRISPR | 2.0% | 0.1% | **20x** |
| Conjugation / ICE | 1.6% | 0.2% | **8x** |
| Toxin-antitoxin | 1.0% | 0.3% | **3.3x** |

Categories depleted in the novel gene set:
| Category | % in novel | % in core | Depletion |
|----------|-----------|----------|-----------|
| Metabolism / enzymes | 12.0% | 25.6% | 0.47x |
| Transport | 5.0% | 10.5% | 0.48x |
| Ribosomal / translation | 0.7% | 2.6% | 0.27x |
| Sporulation | 0.8% | 2.7% | 0.30x |

The enrichment pattern is consistent with mobile element acquisition: IS transposases, phage genes, ICE conjugation machinery, CRISPR defence systems, and toxin-antitoxin maintenance modules are all hallmarks of horizontally acquired DNA. The depletion of metabolic and ribosomal genes confirms that the novel content is not a chromosomal duplication but rather foreign genetic material integrated at specific sites.

---

## 5. Relevance to Spectinomycin Resistance

None of the 736 novel genes encode spectinomycin resistance determinants:

- **No aadA** (spectinomycin adenylyltransferase) — the most common spectinomycin resistance gene
- **No ANT(9)** family genes
- **No novel aminoglycoside-modifying enzymes** of any kind
- **Only 2 drug efflux genes** — neither spectinomycin-specific
- **No mutations in novel ribosomal protein copies** — the extra L33 and S14 copies are standard

The spectinomycin resistance mechanism remains the intrinsic rpsE Δ(A21-K22-V23) deletion (H4), which is a core genome variant, not part of the novel gene set. The mobile genetic elements (ICEs, phage, IS elements) represent active genome plasticity but did not deliver any resistance cargo during this experiment.

---

## 6. Figures

The following figures accompany this analysis:

| Figure | File | Description |
|--------|------|-------------|
| Fig 1 | `results/figures/genomics/fig1_rpsE_msa.png` | rpsE MAFFT MSA showing Δ(A21-K22-V23) deletion |
| Fig 2 | `results/figures/genomics/fig2_paenibacillaceae_RA.png` | Paenibacillaceae relative abundance dynamics |
| Fig 3 | `results/figures/genomics/fig3_genome_comparison.png` | Linear genome comparison M5 vs I6 |
| Fig 4 | `results/figures/genomics/fig4_panaroo_heatmap.png` | Gene presence/absence heatmap |
| Fig 5 | `results/figures/genomics/fig5_rpsE_pileup.png` | rpsE read-level coverage and CIGAR |
| Fig 6 | `results/figures/genomics/fig6_pangenome_functional.png` | Pan-genome functional breakdown |
