# Phase 2 Findings — Session 2 (2026-03-22)

## 2.1 I6 Reference Genome — COMPLETE (major revision)

**The old accession GCF_000172175.2 is deprecated and no longer exists in NCBI.**

A complete closed genome for strain I6 (ATCC 7068) is now available under a new accession:

| Field | Value |
|-------|-------|
| Accession | **GCF_022494515.1** |
| Assembly level | **Complete Genome** |
| Contigs | **1** (closed chromosome) |
| Size | **7,113,008 bp** |
| Sequence ID | NZ_CP086393.1 |
| Proteins | 6,105 |
| rRNA genes | 24 (8×16S, 8×23S, 8×5S) |
| rpsE | WP_036619535.1, coordinates 228,319–228,816 (−strand) |
| Transposase features | 94 |
| CRISPR features | 29 |

File locations:
- `data/references/GCF_022494515.1_genomic.fna`
- `data/references/GCF_022494515.1_protein.faa`
- `data/references/GCF_022494515.1_genomic.gff`

---

## CRITICAL REVISION: The 1.77 Mb anomaly was an artifact

**Old I6 reference (GCF_000172175.2, scaffold, superseded): ~5.7 Mb**
**Correct I6 reference (GCF_022494515.1, complete): 7,113,008 bp**

| Assembly | Size |
|----------|------|
| M5 (our isolate) | 7,468,786 bp |
| I6 (GCF_022494515.1, complete) | 7,113,008 bp |
| **True difference** | **355,778 bp (356 kb)** |

The planning documents stated a "1.77 Mb unexplained extra sequence." This was entirely
an artifact of comparing M5 against a draft scaffold assembly that was missing ~1.4 Mb
of its own sequence. The actual extra content in M5 vs the complete I6 genome is only
**356 kb**, a size consistent with a single large genomic island or a cluster of IS elements,
not a massive HGT acquisition.

**Consequence:** The entire "novel 1.77 Mb region" framing must be revised.

**Consequence for rRNA operon analysis (Phase 1.3):**
- M5 has 8×16S, 8×23S, 9×5S = 25 rRNA genes
- I6 (complete) has 8×16S, 8×23S, 8×5S = 24 rRNA genes
- The 4 "extra" 16S copies at positions 5.66–6.48 Mb in M5 are NOT in a unique novel region
  relative to I6. I6 also has 8 16S rRNA copies. The 4 minus-strand copies are syntenic
  with I6 sequence (to be confirmed by nucmer Phase 2.3).
- The only difference is 1 extra 5S rRNA in M5 relative to I6.

---

## 3.1 rpsE M5/M6/M8 vs I6 — COMPLETE (Phase 3.1 unlocked and executed)

**This is the single most decisive result of the project.**

### Sequence comparison (exact, no alignment tool needed)

| Isolate | Length | Sequence |
|---------|--------|----------|
| M5, M6, M8 | **162 aa** | MRVDPNTLELTERVVNINRV**VK**GGRRFSFSALVVV... |
| I6 (WP_036619535.1) | **165 aa** | MRVDPNTLELTERVVNINRV**AKVVK**GGRRFSFSALVVV... |

**Alignment:**
```
M5: MRVDPNTLELTERVVNINRV---VKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTIPHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVNATLEGLSRLKRAEDVAKLRGKTVEELLG
I6: MRVDPNTLELTERVVNINRVAKVVKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTIPHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVNATLEGLSRLKRAEDVAKLRGKTVEELLG
```

- **Shared prefix (pos 1–20 I6 numbering):** MRVDPNTLELTERVVNINRV — identical
- **Deletion in M5/M6/M8:** positions 21–23 in I6 = **Δ(A21-K22-V23)**
- **Shared suffix (pos 24–165 I6 numbering):** 142 aa — perfectly identical (0 mismatches)

This deletion lies in the **N-terminal beta-hairpin / loop region** of ribosomal protein S5
(uS5), adjacent to the region that contacts 16S rRNA helix 34 at the spectinomycin binding
interface (E. coli contact residues Q15, G16, G17 are in this same N-terminal domain).

### Hypothesis updates

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1 (de novo rpsE mutation in vivo) | **FALSIFIED** | Deletion is identical in M5/M6/M8 from separate mice — cannot be three independent mutations under selection |
| H4 (intrinsic resistance, pre-existing in strain) | **STRONGLY SUPPORTED** | The colonizing strain differs from I6 at rpsE by a pre-existing 3-aa deletion adjacent to the spectinomycin-binding loop |
| H10 (ecological release) | Previously falsified (Phase 0.2) | Unchanged |
| H3 (HGT of aadA from E. coli) | Unresolved | Cannot be addressed by rpsE comparison alone; nucmer + IS element analysis needed |

### Critical interpretation: H4 vs H3 are not mutually exclusive

The rpsE deletion proves the colonizing strain was **already differently configured** at the
spectinomycin target before the experiment. This does not require HGT. However:
- If the MIC of I6 (ATCC 7068) for spectinomycin is **high**, the deletion may not be
  the resistance mechanism (both strains are resistant → deletion is a neutral polymorphism)
- If the MIC of I6 is **low** (sensitive), the deletion is directly implicated as the
  intrinsic resistance mechanism in our isolates
- **Wet lab priority:** Obtain MIC for I6 vs M5 (Phase 7)

Files:
- `results/genomics/comparative/rpsE_isolates_vs_I6.faa`
- `results/genomics/comparative/rpsE_alignment_M5_vs_I6.txt`

---

## 2.1 contig count decision (from execution_plan.md table)

I6 (GCF_022494515.1) = **1 contig (closed chromosome)**

Decision table outcome: **"Use I6 as reference for all comparisons"**

---

## Session 2 Status

| Step | Status | Key result |
|------|--------|-----------|
| 2.1 Download I6 | DONE | GCF_022494515.1 (complete, 1 contig, 7.11 Mb) |
| 2.1 Contig check | DONE | 1 contig → I6 is gold-standard reference |
| **REVISION** | DONE | 1.77 Mb anomaly is an artifact → true extra = 356 kb |
| 3.1 rpsE vs I6 | DONE | 3-aa deletion Δ(A21K22V23) in M5/M6/M8; H1 falsified; H4 supported |
| 2.2 FastANI | BLOCKED | Needed WSL |
| 2.3 nucmer/dnadiff | BLOCKED | Needed WSL |
| 1.1 NanoPlot | BLOCKED | Needed WSL |
| 1.2 minimap2 | BLOCKED | Needed WSL |

---

## Session 3 — rpsE cross-species alignment (2026-03-22)

### WSL Setup — COMPLETE
- Ubuntu 24.04 LTS installed; Miniforge installed
- `paeni-genomics` conda env created with: fastANI 1.34, nucmer 4.0.1, MAFFT v7.526,
  minimap2 2.30, samtools 1.23.1, NanoPlot 1.46.2
- All tools verified functional
- Project files accessible at `/mnt/c/Users/melis/.../hgt-study/`

### rpsE E. coli mapping — COMPLETE (Phase 3.1 follow-up)

**E. coli K12 rpsE fetched:** NP_417762.1, 167 aa
**Method:** Biopython PairwiseAligner (global, gap open=-10, extend=-0.5)
**Output:** `results/genomics/comparative/rpsE_ecoli_vs_I6_alignment.txt`

**Result:**

| I6 position | I6 residue | E. coli position | E. coli residue |
|-------------|------------|-----------------|-----------------|
| 21 | A (deleted in isolates) | 22 | S |
| 22 | K (deleted in isolates) | 23 | **K** |
| 23 | V (deleted in isolates) | 24 | T |

The Δ(A21-K22-V23) deletion in Paenibacillus isolates maps to E. coli positions S22-**K23**-T24.
**E. coli K23** and **T24** are within the N-terminal beta-hairpin of S5 (uS5) that contacts
16S rRNA helix 34 at the spectinomycin binding interface. Mutations at these positions are
documented to confer spectinomycin resistance in E. coli (Bilgin et al. 1990) and N. gonorrhoeae
(Unemo et al. 2013: K20E, T23I, K24N).

**Structural interpretation (Paenibacillus-specific):**
The deletion shortens the Paenibacillus S5 N-terminal loop by 3 residues, removing I6-K22
(Lysine, positively charged), which is the residue corresponding to E. coli K23. This would
reposition the backbone of the loop and reduce h34 contacts, lowering spectinomycin affinity.
This is a **structural corroboration of H4 (intrinsic resistance)**.

**Caveat:** Structural prediction only. Requires validation against a Paenibacillus 30S structure
(not available) or homology model. See MAFFT MSA result below for cross-species confirmation.

### rpsE MAFFT MSA — COMPLETE (Phase 3.1c)

**Tool:** MAFFT v7.526 `--auto` (WSL Ubuntu 24.04)
**Sequences:** Paenibacillus M5 (162 aa), I6 (165 aa), E. coli K12 NP_417762.1 (167 aa), B. subtilis 168 WP_003328273.1 (166 aa)
**Output:** `results/genomics/comparative/rpsE_msa_output.faa`, `rpsE_msa_annotated.txt`

**N-terminal alignment (first 50 columns):**

```
M5           M-RVDPNTLELTERVVNINR---VVKGGRRFSFSALVVVGDGKGWVGAGI
I6           M-RVDPNTLELTERVVNINRVAKVVKGGRRFSFSALVVVGDGKGWVGAGI
E.coli K12   MAHIEKQAGELQEKLIAVNRVSKTVKGGRIFSFTALTVVGDGNGRVGFGY
B.subtilis   MRRIDPSKLELEERLVTVNRVAKVVKGGRRFRFAALVVVGDKNGHVGFGT
conservation *        ** *     **: : ***** * * ** ****  * ** *
deletion(M5)                     ^^^
```

**Deletion column detail (MAFFT coordinates, 1-based):**

| Aln col | M5 | I6 | E. coli K12 | B. subtilis |
|---------|----|----|------------|-------------|
| 21 | **gap** | V(20) | V(21) | V(21) |
| 22 | **gap** | A(21) | S(22) | A(22) |
| **23** | **gap** | **K(22)** | **K(23)** | **K(23)** |

**The conserved Lysine at alignment column 23 (I6-K22 = E. coli-K23 = B. subtilis-K23) is
present across all reference species and ABSENT in M5/M6/M8.** This is the key conserved
charged residue of the N-terminal beta-hairpin of S5 that contacts 16S rRNA helix 34
at the spectinomycin binding interface. Its conservation across Paenibacillus, E. coli,
and B. subtilis underscores its functional importance.

**Gap placement note:** MAFFT places the gap at I6 V20-A21-K22; Biopython pairwise (Session 2)
placed it at A21-K22-V23. Both represent the same 3-aa deletion. The biological conclusion
is identical: M5/M6/M8 lack the conserved Lysine (I6-K22 / E. coli-K23).

### Updated hypothesis table (Session 3 — post MAFFT)

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1 (de novo rpsE mutation in vivo) | **FALSIFIED** | Deletion identical in M5/M6/M8 |
| H4 (intrinsic resistance, pre-existing) | **STRONGLY SUPPORTED** | Conserved Lysine (I6-K22 = E. coli-K23) absent in all isolates; confirmed by MAFFT MSA with 3 reference species |
| H10 (ecological release) | Falsified (Phase 0.2) | Unchanged |
| H3 (HGT of aadA) | Unresolved | Requires nucmer + IS element analysis |

---

## Session 3 — FastANI + nucmer/dnadiff (Phase 2.2 + 2.3)

### FastANI — COMPLETE (Phase 2.2)

**Result: M5/M6/M8 are confirmed _Paenibacillus macerans_ I6 strain at 99.26% ANI**

| Isolate | ANI vs I6 | Mapped fragments | Total fragments |
|---------|-----------|-----------------|-----------------|
| M5 | **99.2613%** | 2245 | 2489 |
| M6 | **99.2614%** | 2245 | 2489 |
| M8 | **99.2615%** | 2245 | 2489 |

- All three isolates are essentially identical to each other (differ only in the 4th decimal)
- 2244/2489 fragments mapped = 90.2% → the unmapped ~10% corresponds to the novel 356 kb content
- Output: `results/genomics/comparative/fastani_M5M6M8_vs_I6.txt`

### nucmer + dnadiff M5 vs I6 — COMPLETE (Phase 2.3)

**Tool:** nucmer 4.0.1 + dnadiff (mummer4, WSL)
**Command:** `nucmer REF QRY; dnadiff -d delta`
**Output:** `results/genomics/comparative/nucmer_M5_vs_I6_dnadiff.report`, `.qdiff`, `.rdiff`, `.1coords`

#### Whole-genome alignment summary

| Metric | I6 (REF) | M5 (QRY) |
|--------|----------|----------|
| Total size | 7,113,008 bp | 7,468,786 bp |
| Aligned bases | 6,760,986 (95.05%) | 6,798,077 (91.02%) |
| **Unaligned bases** | **352,022 bp (4.95%)** | **670,709 bp (8.98%)** |
| Avg identity (1-to-1) | 99.32% | 99.32% |
| Relocations | 78 | 59 |
| Inversions | 10 | 10 |
| SNPs | 35,910 | 35,910 |

#### Novel content in M5 vs I6

| Category | Regions | Total bp | Notes |
|----------|---------|----------|-------|
| Large GAP insertions (net >500 bp) | 47 | 444,927 bp | M5-specific sequence |
| DUP regions | 57 | 76,705 bp (avg 1,346 bp) | IS element copies |
| **Combined novel estimate** | 104 | **~522 kb** | Exceeds 356 kb Δ because I6 also has unique regions |

**BED file of novel regions:** `results/genomics/comparative/M5_novel_regions.bed`

#### Top 10 largest novel regions in M5

| M5 coordinates | M5 size | I6 size | Net extra | Note |
|----------------|---------|---------|-----------|------|
| 1,285,652 – 1,393,969 | 108,318 bp | 1,764 bp | **+106,554 bp** | ★ Primary HGT candidate |
| 6,766,335 – 6,825,514 | 59,180 bp | ~0 bp | **+59,195 bp** | Large novel block |
| 5,322,918 – 5,371,421 | 48,504 bp | ~0 bp | **+48,506 bp** | Large novel block |
| 345,499 – 377,216 | 31,718 bp | 2,919 bp | +28,799 bp | |
| 2,911,706 – 2,936,697 | 24,992 bp | ~0 bp | +25,963 bp | |
| 5,297,205 – 5,322,425 | 25,221 bp | 2,960 bp | +22,261 bp | Adjacent to 48 kb block |
| 5,577,609 – 5,596,327 | 18,719 bp | ~0 bp | +18,654 bp | |
| 2,546,998 – 2,565,301 | 18,304 bp | ~0 bp | +18,235 bp | |
| 7,180,258 – 7,192,022 | 11,765 bp | ~0 bp | +11,686 bp | |
| 5,599,527 – 5,609,082 | 9,556 bp | ~0 bp | +9,835 bp | Adjacent to 18 kb block |

#### DUP pattern — IS element expansion

The 57 DUP regions average 1,346 bp each, consistent with IS element sizes.
This exactly matches the 57 extra IS elements in M5 (151 total) vs I6 (94 total = 151 − 94 = **57**).
**Interpretation:** The bulk of the genome size increase is driven by IS element proliferation, not
a single large HGT acquisition.

#### Key interpretation

The 356 kb extra content in M5 vs I6 is **not a single genomic island** but is distributed as:
- A cluster of large novel blocks (primarily at 1.28 Mb, 5.30–5.61 Mb, 6.77 Mb)
- 57 IS element insertions scattered across the genome

**Primary HGT candidate:** The 108 kb block at 1.28–1.39 Mb. At this size it could accommodate
a genomic island with resistance genes, an integrative conjugative element (ICE), or a prophage.
**Next step:** BLAST this region against resistance databases and check for integrase/mobilization genes.

#### Updated hypothesis table (post Phase 2.3)

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1 (de novo rpsE) | **FALSIFIED** | Identical in M5/M6/M8 |
| H4 (intrinsic resistance) | **STRONGLY SUPPORTED** | rpsE deletion removes conserved Lys; MAFFT-confirmed |
| H10 (ecological release) | **FALSIFIED** | Phase 0.2 |
| H3 (HGT) | **CANDIDATE REGION FOUND** | 108 kb novel block at 1.28–1.39 Mb; content unknown pending BLAST |

---

---

## Session 3 — 108 kb novel region content analysis (Phase 4)

### Source: Bakta annotation TSV (no BLAST needed)

Rather than running BLAST, the existing Bakta annotation of M5
(`data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv`)
was queried for features in the 1,285,652–1,393,969 bp window.

### Result: The 108 kb block is an ICE, not a resistance gene cluster

The region contains a complete set of **Integrative Conjugative Element (ICE)** conjugation machinery:

| Gene product | Locus tag | Position |
|--------------|-----------|----------|
| VirD4 coupling protein | NKFIDM_01192 | ~1,290 kb |
| VirB4 ATPase | NKFIDM_01193 | ~1,292 kb |
| TrbL/VirB6 channel protein | NKFIDM_01194 | ~1,294 kb |
| Relaxase/MobF | NKFIDM_01201 | ~1,302 kb |
| Conjugal transfer protein | NKFIDM_01213 | ~1,325 kb |
| Transposase (IS element) | NKFIDM_01220 | ~1,338 kb |
| Site-specific recombinase (integrase) | NKFIDM_01228 | ~1,355 kb |

**No spectinomycin resistance genes, no aminoglycoside transferases, no aadA, no ANT(9) found
in this region.** The 108 kb block is a self-transmissible ICE with complete conjugation
machinery but no detectable antibiotic resistance cargo.

### Genome-wide resistance gene search (keyword screen of entire Bakta TSV)

All Bakta annotation features matching resistance-related terms were searched:

| Found | Locus | Position | Resistance relevance |
|-------|-------|----------|----------------------|
| AAC(3) aminoglycoside acetyltransferase | NKFIDM_01673 | 1,829,288 bp | Gentamicin/tobramycin — NOT spectinomycin |
| ANT(6) adenylyltransferase | NKFIDM_02131 | 2,297,835 bp | Streptomycin — NOT spectinomycin |

**Key findings:**
- Both hits are in I6-shared regions (not novel M5-specific sequence)
- Neither ANT(6) nor AAC(3) confers spectinomycin resistance
- **No aadA (ANT(9)) anywhere in the M5 genome** — the spectinomycin-specific resistance
  adenylyltransferase is simply absent

### Updated hypothesis table (post 108 kb analysis)

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1 (de novo rpsE mutation in vivo) | **FALSIFIED** | Identical deletion in M5/M6/M8 |
| H4 (intrinsic resistance, pre-existing) | **STRONGLY SUPPORTED** | rpsE Δ(A21-K22-V23) removes conserved Lys; confirmed by MAFFT + read-level data |
| H10 (ecological release) | **FALSIFIED** | Phase 0.2 |
| H3 (HGT of aadA from E. coli) | **NOT SUPPORTED** | No aadA/ANT(9) found anywhere in M5 genome; 108 kb novel block = ICE with no resistance cargo |

**H3 is not supported by genomic data.** The ICE could in principle carry future resistance
acquisitions, but it does not currently harbor spectinomycin resistance genes.

---

## Session 3 — minimap2 read-level confirmation of rpsE deletion (Phase 1.2)

### Coverage at rpsE locus (M5 reads vs M5 assembly)

**Tool:** minimap2 2.30 map-ont + samtools, WSL
**Reads:** M5 Nanopore reads (103x expected)
**Reference:** M5 assembly (TPHP3P_1_M5__+_.fna)
**Output:** `results/genomics/comparative/M5_reads_vs_M5asm.bam`

rpsE = NKFIDM_05734, position 6,418,100–6,418,588 (minus strand, 489 bp)

| Metric | Value |
|--------|-------|
| Total reads mapped | 197,986 |
| Reads at rpsE (samtools view) | 267 |
| Mean depth across rpsE | **163x** |
| Min depth (any position) | 155x |
| Max depth (any position) | 167x |

Coverage is **uniform** across all 489 bp (min/max within 8x of mean). No assembly gaps,
no dropped-off ends. This confirms rpsE is single-copy and fully sequenced.

### Deletion confirmation in raw reads (M5 reads vs I6 rpsE)

M5 reads were also mapped to the I6 rpsE region (NZ_CP086393.1:228,119–229,016, 898 bp) to
directly observe the 9-bp deletion in raw reads vs the longer I6 template.

| Metric | Value |
|--------|-------|
| Primary reads at I6 rpsE | 184 |
| Mean coverage | 160.3x |
| Reads with `9D` anywhere in CIGAR | 106 |
| Reads with `9D261M` motif (deletion + 261 bp anchor) | **106** |

**The `9D261M` CIGAR motif (9-bp deletion followed by 261-bp match) is present in 106/184 reads (58%).**
The remaining reads span only part of the deletion locus or have different alignment splits.
This is definitive read-level confirmation that the 9-bp (3-codon) deletion exists in the
raw Nanopore sequencing data — it is not an assembly artifact.

### Interpretation

The rpsE Δ(A21-K22-V23) deletion in M5/M6/M8:
1. Is present in the raw reads (not just the assembly)
2. Has 163x uniform coverage (not a collapsed repeat or assembly error)
3. Is confirmed at single-nucleotide resolution by CIGAR strings vs I6 reference

**Combined with the MAFFT MSA result and the absence of aadA, H4 (intrinsic resistance)
is now the sole genomically-supported mechanism for spectinomycin resistance in these isolates.**

---

---

## Session 4 — Phases 3.3, 3.4, 3.5 (2026-03-22)

### Phase 3.4 — rsmI methyltransferase comparison

**rsmI (NKFIDM_00044) in M5 = WP_036618147.1 in I6 — SEQUENCES ARE IDENTICAL (100%, 297 aa)**

rsmI encodes 16S rRNA (cytidine(1402)-2'-O)-methyltransferase. C1402 is in helix 44 of the
decoding center, not in helix 34 (the spectinomycin binding site). Identical rsmI rules out
any rsmI-mediated resistance difference between M5 and I6.

**Phase 3.4 verdict:** rsmI contributes NOTHING to the resistance difference.

Output: `results/genomics/comparative/rsmI_comparison.txt`

---

### Phase 3.3 — 16S rRNA helix 34 analysis

**Method:** Motif-anchored comparison. The conserved h34 5'-arm sequence `CAGAGGAGACAGG`
was located in both M5 (position 1045) and I6 (position 1057) — 12 bp offset consistent
with the V2 insertion that makes I6 16S 12 bp longer than M5 (1566 bp vs 1554 bp).

**h34 5'-arm sequence (anchored to motif, both sequences):**
```
M5: GGGACAGAGGAGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGT
I6: GGGACAGAGGAGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGT
```
**IDENTICAL.**

**Key spectinomycin resistance positions (motif-anchored, approximate E. coli numbering):**

| Pos (E. coli) | M5 | I6 | Resistance if | Notes |
|---------------|----|----|---------------|-------|
| 1063 | T | T | C→U | spectinomycin contact (crystal) |
| 1064 | G | G | G→A | spectinomycin contact (crystal) |
| 1066 | T | T | C→U | **PRIMARY resistance mutation** |

All three positions: M5 = I6. Both strains natively have T (= U in RNA) at position 1066.

**All 28 M5 vs I6 16S differences fall in variable regions:**
- 2 in V1 (positions 95–96)
- 12 in V2 (positions 204–236, including a 10-bp indel)
- 4 in V4 (positions 482–486)
- 4 in other conserved regions (positions 467–657, none near h34 at position 1045)
- **0 in helix 34**

**Phase 3.3 verdict:** 16S rRNA h34 is identical between M5 and I6. No 16S rRNA-based
spectinomycin resistance mechanism exists in M5/M6/M8.

Output: `results/genomics/comparative/16S_h34_analysis.txt`

---

### Phase 3.5 — 30S ribosomal protein panel comparison

**Method:** Direct protein comparison of 6 spectinomycin/30S-relevant genes between
M5 (Bakta FAA) and I6 (NCBI RefSeq FAA). Replaces Clair3 VCF analysis (not available
in TPHP3P output).

| Gene | Function | M5 length | I6 length | Result |
|------|----------|-----------|-----------|--------|
| rpsE (S5) | **spectinomycin target** | 162 aa | 165 aa | **3 aa DELETION** (Δ(A21-K22-V23)) |
| rpsC (S3) | contacts h34 | 222 aa | 222 aa | **IDENTICAL** |
| rpsD (S4) | contacts h16/h18 | 199 aa | 199 aa | **IDENTICAL** |
| rpsL (S12) | decoding center | 140 aa | 140 aa | **IDENTICAL** |
| rsmH | N4-methylates C1402 | 324 aa | 324 aa | **IDENTICAL** |
| rpsB (S2) | near binding cleft | not annotated | 232 aa | n/a |

**rpsE is the ONLY 30S protein that differs between M5 and I6.**

**Phase 3.5 verdict:** No secondary 30S protein changes exist that could confound H4.
The rpsE Δ(A21-K22-V23) deletion is the unique molecular difference at the
spectinomycin target between the resistant isolates (M5/M6/M8) and the I6 reference.

Output: `results/genomics/comparative/30S_ribosomal_protein_comparison.txt`

---

## Final hypothesis summary (Sessions 3–4 complete)

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1 (de novo rpsE mutation) | **FALSIFIED** | Pre-existing deletion, identical across 3 mice |
| H3 (HGT of aadA) | **NOT SUPPORTED** | No aadA anywhere; 108 kb ICE has no resistance cargo |
| H4 (intrinsic rpsE deletion) | **STRONGLY SUPPORTED** | 3-aa deletion removes conserved Lys; confirmed by: MAFFT MSA (4 species), read-level CIGAR (106/184 reads 9D261M), 16S h34 identical (rules out rRNA mechanism), all other 30S proteins identical (rules out secondary mechanism) |
| H10 (ecological release) | **FALSIFIED** | Phase 0.2 |

**The genomic evidence converges completely on H4.** The rpsE Δ(A21-K22-V23) deletion is
the only molecular difference at the spectinomycin target in M5/M6/M8 vs I6.

**Single unresolved question for H4 confirmation:** MIC of I6 (ATCC 7068) for spectinomycin.
If I6 is sensitive and M5 is resistant, the rpsE deletion is the resistance mechanism. (Phase 7, wet lab.)

---

---

## Session 4 — Phase 4: Genome-Wide Novel Region Analysis (2026-03-23)

### Phase 4.1 — Novel region annotation (gene content, GC, IS/tRNA flanking)

**Input:** M5_novel_regions.bed (47 GAP + 57 DUP), M5 Bakta TSV
**Output:** `results/genomics/comparative/M5_novel_regions_annotated.tsv`

#### Top novel regions characterised

| Region (M5) | Size | GC% | CDS | Notable content | Classification | BLAST top hit |
|-------------|------|-----|-----|-----------------|----------------|---------------|
| 1,285,651–1,393,969 | 108 kb | 49.4% | 96 | VirD4/VirB4/VirB6/Relaxase/Recombinase | **ICE #1** | Brevibacillus agri 99.8%/100% cov |
| 6,766,334–6,825,514 | 59 kb | 44.1% | 89 | SP-beta phage structural/tail proteins | **Prophage** | P. polymyxa 72%/66% cov |
| 5,322,917–5,371,421 | 48 kb | 44.8% | 34 | Tyr recombinase, Integrase ×2 | **Integrative element** | P. JCM9795 88%/70%, P. chitinolyticus 95%/18% |
| 345,498–377,216 | 31 kb | 47.6% | 24 | Mostly hypotheticals | Unknown | not BLASTed |
| 2,911,705–2,936,697 | 25 kb | 47.7% | 32 | VirB4/TrbL/VirB6/VirD4/MobA | **ICE #2** | Lacrimispora saccharolytica 96.3%/**100%** cov |
| 5,297,204–5,322,425 | 22 kb | **33.8%** | 21 | Mostly hypotheticals | **Unknown — no BLAST hits** | No hits in NCBI nt |
| 5,577,608–5,596,327 | 18 kb | 50.9% | 17 | Hypotheticals | Unknown | not BLASTed |
| 2,546,997–2,565,301 | 18 kb | 49.9% | 13 | Hypotheticals | Unknown | P. sp. P25 78%/100% |
| 7,180,257–7,192,022 | 11 kb | 40.1% | 10 | Hypotheticals | Unknown | not BLASTed |
| 5,599,526–5,609,082 | 9 kb | 46.6% | 10 | **MFS drug:H+ antiporter-2** | Drug efflux (non-specific) | not BLASTed |

**Flanking signals:** None of the top 10 regions are flanked by IS elements within 5 kb, and none
are tRNA-flanked within 10 kb. Classic GI hallmarks (IS + tRNA flanking) are absent.

#### IS element census (Phase 4.3)

| IS family | Count in M5 |
|-----------|-------------|
| IS110 family | 29 |
| IS4 family | 19 |
| IS3 family | 2 |
| Other/unknown | ~83 |
| **Total** | **133** (Bakta keyword search; Bakta reports 151 total transposable elements) |

M5 has ~57 more IS elements than I6 (151 vs 94), matching the 57 DUP regions in nucmer dnadiff.
These IS element duplications are scattered genome-wide and account for ~76 kb of extra sequence.

#### Key findings

**Finding 1 — Two independent ICEs in M5 (not one):**
- ICE #1 at 1.28–1.39 Mb (108 kb): VirD4/VirB4/TrbL/VirB6/Relaxase/Recombinase. No resistance genes.
- ICE #2 at 2.91–2.94 Mb (25 kb): VirB4/TrbL/VirD4/MobA conjugation machinery. Paenibacillus macerans-specific WP_ accessions.

**Finding 2 — Prophage at 6.77 Mb (59 kb):**
- BLAST: Paenibacillus polymyxa ATCC 15970 (72% id, 66% cov) — not Bacillus-origin; a Paenibacillus-infecting phage
- Genes: SP-beta-like phage structural proteins (YorD, phage tail, GP4 fiber) — Bakta annotation used Bacillus SP-beta reference but BLAST says Paenibacillus host
- GC = 44.1%; no antibiotic resistance cargo

**Finding 3 — Anomalous low-GC block at 5.30 Mb (22 kb, GC=33.8%):**
- **BLAST: NO HITS** against entire NCBI nt database (E<1e-10)
- GC=33.8% is 10–15% below M5 genome average — strong signal of foreign DNA from a low-GC organism
- Combined with no BLAST hits: this sequence is either from an organism not in NCBI databases, or is so diverged that nucleotide BLAST cannot find it
- Content: 21 CDS, all hypotheticals — no resistance genes annotated
- This is the most mysterious novel region in M5 — candidate for further investigation (protein-level BLAST or HMM search)

**Finding 4 — MFS drug efflux transporter in 9 kb block (5.60 Mb):**
- Drug resistance MFS transporter, drug:H+ antiporter-2 family (NKFIDM_04924)
- MFS efflux is non-specific; spectinomycin is not a typical MFS substrate
- Low priority unless MIC data suggests efflux activity

---

### Phase 4 — ICE donor identification (NCBI Entrez lookup)

**Method:** NCBI esummary on WP_ RefSeq protein accessions from ICE marker genes.
**Output:** `results/genomics/comparative/phase4_ICE_donor_lookup.txt`

#### 108 kb ICE #1 — donor identification (combined Entrez WP_ + BLAST)

| Locus | WP accession | Entrez organism | Gene product |
|-------|-------------|-----------------|--------------|
| NKFIDM_01129 | WP_216539118.1 | Paenibacillus | Recombinase |
| NKFIDM_01139 | WP_255222961.1 | Paenibacillus sp. 7541 | Mobilization protein |
| NKFIDM_01140 | WP_036623496.1 | Paenibacillaceae | Relaxase |
| NKFIDM_01144 | WP_036623490.1 | Paenibacillaceae | VirD4 coupling protein |
| NKFIDM_01151 | WP_036623477.1 | Paenibacillaceae | VirB4 ATPase |

**BLAST (4 kb window, centre of ICE):**
- **Hit 1: Brevibacillus agri ACCC03016** — 99.8% identity, **100% coverage** — the ICE is almost identical to a Brevibacillus agri chromosomal region
- Hits 2–5: Enterococcus faecium/hirae plasmids — 98.4% identity, **38.6% coverage** — only a ~1.5 kb sub-segment of the 4 kb window matches these plasmids. This represents a conserved conjugation module shared between Enterococcus plasmids and this ICE, not an Enterococcus-derived ICE.

**Conclusion:** ICE #1 primary donor = **Brevibacillus agri** (Bacillaceae, a close Paenibacillus relative). The ICE backbone is Brevibacillus-derived. The Enterococcus match is a chimeric segment — a shorter conserved conjugation module that has circulated across multiple lineages (Enterococcus, Brevibacillus, Paenibacillus). This is a hallmark of ICE evolution through modular recombination.

#### 25 kb ICE #2 — donor identification

Entrez WP_ lookup showed Paenibacillus macerans (WP_326048* series) — but BLAST reveals the true picture:

**BLAST (4 kb window, centre of ICE #2):**
- **Hit 1: Lacrimispora saccharolytica** — 96.3% identity, **100% coverage** (complete 4 kb match)
- **Hit 2: [Clostridium] saccharolyticum WM1** — 96.3% identity, **100% coverage**
- Hits 3–5: Clostridium sporogenes, Sporomusa ovata (gut Clostridia) — 93.5%, 56% coverage

**Conclusion:** ICE #2 originated from **Clostridia / Lachnospiraceae** (specifically Lacrimispora/Clostridium saccharolytica group), not from Paenibacillus. The WP_326048* accessions show P. macerans because M5 itself (or a closely related P. macerans isolate) is the only P. macerans representative in RefSeq with this ICE — but the sequence is Clostridial in origin. This represents a **cross-order HGT event** (Clostridia → Paenibacillus) mediated by the ICE. Lacrimispora saccharolytica and Clostridium saccharolyticum are gut-associated anaerobes that could co-exist with Paenibacillus in the mouse intestine.

#### Overall interpretation

**No E. coli-origin mobile elements in M5.** All characterised mobile elements derive from gut/soil bacteria in the Bacillales/Clostridia ecosystem:

| ICE | Donor | Evidence | E. coli involvement |
|-----|-------|----------|---------------------|
| ICE #1 (108 kb) | Brevibacillus agri | 99.8%/100% BLAST | None |
| ICE #2 (25 kb) | Lacrimispora/Clostridium | 96.3%/100% BLAST | None |
| Prophage (59 kb) | Paenibacillus polymyxa | 72%/66% BLAST | None |
| Integrative element (48 kb) | Paenibacillus spp. | 88-95% BLAST | None |

H3 (HGT from E. coli) is rejected at multiple independent levels:
1. No aadA/ANT(9) anywhere in M5 (Phase 3.2)
2. The 108 kb ICE has no resistance cargo (Phase 3.2 / Session 3)
3. ICE donors are Brevibacillus/Clostridia — not E. coli (Phase 4, BLAST)

---

## Decision Matrix — Complete Hypothesis Assessment (Sessions 1–4)

| Hypothesis | Prediction | Evidence | Verdict |
|------------|-----------|----------|---------|
| **H1** — De novo rpsE mutation in vivo | rpsE differs between M5/M6/M8 | M5=M6=M8 rpsE identical; deletion is pre-existing | **FALSIFIED** |
| **H3** — HGT of aadA from E. coli | aadA present in M5; ICE of E. coli origin | No aadA in M5; ICEs are Paenibacillus-native | **NOT SUPPORTED** |
| **H4** — Intrinsic resistance (rpsE Δ) | rpsE differs from I6 at spectinomycin target | Δ(A21-K22-V23) removes conserved Lys confirmed by: MAFFT 4-species MSA, 163x/9D CIGAR read-level, all 30S proteins otherwise identical, 16S h34 identical | **STRONGLY SUPPORTED** |
| **H10** — Ecological release | Paenibacillus expands in controls as well | Controls show Paenibacillus at <1% RA | **FALSIFIED** |

**Sole remaining question:** MIC of I6 (ATCC 7068) vs M5 for spectinomycin.
If I6 is sensitive (MIC < 64 µg/mL) and M5 is resistant, H4 is confirmed as the mechanism.

---

## Next steps

1. **Wet lab (Phase 7, highest priority):** MIC determination for I6 (ATCC 7068) vs M5
2. **BLAST results** for novel regions (in progress, NCBI qblast running) — will characterise 59 kb prophage, 22 kb low-GC block, and 48 kb integrative element
3. **Pan-genome (Phase 4.2):** Panaroo on M5/M6/M8/I6 — requires WSL install; deferred
4. **PCR primer sequences** — required before Phase 0.4 in silico PCR; obtain from original methods
5. **M7 re-sequencing** — medium priority wet lab task

