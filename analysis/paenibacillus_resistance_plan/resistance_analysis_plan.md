# Paenibacillus macerans Spectinomycin Resistance — Exhaustive Analysis Plan

**Author:** Melis Gencel
**Date:** 2026-03-22
**Status:** PLANNING ONLY — no code written

---

## 0. What We Have: Data Inventory

### Sequencing data (4 isolates)

| Sample | Mouse | FASTQ | Coverage (raw) | Coverage (assembled) | Contigs | Genome size | Genes | rRNAs |
|--------|-------|-------|---------------|----------------------|---------|-------------|-------|-------|
| TPHP3P_1_M5 | m5 | 1.17 Gb, 191,702 reads | 156x | 103x | **1** | 7.47 Mb | 6,695 | 25 |
| TPHP3P_2_M6 | m6 | 1.03 Gb, 145,608 reads | 138x | 101x | **1** | 7.47 Mb | 6,695 | 25 |
| TPHP3P_3_M7 | m7 | 69.9 Mb, 8,693 reads | **9x** | **9x** | **7** | 7.46 Mb | 6,772 | — |
| TPHP3P_4_M8 | m8 | 637 Mb, 87,841 reads | 85x | 84x | **1** | 7.47 Mb | 6,695 | 25 |

Assembly pipeline: **PlasmiSaurus** (Snakemake), Bakta v1.11.3 (database v6.0 full) for annotation.

### Critical observation about M7
M7 has only 9x coverage and assembled into 7 completely disconnected contigs (100% dead ends, N50 ~1.8 Mb). The mash species match shows each contig hitting a completely different Paenibacillaceae species (contig_1 → P. azoreducens; contig_2 → P. stellifer; contig_3 → Brevibacillus composti; contig_4 → P. yonginensis; contig_5 → P. cellulositrophicus; contig_6 → Clostridium scatologenes; contig_7 → nomatch/Thermobacillus). This is highly suspicious: either the isolate is a mixed culture, or the low-coverage assembly produced chimeric contigs, or the mash reference database is inadequate. **M7 must be treated as unreliable until validated.**

### Species identification discrepancy
- **Sourmash** (k-mer genus-level): All isolates → *Paenibacillus macerans*
- **Mash** (against chromosome database): M5/M6/M8 → *P. woosongensis B2_4* at only 0.79 identity, 38/5000 shared hashes — very low match
- **FastANI** (from prior analysis): 99% identity to *P. macerans* I6 (NCBI ATCC 778052)
- The sourmash/FastANI agreement vs low mash identity reflects a database gap: P. macerans may be poorly represented in the mash chromosome database. FastANI is the most authoritative whole-genome identity metric.

### CheckM assembly quality
All four isolates: 99.85% completeness, 0.06–0.12% contamination, strain heterogeneity = 0. This indicates pure, high-quality assemblies (except M7 which must be re-evaluated for coverage and chimeras).

### Key genes found in M5 annotation (Bakta)
Directly relevant to spectinomycin resistance:
- **rpsE** (NKFIDM_05734, pos 6,418,100–6,418,588, –strand): 30S ribosomal protein S5 — primary spectinomycin target. Notably: no RefSeq/WP_* accession assigned by Bakta, only UniRef hits (`UniRef50_P21467`, `UniRef90_E3E700`), suggesting the protein may be divergent from reference databases.
- **[uS5]-alanine N-acetyltransferase** (NKFIDM_05151, pos 5,830,084–5,830,644, –strand): acetyltransferase that modifies the N-terminal alanine of ribosomal protein S5/uS5. Acetylation of S5 is known to affect ribosome structure; this enzyme is directly adjacent functionally to the spectinomycin resistance locus.
- **rsmI** (NKFIDM_00044, pos 51,713–52,606, +strand): 16S rRNA (cytidine-1402, 2'-O)-methyltransferase. C1402 is in helix 34 of 16S rRNA — the primary spectinomycin binding site. The presence and activity of this enzyme is directly relevant to resistance.
- **25 rRNA genes**: 5 complete ribosomal RNA operons (5 × 16S + 5 × 23S + 5 × 5S). Any 16S mutation must be evaluated across all copies.
- **rsmA** (NKFIDM_00050): 16S rRNA adenine dimethylase (A1518/A1519) — normal ribosome maturation.
- **151 annotated transposable elements** across the M5 chromosome — this genome is highly plastic.
- **11 CRISPR arrays** — likely explaining the IS/transposon abundance (defense vs. attack arms race).
- **marR, merR family regulators**, multiple ABC transporters, MFS transporters — putative drug efflux systems.

### Genome size anomaly
P. macerans ATCC 7068 (reference I6 genome, GCF_000172175.2) is ~5.7 Mb. Our isolates are ~7.47 Mb — approximately **1.77 Mb larger**. This is substantial and must be explained. Possible causes: large plasmid integrated into chromosome, HGT from other gut bacteria, or our isolate is a different strain than labeled.

---

## 1. Phase 1: QC and Assembly Validation

**Purpose:** Every downstream interpretation depends on having clean, accurate assemblies. Nothing proceeds to interpretation until these pass.

### 1.1 Raw read QC
- **Tool:** NanoPlot or NanoStat on each FASTQ
- **Metrics to assess:** read length distribution, quality score distribution (Q-score), N50, longest read, total yield
- **Pass criteria:** median Q ≥ 10 (Q15 preferred for Guppy HAC or Dorado), rN50 > 10 kb (good for assembly), no unexpected bimodal length distributions suggesting barcoding errors
- **Nanopore advantage:** long reads allow direct assessment of read length; long reads (rN50 = 9.7–39.2 kb across our samples) will directly span repeat regions including rRNA operons
- **Nanopore disadvantage:** systematic homopolymer errors (especially in A/T runs); errors in rpsE or 16S rRNA coding regions must be verified at high coverage

### 1.2 Assembly graph inspection
- **Tool:** Bandage on the assembly graph (GFA format, if available from Flye/Unicycler/Raven within PlasmiSaurus)
- **For M5/M6/M8:** One circular contig with self-edge (0 bp overlap, 0 dead ends as reported) — confirms circular chromosome. Inspect whether the single contig is correctly circularized at dnaA (chromosomal replication origin marker). Check that dnaA is at the start.
- **For M7:** 7 disconnected linear fragments (100% dead ends). This assembly is fragmented at 9x coverage — below the minimum recommended coverage for Flye (typically 20–30x for reliable assembly). The mash hits to 7 completely different species strongly suggest chimeric assembly or contamination. **M7 should not be used for mutation analysis without re-sequencing.**
- **Positive result:** M5/M6/M8 assemble as single circular chromosomes, graph confirms circularity
- **Negative result:** Assembly graph shows multiple components suggesting plasmids or contaminating chromosomes

### 1.3 Assembly coverage validation
- **Tool:** minimap2 (map original reads back to assembly) → samtools depth → custom coverage plots
- **What to check:**
  - Uniform coverage across the chromosome for M5/M6/M8 (expected ~100x everywhere)
  - Look for regions of elevated coverage (2x median = potential tandem duplication) or dropped coverage (misassembly)
  - Specifically check the rpsE locus, all 16S rRNA operons, and any annotated transposase regions
  - rRNA operons will show collapsed coverage because all 5 copies map to the same reference location — this is expected but needs to be noted
- **Nanopore advantage:** Long reads can uniquely span between and through rRNA operons and resolve which of the 5 rRNA copies carries which variant
- **Nanopore disadvantage:** At rRNA repeats, assembly may have collapsed all 5 copies into one — check read depth at rRNA regions for expected 5x multiplier over chromosomal coverage

### 1.4 Contamination check
- **CheckM results:** Already run — 99.85% completeness, 0.06–0.12% contamination for all isolates. This is excellent.
- **Additional check:** Run Kraken2 or Centrifuge on the raw reads to assess what fraction map to non-Paenibacillus organisms. This is important given the concern about E. coli DNA in cultures (relevant to PCR false positive hypothesis).
- **Tool:** Kraken2 with PlusPF database (includes E. coli)
- **Positive result (clean):** >99% reads classified as Paenibacillus/Firmicutes; <0.1% E. coli
- **Concerning result:** >0.1% reads classified as E. coli → raises probability of PCR contamination being the source of the "positive" spectinomycin resistance PCR result

### 1.5 Re-annotation and annotation validation
- The Bakta annotation is already available and is high quality (v1.11.3, full database). However:
  - Verify that rpsE is in the annotation and correctly identified (confirmed: NKFIDM_05734)
  - Verify that the [uS5]-alanine acetyltransferase (NKFIDM_05151) is correctly identified
  - Check whether any hypothetical proteins flank rpsE that could be resistance cassettes
  - Count and catalogue all 25 rRNA genes — verify they appear as 5 complete operons (5 × 16S + 23S + 5S)
  - Count all pseudogenes (23 reported in M5) — determine if any are in resistance-related loci
- **Pass criterion:** rpsE, rsmI, all rRNA operons correctly annotated with confident coordinates

---

## 2. Phase 2: Reference Genome Validation

**Purpose:** Every comparative analysis depends on the right reference. P. macerans I6 (NCBI ATCC 778052, accession GCF_000172175.2) is claimed to be 99% ANI identical to our isolate, but the PlasmiSaurus mash analysis only found P. woosongensis at 79% identity. This discrepancy must be resolved first.

### 2.1 Confirm the reference genome
- **Tool:** FastANI, comparing M5/M6/M8 assemblies to:
  1. P. macerans I6 (GCF_000172175.2 — the claimed reference)
  2. P. woosongensis B2_4 (GCF_030122845.1 — the mash hit)
  3. All available P. macerans genomes in NCBI (as of 2025, there are ~20)
  4. The broader Paenibacillaceae family (download all complete genomes)
- **Expected result:** If FastANI truly shows 99% to I6, the mash discrepancy is a database issue (I6 not in the mash chromosome database). Both tools agree on species identity.
- **Concerning result:** If FastANI shows <95% to I6, a different reference may be more appropriate
- **Why it matters:** Our SNP calls, structural variant calls, and insertion/deletion calls are all relative to the reference. Using the wrong reference inflates the number of apparent "mutations."

### 2.2 Check genome size discrepancy
- P. macerans I6 is ~5.7 Mb; our isolates are ~7.47 Mb — a 1.77 Mb difference.
- **Tool:** Progressive Mauve or nucmer (MUMmer) whole-genome alignment of M5 vs I6
- **Questions to answer:**
  1. Where are the extra 1.77 Mb? Are they one large region or many small insertions?
  2. Are the extra regions plasmid-like (circular signature, origin of transfer)?
  3. Do they contain recognizable mobile genetic element cargo?
  4. Is the extra sequence enriched for known HGT features (different GC content, codon usage bias, flanking repeats)?
- **Nanopore advantage:** Long reads can span the junctions between ancestral chromosome and inserted sequence; junction sequences are key evidence for HGT mechanism
- **Positive HGT signal:** Extra regions flanked by direct repeats, att sites, or IS elements; GC content deviating ±5% from genome average
- **Negative/ambiguous:** Extra regions are seamlessly integrated with no flanking mobile element signatures

### 2.3 Assess I6 reference quality
- Download GCF_000172175.2 and inspect: is it a complete closed chromosome or a draft with many contigs?
- If I6 is a draft assembly (many contigs), SNP calling will be compromised by reference gaps
- **Action:** If I6 is a draft, use M5 or M6 assembly as the reference instead (since M5/M6 are single complete circular chromosomes at high coverage, they make better references for within-isolate comparisons)

---

## 3. Phase 3: Resistance Mechanism — All Hypotheses in Parallel

The following hypotheses are **not mutually exclusive** and should all be tested simultaneously.

---

### H1: De novo point mutation in rpsE (protein S5)

**Mechanism:** Spectinomycin binds the 30S ribosomal subunit and blocks translocation by stiffening the ratchet movement between 30S head and body. The ribosomal protein S5 (encoded by rpsE) directly contacts the spectinomycin binding site on 16S rRNA. Mutations in S5 — particularly in the region that contacts helix 34 (the spectinomycin binding site) — can sterically or electrostatically prevent drug binding. Known examples: K61N, Q16K in E. coli S5; G64D, P90L in other organisms.

**Evidence for (prior analysis):** Mutations reported in rpsE.
**Evidence against:** AMRFinderPlus/RGI/ResFinder found nothing — but these tools focus on known resistance genes, not point mutations in core chromosomal genes.

**Tests:**
1. **Pairwise alignment of rpsE amino acid sequence (M5/M6/M8) vs I6:** Extract rpsE protein sequences from all assemblies and align to I6 rpsE. Identify all amino acid substitutions.
   - Tool: muscle or MAFFT for alignment; any alignment viewer
   - What to look for: substitutions in the S5 regions that contact 16S rRNA helix 34, particularly residues that are conserved in spectinomycin-sensitive organisms but altered here
   - **Positive:** One or more substitutions in the known drug-contact residues (compare against published cryo-EM structures of spectinomycin-bound 30S, e.g., PDB 4WFA, 1FJG)
   - **Negative:** rpsE sequence identical to I6, or mutations in positions not contacting the drug

2. **Structural mapping of rpsE mutations:** Map any amino acid substitutions onto the 3D structure of S5 (PDB: 1FJG for T. thermophilus 30S, or use AlphaFold2 prediction for P. macerans S5)
   - Assess whether the substituted residue is within the spectinomycin-S5 interface
   - This requires no additional sequencing — purely computational

3. **Validate with read-level pileup at rpsE:** Even with clean assembly, Nanopore errors could create false-positive or false-negative variant calls at rpsE.
   - Tool: minimap2 + samtools mpileup or medaka variant calling
   - With 100x coverage, genuine SNPs are reliable; but indels in homopolymer runs within rpsE require manual inspection
   - **Nanopore disadvantage:** If the resistance mutation is in a homopolymer run within rpsE (e.g., a poly-A or poly-T stretch), Nanopore error rates of ~5-10% at homopolymers mean the mutation may be miscalled. Illumina short-read polishing (Medaka + Pilon) would resolve this.

4. **Check whether rpsE is present in one or multiple copies:** Most bacteria have a single rpsE, but given the enlarged genome, there may be a duplication. Count occurrences in the annotation.

5. **Compare rpsE across all three good isolates (M5, M6, M8):** If the same mutation is in all three mice, it is likely ancestral (in I6 or in the in vivo ancestor). If different mice have different mutations, they evolved independently in the mouse gut.
   - **Positive (shared):** Same mutation in M5/M6/M8 → strong evidence mutation was present before or at the start of the experiment; compare to I6
   - **Positive (different):** Different mutations per mouse → convergent evolution, very strong evidence for positive selection by spectinomycin
   - **Negative:** No difference from I6 at rpsE → rpsE is not the explanation

**Limitations:**
- Single isolate per mouse (M5, M6, M8) — we don't know if this represents the dominant strain or a variant
- No I6 at the exact time of colonization — the reference may itself have variants relative to the stock culture
- Nanopore homopolymer errors could mask or create mutations in rpsE

---

### H2: De novo mutation in 16S rRNA helix 34

**Mechanism:** Spectinomycin intercalates into helix 34 of 16S rRNA (nucleotides C1063, A1064, A1192, A1193 in E. coli 16S numbering). Mutations at these positions reduce binding affinity. This is the same mechanism as reported in the prior analysis of this isolate.

**Special complication:** P. macerans has 5 ribosomal operons (evidenced by 25 rRNA genes = 5 × 16S + 5 × 23S + 5 × 5S). A mutation in only ONE copy produces heteroplasmy — the remaining 4 copies still bind spectinomycin. For functional resistance, the mutation likely needs to be in all or most copies (or the mutant copy confers a dominant-negative effect — unlikely).

**Evidence for:** Mutations reported in 16S rRNA helix 34 in prior analysis.

**Tests:**
1. **Extract all 16S rRNA sequences from each assembly:** The 25 annotated rRNAs include 5 × 16S. In a collapsed assembly (all copies mapped to same contig region), variants may be visible in the pileup.
   - Tool: barrnap or RNAMMER to re-identify all rRNA positions; or use Bakta GFF3 coordinates
   - Extract the 5 16S sequences and align them to each other and to I6 16S

2. **Haplotype phasing of 16S copies using long reads:**
   - This is the key Nanopore advantage. Individual reads are long enough to span an entire 16S rRNA gene (~1,550 bp) and into flanking sequence. By clustering reads that span a 16S gene by their flanking context, we can determine whether the resistance mutation is in all 5 copies, some copies, or just one.
   - Tool: Extract reads spanning 16S regions (minimap2 → filter by alignment position), perform read clustering by flanking sequence
   - **Positive (all copies):** All 5 16S copies carry the same helix 34 mutation → full resistance expected
   - **Positive (partial):** Only 1-2 copies carry the mutation → intermediate resistance; implies mutation is recent and has not yet fixed
   - **Negative:** No 16S copies differ from I6 → 16S mutation is not the mechanism
   - **Nanopore disadvantage:** Systematic errors in homopolymer runs within 16S rRNA (particularly A/U-rich regions) may produce false variant calls. The Nanopore error rate in rRNA is typically higher than for protein-coding genes because rRNA has complex secondary structure that stalls the motor enzyme.

3. **Verify helix 34 variant with read pileup:** At the exact spectinomycin-binding positions in each of the 5 annotated 16S genes, assess the per-read base call and quality score
   - Look for consistent variant at ≥80% of reads (genuine mutation) vs scattered noise (<20% = likely sequencing error)
   - If variant is at ~20% frequency: this is the signature of 1/5 copies being mutated (1 mutant out of 5 operons = 20% reads from that locus)

4. **Compare 16S helix 34 variants across M5/M6/M8:** Same logic as rpsE — shared mutations suggest ancestral origin; independent mutations suggest convergent evolution.

**Limitations:**
- 16S rRNA is transcribed from 5 operons — Nanopore assembly typically collapses these into a consensus. The "assembly" sequence may represent a majority vote across copies, masking a minority variant.
- The helix 34 mutation may be in the I6 reference itself (natural Paenibacillus divergence from E. coli 16S numbering) — we cannot interpret this without the I6 16S sequence at these exact positions.
- Homopolymer errors are most likely in rRNA genes.

---

### H3: Horizontal acquisition of E. coli spectinomycin resistance gene (aadA or speR)

**Mechanism:** E. coli in the gut carries a chromosomally integrated, spectinomycin-resistance spectinomycin adenylyltransferase (aadA/speR) as part of the original inoculum (the colonizing strain was "spectinomycin-resistant E. coli K12"). HGT from E. coli to Paenibacillus could transfer this gene. This hypothesis was the original one (hence the PCR result), but AMR databases found nothing.

**Prior evidence for:** PCR with E. coli-specific spectinomycin resistance primers gave a positive band.
**Prior evidence against:** AMRFinderPlus, RGI, ResFinder found no known spectinomycin resistance genes in the assembly.

**This is a critical contradiction that must be resolved.**

**Tests:**

1. **Map the original PCR amplicon to the assemblies:**
   - If you have the primer sequences, design an in-silico PCR using the M5/M6/M8 assemblies as templates (tool: isPCR, primersearch, or BLAST of primer sequences)
   - **Positive:** Primers amplify a region in the Paenibacillus assembly → a gene is present that the databases missed
   - **Negative:** Primers do not amplify in Paenibacillus assembly → the PCR band was from E. coli contamination

2. **BLAST the spectinomycin resistance gene sequence against M5 assembly:**
   - Download the aadA gene from the E. coli K12 inoculum strain's genome sequence (if available)
   - BLAST against M5 assembly at high sensitivity (blastn with word_size 7, or tblastn with the protein)
   - **Positive:** Hit with >80% identity over >100 bp → gene present, possibly divergent enough to escape AMR databases
   - **Negative:** No significant hit → gene truly absent

3. **Re-run AMR databases at lower thresholds:**
   - RGI: run in "nudge" or "loose" mode to detect partial matches
   - AMRFinderPlus: --plus flag to include stress/virulence genes
   - ResFinder with lower identity threshold (e.g., 70% instead of 90%)
   - abricate with CARD, NCBI, ARG-ANNOT, RESFINDER databases
   - **Positive:** Hit at 70-80% identity → distant homolog present

4. **Search for novel adenylyltransferase-family proteins in M5:**
   - Spectinomycin adenylyltransferases are in the ANT(9) family (aminoglycoside nucleotidyltransferases)
   - BLAST the M5 protein set (NKFIDM_*.faa) against a curated ANT database or CARD
   - Tool: diamond blastp against CARD protein sequences
   - **Positive:** A hypothetical protein in M5 has structural homology to ANT(9) family

5. **Map raw reads to E. coli aadA gene:**
   - This is critical: raw reads are more sensitive than assembly for detecting minority species or integrated genes that may have been fragmented or misassembled
   - Tool: minimap2, map FASTQ reads to E. coli K12 aadA/speR sequence
   - **Positive:** Any reads map to aadA with >80% identity → gene sequence is present in the raw data, even if not in the final assembly
   - **Negative:** Zero reads map to aadA → gene is truly absent in the sequenced isolate

6. **Assess E. coli contamination of culture:**
   - Map raw reads to the complete E. coli K12 inoculum genome
   - Assess what fraction of reads are E. coli vs Paenibacillus
   - **Positive (concerning):** >0.1% reads map to E. coli → culture likely contaminated; PCR result is probably E. coli DNA in the template
   - **Positive (supporting HGT):** Very few or no reads map globally to E. coli, but reads do map specifically to the aadA gene region → the aadA gene has been acquired by Paenibacillus but the rest of the E. coli genome is absent

**Nanopore advantage:** Long reads can show whether an E. coli-origin gene is flanked by Paenibacillus sequence (proving integration) vs floating in a metagenomics-style read cloud (suggesting contamination).
**Nanopore disadvantage:** Low coverage contaminating sequences may not assemble; they're detectable only in read-level mapping.

**Limitation:** We do not have the sequence of the exact E. coli K12 inoculum strain. If the aadA gene is in a non-standard location or has been modified, BLAST searches may miss it.

---

### H4: Native intrinsic resistance — rpsE/16S naturally diverged from E. coli homologs

**Mechanism:** Spectinomycin was developed based on E. coli ribosome biochemistry. Paenibacillus macerans is phylogenetically distant from E. coli (Firmicutes vs Proteobacteria). The natural P. macerans rpsE and 16S helix 34 may already have sequence variants at the drug-binding interface that confer partial or full resistance, independent of any mutation during the experiment.

Under this hypothesis, Paenibacillus was always resistant to spectinomycin at the concentration used in drinking water. Its expansion under spectinomycin selection was not because it acquired resistance, but because it was released from competition (other bacteria suppressed, leaving an ecological niche).

**This is the null hypothesis and must be tested before any mechanism can be claimed.**

**Tests:**

1. **Compare rpsE protein sequence to spectinomycin-sensitive organisms:**
   - Align P. macerans rpsE (M5) to rpsE from E. coli K12, S. pneumoniae, H. influenzae (organisms where spectinomycin MIC data is known)
   - Specifically at the residues that contact the drug (from cryo-EM structures)
   - If P. macerans rpsE naturally differs at these positions in the same way as known resistance mutations → intrinsic resistance is likely

2. **Compare M5/M6/M8 rpsE to I6 rpsE:** If the mutations found are also in I6, they are ancestral and not acquired during the experiment. This is the single most important comparison.
   - **Positive (intrinsic):** M5 rpsE = I6 rpsE (no difference → P. macerans is naturally resistant and our isolate did not change)
   - **Positive (acquired):** M5 rpsE ≠ I6 rpsE at known resistance positions → mutation occurred after colonization

3. **MIC testing (wet lab — Phase 7, but flagged here):** Actually measure the spectinomycin MIC for P. macerans I6 (ATCC culture) vs our isolates. If I6 is already resistant, the expansion is ecological not resistance-driven.

4. **Ecological release check:** Examine the 16S microbiome data (already generated in module 01_16S) to assess whether Paenibacillus expanded in the CONTROL mice (antibiotic-only, no E. coli inoculation) as well as in colonized mice. If Paenibacillus expanded only in colonized mice (not in antibiotic-only controls), then the E. coli colonization may have initially suppressed Paenibacillus, and its expansion later in colonized mice reflects E. coli dying off or ecological release, not Paenibacillus acquiring resistance.
   - **This analysis can be done RIGHT NOW with existing 16S data.**

**Limitation:** We cannot measure the MIC of the in-vivo ancestral Paenibacillus population because we have no early-timepoint samples. We can only compare the late-timepoint isolate to the I6 reference.

---

### H5: 16S rRNA methylation at helix 34 (enzymatic resistance)

**Mechanism:** rsmI (NKFIDM_00044) catalyzes 2'-O-methylation of cytidine-1402 of 16S rRNA. C1402 is in helix 34, adjacent to the spectinomycin binding site. If this enzyme is hyperactive, overexpressed, or if a novel enzyme methylates positions directly contacted by spectinomycin (C1190, A1064 in some numbering schemes), the drug cannot bind even without mutations in the rRNA sequence.

This is distinct from H2 (sequence mutation in 16S) — the rRNA nucleotides are unmutated, but their chemical modification blocks drug access.

**Tests:**

1. **Verify rsmI is present in all isolates and I6:**
   - BLAST rsmI (NKFIDM_00044) protein against I6 genome
   - **If rsmI is absent from I6 but present in our isolates** → acquired enzyme, potential HGT of a methyltransferase rather than adenylyltransferase resistance
   - **If rsmI is present in I6** → native gene, not acquired; however, its promoter or regulatory region may differ

2. **Check rsmI promoter and regulatory context vs I6:**
   - Compare the 300 bp upstream of rsmI in M5 vs I6
   - Any promoter mutation or regulatory change could increase rsmI expression → more methylated 16S → reduced spectinomycin binding
   - Tool: align promoter regions by pairwise alignment

3. **Search for other 16S rRNA methyltransferases:**
   - RlmAII (methylates G1405, known aminoglycoside resistance in mycobacteria)
   - Clf (methylates A1408)
   - Search the M5 proteome for SAM-dependent methyltransferase domain (PFAM PF05175, PF08241) with a BLAST approach
   - Any novel methyltransferase near an IS element or with unusual phylogenetic distribution is a candidate

4. **Compare rsmI gene copy number vs I6:** In the M5 annotation, verify whether rsmI is single-copy or duplicated. A tandem duplication would increase expression without promoter mutation.

**Limitations:** Without transcriptomic data (RNA-seq), we cannot directly measure rsmI expression. The best we can do is identify sequence differences in the coding or regulatory region vs I6 and infer expression changes.

---

### H6: Gene duplication of rpsE or other resistance locus

**Mechanism:** Tandem gene duplication of rpsE could allow one copy to accumulate resistance mutations (drug-binding-site alterations) while the other copy maintains normal ribosome function. This is a recognized mechanism for evolving resistance to targets where the wild-type function is essential.

**Tests:**

1. **Count rpsE copies in each assembly:** Search Bakta annotation and GFF3 for all annotations matching rpsE/S5/ribosomal protein S5.
2. **Check read coverage at the rpsE locus:** Elevated coverage (e.g., 2x median) at rpsE in M5 assembly would indicate a tandem duplication that the assembler collapsed.
3. **Check whether the [uS5]-alanine N-acetyltransferase (NKFIDM_05151) is adjacent to rpsE:** In many bacteria, the S5 acetyltransferase is encoded in the same operon as rpsE (e.g., in B. subtilis, rimJ/yaaA is near rpsE). If the acetyltransferase has been duplicated alongside rpsE, this is a candidate functional duplication.
   - In the M5 annotation, rpsE is at position 6,418,100 and the acetyltransferase is at 5,830,084 — these are 588 kb apart, so they are NOT in the same operon. This reduces the tandem duplication hypothesis for these two genes.
4. **Structural variant analysis:** Use tools like Sniffles2 or SVIM on the read-to-assembly alignment to detect any tandem duplications at or near rpsE.

**Limitations:** Assembly may collapse tandem duplicates. Only raw read evidence (coverage + structural variant tools) can resolve this.

---

### H7: Efflux pump upregulation

**Mechanism:** Multiple antibiotic resistance (mar/mex) efflux pumps can export spectinomycin from the cell, lowering the intracellular drug concentration below the inhibitory level. Upregulation can occur via mutation in the pump promoter or in a repressor gene (marR family). P. macerans M5 has annotated marR (NKFIDM_03886) and merR (NKFIDM_03252) family regulators.

**Tests:**

1. **Catalogue all efflux pumps in M5:**
   - Search M5 proteome for RND (resistance-nodulation-division), ABC, and MFS transporter families known to export antibiotics
   - Tool: BLAST against the CARD efflux pump category; or use TransportDB
   - Count and identify candidate transporters

2. **Compare efflux pump copy number vs I6:** Are there more efflux pumps in M5 than in I6? HGT of an additional efflux pump?

3. **Inspect marR promoter vs I6:** marR is a repressor of efflux pump genes. Mutations in marR that impair its repressor function lead to constitutive efflux upregulation. Compare the marR coding sequence and its binding sites in M5 vs I6.

4. **Assess whether any IS elements are inserted near efflux pump promoters:** IS element insertion into a repressor gene or upstream of a pump gene is a classic mechanism for resistance upregulation.
   - Nanopore advantage: Long reads can span an IS element and its flanking sequence in a single read, definitively placing the insertion.

**Limitations:** Efflux typically confers lower-level or broader-spectrum resistance. Spectinomycin is not a great efflux substrate. This hypothesis is considered less likely for high-level spectinomycin resistance but cannot be excluded.

---

### H8: Novel resistance gene — no database match

**Mechanism:** Paenibacillus macerans may carry a spectinomycin resistance gene that is not yet in any AMR database. This could be:
- A novel adenylyltransferase family enzyme that has diverged beyond the database similarity threshold
- A novel acetyltransferase that modifies the drug
- A ribosome protection protein not previously characterized

The AMR databases (CARD, ResFinder, ARG-ANNOT) are biased toward clinical isolates and known Gram-negatives. Novel resistance mechanisms in Firmicutes (especially non-pathogenic environmental organisms) are systematically underrepresented.

**Tests:**

1. **Focused search for ANT/AAD family proteins:** Search M5 proteome with HMM profiles for aminoglycoside nucleotidyltransferases (Pfam PF05347, PF13612) using HMMER
   - A hit with E-value < 1e-5 but <70% identity to any known resistance gene is the signature of a novel homolog

2. **Search for aminoglycoside acetyltransferase domain:** AAC(3), AAC(6') family (Pfam PF02110, PF13508) — these sometimes confer spectinomycin resistance

3. **Phylogenetic analysis of any ANT/AAC hits:** Build a phylogenetic tree with known spectinomycin resistance enzymes + any candidates from M5. If M5 candidate branches outside all named resistance genes but with a reasonable root position, it could be a novel resistance enzyme.

4. **Check for novel ribosome-binding proteins:** Proteins with ribosome-binding domains (S1 domain, KH domain, RRM domain) that are unique to M5 vs I6 and have no database hits could be ribosome protection proteins.

5. **Screen all hypothetical proteins in M5 (n=23 pseudogenes + many hypotheticals):** For each hypothetical protein, run PHYRE2 or AlphaFold2 structure prediction and assess whether the predicted structure resembles any known AMR enzyme class. This is computationally intensive but feasible for a single genome.

**Limitations:** Structural prediction for novel resistance mechanisms is speculative. Wet lab validation (expression in sensitive host, MIC testing) is ultimately required.

---

### H9: Conjugative or transductive HGT — mechanism of acquisition

**This is a separate question from WHAT resistance mechanism, but critically important for the HGT claim in the paper.**

Even if H3 (E. coli aadA) is supported, we need to determine HOW it was transferred.

**Mechanism options:**
- **Conjugation:** E. coli transfers a plasmid or chromosomal region via a conjugative pilus to Paenibacillus. Requires: (1) E. coli has a conjugative plasmid or ICE; (2) conjugation can occur across the Firmicutes/Proteobacteria barrier (unusual); (3) the recipient has recombination machinery to integrate the DNA.
- **Transformation:** Naked DNA released from lysed E. coli cells is taken up by naturally competent Paenibacillus. Paenibacillus species can be naturally competent.
- **Transduction:** A phage infects E. coli, packages the resistance gene, and transfers it to Paenibacillus. Requires phage with broad host range.

**Tests:**

1. **Search for oriT (origin of transfer) in M5:** oriT sequences mark conjugative elements. The Bakta annotation reported 0 oriT — but this uses a database of known oriTs. Search with a broader tool: CONJscan (Abby et al.) or ICEfinder.

2. **Search for ICE (integrative conjugative element) in M5:** ICEs are chromosomally integrated conjugative elements. Tools: ICEfinder, ICEPAN (web server)
   - If an ICE is found that carries a resistance gene, this is the HGT mechanism

3. **Search for prophages:** Map M5 sequence against PHROG (phage protein database) or use PHASTER. If a prophage is present, and it matches a phage known to infect E. coli (transducing phage), this supports transduction.
   - The M5 annotation shows 11 CRISPR arrays — CRISPR is a phage defense system. An organism with active CRISPRs has likely been under phage attack. Check if any CRISPR spacers match E. coli phages.

4. **Inspect mobile element context of any resistance gene:** If a resistance gene is found (by any of the above tests), assess its immediate genomic context:
   - Flanking direct repeats (= IS element footprint)
   - Flanking IS elements (= composite transposon)
   - Flanking attL/attR sites (= prophage or integron integration)
   - GC content deviation from genome average (= recently acquired sequence)
   - Codon usage analysis (= donor organism identity)

5. **Codon usage analysis of any candidate resistance gene:** If a novel gene is found, compare its codon usage (via EMBOSS cusp or CodonW) to the P. macerans genome average vs E. coli genome average. A gene with E. coli-like codon usage in a Paenibacillus genome is evidence of recent HGT from E. coli.
   - **Nanopore advantage:** Long reads capture the full mobile element context in single reads — the read can show the resistance gene AND its flanking IS elements AND the chromosomal integration point all at once.

---

### H10: The expansion was ecological, not resistance-mediated

**Mechanism:** This is the alternative to all resistance hypotheses. Spectinomycin in drinking water at the concentrations used may not achieve sufficient gut luminal concentrations to inhibit a spore-forming bacterium like Paenibacillus macerans that spends much of its time in a metabolically inactive spore state. The expansion could instead reflect:
- Ecological release: E. coli and its spectinomycin resistance gene are suppressing other bacteria via quorum sensing, bacteriocins, or competitive exclusion. As E. coli declines (or its population stabilizes), Paenibacillus fills the niche.
- Spore protection: Paenibacillus spores are completely resistant to spectinomycin (drug cannot inhibit translation in a non-translating spore). If the in vivo spectinomycin concentration is below the vegetative cell MIC but the bacterium can cycle through spore states, it may persist and grow intermittently.
- Nutrient niche: Paenibacillus can fix nitrogen (nifH annotated in M5!) and degrade complex polysaccharides — it may occupy a metabolic niche that E. coli does not.

**Tests:**

1. **16S data (already available):** Plot Paenibacillus relative abundance over time in control mice (c_m1–c_m4, antibiotic-only, no E. coli) vs colonized mice (m1–m8). If Paenibacillus expands only in colonized mice with E. coli, it suggests E. coli/spectinomycin interaction is driving the expansion.

2. **Timing analysis from 16S data:** Does Paenibacillus begin expanding before or after the peak of E. coli colonization? If Paenibacillus expands only as E. coli declines, this is ecological release not resistance.

3. **Barcode data (already available):** Within colonized mice, does Paenibacillus abundance correlate with specific E. coli barcode lineages? If a specific E. coli lineage suppresses Paenibacillus and its decline allows Paenibacillus expansion, this points to ecological mechanism.

---

## 4. Phase 4: Comparative Genomics

**Purpose:** Understand the full genomic context — what else differs between our isolates and I6, and what is shared/divergent among M5, M6, M8.

### 4.1 Pan-genome analysis across three good isolates
- **Tool:** Roary, PIRATE, or Panaroo on M5, M6, M8 (+ I6 if available)
- **Output:** Core genome (genes in all isolates), accessory genome (genes in only some), unique genes per isolate
- **What to look for:**
  - Genes unique to M5/M6/M8 but absent from I6 → acquired during the experiment or present in the mouse gut ancestor but not the reference
  - Genes present in some but not all of M5/M6/M8 → heterogeneity among isolates from different mice
  - Resistance genes in the accessory genome

### 4.2 Whole-genome SNP analysis vs I6
- **Tool:** minimap2 + Medaka (Nanopore-native) or Snippy (BWA-MEM based) for SNP calling
- **For M5/M6/M8 vs I6:** Identify all SNPs, indels, and structural variants
- **Of special interest:**
  - SNPs in rpsE (see H1)
  - SNPs in 16S rRNA (see H2)
  - SNPs in rsmI or its promoter (see H5)
  - SNPs in marR, marA, or efflux pump genes (see H7)
  - Large insertions absent from I6 (= HGT candidate)
  - IS element insertions (disrupting genes = potential resistance mechanism)
- **Nanopore advantage:** Structural variants (inversions, translocations, large insertions) are detectable in single reads; short reads miss most SVs
- **Tool for SVs:** Sniffles2 (for long reads), Syri (for genome-level rearrangements)

### 4.3 Synteny analysis vs I6 and vs P. macerans type strain
- **Tool:** Progressive Mauve or Sibelia for synteny blocks
- Map which regions of M5 have no synteny to I6 — these are insertions/novel regions
- **Expected:** Given the 1.77 Mb size difference, large novel regions should be identified

### 4.4 GC content and codon usage of novel regions
- For every region in M5 with no I6 synteny: plot GC content in 1 kb windows
- Regions with GC content deviating >5% from the genome average (expected ~45% for P. macerans) are likely horizontally acquired
- **Tool:** seqkit sliding window GC calculation; or bedtools with custom script
- **Codon usage:** CodonW or EMBOSS cusp on genes within novel regions

### 4.5 Mobile element census
- 151 transposase-containing genes were detected in M5. Catalogue:
  - IS families present (IS4 confirmed, others likely)
  - Insertion sites (which genes are disrupted or flanked)
  - Whether IS elements are in the same positions in M5/M6/M8 vs I6 (transposon activity during experiment)
- **Nanopore advantage:** Long reads uniquely resolve IS element insertions because short reads cannot span the repeat sequences at IS boundaries

---

## 5. Phase 5: HGT Signal — Genome-Wide

**Purpose:** Distinguish whether HGT is localized (one gene) or a broader phenomenon (multiple transfers, mobile island).

### 5.1 Alien Hunter or similar HGT prediction tools
- **Tools:** Alien Hunter, IslandViewer4 (GI prediction), SIGI-HMM, PAI-Finder
- Run on M5 assembly to identify all putative genomic islands (GIs) based on:
  - Unusual GC content
  - Unusual dinucleotide composition
  - Flanked by tRNA genes (classic GI integration sites)
  - Flanked by IS elements or direct repeats
- **Output:** List of GI candidates with coordinates
- **What to look for:** A GI that contains a resistance gene candidate (even if weak database match)

### 5.2 CRISPR spacer analysis
- M5 has 11 CRISPR arrays. CRISPR spacers record the history of phage/plasmid invasions.
- Extract all spacer sequences from the Bakta annotation (CRISPR arrays are annotated)
- BLAST spacers against:
  - E. coli phage databases (λ, M13, P1, etc.)
  - E. coli plasmid sequences
  - If spacers match E. coli plasmids → Paenibacillus has previously encountered E. coli-origin DNA and mounted an immune response → implies prior DNA transfer events (consistent with HGT history)
- **Tool:** CRISPRdetect + BLAST

### 5.3 Test for recent HGT between E. coli and Paenibacillus in gene trees
- For any gene found in M5 that is also in E. coli: build a phylogenetic tree of that gene including:
  - M5 sequence
  - Other Paenibacillus spp. sequences
  - E. coli sequences
- A gene that phylogenetically groups WITH E. coli rather than with other Paenibacillus spp. is evidence of HGT from E. coli to Paenibacillus
- **Tool:** MAFFT + IQ-TREE + iTOL for visualization

### 5.4 Assess whether the extra 1.77 Mb is a self-transferable element
- If the extra 1.77 Mb in M5 vs I6 contains: integrase + conjugation machinery + resistance genes → this is an integrative conjugative element (ICE) or a mega-transposon
- The Bakta annotation shows 0 oriT — but check the novel regions with ICEfinder
- If oriT is found in the novel 1.77 Mb region → this element could transfer to other bacteria in the gut

---

## 6. Phase 6: Decision Logic

Given all results from Phases 1–5, the following decision tree determines what conclusion is warranted:

### 6.1 Sufficient evidence for HGT claim
**Required conditions (all must be met):**
1. An E. coli-origin sequence is found in the M5/M6/M8 assembly (H3)
2. The sequence is absent from I6 (not ancestral)
3. The sequence phylogenetically groups with E. coli (not other Paenibacillus spp.)
4. The sequence has codon usage consistent with E. coli origin
5. The sequence is flanked by mobile element signatures (IS elements, att sites, direct repeats)
6. The PCR result is explained by this finding (in silico PCR with original primers amplifies this sequence in M5)
7. M5/M6/M8 all carry the same sequence (monophyletic acquisition, or independent parallel acquisitions)

**If all 7 met:** Strong evidence for HGT from E. coli to Paenibacillus as the resistance mechanism.

### 6.2 Sufficient evidence for de novo mutation
**Required conditions (all must be met):**
1. rpsE or 16S helix 34 sequence differs from I6 at known drug-contact positions
2. The mutation is in all three good isolates (M5/M6/M8) OR different mutations in each (convergent)
3. No E. coli-origin sequence is found in the assemblies
4. The mutation is absent from I6 reference
5. The mutation maps structurally to the drug-binding interface

**If all 5 met:** Evidence for convergent de novo mutation as the mechanism. If three mice have independent mutations → very strong evidence for selection in vivo.

### 6.3 Intrinsic resistance (null hypothesis not rejected)
**Conditions:**
1. rpsE sequence in M5 = rpsE sequence in I6 (no mutation)
2. 16S helix 34 in M5 = 16S in I6 (no mutation)
3. No E. coli-origin gene found
4. 16S data shows Paenibacillus expands in control mice too (H10)

**Conclusion:** Paenibacillus expansion may not be resistance-mediated. The isolates grew on spectinomycin plates due to naturally low spectinomycin sensitivity. The PCR positive is E. coli contamination.

### 6.4 Ambiguous outcome
Likely when:
- rpsE mutations exist but are also in I6 (ancestral)
- The 1.77 Mb extra genomic content is present but contains no identifiable resistance gene
- PCR in silico fails but raw reads map weakly to aadA

**In this case:** The study can conclude: "mechanism remains unknown, but we definitively exclude the following: [list all negatives]." Wet lab validation required.

### 6.5 Cross-isolate consistency check
- If M5, M6, M8 all have the same resistance mutation/gene: the mutation was acquired once (in the common gut ancestor or before colonization) → supports a clonal expansion hypothesis
- If M5, M6, M8 have different mutations at the same locus: convergent evolution under selection → strongest possible evidence that the mutation is under positive selection by spectinomycin
- If M5, M6, M8 have completely different genotypes: independent acquisitions → HGT is occurring reproducibly and the transfer rate is high

---

## 7. Phase 7: What Wet Lab Experiments Resolve Remaining Ambiguity

These are ordered by importance and feasibility:

### 7.1 MIC determination for P. macerans I6 ATCC stock
- **Purpose:** Tests H4 (intrinsic resistance) directly
- **Method:** Standard broth microdilution (EUCAST/CLSI) for spectinomycin against I6 ATCC culture
- **Expected outcome:** If I6 MIC > clinical breakpoint → intrinsic resistance; if I6 MIC < breakpoint → our isolates have acquired resistance
- **Priority: HIGH** — this is the single most important experiment to do

### 7.2 Sanger sequencing of rpsE and 16S helix 34 from multiple colonies
- **Purpose:** Validates Nanopore-called mutations; independent of assembly errors
- **Method:** PCR amplification of rpsE and 16S helix 34 from M5/M6/M8 colonies, Sanger sequencing, compare to I6
- **For 16S:** sequence multiple colonies and use genus-specific primers to amplify individual 16S copies rather than the full set
- **Priority: HIGH**

### 7.3 Southern blot or specific PCR for the E. coli aadA gene
- **Purpose:** Validates the original PCR result and resolves H3
- **Method:** Re-do the PCR on:
  1. Freshly grown M5/M6/M8 colonies
  2. A P. macerans I6 culture (negative control)
  3. The E. coli K12 inoculum (positive control)
  4. A DNA-free blank (contamination control)
- **Priority: HIGH** — the PCR result is the foundation of the entire HGT hypothesis

### 7.4 Re-sequencing M7 with more DNA
- **Purpose:** M7's fragmented assembly at 9x is unusable. Re-extract DNA and re-sequence.
- If M7 isolate is still available, re-sequence with R10.4.1 (higher accuracy Nanopore chemistry) and/or add Illumina short reads for polishing
- **Priority: MEDIUM**

### 7.5 Complementation / expression of M5 rpsE in E. coli spectinomycin-sensitive strain
- **Purpose:** Tests whether the M5 rpsE mutation is sufficient for resistance (H1)
- **Method:** Clone M5 rpsE under an inducible promoter in E. coli; test MIC
- If E. coli expressing M5 rpsE has elevated spectinomycin MIC → direct evidence that the M5 mutation confers resistance
- **Priority: MEDIUM**

### 7.6 Early-timepoint Paenibacillus isolation (if experiment can be repeated)
- **Purpose:** Directly tests whether resistance is ancestral or acquired
- **Method:** In a new experiment with fresh I6 inoculation and spectinomycin, take Paenibacillus isolates at day 1, day 5, day 12 — if early isolates are susceptible but late isolates resistant, resistance evolved in vivo
- **Priority: LOW** (requires new animal experiment)

### 7.7 Transfer experiment
- **Purpose:** Directly tests whether the resistance can be transferred from Paenibacillus back to E. coli
- **Method:** Mix cultured Paenibacillus isolates with spectinomycin-sensitive E. coli under conditions that favor conjugation; plate on spectinomycin to select for resistant E. coli transconjugants
- **Priority: LOW** (speculative — only relevant if a transferable element is found in Phase 5)

### 7.8 RNA-seq to measure rsmI and efflux pump expression
- **Purpose:** Tests H5 (methyltransferase overexpression) and H7 (efflux upregulation)
- **Method:** Grow M5/M6/M8 and I6 in identical conditions; extract RNA; RNA-seq; compare transcript levels of rsmI, efflux pump genes
- **Priority: LOW** (expensive; pursue only if rsmI promoter mutations found in Phase 3)

---

## 8. Critical Gaps and Caveats

### 8.1 The M7 problem
M7 assembled into 7 fragments at only 9x coverage. The mash hits suggest it may be a contaminated culture (multiple species assembled together). Until M7 is re-sequenced or validated by culture purity testing, it should be excluded from any resistance mechanism conclusions.

### 8.2 No early-timepoint samples
We only have day 12 isolates. We cannot track the evolutionary trajectory within the gut. All we can say is the end-state genomic sequence. If Paenibacillus was already resistant before colonization (H4), we cannot detect this without MIC data for the pre-experiment ancestral population.

### 8.3 Single isolate per mouse (M5/M6/M8)
Each mouse contributed one isolate. This does not capture population-level diversity. There may be co-existing resistant and susceptible Paenibacillus subpopulations in the gut that we miss. A single colony pick is a single lineage.

### 8.4 Nanopore-specific error issues for key loci

| Locus | Nanopore concern | Mitigation |
|-------|-----------------|------------|
| rpsE | Homopolymer errors may create or mask SNPs | Verify with Sanger sequencing; check base quality at suspect positions |
| 16S rRNA helix 34 | rRNA genes have complex secondary structure → stalling errors; all 5 copies map to same assembly region | Phase individual reads by flanking sequence; validate with Sanger on isolated copies |
| IS element boundaries | AT-rich direct repeats are hard to call | High coverage mitigates; use long reads spanning repeat entirely |
| aadA/novel gene | If AT-rich, may be misassembled | Map raw reads directly rather than relying on assembly |

### 8.5 The genome size anomaly must be explained
The 1.77 Mb excess vs I6 is too large to ignore. If this extra DNA contains the resistance mechanism (even partially), not explaining it would be a major gap in the analysis. The expanded genome could represent:
- A large integrated element (ICE, genomic island)
- Tandem repeat expansions
- Gene family expansions (the 151 transposons alone account for perhaps 150 × 1.5 kb = ~225 kb)
- A plasmid that circularized into the main assembly

### 8.6 The contradiction between PCR positive and database negative
This is the central unresolved tension. The PCR with E. coli-specific primers should not amplify Paenibacillus DNA — unless: (a) there is genuine HGT, (b) there is E. coli DNA contamination in the culture, or (c) the primers amplify a distant Paenibacillus homolog. This must be explicitly resolved. The analysis plan must answer: **was the PCR result real or an artifact?**

---

## Appendix A: File Structure for This Analysis

```
analysis/paenibacillus_resistance_plan/
└── resistance_analysis_plan.md   ← this document

data/sequence_data/
├── TPHP3P_fastq/
│   ├── TPHP3P_1_M5__+_.fastq.gz
│   ├── TPHP3P_2_M6__+_.fastq.gz
│   ├── TPHP3P_3_M7__+_.fastq.gz
│   └── TPHP3P_4_M8__+_.fastq.gz
└── TPHP3P_results/
    ├── TPHP3P_1_M5__+_/           ← 1 contig, 7.47 Mb, 103x ★ PRIMARY
    ├── TPHP3P_2_M6__+_/           ← 1 contig, 7.47 Mb, 101x ★ PRIMARY
    ├── TPHP3P_3_M7__+_/           ← 7 contigs, 7.46 Mb, 9x ⚠ UNRELIABLE
    └── TPHP3P_4_M8__+_/           ← 1 contig, 7.47 Mb, 84x ★ PRIMARY
```

Each sample directory contains:
- `*_summary.txt` — assembly QC overview
- `*_stats.tsv` — read statistics
- `*_contigs.txt` — assembly graph stats
- `*_mash-species.txt` — closest genome match (mash)
- `*_sourmash-species.txt` — species classification (sourmash)
- `*_checkm-results.tsv` — assembly completeness/contamination
- `annotation/` — Bakta annotation (.gff3, .tsv, .fna, .faa, .gbff)

## Appendix B: Tool Summary

| Analysis | Tool(s) | Nanopore-specific |
|----------|---------|-------------------|
| Read QC | NanoPlot, NanoStat | Yes |
| Assembly graph | Bandage | Neutral |
| Read mapping | minimap2 (map-ont preset) | Yes |
| Variant calling | Medaka, Clair3 | Yes (Nanopore-native) |
| SNP calling vs reference | Snippy (or minimap2+bcftools) | Hybrid |
| Structural variants | Sniffles2, SVIM | Yes |
| AMR gene screening | AMRFinderPlus, RGI (CARD), abricate | No |
| Methyltransferase HMM | HMMER3 (pfam models) | No |
| Species ID | FastANI, sourmash, MASH | No |
| Whole-genome alignment | MUMmer/nucmer, Progressive Mauve | No |
| Pan-genome | Roary, Panaroo | No |
| ICE detection | ICEfinder, CONJscan | No |
| Prophage detection | PHASTER, PhiSpy | No |
| Genomic islands | IslandViewer4, Alien Hunter | No |
| HGT signals | Alien Hunter, SIGI-HMM | No |
| CRISPR analysis | CRISPRdetect + BLAST | No |
| Structure prediction | AlphaFold2, PHYRE2 | No |
| Phylogenetics | MAFFT + IQ-TREE | No |
| Codon usage | CodonW, EMBOSS cusp | No |
