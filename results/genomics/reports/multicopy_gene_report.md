# Multi-Copy Gene Comparison: Unified Bakta Annotation

**Date:** 2026-03-23
**Method:** All genomes (M5, M6, M8, I6) annotated with Bakta v1.12.0, db-light v6.0
**Threshold:** >2 copies of the same product name per genome

---

## Summary

| Genome | Total CDS | Unique products | Multi-copy families (>2) |
|--------|-----------|-----------------|-------------------------|
| M5 | 6576 | 3361 | 358 |
| M6 | 6576 | 3361 | 358 |
| M8 | 6576 | 3361 | 358 |
| I6 | 6154 | 3145 | 256 |

M5, M6, and M8 are identical in gene content (clonal population).

## Classification

| Classification | Count | Description |
|---------------|-------|-------------|
| Shared (same copy #) | 50 | Ancestral paralogs, same count |
| Isolate expanded | 31 | Both have >2, isolates have more |
| I6 expanded | 25 | Both have >2, I6 has more |
| Isolate-only multi-copy | 252 | >2 in isolates, ≤2 in I6 |
| I6-only multi-copy | 150 | >2 in I6, ≤2 in isolates |

## Top isolate-expanded genes

| Gene | Product | M5 | M6 | M8 | I6 | Δ | Category |
|------|---------|----|----|----|----|---|----------|
| ABC transporter permease subun | ABC transporter permease subunit | 110 | 110 | 110 | 9 | +101 | ABC transporter |
| Helix-turn-helix | Helix-turn-helix domain protein | 81 | 81 | 81 | 7 | +74 | Other |
| Binding--dependent transport s | Binding--dependent transport system inne | 73 | 73 | 73 | 5 | +68 | Other |
| Bacterial extracellular solute | Bacterial extracellular solute-binding f | 71 | 71 | 71 | 5 | +66 | Other |
| Response regulator | Response regulator | 43 | 43 | 43 | 5 | +38 | Other |
| Helix-turn-helix | Helix-turn-helix domain-containing prote | 25 | 25 | 25 | 4 | +21 | Other |
| araJ | MFS transporter | 30 | 30 | 30 | 9 | +21 | MFS transporter |
| RNA polymerase sigma factor, s | RNA polymerase sigma factor, sigma-70 fa | 24 | 24 | 24 | 4 | +20 | Other |
| Helix-turn-helix | Helix-turn-helix family protein | 22 | 22 | 22 | 3 | +19 | Other |
| Lipoprotein | Lipoprotein | 35 | 35 | 35 | 18 | +17 | Other |
| Gfo/Idh/MocA family oxidoreduc | Gfo/Idh/MocA family oxidoreductase | 17 | 17 | 17 | 4 | +13 | Metabolism |
| Glycosyl transferases group 1 | Glycosyl transferases group 1 family pro | 16 | 16 | 16 | 3 | +13 | Metabolism |
| merR | MerR family transcriptional regulator | 14 | 14 | 14 | 4 | +10 | Regulator |
| ABC-2 transporter | ABC-2 transporter family protein | 12 | 12 | 12 | 3 | +9 | Other |
| Extracellular solute-binding p | Extracellular solute-binding protein | 20 | 20 | 20 | 12 | +8 | ABC transporter |
| MBL fold metallo-hydrolase | MBL fold metallo-hydrolase | 10 | 10 | 10 | 3 | +7 | Metabolism |
| SDR family NAD(P)-dependent ox | SDR family NAD(P)-dependent oxidoreducta | 10 | 10 | 10 | 3 | +7 | Metabolism |
| tetR | TetR family transcriptional regulator | 10 | 10 | 10 | 4 | +6 | Regulator |
| GNAT family N-acetyltransferas | GNAT family N-acetyltransferase | 9 | 9 | 9 | 4 | +5 | Metabolism |
| DUF2179 | DUF2179 domain-containing protein | 10 | 10 | 10 | 6 | +4 | DUF/domain |

## Top I6-expanded genes

| Gene | Product | I6 | M5 | M6 | M8 | Δ | Category |
|------|---------|----|----|----|----|---|----------|
| ABC transmembrane type-1 | ABC transmembrane type-1 domain-containi | 104 | 7 | 7 | 7 | -97 | Other |
| ABC transporter | ABC transporter domain-containing protei | 45 | 8 | 8 | 8 | -37 | ABC transporter |
| N-acetyltransferase | N-acetyltransferase domain-containing pr | 33 | 3 | 3 | 3 | -30 | Metabolism |
| ABC transporter substrate-bind | ABC transporter substrate-binding protei | 47 | 21 | 21 | 21 | -26 | ABC transporter |
| yozG | HTH cro/C1-type domain-containing protei | 27 | 6 | 6 | 6 | -21 | Other |
| Transcriptional regulator | Transcriptional regulator | 33 | 18 | 18 | 18 | -15 | Regulator |
| gerE | DNA-binding response regulator | 16 | 4 | 4 | 4 | -12 | Regulator |
| HTH arsR-type | HTH arsR-type domain-containing protein | 13 | 3 | 3 | 3 | -10 | Other |
| Sugar ABC transporter permease | Sugar ABC transporter permease | 13 | 3 | 3 | 3 | -10 | ABC transporter |
| araC | AraC family transcriptional regulator | 18 | 9 | 9 | 9 | -9 | Regulator |
| Transposase zinc-ribbon | Transposase zinc-ribbon domain-containin | 43 | 35 | 35 | 35 | -8 | IS transposase |
| ABC transporter permease | ABC transporter permease | 29 | 22 | 22 | 22 | -7 | ABC transporter |
| yobV | HTH deoR-type domain-containing protein | 10 | 4 | 4 | 4 | -6 | Other |
| Spore germination protein | Spore germination protein | 10 | 4 | 4 | 4 | -6 | Sporulation |
| histidine kinase | histidine kinase | 81 | 75 | 75 | 75 | -6 | Metabolism |
| Permease | Permease | 8 | 3 | 3 | 3 | -5 | ABC transporter |
| DNA 3'-5' helicase | DNA 3'-5' helicase | 7 | 3 | 3 | 3 | -4 | Other |
| SAM-dependent methyltransferas | SAM-dependent methyltransferase | 7 | 3 | 3 | 3 | -4 | Metabolism |
| VOC | VOC domain-containing protein | 7 | 3 | 3 | 3 | -4 | Other |
| ABC transporter ATP-binding pr | ABC transporter ATP-binding protein | 11 | 8 | 8 | 8 | -3 | ABC transporter |

## Shared ancestral paralogs (same copy number)

| Gene | Product | Copies | Category |
|------|---------|--------|----------|
| SLH | SLH domain-containing protein | 14 | Other |
| Copper amine oxidase-like N-te | Copper amine oxidase-like N-terminal domain-c | 12 | Other |
| DUF4367 | DUF4367 domain-containing protein | 9 | DUF/domain |
| 6-phospho-beta-glucosidase | 6-phospho-beta-glucosidase | 6 | Other |
| Glycosyl hydrolase | Glycosyl hydrolase | 6 | Metabolism |
| ubiE | Methyltransferase domain-containing protein | 6 | Metabolism |
| Pseudouridine synthase | Pseudouridine synthase | 6 | Metabolism |
| RNA polymerase sigma factor | RNA polymerase sigma factor | 6 | Other |
| Aminopeptidase | Aminopeptidase | 5 | Other |
| arsR | ArsR family transcriptional regulator | 5 | Regulator |
| Phage tail protein | Phage tail protein | 5 | Phage |
| Sigma-70 family RNA polymerase | Sigma-70 family RNA polymerase sigma factor | 5 | Other |
| Tetratricopeptide repeat prote | Tetratricopeptide repeat protein | 5 | Regulator |
| Beta-galactosidase | Beta-galactosidase | 4 | Other |
| Beta-xylanase | Beta-xylanase | 4 | Other |
| DUF58 | DUF58 domain-containing protein | 4 | DUF/domain |
| Efflux RND transporter peripla | Efflux RND transporter periplasmic adaptor su | 4 | Drug efflux |
| Gp5/Type VI secretion system V | Gp5/Type VI secretion system Vgr protein OB-f | 4 | Other |
| IDEAL | IDEAL domain-containing protein | 4 | Other |
| OmpR/PhoB-type | OmpR/PhoB-type domain-containing protein | 4 | Other |
| pucR | PucR family transcriptional regulator | 4 | Regulator |
| Sugar ABC transporter ATP-bind | Sugar ABC transporter ATP-binding protein | 4 | ABC transporter |
| VWA | VWA domain-containing protein | 4 | Other |
| yitT | YitT family protein | 4 | Other |
| nagZ | beta-N-acetylhexosaminidase | 4 | Other |
| norM | putative multidrug resistance protein NorM | 4 | Drug efflux |
| dapA | 4-hydroxy-tetrahydrodipicolinate synthase | 3 | Regulator |
| ABC transporter, permease prot | ABC transporter, permease protein | 3 | ABC transporter |
| ABC-2 type transporter transme | ABC-2 type transporter transmembrane domain-c | 3 | Other |
| ATPase | ATPase | 3 | Other |
| Alpha/beta hydrolase | Alpha/beta hydrolase | 3 | Metabolism |
| Binding-protein-dependent tran | Binding-protein-dependent transport systems i | 3 | Other |
| CBS | CBS domain-containing protein | 3 | Other |
| Cellobiose phosphorylase | Cellobiose phosphorylase | 3 | Other |
| cheA | Chemotaxis protein CheA | 3 | Motility |
| DNA-binding protein | DNA-binding protein | 3 | ABC transporter |
| DUF5643 | DUF5643 domain-containing protein | 3 | DUF/domain |
| Endospore germination permease | Endospore germination permease | 3 | ABC transporter |
| Glycoside hydrolase family 2 | Glycoside hydrolase family 2 | 3 | Metabolism |
| Glycosyl transferase | Glycosyl transferase | 3 | Metabolism |
| gntR | GntR family transcriptional regulator | 3 | Regulator |
| Permease IIC component | Permease IIC component | 3 | ABC transporter |
| Phosphoenolpyruvate-protein ph | Phosphoenolpyruvate-protein phosphotransferas | 3 | Metabolism |
| Sensor histidine kinase | Sensor histidine kinase | 3 | Metabolism |
| Spore coat protein | Spore coat protein | 3 | Sporulation |
| Sucrose-6-phosphate hydrolase | Sucrose-6-phosphate hydrolase | 3 | Metabolism |
| Transporter associated | Transporter associated domain protein | 3 | Other |
| beta-fructofuranosidase | beta-fructofuranosidase | 3 | Other |
| fructokinase | fructokinase | 3 | Metabolism |
| peptidoglycan glycosyltransfer | peptidoglycan glycosyltransferase | 3 | Metabolism |

## Is it normal to have multiple copies?

Yes. The following are **expected** multi-copy families:

- **IS transposases**: IS110 (29 in M5 vs 0 in I6), IS4 (19 vs 0). Mobile elements that self-replicate.
- **ABC transporters/permeases**: Largest paralogous family in bacteria. Normal for 7+ Mb genomes.
- **Transcriptional regulators**: Scale with genome size. TetR, MarR, MerR, LysR families.
- **rRNA operons**: 8× 16S in M5 — normal for fast-growing Firmicutes.
- **Phage genes**: SP-beta prophage contributes ~30 structural genes.

## Relevance to spectinomycin resistance

**None.** rpsE is single-copy in all genomes. No spectinomycin resistance genes are among the multi-copy families. The 2 novel drug efflux genes (MFS drug:H+ antiporter, MATE efflux) are not spectinomycin-specific.
