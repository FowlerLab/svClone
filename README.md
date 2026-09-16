# svClone

Design of Golden Gate variant libraries built from **mutagenic oligo tiles** plus
**synthesized gene blocks**.

Given one or more codon-optimized CDSs, the notebooks here split each gene into
Golden Gate blocks, pick BsaI breakpoints whose 4-nt overhangs assemble with high
fidelity, write one 250-nt oligo per variant covering the mutagenic window of a
block, emit the remainder of each block as a gBlock/PCR amplicon, and assign an
orthogonal PCR primer pair to every block so each block's variants can be
amplified as an independent subpool. Outputs are order-ready sheets for an oligo
array (Twist/Agilent), gBlocks, and IDT primer plates.

## Contents

```
find_breakpoints_BsaI_pilotcircRNA_010324.ipynb   notebook 1 — Jan 2024 pilot library
design_MAPK_library_validatedprimers_V2_053024.ipynb   notebook 2 — production MAPK library

RR_orthogonalFv1_plate.xlsx / RR_orthogonalRv1_plate.xlsx   82-well orthogonal primer plates
bsaI_empirical.csv                                BsaI overhang fidelity matrix
MANE.GRCh38.v1.3.summary.txt / .refseq_rna.gbff    RefSeq MANE annotations (notebook 2)
Validated Gene Primer Pair Combinations.csv        bench-validated primer pairs (notebook 2)

pilot_order_011324/     notebook 1 output — the pilot array as submitted
MAPK_ordering/          notebook 2 output — the main MAPK library sheets
final_order/            what was actually ordered (MAPK + PTEN + active-learning, merged by hand)
orthogonal_primer_random_combos_used_SP_060524_v2.csv   notebook 2 output — record of the primer draw
```

Everything else at the top level is reference data that neither notebook reads —
see [Files not on the critical path](#files-not-on-the-critical-path).

---

## The design logic

Every gene is processed the same way:

1. **Cap the CDS.** The gene is flanked with a BsaI overhang (`CGTC` … `GCAT`) and
   `GCTCTTCC` SapI sites, so the finished insert can be dropped into the vector
   by a SapI cut and the internal joints are BsaI.
2. **Pick blocks and breakpoints.** The CDS is divided into blocks of roughly
   `block_size_range` nt. Each candidate breakpoint is allowed to slide ±`slack`
   nt; every combination is scored against an **empirical overhang-ligation
   matrix** (`bsaI_empirical.csv`, the 256×256 overhang fidelity table). The
   score is the product of on-target ligation fractions for the overhangs in
   play, and a breakpoint set scores 0 if any two overhangs are duplicated or if
   an overhang is on the blacklist (palindromic, or < 300 reads on-target).
   `optimize_gene()` reports whether all regions came out high (≥ 0.95), medium
   (≥ 0.9) or low fidelity.
3. **Split each block into oligo + amplicon.** One sub-fragment of each block is
   the *mutagenic window* — short enough to fit on a synthesized oligo together
   with its primer arms and BsaI sites (`max_oligo_size = 250`). The rest of the
   block is ordered as a gBlock and/or amplified with designed PCR primers.
   Windows are checked to overlap so that no codon falls between two blocks.
4. **Enumerate variants.** For every codon in the window: all 19 missense
   substitutions (most-used codon for the target amino acid), a synonymous
   variant, a stop (`TAA`), and all three 3-nt deletions. Any of those categories
   can be switched off, or replaced with an explicit variant list.
5. **Assign orthogonal primers.** Each block gets one forward/reverse pair from
   the 82-well orthogonal primer plates, so blocks are separately amplifiable.
   Primers that contain a BsaI site, or that cross-prime anywhere in the gene
   set, are dropped first (`check_nonspecific()`, adapted from DIMPLE,
   doi:10.1186/s13059-023-02880-6).
6. **QC and write.** `post_qc()` re-checks every chosen *and* unused primer
   against every WT oligo for off-target annealing; oligos carrying stray
   BsaI/SapI sites or duplicate sequences are deleted; gBlocks under 300 bp are
   padded with a fixed random spacer for the vendor minimum. Everything is then
   written to the order sheets.

Notebook 1 is the reference implementation of the above. Notebook 2 is the same
algorithm with the extensions listed in its section.

---

## Requirements

Python 3 with:

```
biopython  numpy  pandas  scipy  openpyxl  primer3-py  matplotlib
```

`matplotlib` is only needed by `design_MAPK_library_validatedprimers_V2_053024.ipynb`.
The first cell of each notebook has `pip install` lines (commented out in the
MAPK notebook); note that neither `pip` cell lists `primer3-py`, which both
notebooks import, so install it yourself.

Run the notebooks from the repository root — every path in them is relative
(`./…`).

---

## Notebook 1 — `find_breakpoints_BsaI_pilotcircRNA_010324.ipynb`

**Purpose:** the January 2024 pilot library. Small and self-contained — gene
sequences are pasted into the notebook, and all the design functions live in one
cell — so this is the one to read first.

**Inputs** (both present in the repo):

| File | Used for |
| --- | --- |
| `RR_orthogonalFv1_plate.xlsx`, `RR_orthogonalRv1_plate.xlsx` | 82-well orthogonal primer plates; the last 20 nt of each sequence is the variable priming region |
| `bsaI_empirical.csv` | 256×256 BsaI overhang ligation-fidelity matrix |

**Walkthrough**

| Cells | What happens |
| --- | --- |
| 0 | Imports / `pip install` |
| 1–3 | Read both primer plates, take the 20-nt `PrimerEnd`, flag primers containing a BsaI site |
| 4 | `check_nonspecific()` — Tm-based off-target annealing scan (BioPython `Tm_NN`, falling back to `primer3.calcHeterodimer` where the NN table is insufficient) |
| 5 | **Gene definitions.** 21 codon-optimized CDS fragments pasted as `Seq` objects (BRAF, KRAS, MRAS, EGFR, ERBB2, SHP2, SOS2, ARAF, CRAF, KSR1/2, MEK1/2), split into n-term/c-term halves where the protein is too long for one superblock |
| 6 | Cross-prime every plate primer against every gene → `Num_Nonspecific_Binding_Sites`. **This is the slow cell — hours, not minutes.** |
| 7 | Keep only primers with no BsaI site and zero off-target hits, drop one hand-blacklisted sequence, pair F with R into `orthogonal_primers_touse` |
| 8–11 | Load the BsaI matrix and the codon-usage table, build the overhang blacklist, set defaults (`block_size_range = [155, 175]`, `max_oligo_size = 250`, `slack = 5`, gBlock padding spacer) |
| 12 | All design functions: `post_qc`, `compute_overlaps`, `score_breakpoints`, `optimize_breakpoints`, `optimize_gene`, `generate_primer`, `make_all_mutations`, `write_oligo_library` |
| 13 | Expand the primer pool by cyclically permuting the reverse primers by 18/36/54 → `all_combos` (4× as many usable pairs) |
| **14** | **The pilot order.** Full variant libraries for `BRAFnterm`, `KRAS`, `CRAFnterm` → `pilot_order_011324/circRNApilot_*_v3.*` |
| 15–16 | WT-only "test amplicon" runs over the other 18 genes at two different primer-set offsets → `circRNApilot_testamp*_v3.*` / `_v4.*`, used to check that each primer pair amplifies only its own subpool |
| 17 | `# TOTAL LIBRARY` — the full 18-gene variant library, written to `./bsaI_test/…`. **Create `bsaI_test/` first or this cell fails**; it was superseded by notebook 2. |

**Outputs in the repo** (`pilot_order_011324/`): `circRNApilot_oligos_v3.csv`,
`_primers_v3.tsv`, `_gblocks_v3.tsv`, `_ampkey_v3.tsv` from cell 14, plus
`circRNApilot_FINALoligos.csv` — the array actually submitted, which is
`oligos_v3` plus the 960 WT test-amplicon oligos from cells 15–16.

The `circRNApilot_testamp*` files from cells 15–16 are not in the repo, and
`bsaI_test/` (cell 17) was never committed. Re-running those cells recreates them
— but see the note on unseeded primer draws below: the regenerated sheets will
not match what `circRNApilot_FINALoligos.csv` was built from.

---

## Notebook 2 — `design_MAPK_library_validatedprimers_V2_053024.ipynb`

**Purpose:** the production MAPK-pathway library (May/June 2024). Same core
algorithm as notebook 1, extended with:

- **Gene sequences pulled from RefSeq MANE** instead of pasted (alongside pasted
  codon-optimized versions of the fragments being synthesized).
- **Pre-validated primer pairs** — combinations already shown to work at the
  bench are used first, and the rest are drawn randomly from the surviving pool.
- **`aa_start`** so variants in a c-terminal fragment are numbered by their real
  position in the full-length protein.
- **`tile_boundaries`** to hand-specify breakpoints instead of optimizing them
  (used to re-order PTEN against previously cloned tiles).
- **`mutations_to_use`** to build a library from an explicit variant list rather
  than all singles (used for the active-learning libraries).
- **`breakpoint_file`** output recording each block's mutagenic window, plus a
  check that consecutive windows overlap.
- gBlocks split by length into a normal and a `>1000 bp` order sheet.

**Inputs**

| File | Present? | Used for |
| --- | --- | --- |
| `RR_orthogonalFv1_plate.xlsx`, `RR_orthogonalRv1_plate.xlsx` | yes | orthogonal primer plates |
| `MANE.GRCh38.v1.3.summary.txt` | yes | MANE Select transcript per gene symbol |
| `MANE.GRCh38.v1.3.refseq_rna.gbff` | yes (290 MB) | GenBank stream the CDS features are extracted from |
| `Validated Gene Primer Pair Combinations.csv` | yes | bench-validated F/R well pairs |
| `bsaI_empirical.csv` | yes | overhang fidelity matrix |
| `RAF1_variants_v2.json`, `BRAF_variants_v2.json`, `HOXA11_variants_v2.json` | **no** | variant lists for the active-learning libraries (cells 28–33) |
| `DH_Library_20240507_fixed.csv` | **no** | cross-checking our primers against a separate oligo array (cell 26, QC only) |
| `MAPK_ordering/circRNA_oligos_MAPKminus.csv` | **no** | cross-checking AAV HDR primers against the WT oligos (cell 25, QC only) |

**Walkthrough**

| Cells | What happens |
| --- | --- |
| 0 | Imports |
| 1–4 | Read the primer plates; take the 20-nt `PrimerEnd`; mark wells empirically known to be bad (`Exclude`) or weak (`Worse`); flag BsaI sites |
| 5 | `check_nonspecific()` |
| 6–7 | Read the MANE summary, stream the `.gbff`, pull the CDS for the 20 genes in `genes_of_interest`, print CDS + translation to eyeball them |
| 8 | **Gene definitions** — codon-optimized CDS fragments as `Seq` objects, named `<GENE>nterm` / `<GENE>cterm` |
| 9 | Cross-prime every plate primer against every gene fragment. **Slow cell — hours.** |
| 10 | Filter primers: no BsaI site, no off-target hit, not `Exclude`, not hand-blacklisted |
| 11–12 | Read `Validated Gene Primer Pair Combinations.csv`, drop validated pairs whose wells did not survive filtering, attach sequences |
| 13 | Draw 3 random F×R pairings of the surviving primers, drop any duplicating a validated pair, and **save the draw** to `orthogonal_primer_random_combos_used_SP_060524_v2.csv` |
| 14–17 | BsaI matrix, codon table, overhang blacklist, defaults (`block_size_range = [153, 172]`, `first_last_block_reduction = 8`, `slack = 5`, padding spacer) |
| 18 | All design functions, including the extended `write_oligo_library()` |
| 19 | **PTEN re-order** — fixed `tile_boundaries`, blocks 1, 2, 4–8 only, no PCR primers → `./PTEN/` |
| **20** | **The main MAPK library** — 21 gene fragments, validated primers first, `aa_start` per fragment → `MAPK_ordering/circRNA_*` |
| 21–22 | Histograms of oligo lengths and of mutagenic-window overlaps (sanity checks) |
| 23–24 | The gene pieces bought as cloned WT constructs rather than tiled, SapI-capped → `MAPK_ordering/circRNA_WT_gblocks_cloned.tsv` |
| 25–27 | QC cross-checks against two external oligo arrays (needs the two missing CSVs) |
| 28–30 | Load the active-learning variant JSONs, shift HOXA11 numbering to start at 243, plot variants per position |
| 31–33 | **Active-learning libraries** — single-block designs driven by `mutations_to_use` (4160 variants per gene), plus an all-singles HOXA11 library → `./Active_Learning/` |

**Outputs**

- `MAPK_ordering/` — cell 20's order sheets (`circRNA_oligos.csv`,
  `circRNA_primers.tsv`, `circRNA_gblocks.tsv`, `circRNA_gblocks_large.tsv`,
  `circRNA_ampkey.tsv`, `circRNA_breakpoints.tsv`).
- `final_order/` — **what was actually ordered**: the same sheets with cell 19's
  PTEN blocks and cells 32–33's active-learning blocks merged in (187 vs. 176
  ampkey rows), renamed `MAPK_*_0605*`. Not written by any notebook; assembled by
  hand from the three runs. `MAPK_WT_gblocks_cloned_060324_v1.tsv` is cell 24's
  output.
- `orthogonal_primer_random_combos_used_SP_060524_v2.csv` — the record of cell
  13's random primer draw.
- `./PTEN/` and `./Active_Learning/` were never committed. **Create those
  directories before running cells 19 and 31–33**, or the writes fail.

> **The random primer draw is not seeded.** Cell 13 calls `.sample()` with no
> `random_state`, so re-running the notebook assigns *different* primer pairs to
> blocks than the committed order sheets used. Treat
> `orthogonal_primer_random_combos_used_SP_060524_v2.csv` and
> `MAPK_ordering/circRNA_ampkey.tsv` as the authoritative record of which primer
> amplifies which block; do not assume a fresh run reproduces them.

---

## Output file formats

None of the order sheets carry a header except the amp keys and breakpoints.

| File | Format |
| --- | --- |
| `*oligos*.csv` | `oligo_name,sequence` — the synthesis array |
| `*primers*.tsv` | `name`, `sequence`, `25nm`, `STD` — IDT plate order format |
| `*gblocks*.tsv`, `*gblocks_large*.tsv` | `name`, `sequence` — split at 1000 bp; anything under 300 bp is 5′-padded with the fixed random spacer |
| `*ampkey*.tsv` | `Gene`, `Block`, `Forward Primer Well`, `Forward Primer`, `Reverse Primer Well`, `Reverse Primer`[, `Validated`] — which primer pair amplifies which block's subpool |
| `*breakpoints*.tsv` | unnamed index column, then `Gene`, `Block`, `Mutagenesis Start`, `Mutagenesis End` — the mutagenized nt window of each block (notebook 2 only) |

**Oligo naming:** `<gene>_block<N>_<variant>`, where `<variant>` is

- `WT` — the unmutated block,
- `K427A` — missense, `<wt aa><position><new aa>`; position is 1-based in the
  fragment, or in the full-length protein if `aa_start` was given,
- `S429S` — synonymous (most-used codon that isn't the one in the gene),
- `R509X` — stop (`TAA`),
- `del1287` — 3-nt deletion starting at that nucleotide,
- `variant7` — when built from an explicit `mutations_to_use` list; the number
  indexes into that list, so keep the source JSON to decode it.

**Amplicon/gBlock naming:** `<gene>_block<N>_s<k>` for the non-oligo
sub-fragments, and `<gene>_gene_ampF` / `_ampR` for the primers that amplify the
whole fragment.

---

## Files not on the critical path

Neither notebook reads any of the following. Everything the two notebooks *do*
need is present.

**Downstream sequencing**

- `P5_information.tsv`, `P5_order_IDT.tsv`, `P7_information.tsv`,
  `P7_order_IDT.tsv` — indexed Illumina P5/P7 primers for the NGS readout, not
  for library construction.

**Housekeeping**

- `.DS_Store` at the top level and in `pilot_order_011324/` — Finder artifacts,
  both currently tracked by git. Safe to remove; there is no `.gitignore` in the
  repo, so they will come back.

---

## Known rough edges

- **Output directories are not created for you.** Before running, `mkdir` any of
  `bsaI_test/`, `PTEN/`, `Active_Learning/` you intend to write to.
- **Missing inputs.** The active-learning cells (notebook 2, 28–33) need the
  three `*_variants_v2.json` files; the QC cells (25–27) need
  `DH_Library_20240507_fixed.csv` and `MAPK_ordering/circRNA_oligos_MAPKminus.csv`.
  None are in the repo. The main library cells (19, 20) run without them.
- **The cross-priming cells are the bottleneck** (notebook 1 cell 6, notebook 2
  cell 9): `check_nonspecific()` is a pure-Python scan of 82 primers × every
  position of every gene, in both orientations. Expect hours, and run them once
  per gene set rather than once per design iteration.
- **Genes must be clean and in frame.** `write_oligo_library()` silently skips
  any gene whose length isn't divisible by 3, or that contains a BsaI or SapI
  site — it prints the reason and moves on, so read the per-gene log rather than
  assuming every gene you passed made it into the output.
- **Proteins over ~620 aa** are rejected by `optimize_gene()`; split them into
  n-term/c-term fragments by hand, as the gene definition cells do.
