# Phase 1 — Convergence metric + viability decomposition (Aim 1 + H0)

This is the first executable phase of the in silico research protocol
(*Neuronal_Convergence_InSilico_Research_Protocol.md*). It does two things, on a
real Replogle-format Perturb-seq AnnData or on a built-in synthetic dataset:

1. **Aim 1 — quantify cross-annotation convergence.** Represent each
   perturbation by its per-gene z-normalized pseudobulk profile, measure
   pairwise profile correlation, and report convergence **relative to a
   non-targeting-control null**, split into within-annotation and
   cross-annotation (functionally unrelated) pairs.
2. **H0 — remove the generic distress / viability axis.** Score each
   perturbation on a cell-cycle + death/stress signature, regress it out of
   every gene, and re-measure convergence on the **same features**. The drop is
   the fraction of apparent convergence that is generic distress; the residual
   is the mechanism-specific signal that H1–H5 then explain.

## Validation status

The pipeline is validated on a synthetic dataset engineered to match the
Replogle AnnData schema, with a planted viability axis and a *separate* planted
"genuine" cross-annotation program (the ground truth). On the latest run:

- cross-annotation convergence above null **before** = +0.116 (z ≈ 3.6)
- **after** viability regression = +0.024 → **~80% of convergence was the
  viability axis, ~20% mechanism-specific residual**
- the viability-similarity ↔ convergence relationship collapses (Spearman
  rho 0.45 → 0.06)
- the planted genuine program survives and is cleanly separated from
  non-program perturbations (Mann-Whitney p ≈ 0), so the **PASS** check is True.

In short: the decomposition provably separates generic distress from genuine
convergence on data where we know the answer, *before* it is applied to real
data where we do not.

## Important: sandbox vs. your environment

This pipeline was developed in an environment whose network reaches package
registries and GitHub but **not** the data hosts (Figshare+, Zenodo, GEO/SRA).
So it was validated on synthetic data here. The **identical code** runs on the
real Replogle file in any environment with normal internet access — see below.

## Install

```bash
pip install -r requirements.txt
```

## Run on the synthetic validation dataset

```bash
cd src
python run_phase1.py --outdir ../outputs
```

## Run on the real Replogle atlas

1. Download a processed AnnData from Figshare+ (DOI `10.25452/figshare.plus.20029387`).
   The recommended starting file is the K562 *essential* z-normalized pseudobulk,
   `K562_essential_normalized_bulk_01.h5ad` (smaller than the genome-wide file).
   (`pertpy`/`scPerturb` loaders are an alternative source.)
2. Prepare two `obs` columns:
   - **control flag** (`--control-key`): boolean, `True` for non-targeting
     controls. In the Replogle obs this is the non-targeting guide class; create
     a boolean column from it.
   - **annotation group** (`--group-key`): the grouping whose *between-group*
     pairs are treated as cross-annotation. For Phase 1 a coarse placeholder
     grouping is fine — the strict GO/CORUM/STRING **orthogonality filter** is a
     later module. Until then, "cross-annotation" = "different placeholder
     group," which you should state as a limitation.
3. Run:

```bash
cd src
python run_phase1.py \
  --input /path/to/K562_essential_normalized_bulk_01.h5ad \
  --group-key annotation_group \
  --control-key is_control \
  --n-top-genes 2000 \
  --outdir ../outputs
```

Notes for real data:
- The cell-cycle gene lists in `gene_sets.py` use a few legacy symbols
  (`MLF1IP`=CENPU, `FAM64A`=PIMREG, `HN1`=JPT1). Add an alias-mapping step if your
  `var_names` use current symbols, so the viability signature matches fully.
- `--n-top-genes` selects highly variable genes by variance; tune to the atlas.
- Outputs: a 3-panel figure (`phase1_convergence.png`), the residual AnnData
  (`residual_after_viability.h5ad`, the input to H1/H3 signature scoring), and a
  printed report with the decomposition headline.

## Layout

```
neuro_convergence_phase1/
  README.md
  requirements.txt          pinned, verified environment
  REGISTRY.yaml             data provenance (accessions/DOIs; fill in as you download)
  src/
    gene_sets.py            cell-cycle + death/stress (viability signature); ISR set
    convergence.py          Aim 1: feature selection, profile correlation, null-calibrated summary
    viability.py            H0: signature scoring, regress-out, relationship test
    attribution.py          Phase 2: generic signature-axis attribution + dominant-axis alignment
    simulate.py             synthetic Replogle-schema generator (doubles as a unit test)
    run_phase1.py           Phase 1 runner (Aim 1 + H0; synthetic by default, --input for real)
    run_phase2.py           Phase 2 runner (H1 ISR attribution; synthetic by default, --input for real)
  outputs/
    phase1_convergence.png
    residual_after_viability.h5ad
    phase2_isr_attribution.png
    residual_after_viability_and_isr.h5ad
  neuro_convergence_phase1_colab.ipynb   self-contained Colab notebook (Aim 1 + H0)
  neuro_convergence_phase2_colab.ipynb   self-contained Colab notebook (H1 ISR attribution)
```

## Phase 2 - mechanism attribution (Aim 2 / H1)

Phase 2 takes the post-viability residual and asks whether it is the integrated
stress response (ISR/ATF4). The attribution engine (`attribution.py`) is generic:
score the residual on a signature, regress it out, measure the drop
(`frac_explained`, with a bootstrap CI), and test whether the dominant convergent
axis (top principal component) aligns with the signature. The identical call later
serves the NMD axis (H3).

Run it: `python run_phase2.py` (synthetic) or `python run_phase2.py --input ...`.

On the synthetic ground truth the engine recovers the planted ISR program: ISR
loading cleanly separates true carriers (MWU p~1e-42), the dominant convergent axis
aligns with ISR well above chance (|cosine| 0.51 vs ~0.15), and among
ISR-engaging perturbations the convergence collapses to null when ISR is removed
(0.082 -> -0.006); PASS = True.

**Read the headline fraction as a lower bound.** A single global "ISR explains X%
of the residual" understates ISR's role, because cross-annotation convergence is
diluted by perturbations that do not engage ISR and because a marker set only
partially proxies a broad program. The robust evidence is the dominant-axis
alignment and the among-engagers collapse, not the percentage.

## What Phase 1 deliberately does NOT do yet

- **No annotation-orthogonality filter.** Cross-annotation = between placeholder
  group. The GO-semantic-similarity / CORUM / STRING filter is the next module;
  until it is in place, convergence among "unrelated" genes is only as strict as
  your grouping.
- **No mechanism attribution beyond H0.** H1 (ISR/ATF4) and H3 (NMD) scoring of
  the residual is the immediate next step (decoupler signature scoring), using
  the saved `residual_after_viability.h5ad`.
- **Pseudobulk correlation, not E-distance.** E-distance (pertpy) is the
  single-cell complement and a planned sensitivity analysis; Phase 1 uses the
  correlation convention on pseudobulk for speed and comparability with Replogle.
- **No causal claims.** Per the protocol's epistemic scope, this is consistency/
  decomposition, not intervention.

## Next

- **H3 (NMD / transcriptional adaptation)** - reuse `attribute_to_signature` with
  an NMD-response gene set, then the cross-modality natural experiment (CRISPR-KO
  vs CRISPRi in scPerturb) and the gnomAD/GTEx pLoF NMD-trigger-vs-escape contrast.
- **H2 (network pleiotropy)** - predict pairwise convergence from trans-network
  distance (eQTLGen) vs GO semantic similarity; needs external downloads.
- **Empirical ISR axis** - on real data, define the ISR direction from the atlas's
  own ATF4 / eIF2-pathway perturbations and project onto it (stronger than the
  curated marker set used here).
- **Annotation-orthogonality filter** - turns the placeholder grouping into a
  defensible "functionally unrelated" definition (GO semantic similarity + CORUM +
  STRING); the key infrastructure step for the cross-annotation claim on real data.
- Draft the **OSF pre-registration** from the protocol before running on real data,
  to lock the metric, signatures, and thresholds.
