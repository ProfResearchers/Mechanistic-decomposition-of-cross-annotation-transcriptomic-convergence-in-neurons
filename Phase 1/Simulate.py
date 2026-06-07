"""
simulate.py
===========
Generate a synthetic perturbation x gene AnnData that matches the schema of the
Replogle et al. 2022 z-normalized pseudobulk file, so the SAME pipeline runs on
synthetic and real data.

Planted structure (ground truth for validation)
------------------------------------------------
1) VIABILITY axis  : a gene-weight vector loading on cell-cycle + death/stress
   genes. A random subset of perturbations carries a (signed) dose of it. This
   is the GENERIC DISTRESS confound (H0) and is spread across annotation groups.
2) ISR-like axis   : a separate gene-weight vector loading on ISR/ATF4 genes. A
   *cross-annotation* subset of perturbations carries a positive dose. This is
   the GENUINE mechanism-specific convergence that should SURVIVE viability
   regression.
3) Idiosyncratic   : sparse gene-specific effects unique to each perturbation
   (the non-convergent part).
4) Noise           : Gaussian.

Controls (non-targeting) carry noise only. Finally every gene is z-normalized
across observations (matching the real file).

Expected pipeline behaviour
----------------------------
- Before regression: high cross-annotation convergence (viability + ISR).
- After regressing out the viability signature: convergence drops by the
  viability contribution but a positive residual (ISR) remains.
- Viability-pair-similarity should correlate with pairwise convergence before
  regression and much less after.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import anndata as ad

from gene_sets import VIABILITY_SIGNATURE, ISR_ATF4


def simulate(
    n_perturbations=300,
    n_controls=60,
    n_generic_genes=800,
    n_annotation_groups=6,
    frac_viability=0.60,     # fraction of perturbations carrying the distress axis
    frac_isr=0.45,           # fraction carrying the genuine ISR convergence
    viability_strength=1.5,
    isr_strength=1.8,
    idiosyncratic_strength=0.7,
    noise=1.0,
    seed=0,
):
    rng = np.random.default_rng(seed)

    # ---- gene universe: signature genes + generic genes ----
    sig_genes = list(dict.fromkeys(VIABILITY_SIGNATURE + ISR_ATF4))  # unique, ordered
    generic_genes = [f"GENE{i:04d}" for i in range(n_generic_genes)]
    var_names = sig_genes + generic_genes
    G = len(var_names)
    gidx = {g: i for i, g in enumerate(var_names)}

    # ---- planted axes (gene-weight vectors) ----
    viability_axis = np.zeros(G)
    for g in VIABILITY_SIGNATURE:
        viability_axis[gidx[g]] = rng.normal(1.0, 0.2)
    # broaden the program onto a halo of generic genes (real stress programs
    # involve hundreds of genes); coherent sign so distressed cells share it
    generic_pool = np.arange(len(sig_genes), G)
    halo_v = rng.choice(generic_pool, size=120, replace=False)
    viability_axis[halo_v] += np.abs(rng.normal(0.6, 0.2, size=halo_v.size))

    isr_axis = np.zeros(G)
    for g in ISR_ATF4:
        isr_axis[gidx[g]] = rng.normal(1.0, 0.2)
    # disjoint halo so the ISR program is well separated from the viability axis
    remaining_pool = np.setdiff1d(generic_pool, halo_v)
    halo_i = rng.choice(remaining_pool, size=120, replace=False)
    isr_axis[halo_i] += np.abs(rng.normal(0.6, 0.2, size=halo_i.size))

    n = n_perturbations
    # annotation groups assigned at random (independent of the planted axes)
    groups = rng.integers(0, n_annotation_groups, size=n)

    # which perturbations carry each axis
    carries_via = rng.random(n) < frac_viability
    via_dose = np.where(carries_via, np.abs(rng.normal(1.0, 0.4, size=n)), 0.0) * viability_strength

    carries_isr = rng.random(n) < frac_isr
    isr_dose = np.where(carries_isr, np.abs(rng.normal(1.0, 0.4, size=n)), 0.0) * isr_strength
    # ensure the ISR carriers span multiple annotation groups (cross-annotation)
    # (they already do, since group assignment is independent; assert later)

    # ---- build perturbation matrix ----
    Xp = (
        np.outer(via_dose, viability_axis)
        + np.outer(isr_dose, isr_axis)
        + idiosyncratic_strength * _sparse_idiosyncratic(rng, n, G, density=0.02)
        + noise * rng.normal(0, 1.0, size=(n, G))
    )

    # ---- controls: noise only ----
    Xc = noise * rng.normal(0, 1.0, size=(n_controls, G))

    X = np.vstack([Xp, Xc])

    # ---- z-normalize each gene RELATIVE TO THE CONTROLS ----
    # Perturbation effects are deviations from the non-targeting baseline; this
    # keeps the control-control correlation null near zero, exactly as in the
    # real Replogle z-normalized file (where NTCs sit at the population centre).
    ctrl_slice = X[n:]                      # the n_controls control rows
    mu = ctrl_slice.mean(axis=0, keepdims=True)
    sd = ctrl_slice.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    X = (X - mu) / sd

    # ---- obs metadata ----
    pert_names = [f"PERT{i:04d}" for i in range(n)]
    ctrl_names = [f"NTC{i:03d}" for i in range(n_controls)]
    obs = pd.DataFrame(
        {
            "perturbation": pert_names + ctrl_names,
            "annotation_group": [f"grp{g}" for g in groups] + ["control"] * n_controls,
            "is_control": [False] * n + [True] * n_controls,
            # ground-truth flags (for validation only; not used by the pipeline)
            "_truth_viability": list(carries_via) + [False] * n_controls,
            "_truth_isr": list(carries_isr) + [False] * n_controls,
        }
    )
    obs.index = obs["perturbation"].values

    var = pd.DataFrame(index=var_names)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["sim_params"] = dict(
        n_perturbations=n, n_controls=n_controls, frac_viability=frac_viability,
        frac_isr=frac_isr, viability_strength=viability_strength,
        isr_strength=isr_strength, seed=seed,
    )
    return adata


def _sparse_idiosyncratic(rng, n, G, density=0.02):
    M = np.zeros((n, G))
    k = max(1, int(density * G))
    for i in range(n):
        cols = rng.choice(G, size=k, replace=False)
        M[i, cols] = rng.normal(0, 1.0, size=k)
    return M
