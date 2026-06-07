"""
viability.py
============
H0: the generic distress / viability decomposition.

Rationale (Szalai et al. 2019): a cell-death/proliferation signature is a major
component of perturbation transcriptomes and makes unrelated perturbations look
alike for non-mechanistic reasons. Before attributing convergence to any
specific mechanism (H1-H5), the viability axis must be removed and the residual
convergence re-measured. The residual is the scientific signal.

This module:
  1) scores each observation on the viability signature,
  2) regresses that score out of every gene (linear residualization),
  3) provides the relationship test: do perturbations with *similar* viability
     scores converge more? (sick cells looking alike).
"""

from __future__ import annotations
import numpy as np


def score_signature(adata, genes):
    """Per-observation signature score = mean of z-normalized expression over
    the signature genes that are present in var_names.

    Because .X is already per-gene z-normalized, this score is a signed loading
    on the signature axis (positive = signature up).
    """
    present = [g for g in genes if g in set(adata.var_names)]
    if not present:
        raise ValueError("None of the signature genes are in adata.var_names.")
    cols = [adata.var_names.get_loc(g) for g in present]
    X = np.asarray(adata.X, dtype=float)
    return X[:, cols].mean(axis=1), present


def regress_out_score(adata, score):
    """Return a copy of adata with the linear effect of `score` removed from
    every gene (column-wise OLS residual across observations).

    For gene column x and centered score s_c:
        beta = (x . s_c) / (s_c . s_c);  residual = x - beta * s_c
    Vectorized across all genes simultaneously.
    """
    X = np.asarray(adata.X, dtype=float).copy()
    s = np.asarray(score, dtype=float)
    # `score` is expected to be centered by the caller so controls sit at ~0
    # (then regression leaves controls untouched and the null is stable).
    denom = float(s @ s)
    if denom <= 0:
        return adata.copy()
    betas = (X.T @ s) / denom              # length n_genes
    X_resid = X - np.outer(s, betas)       # remove the score-aligned component
    out = adata.copy()
    out.X = X_resid
    return out


def viability_pair_similarity(score, is_ctrl, groups):
    """For each cross-annotation perturbation pair, return (similarity, ) arrays
    aligned with the upper-triangle of the perturbation correlation matrix.

    similarity = product of (mean-centered) viability scores: large positive
    when both perturbations are on the same side of the viability axis
    (e.g. both 'sick'). Used to test whether shared viability explains
    pairwise convergence.
    """
    s = np.asarray(score, dtype=float)[~is_ctrl]
    s_c = s - s.mean()
    g = np.asarray(groups)[~is_ctrl]
    n = s_c.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    between = g[iu] != g[ju]
    sim = s_c[iu] * s_c[ju]
    return sim[between], (iu[between], ju[between])
