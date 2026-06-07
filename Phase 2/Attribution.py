"""
attribution.py
==============
Aim 2 / H1: attribute the *residual* convergence (what survives the H0 viability
regression) to a mechanistic signature axis.

The engine is generic: `attribute_to_signature` takes any gene set and answers
two questions about the residual cross-annotation convergence:

  1. **How much of it does this signature explain?** Score each perturbation on
     the signature, regress that score out (same primitive as H0), and measure
     the drop in cross-annotation convergence -> `frac_explained` (with a
     bootstrap CI).
  2. **Is the dominant convergent axis this signature?** Take the top principal
     component of the residual perturbation profiles (the leading axis of shared
     variation) and measure its alignment (|cosine|) with the signature ->
     `pc1_align`, alongside the variance it explains.

For H1 the signature is ISR/ATF4. The identical call with an NMD gene set serves
H3 later. This reuses `score_signature` / `regress_out_score` from viability.py
and `convergence_summary` from convergence.py, so H0 and H1 share one mechanism.
"""

from __future__ import annotations
import numpy as np

from convergence import convergence_summary
from viability import score_signature, regress_out_score


def dominant_axis(X):
    """Top principal component of perturbation profiles (rows=perturbations,
    cols=features). Returns (loadings_over_features, fraction_variance_explained).
    """
    Xc = X - X.mean(axis=0, keepdims=True)
    # economy SVD; first right-singular vector = top PC loadings over genes
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Vt[0]
    var_explained = float((S[0] ** 2) / (S ** 2).sum()) if S.size else float("nan")
    return pc1, var_explained


def axis_alignment(loadings, genes, feature_names):
    """|cosine| between a loading vector and the signature indicator over the
    same features. PC sign is arbitrary, so the absolute value is taken.
    """
    gs = set(genes)
    ind = np.array([1.0 if g in gs else 0.0 for g in feature_names])
    den = np.linalg.norm(loadings) * np.linalg.norm(ind)
    return abs(float(loadings @ ind)) / den if den > 0 else float("nan")


def _frac_once(adata_resid, genes, group_key, control_key, feature_index):
    """One pass: residual convergence, regress signature out, recompute. Returns
    (frac_explained, base_summary, after_summary, resid2, score, present)."""
    is_ctrl = adata_resid.obs[control_key].to_numpy().astype(bool)
    base = convergence_summary(adata_resid, group_key=group_key,
                               control_key=control_key, feature_index=feature_index)
    score, present = score_signature(adata_resid, genes)
    score = score - score[is_ctrl].mean()          # controls -> 0, untouched
    resid2 = regress_out_score(adata_resid, score)
    after = convergence_summary(resid2, group_key=group_key,
                                control_key=control_key, feature_index=feature_index)
    b0, b1 = base["between_above_null"], after["between_above_null"]
    frac = (b0 - b1) / b0 if (b0 == b0 and b0 != 0) else float("nan")
    return frac, base, after, resid2, score, present


def attribute_to_signature(adata_resid, genes, group_key="annotation_group",
                           control_key="is_control", feature_index=None,
                           label="ISR", n_boot=100, seed=0):
    """Attribute residual cross-annotation convergence to a signature axis.

    `adata_resid` should be the post-viability residual (from Phase 1). If
    `feature_index` is None it is selected inside `convergence_summary`; pass the
    Phase 1 feature index to keep the comparison on identical genes.
    """
    frac, base, after, resid2, score, present = _frac_once(
        adata_resid, genes, group_key, control_key, feature_index)

    # dominant convergent axis (perturbations only) and its signature alignment
    is_ctrl = adata_resid.obs[control_key].to_numpy().astype(bool)
    fidx = base["feature_index"]
    X = np.asarray(adata_resid.X, dtype=float)[~is_ctrl][:, fidx]
    pc1, ve = dominant_axis(X)
    feat_names = [str(adata_resid.var_names[i]) for i in fidx]
    align = axis_alignment(pc1, genes, feat_names)

    # bootstrap CI for frac_explained over perturbations
    ci = (float("nan"), float("nan"))
    if n_boot and n_boot > 0:
        rng = np.random.default_rng(seed)
        pert_idx = np.where(~is_ctrl)[0]
        ctrl_idx = np.where(is_ctrl)[0]
        boots = []
        for _ in range(int(n_boot)):
            samp = np.concatenate([rng.choice(pert_idx, size=pert_idx.size, replace=True),
                                   ctrl_idx])
            sub = adata_resid[samp].copy()
            sub.obs_names_make_unique()
            f, *_ = _frac_once(sub, genes, group_key, control_key, fidx)
            if f == f:
                boots.append(f)
        if boots:
            ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))

    return {
        "label": label,
        "present": present,
        "n_signature_found": len(present),
        "base": base,                      # residual convergence (after viability)
        "after": after,                    # after also removing this signature
        "frac_explained": frac,
        "frac_ci": ci,
        "pc1_var_explained": ve,
        "pc1_align": align,
        "pc1": pc1,
        "feat_names": feat_names,
        "score": score,
        "residual2": resid2,
    }
