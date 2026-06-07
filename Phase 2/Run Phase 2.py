"""
run_phase2.py
=============
Phase 2 end-to-end: take the post-viability residual (H0) and attribute the
surviving cross-annotation convergence to the ISR/ATF4 axis (H1).

    python run_phase2.py                       # synthetic validation
    python run_phase2.py --input path.h5ad \\
        --group-key annotation_group --control-key is_control
"""

from __future__ import annotations
import argparse
import os
import sys
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

sys.path.insert(0, os.path.dirname(__file__))
import gene_sets
from convergence import convergence_summary
from viability import score_signature, regress_out_score, viability_pair_similarity
from attribution import attribute_to_signature


def _fmt(x):
    return "nan" if x != x else f"{x:+.4f}"


def get_residual(adata, group_key, control_key, n_top_genes):
    """H0 step: regress the viability axis out, return (before, after_v, residual, features)."""
    is_ctrl = adata.obs[control_key].to_numpy().astype(bool)
    before = convergence_summary(adata, group_key=group_key, control_key=control_key,
                                 n_top_genes=n_top_genes)
    feat = before["feature_index"]
    score, _ = score_signature(adata, gene_sets.VIABILITY_SIGNATURE)
    score = score - score[is_ctrl].mean()
    resid = regress_out_score(adata, score)
    after_v = convergence_summary(resid, group_key=group_key, control_key=control_key,
                                  feature_index=feat)
    return before, after_v, resid, feat


def _per_pert_cross_mean(C, groups):
    """Mean cross-annotation correlation for each perturbation (row)."""
    n = C.shape[0]
    out = np.full(n, np.nan)
    for i in range(n):
        diff = groups != groups[i]
        diff[i] = False
        if diff.any():
            out[i] = C[i, diff].mean()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None)
    ap.add_argument("--group-key", default="annotation_group")
    ap.add_argument("--control-key", default="is_control")
    ap.add_argument("--n-top-genes", type=int, default=2000)
    ap.add_argument("--n-boot", type=int, default=100)
    ap.add_argument("--outdir", default="../outputs")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.input:
        import anndata as ad
        adata = ad.read_h5ad(args.input)
        source = f"real file: {args.input}"
    else:
        from simulate import simulate
        adata = simulate(seed=args.seed)
        source = "synthetic (Replogle-schema) validation dataset"

    gk, ck = args.group_key, args.control_key
    is_ctrl = adata.obs[ck].to_numpy().astype(bool)
    groups = adata.obs[gk].astype(str).to_numpy()

    print("=" * 72)
    print("PHASE 2  -  Aim 2 / H1: attribute residual convergence to the ISR axis")
    print("=" * 72)
    print(f"Source                : {source}")

    # ---- H0 recap: get the residual ----
    before, after_v, resid, feat = get_residual(adata, gk, ck, args.n_top_genes)
    tot = before["between_above_null"]
    res_v = after_v["between_above_null"]
    frac_v = (tot - res_v) / tot if (tot == tot and tot != 0) else float("nan")
    print("\n--- H0 recap (from Phase 1) ---")
    print(f"Total cross-annotation convergence above null   : {_fmt(tot)} (z {before['between_z']:.2f})")
    print(f"Residual after removing viability               : {_fmt(res_v)}")
    if frac_v == frac_v:
        print(f"-> viability explains {frac_v*100:.1f}% of total convergence")

    # ---- H1: ISR attribution of the residual ----
    R = attribute_to_signature(resid, gene_sets.ISR_ATF4, group_key=gk, control_key=ck,
                               feature_index=feat, label="ISR", n_boot=args.n_boot,
                               seed=args.seed)
    print("\n--- H1: ISR/ATF4 attribution of the residual ---")
    print(f"ISR signature genes found                       : {R['n_signature_found']} / {len(gene_sets.ISR_ATF4)}")
    print(f"Residual convergence (after viability)          : {_fmt(R['base']['between_above_null'])}")
    print(f"Residual convergence after ALSO removing ISR    : {_fmt(R['after']['between_above_null'])}")
    ci = R["frac_ci"]
    print(f"-> ISR explains {R['frac_explained']*100:.1f}% of the residual "
          f"(95% CI {ci[0]*100:.1f}-{ci[1]*100:.1f}%)")
    print("   (lower bound: diluted by non-engaging perturbations + marker-vs-program;")
    print("    the axis alignment and among-engagers collapse below are the primary evidence)")
    print(f"Dominant convergent axis (top PC)               : {R['pc1_var_explained']*100:.1f}% of residual variance")
    print(f"  alignment of that axis with the ISR signature : {R['pc1_align']:.3f}  (|cosine|)")

    # relationship: do high-ISR perturbations converge more? (before ISR removal)
    score_isr, _ = score_signature(resid, gene_sets.ISR_ATF4)
    score_isr = score_isr - score_isr[is_ctrl].mean()
    sim_isr, _ = viability_pair_similarity(score_isr, is_ctrl, groups)
    rho, p = spearmanr(sim_isr, R["base"]["_between"])
    print(f"Spearman(ISR-pair-similarity, residual cross-r) : rho={rho:+.3f} (p={p:.1e})")

    # ---- full decomposition chain ----
    print("\n--- Decomposition chain (cross-annotation convergence) ---")
    if frac_v == frac_v and R["frac_explained"] == R["frac_explained"]:
        isr_of_total = (1 - frac_v) * R["frac_explained"]
        remaining = (1 - frac_v) * (1 - R["frac_explained"])
        print(f"  viability (H0)            : {frac_v*100:5.1f}% of total")
        print(f"  ISR/ATF4 (H1)             : {isr_of_total*100:5.1f}% of total")
        print(f"  unexplained residual      : {remaining*100:5.1f}% of total  (-> H2/H3/H4/H5)")

    # ---- synthetic ground-truth validation ----
    if "_truth_isr" in adata.obs.columns:
        _validate(adata, resid, R, is_ctrl, groups, score_isr)

    # ---- figure ----
    try:
        _figure(R, resid, groups, is_ctrl, score_isr, args.outdir)
        print(f"\nFigure saved          : {os.path.join(args.outdir, 'phase2_isr_attribution.png')}")
    except Exception as e:  # pragma: no cover
        print(f"\n[figure skipped: {e}]")

    try:
        out_h5 = os.path.join(args.outdir, "residual_after_viability_and_isr.h5ad")
        R["residual2"].write_h5ad(out_h5)
        print(f"Residual AnnData saved: {out_h5}")
    except Exception as e:  # pragma: no cover
        print(f"[residual save skipped: {e}]")

    print("=" * 72)


def _validate(adata, resid, R, is_ctrl, groups, score_isr):
    """Robust ground-truth checks that the residual convergence IS the ISR axis.
    A single global 'fraction explained' is a lower bound (cross-annotation
    convergence is diluted by perturbations that do not engage ISR, and a marker
    set only partially proxies a broad program), so PASS rests on three robust
    signals instead: loading separation, the among-engagers collapse, and the
    dominant-axis alignment being above chance."""
    truth = adata.obs["_truth_isr"].to_numpy().astype(bool)[~is_ctrl]
    g = groups[~is_ctrl]
    loadings = score_isr[~is_ctrl]

    # 1. ISR loading higher in true-ISR carriers
    p_load = mannwhitneyu(loadings[truth], loadings[~truth], alternative="greater").pvalue

    # 2. among ISR-engaging perturbations, removing ISR collapses their convergence to ~null
    def cross_mean(C, mask):
        idx = np.where(mask)[0]
        sub = C[np.ix_(idx, idx)]
        n = sub.shape[0]
        iu, ju = np.triu_indices(n, k=1)
        cross = g[idx][iu] != g[idx][ju]
        return float(np.mean(sub[iu[cross], ju[cross]])) if cross.any() else float("nan")

    before_isr = cross_mean(R["base"]["_C"], truth)
    after_isr = cross_mean(R["after"]["_C"], truth)
    null_m, null_s = R["after"]["null_mean"], R["after"]["null_std"]
    collapse = (after_isr < before_isr) and (after_isr <= null_m + 3 * null_s)

    # 3. the dominant convergent axis aligns with ISR above chance
    feat_names = R["feat_names"]
    nsig = sum(x in set(gene_sets.ISR_ATF4) for x in feat_names)
    chance = (nsig / len(feat_names)) ** 0.5 if feat_names else float("nan")
    align_ok = R["pc1_align"] > 2 * chance

    print("\n--- Synthetic ground-truth validation ---")
    print(f"ISR loading: true-ISR vs others (MWU one-sided) : p = {p_load:.1e}")
    print(f"Cross-r among TRUE-ISR perturbations            : {_fmt(before_isr)} -> {_fmt(after_isr)} after ISR removal")
    print(f"   collapses to ~null ({_fmt(null_m)} +/- {null_s:.3f})    : {collapse}")
    print(f"Dominant-axis alignment {R['pc1_align']:.2f} vs chance ~{chance:.2f}  : {'above chance' if align_ok else 'NOT above chance'}")
    ok = (p_load < 1e-3) and collapse and align_ok
    print(f"PASS: residual convergence identified as ISR    : {ok}")


def _figure(R, resid, groups, is_ctrl, score_isr, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    bins = np.linspace(-0.6, 1.0, 60)

    ax = axes[0]
    ax.hist(R["base"]["_null"], bins=bins, alpha=.5, color="grey", density=True, label="null")
    ax.hist(R["base"]["_between"], bins=bins, alpha=.5, color="#d62728", density=True,
            label="residual (after viability)")
    ax.hist(R["after"]["_between"], bins=bins, alpha=.5, color="#2ca02c", density=True,
            label="after also removing ISR")
    ax.set_title("Residual convergence, ISR removed")
    ax.set_xlabel("cross-annotation r"); ax.set_ylabel("density"); ax.legend(fontsize=8)

    ax = axes[1]
    ppm = _per_pert_cross_mean(R["base"]["_C"], groups[~is_ctrl])
    load = score_isr[~is_ctrl]
    truth = resid.obs["_truth_isr"].to_numpy().astype(bool)[~is_ctrl] if "_truth_isr" in resid.obs.columns else None
    if truth is not None:
        ax.scatter(load[~truth], ppm[~truth], s=8, alpha=.4, color="grey", label="non-ISR (truth)")
        ax.scatter(load[truth], ppm[truth], s=8, alpha=.5, color="#d62728", label="ISR (truth)")
        ax.legend(fontsize=8)
    else:
        ax.scatter(load, ppm, s=8, alpha=.4)
    ax.set_title("ISR loading vs convergence")
    ax.set_xlabel("per-perturbation ISR loading"); ax.set_ylabel("mean cross-annotation r")

    ax = axes[2]
    pc1 = R["pc1"]; feat_names = R["feat_names"]
    isr_set = set(gene_sets.ISR_ATF4)
    is_isr = np.array([g in isr_set for g in feat_names])
    x = np.arange(len(pc1))
    ax.scatter(x[~is_isr], pc1[~is_isr], s=4, alpha=.3, color="grey", label="other genes")
    ax.scatter(x[is_isr], pc1[is_isr], s=18, alpha=.9, color="#d62728", label="ISR genes")
    ax.set_title(f"Dominant convergent axis (align {R['pc1_align']:.2f})")
    ax.set_xlabel("feature gene index"); ax.set_ylabel("top-PC loading"); ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "phase2_isr_attribution.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
