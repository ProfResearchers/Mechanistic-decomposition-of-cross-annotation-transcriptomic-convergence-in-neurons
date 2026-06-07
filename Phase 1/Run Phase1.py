"""
run_phase1.py
=============
Phase 1 end-to-end: Aim 1 (quantify cross-annotation convergence) + H0 (remove
the generic viability axis and re-measure). Runs on synthetic data by default;
point --input at a real Replogle-format .h5ad to run the identical analysis.

Usage
-----
    python run_phase1.py                       # synthetic validation
    python run_phase1.py --input path.h5ad \\
        --group-key annotation_group --control-key is_control

Real-data note: the Replogle z-normalized pseudobulk file has perturbations as
obs and genes as var. You must supply (or precompute) two obs columns:
  * a boolean control flag (True for non-targeting controls), and
  * an annotation-group column. In Phase 1 a placeholder grouping is acceptable;
    the strict GO/CORUM/STRING orthogonality grouping is a later module.
"""

from __future__ import annotations
import argparse
import os
import sys
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
import gene_sets
from convergence import convergence_summary
from viability import score_signature, regress_out_score, viability_pair_similarity


def _fmt(x):
    return "nan" if x != x else f"{x:+.4f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=None, help="Path to a .h5ad; omit for synthetic.")
    ap.add_argument("--group-key", default="annotation_group")
    ap.add_argument("--control-key", default="is_control")
    ap.add_argument("--n-top-genes", type=int, default=2000)
    ap.add_argument("--outdir", default="../outputs")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- load ----
    if args.input:
        import anndata as ad
        adata = ad.read_h5ad(args.input)
        source = f"real file: {args.input}"
    else:
        from simulate import simulate
        adata = simulate(seed=args.seed)
        source = "synthetic (Replogle-schema) validation dataset"

    is_ctrl = adata.obs[args.control_key].to_numpy().astype(bool)
    groups = adata.obs[args.group_key].astype(str).to_numpy()

    print("=" * 72)
    print("PHASE 1  -  Aim 1 (convergence) + H0 (viability decomposition)")
    print("=" * 72)
    print(f"Source                : {source}")
    print(f"Observations          : {adata.n_obs}  "
          f"({int((~is_ctrl).sum())} perturbations, {int(is_ctrl.sum())} controls)")
    print(f"Genes                 : {adata.n_vars}")
    print(f"Annotation groups     : {sorted(set(groups[~is_ctrl]))}")

    # ---- Aim 1: convergence BEFORE ----
    before = convergence_summary(
        adata, group_key=args.group_key, control_key=args.control_key,
        n_top_genes=args.n_top_genes,
    )
    feat = before["feature_index"]  # reuse identical features after regression

    print("\n--- Aim 1: cross-annotation convergence (BEFORE) ---")
    print(f"Feature genes used    : {before['n_features']}")
    print(f"Null (ctrl-ctrl) mean : {_fmt(before['null_mean'])}  (sd {before['null_std']:.4f})")
    print(f"Within-annotation r   : {_fmt(before['within_mean'])}   "
          f"[above null {_fmt(before['within_above_null'])}, z {before['within_z']:.2f}]")
    print(f"CROSS-annotation r    : {_fmt(before['between_mean'])}   "
          f"[above null {_fmt(before['between_above_null'])}, z {before['between_z']:.2f}]")

    # ---- H0: viability signature ----
    score, present = score_signature(adata, gene_sets.VIABILITY_SIGNATURE)
    score = score - score[is_ctrl].mean()   # center on controls -> controls ~0
    print("\n--- H0: viability / generic-distress axis ---")
    print(f"Signature genes found : {len(present)} / {len(gene_sets.VIABILITY_SIGNATURE)}")

    # relationship test (BEFORE): do shared-viability pairs converge more?
    sim_between, _ = viability_pair_similarity(score, is_ctrl, groups)
    rho_b, p_b = spearmanr(sim_between, before["_between"])
    print(f"Spearman(viability-pair-similarity, cross-annotation r)  BEFORE : "
          f"rho={rho_b:+.3f}  (p={p_b:.1e})")

    # regress out and recompute on the SAME features
    adata_resid = regress_out_score(adata, score)
    after = convergence_summary(
        adata_resid, group_key=args.group_key, control_key=args.control_key,
        feature_index=feat,
    )
    score_after, _ = score_signature(adata_resid, present)
    sim_after, _ = viability_pair_similarity(score_after, is_ctrl, groups)
    rho_a, p_a = spearmanr(sim_after, after["_between"])

    print("\n--- Aim 1: cross-annotation convergence (AFTER regressing viability) ---")
    print(f"Null (ctrl-ctrl) mean : {_fmt(after['null_mean'])}  (sd {after['null_std']:.4f})")
    print(f"CROSS-annotation r    : {_fmt(after['between_mean'])}   "
          f"[above null {_fmt(after['between_above_null'])}, z {after['between_z']:.2f}]")
    print(f"Spearman(viability-pair-similarity, cross-annotation r)  AFTER  : "
          f"rho={rho_a:+.3f}  (p={p_a:.1e})")

    # ---- attribution ----
    b0 = before["between_above_null"]
    b1 = after["between_above_null"]
    frac = (b0 - b1) / b0 if (b0 == b0 and b0 != 0) else float("nan")
    print("\n--- Decomposition headline ---")
    print(f"Cross-annotation convergence above null  BEFORE : {_fmt(b0)}")
    print(f"Cross-annotation convergence above null  AFTER  : {_fmt(b1)}")
    print(f"Fraction attributable to the viability axis     : "
          f"{frac*100:5.1f}%" if frac == frac else "  nan")
    print(f"Residual (mechanism-specific) convergence       : {_fmt(b1)} "
          f"({(1-frac)*100:.1f}% of original)" if frac == frac else "")

    # ---- synthetic-only ground-truth check ----
    if "_truth_isr" in adata.obs.columns:
        _validate_synthetic(adata, before, after, is_ctrl, groups)

    # ---- figure ----
    try:
        _figure(before, after, sim_between, args.outdir)
        print(f"\nFigure saved          : {os.path.join(args.outdir, 'phase1_convergence.png')}")
    except Exception as e:  # pragma: no cover
        print(f"\n[figure skipped: {e}]")

    # ---- save residual AnnData for downstream phases ----
    try:
        out_h5 = os.path.join(args.outdir, "residual_after_viability.h5ad")
        adata_resid.write_h5ad(out_h5)
        print(f"Residual AnnData saved: {out_h5}")
    except Exception as e:  # pragma: no cover
        print(f"[residual save skipped: {e}]")

    print("=" * 72)


def _validate_synthetic(adata, before, after, is_ctrl, groups):
    """Confirm the planted ISR convergence survives viability regression and the
    planted viability convergence is largely removed.

    Uses one-sided Mann-Whitney U tests on the distributions of cross-annotation
    pair correlations (the statistically correct comparison: the spread of the
    null is over individual pairs, so we test distributions, not a mean against
    3*SD).
    """
    from scipy.stats import mannwhitneyu

    truth_isr = adata.obs["_truth_isr"].to_numpy().astype(bool)[~is_ctrl]
    g = groups[~is_ctrl]

    def cross_pair_corrs(C, mask):
        idx = np.where(mask)[0]
        if idx.size < 2:
            return np.array([])
        sub = C[np.ix_(idx, idx)]
        n = sub.shape[0]
        iu, ju = np.triu_indices(n, k=1)
        cross = g[idx][iu] != g[idx][ju]
        return sub[iu[cross], ju[cross]]

    isr_b = cross_pair_corrs(before["_C"], truth_isr)
    isr_a = cross_pair_corrs(after["_C"], truth_isr)
    non_a = cross_pair_corrs(after["_C"], ~truth_isr)
    null_a = after["_null"]

    p_vs_null = mannwhitneyu(isr_a, null_a, alternative="greater").pvalue
    p_vs_non = mannwhitneyu(isr_a, non_a, alternative="greater").pvalue

    print("\n--- Synthetic ground-truth validation ---")
    print(f"Cross-annotation r among TRUE ISR perturbations   : "
          f"before {_fmt(float(np.mean(isr_b)))} -> after {_fmt(float(np.mean(isr_a)))}  (should stay > null)")
    print(f"Cross-annotation r among non-ISR perturbations    : "
          f"after  {_fmt(float(np.mean(non_a)))}  (should fall toward null {_fmt(after['null_mean'])})")
    print(f"ISR residual vs null  (Mann-Whitney, one-sided)   : p = {p_vs_null:.1e}")
    print(f"ISR residual vs non-ISR residual (one-sided)      : p = {p_vs_non:.1e}")
    ok = (p_vs_null < 1e-3) and (p_vs_non < 1e-3) and (np.mean(isr_a) > np.mean(non_a))
    print(f"PASS: genuine convergence isolated from distress  : {ok}")


def _figure(before, after, sim_between, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    bins = np.linspace(-0.6, 1.0, 60)
    ax = axes[0]
    ax.hist(before["_null"], bins=bins, alpha=.5, label="null (ctrl-ctrl)", color="grey", density=True)
    ax.hist(before["_within"], bins=bins, alpha=.5, label="within-annotation", color="#1f77b4", density=True)
    ax.hist(before["_between"], bins=bins, alpha=.5, label="cross-annotation", color="#d62728", density=True)
    ax.set_title("BEFORE viability regression")
    ax.set_xlabel("profile correlation"); ax.set_ylabel("density"); ax.legend(fontsize=8)

    ax = axes[1]
    ax.hist(after["_null"], bins=bins, alpha=.5, label="null (ctrl-ctrl)", color="grey", density=True)
    ax.hist(after["_between"], bins=bins, alpha=.5, label="cross-annotation", color="#d62728", density=True)
    ax.set_title("AFTER viability regression")
    ax.set_xlabel("profile correlation"); ax.legend(fontsize=8)

    ax = axes[2]
    ax.scatter(sim_between, before["_between"], s=6, alpha=.3)
    ax.set_title("Shared viability vs convergence (before)")
    ax.set_xlabel("viability-pair similarity"); ax.set_ylabel("cross-annotation r")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "phase1_convergence.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
