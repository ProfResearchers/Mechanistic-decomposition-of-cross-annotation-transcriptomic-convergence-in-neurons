"""
gene_sets.py
============
Curated gene sets used by the Phase 1 pipeline.

- The VIABILITY signature (cell-cycle proliferation + death/stress) is what H0
  regresses out: it captures the generic "sick / dividing cell" axis that makes
  functionally unrelated perturbations look alike for reasons that are NOT
  mechanism-specific.
- The ISR_ATF4 set is NOT regressed out in H0. It is included so the residual
  convergence (the part that survives viability regression) can be interpreted,
  and it is the anchor for H1 in a later phase.

Cell-cycle lists are the widely used Tirosh et al. (2016) / Seurat `cc.genes`
S-phase and G2M lists. A handful of symbols use legacy names (e.g. MLF1IP=CENPU,
FAM64A=PIMREG, HN1=JPT1); update aliases when applying to real data with a
symbol-mapping step. For the synthetic validation the literal symbols are used.
"""

CELL_CYCLE_S = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]

CELL_CYCLE_G2M = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# Generic apoptosis / cellular-stress markers (the "death" half of the
# viability axis). Kept distinct from the ISR set below where possible.
DEATH_STRESS = [
    "GADD45A", "GADD45B", "PMAIP1", "BBC3", "BAX", "BAK1", "CASP3", "CASP7",
    "CDKN1A", "MDM2", "TP53I3", "BTG2", "PLK3", "FAS", "TNFRSF10B", "GDF15",
]

# Integrated stress response / ATF4 program. NOT regressed in H0.
# Used to interpret the residual convergence and as the H1 anchor.
ISR_ATF4 = [
    "ATF4", "DDIT3", "TRIB3", "ASNS", "PSAT1", "PHGDH", "SLC7A11", "ATF3",
    "CEBPB", "CEBPG", "CHAC1", "DDIT4", "STC2", "VEGFA", "MTHFD2", "SHMT2",
    "ALDH18A1", "WARS", "SARS", "CTH",
]

# The viability signature regressed out in H0.
VIABILITY_SIGNATURE = CELL_CYCLE_S + CELL_CYCLE_G2M + DEATH_STRESS


def available(symbols, var_names):
    """Return the subset of `symbols` present in `var_names` (an iterable)."""
    s = set(var_names)
    return [g for g in symbols if g in s]
