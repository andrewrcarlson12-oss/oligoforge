"""Objective profiles for evidence-structured assay ranking.

The profiles deliberately separate hard requirements from ranking preferences.  A
user may select a legitimate design objective, but cannot use a weight to make a
known false-positive product or broken assay geometry disappear.
"""
from copy import deepcopy

RANKER_VERSION = "2.2.0"
PROFILE_VERSION = "2026-07-ranking-truth-3"

_BASE = dict(
    label="Balanced hydrolysis-probe assay",
    description="Balance target recovery, exclusivity, robustness, and practical qPCR behavior.",
    min_target_coverage=1.0,
    min_probe_coverage=1.0,
    require_no_signal_offtargets=True,
    require_no_product_offtargets=False,
    require_junction=False,
    prefer_short_amplicon=110,
    max_degeneracy_fold=64,
    priority=("offtarget", "coverage", "robustness", "triplet", "practical"),
)

OBJECTIVE_PROFILES = {
    "balanced": _BASE,
    "single_target_detection": dict(_BASE,
        label="Single-target detection",
        description="Prioritize clean amplification of one declared target with no signal-generating off-target product.",
        min_target_coverage=1.0,
        min_probe_coverage=1.0,
        priority=("offtarget", "coverage", "triplet", "robustness", "practical")),
    "broad_inclusivity": dict(_BASE,
        label="Broad inclusivity",
        description="Prioritize coherent amplification and probe binding across the declared target set.",
        min_target_coverage=0.95,
        min_probe_coverage=0.95,
        priority=("coverage", "offtarget", "robustness", "triplet", "practical")),
    "discrimination": dict(_BASE,
        label="Species/strain discrimination",
        description="Prioritize exclusion of the supplied near-neighbor set while retaining declared target coverage.",
        min_target_coverage=0.90,
        min_probe_coverage=0.90,
        require_no_product_offtargets=True,
        priority=("offtarget", "coverage", "robustness", "triplet", "practical")),
    "confirmatory": dict(_BASE,
        label="Confirmatory high-exclusivity assay",
        description="Use the strictest off-target hierarchy; a predicted off-target product is disqualifying.",
        min_target_coverage=0.90,
        min_probe_coverage=0.90,
        require_no_product_offtargets=True,
        priority=("offtarget", "coverage", "robustness", "triplet", "practical")),
    "screening": dict(_BASE,
        label="Screening assay",
        description="Favor inclusive target recovery while still forbidding signal-generating false-positive products.",
        min_target_coverage=0.95,
        min_probe_coverage=0.90,
        priority=("coverage", "offtarget", "triplet", "robustness", "practical")),
    "transcript_specific": dict(_BASE,
        label="Transcript-specific RT-qPCR",
        description="Require a declared exon-junction relationship and otherwise balance coverage and exclusivity.",
        require_junction=True,
        priority=("junction", "offtarget", "coverage", "robustness", "triplet", "practical")),
    "degraded_template": dict(_BASE,
        label="Short degraded-template assay",
        description="Prefer short products after validity, target recovery, and exclusivity are satisfied.",
        prefer_short_amplicon=85,
        priority=("offtarget", "coverage", "practical", "robustness", "triplet")),
    "multiplex": dict(_BASE,
        label="Multiplex panel member",
        description="Prioritize panel compatibility after target recovery and off-target validity.",
        priority=("offtarget", "multiplex", "coverage", "robustness", "triplet", "practical")),
    "sybr": dict(_BASE,
        label="SYBR assay",
        description="Primer-pair assay with stricter dimer and product-specificity emphasis.",
        min_probe_coverage=0.0,
        require_no_product_offtargets=True,
        priority=("offtarget", "coverage", "triplet", "robustness", "practical")),
}


def resolve_objective(name=None, *, no_probe=False):
    """Resolve the scientifically appropriate objective key.

    ``balanced`` historically meant the hydrolysis-probe baseline.  Passing it
    explicitly from a generic UI accidentally disabled SYBR's stricter product-
    specificity rule.  Probe-less chemistry now resolves that generic default to
    ``sybr`` everywhere; an explicitly chosen non-default objective is preserved.
    """
    key = str(name or ("sybr" if no_probe else "balanced")).strip().lower()
    if no_probe and key == "balanced":
        key = "sybr"
    if key not in OBJECTIVE_PROFILES:
        key = "sybr" if no_probe else "balanced"
    return key


def get_profile(name=None, *, no_probe=False, overrides=None):
    """Return a defensive copy of one ranking objective.

    ``overrides`` is limited to explicit requirements/preferences; callers should
    validate user values before passing them here.
    """
    key = resolve_objective(name, no_probe=no_probe)
    out = deepcopy(OBJECTIVE_PROFILES[key])
    out["key"] = key
    out["ranker_version"] = RANKER_VERSION
    out["profile_version"] = PROFILE_VERSION
    if overrides:
        allowed = {"min_target_coverage", "min_probe_coverage", "require_no_signal_offtargets",
                   "require_no_product_offtargets", "require_junction", "prefer_short_amplicon",
                   "max_degeneracy_fold"}
        for k, v in overrides.items():
            if k in allowed and v is not None:
                out[k] = v
    return out


def public_profiles():
    return {k: {x: v for x, v in p.items() if x not in {"priority"}}
            for k, p in OBJECTIVE_PROFILES.items()}
