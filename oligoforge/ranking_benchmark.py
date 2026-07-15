"""Independent ranking benchmark utilities and machine-readable metrics."""
import hashlib, json, time

BENCHMARK_SCHEMA = "1.2.0"



def wilson_interval(successes, total, z=1.959963984540054):
    """Wilson score interval for a binomial proportion.

    This keeps tiny synthetic benchmark results from being reported as if 100%
    observed accuracy implied certainty about unseen cases.
    """
    successes = int(successes); total = int(total)
    if total <= 0 or successes < 0 or successes > total:
        return None
    p = successes / total
    den = 1.0 + (z * z) / total
    center = (p + (z * z) / (2.0 * total)) / den
    half = (z / den) * ((p * (1.0 - p) / total + (z * z) / (4.0 * total * total)) ** 0.5)
    return {"estimate": round(p, 6), "lower": round(max(0.0, center - half), 6),
            "upper": round(min(1.0, center + half), 6), "n": total,
            "method": "Wilson 95% interval"}

def corpus_manifest(cases, split="development"):
    canonical = json.dumps(cases, sort_keys=True, separators=(",", ":"), default=str)
    return dict(schema_version=BENCHMARK_SCHEMA, split=split, n_cases=len(cases),
                corpus_sha256=hashlib.sha256(canonical.encode()).hexdigest())


def topk(rankings, accepted_ids, ks=(1, 3, 5, 10)):
    ids = [x.get("id") for x in rankings]
    acc = set(accepted_ids or [])
    return {"top_%d" % k: int(any(x in acc for x in ids[:k])) for k in ks}


def diversity(candidates):
    if not candidates:
        return dict(n=0, regions=0, unique_triplets=0)
    ids = {(x.get("forward"), x.get("reverse"), x.get("probe")) for x in candidates}
    regs = {(x.get("f_span") or x.get("f_xy") or [0])[0] // 250 for x in candidates}
    return dict(n=len(candidates), regions=len(regs), unique_triplets=len(ids))


def compare(legacy, new, accepted_ids):
    return dict(legacy=topk(legacy, accepted_ids), new=topk(new, accepted_ids),
                legacy_diversity=diversity(legacy), new_diversity=diversity(new))



def validate_corpus(corpus):
    """Validate frozen benchmark structure and target-group split isolation.

    Synthetic fixtures may omit ``target_group``; in that case leakage checking is
    reported as not applicable rather than silently passing.  Biological or
    experimental fixtures must supply target_group to support leakage-safe claims.
    """
    errors, warnings = [], []
    cases = list((corpus or {}).get("cases") or [])
    splits = dict((corpus or {}).get("splits") or {})
    ids = [c.get("id") for c in cases]
    if any(not x for x in ids):
        errors.append("every case requires a non-empty id")
    dup = sorted({x for x in ids if ids.count(x) > 1})
    if dup:
        errors.append("duplicate case ids: %s" % ", ".join(dup))
    byid = {c.get("id"): c for c in cases if c.get("id")}
    assigned = {}
    for split, members in splits.items():
        for cid in members or []:
            if cid not in byid:
                errors.append("split %s references unknown case %s" % (split, cid))
            if cid in assigned:
                errors.append("case %s appears in both %s and %s" % (cid, assigned[cid], split))
            assigned[cid] = split
    missing = sorted(set(byid) - set(assigned))
    if missing:
        errors.append("cases missing from splits: %s" % ", ".join(missing))
    for case in cases:
        cids = {c.get("id") for c in case.get("candidates") or []}
        if not cids:
            errors.append("case %s has no candidates" % case.get("id"))
        expected = set(case.get("expected") or [])
        if not expected:
            errors.append("case %s has no accepted candidate" % case.get("id"))
        unknown = sorted(expected - cids)
        if unknown:
            errors.append("case %s expects unknown candidates: %s" % (case.get("id"), ", ".join(unknown)))

    group_splits = {}
    groups_present = 0
    biological_without_group = []
    for case in cases:
        group = case.get("target_group")
        cls = str(case.get("class") or "")
        if group:
            groups_present += 1
            split = assigned.get(case.get("id"))
            prior = group_splits.setdefault(group, split)
            if prior != split:
                errors.append("target group %s leaks across %s and %s" % (group, prior, split))
        elif cls not in {"synthetic_adversarial", "condition_perturbation", "transcript", "multiplex", "multi_isolate", "discrimination", "hydrolysis_probe"}:
            biological_without_group.append(case.get("id"))
    if groups_present == 0:
        warnings.append("target_group leakage check not applicable: no cases declare target_group")
    if biological_without_group:
        warnings.append("non-synthetic cases without target_group: %s" % ", ".join(biological_without_group))
    canonical = json.dumps(corpus, sort_keys=True, separators=(",", ":"), default=str)
    return {
        "schema_version": BENCHMARK_SCHEMA,
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "n_cases": len(cases),
        "n_splits": len(splits),
        "n_target_groups": len(group_splits),
        "target_group_leakage": any("leaks across" in x for x in errors),
        "corpus_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
    }
