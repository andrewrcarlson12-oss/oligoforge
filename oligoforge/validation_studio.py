"""Deterministic, conservative validation-experiment planning for OligoForge.

Validation Studio turns *differences between computational assay predictions* into a
small, auditable comparison experiment.  It does not claim that a plate plan is
globally optimal and it never converts modeled amplification into analytical or
clinical performance.  The same isolate-level product reconstruction used elsewhere
in OligoForge is authoritative here.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import random
import re
from collections import defaultdict

from . import isolates as ISO
from .provenance import sha256_value


VALIDATION_STUDIO_VERSION = "1.0.0"
SCHEMA_VERSION = "oligoforge-validation-plan/v1"
MAX_CANDIDATES = 16
MAX_CASES = 96
MAX_SEQUENCE = 120_000
MAX_TOTAL_BASES = 2_000_000


def _clean_seq(value):
    return re.sub(r"[^ACGTRYSWKMBDHVN]", "", str(value or "").upper().replace("U", "T"))


def _candidate_identity(row):
    assay = row.get("assay") if isinstance(row.get("assay"), dict) else row
    chemistry = str(row.get("chemistry") or row.get("profile") or "unspecified").strip()
    payload = {
        "forward": _clean_seq(assay.get("forward")),
        "reverse": _clean_seq(assay.get("reverse")),
        "probe": _clean_seq(assay.get("probe")),
        "chemistry": chemistry,
    }
    return payload, sha256_value(payload)[:16]


def normalize_candidates(candidates):
    if not isinstance(candidates, list) or len(candidates) < 2:
        raise ValueError("at least two candidate assays are required")
    if len(candidates) > MAX_CANDIDATES:
        raise ValueError("candidate count exceeds %d" % MAX_CANDIDATES)
    out, seen = [], set()
    for i, raw in enumerate(candidates):
        if not isinstance(raw, dict):
            raise ValueError("candidate %d must be an object" % (i + 1))
        ident, digest = _candidate_identity(raw)
        if not ident["forward"] or not ident["reverse"]:
            raise ValueError("candidate %d requires forward and reverse primers" % (i + 1))
        if digest in seen:
            continue
        seen.add(digest)
        cid = str(raw.get("candidate_id") or raw.get("id") or ("candidate-%d" % (i + 1))).strip()
        if not cid:
            cid = "candidate-%d" % (i + 1)
        if any(x["candidate_id"] == cid for x in out):
            raise ValueError("candidate_id values must be unique")
        out.append({
            "candidate_id": cid,
            "name": str(raw.get("name") or cid),
            "identity_sha256": digest,
            "assay": {"forward": ident["forward"], "reverse": ident["reverse"],
                      "probe": ident["probe"] or None},
            "chemistry": ident["chemistry"],
            "source_rank": raw.get("rank") or raw.get("candidate_rank"),
            "rank_explanation": raw.get("rank_explanation"),
        })
    if len(out) < 2:
        raise ValueError("at least two distinct candidate assay/chemistry identities are required")
    return out


def normalize_cases(records):
    if not isinstance(records, list) or not records:
        raise ValueError("at least one target or off-target case is required")
    if len(records) > MAX_CASES:
        raise ValueError("case count exceeds %d" % MAX_CASES)
    out, seen, total = [], set(), 0
    for i, raw in enumerate(records):
        if not isinstance(raw, dict):
            raise ValueError("case %d must be an object" % (i + 1))
        rid = str(raw.get("case_id") or raw.get("record_id") or raw.get("id") or ("case-%d" % (i + 1))).strip()
        if not rid or rid in seen:
            raise ValueError("case identifiers must be non-empty and unique")
        seen.add(rid)
        seq = _clean_seq(raw.get("sequence"))
        if not seq:
            raise ValueError("case %s has no usable sequence" % rid)
        if len(seq) > MAX_SEQUENCE:
            raise ValueError("case %s exceeds %d nt" % (rid, MAX_SEQUENCE))
        total += len(seq)
        if total > MAX_TOTAL_BASES:
            raise ValueError("case panel exceeds %d total bases" % MAX_TOTAL_BASES)
        role = str(raw.get("role") or raw.get("kind") or "target").lower().replace("-", "_")
        if role in {"near_neighbor", "near_neighbour", "offtarget"}:
            role = "off_target"
        if role not in {"target", "off_target"}:
            raise ValueError("case %s role must be target or off_target" % rid)
        source_type = str(raw.get("source_type") or "natural").lower()
        if source_type not in {"natural", "synthetic", "computational"}:
            raise ValueError("case %s source_type is invalid" % rid)
        out.append({
            "case_id": rid,
            "name": str(raw.get("name") or rid),
            "role": role,
            "group": str(raw.get("group") or "ungrouped"),
            "source_type": source_type,
            "sequence": seq,
            "sequence_sha256": hashlib.sha256(seq.encode("ascii")).hexdigest(),
            "metadata": raw.get("metadata") if isinstance(raw.get("metadata"), dict) else {},
        })
    return out


def _prediction_state(result, has_probe):
    if not result.get("amplifies"):
        reason = str(result.get("reason") or "no coherent product").lower()
        subtype = "terminal-primer concern" if "3\u2032" in reason or "mismatch" in reason else "no coherent product"
        return ("no_product", subtype)
    if has_probe and not result.get("probe_binds"):
        return ("product_without_probe", "probe recognition concern")
    return ("signal_product" if has_probe else "product", "modeled product")


def _pair_key(a, b):
    return tuple(sorted((a, b)))


def evaluate_disagreements(candidates, cases, model=None):
    model = dict(model or {})
    kwargs = {
        "max_mm": max(0, min(int(model.get("max_mm", 5)), 10)),
        "clamp_n": max(1, min(int(model.get("clamp_n", 2)), 8)),
        "min_product": max(20, int(model.get("min_product", 40))),
        "max_product": min(20_000, max(40, int(model.get("max_product", 3000)))),
        "min_probe_ident": max(0.0, min(float(model.get("min_probe_ident", 85.0)), 100.0)),
    }
    rows = []
    for case in cases:
        predictions, states = {}, {}
        for cand in candidates:
            assay = cand["assay"]
            res = ISO.amplify(assay["forward"], assay["reverse"], assay.get("probe") or "",
                              case["sequence"], **kwargs)
            predictions[cand["candidate_id"]] = {
                "amplifies": bool(res.get("amplifies")),
                "product": res.get("product"),
                "n_products": res.get("n_products"),
                "forward_identity": res.get("f_ident"),
                "reverse_identity": res.get("r_ident"),
                "probe_identity": res.get("probe_ident"),
                "probe_binds": res.get("probe_binds"),
                "reason": res.get("reason"),
            }
            states[cand["candidate_id"]] = _prediction_state(res, bool(assay.get("probe")))
        distinguishes = []
        for i, left in enumerate(candidates):
            for right in candidates[i + 1:]:
                a, b = left["candidate_id"], right["candidate_id"]
                if states[a] != states[b]:
                    distinguishes.append(list(_pair_key(a, b)))
        state_count = len(set(states.values()))
        rows.append({
            "case_id": case["case_id"], "name": case["name"], "role": case["role"],
            "group": case["group"], "source_type": case["source_type"],
            "sequence_sha256": case["sequence_sha256"],
            "predictions": predictions,
            "modeled_states": {k: {"state": v[0], "basis": v[1]} for k, v in states.items()},
            "distinguishes_candidate_pairs": distinguishes,
            "n_distinct_modeled_states": state_count,
            "informative": bool(distinguishes),
        })
    return rows, kwargs


def select_informative_cases(rows, max_cases=12):
    max_cases = max(1, min(int(max_cases), 48))
    informative = [r for r in rows if r["informative"]]
    selected, covered_pairs, covered_groups = [], set(), set()
    while informative and len(selected) < max_cases:
        def score(row):
            pairs = {_pair_key(*p) for p in row["distinguishes_candidate_pairs"]}
            new_pairs = len(pairs - covered_pairs)
            group_key = (row["role"], row["group"])
            novelty = 1 if group_key not in covered_groups else 0
            priority = 1 if row["role"] == "off_target" else 0
            return (-new_pairs, -row["n_distinct_modeled_states"], -novelty, -priority, row["case_id"])
        informative.sort(key=score)
        best = informative.pop(0)
        selected.append(best)
        covered_pairs.update(_pair_key(*p) for p in best["distinguishes_candidate_pairs"])
        covered_groups.add((best["role"], best["group"]))
    for row in selected:
        bases = []
        if row["role"] == "target":
            bases.append("tests candidate-specific modeled target support")
        else:
            bases.append("tests candidate-specific modeled near-neighbor rejection")
        if any(x["state"] == "product_without_probe" for x in row["modeled_states"].values()):
            bases.append("includes a probe-recognition disagreement")
        if any("terminal" in x["basis"] for x in row["modeled_states"].values()):
            bases.append("includes a terminal-primer mismatch concern")
        row["selection_rationale"] = "; ".join(bases)
    return selected


def _well_names(fmt):
    if fmt == 96:
        rows, cols = "ABCDEFGH", 12
    elif fmt == 384:
        rows, cols = "ABCDEFGHIJKLMNOP", 24
    else:
        raise ValueError("plate_format must be 96 or 384")
    # Column-major order naturally interleaves candidates across rows and avoids
    # placing one candidate in a contiguous corner of the plate.
    return [r + str(c) for c in range(1, cols + 1) for r in rows]


def _is_edge(well, fmt):
    rows = "ABCDEFGH" if fmt == 96 else "ABCDEFGHIJKLMNOP"
    max_col = 12 if fmt == 96 else 24
    return well[0] in {rows[0], rows[-1]} or int(well[1:]) in {1, max_col}


def _safe_cell(value):
    text = "" if value is None else str(value)
    return "'" + text if text[:1] in {"=", "+", "-", "@"} else text


def build_plate_layout(candidates, selected_cases, plate_format=96, replicates=3,
                       controls=None, seed=0, use_edge_wells=True):
    plate_format = int(plate_format)
    replicates = max(1, min(int(replicates), 12))
    controls = dict(controls or {})
    entries = []
    # Interleave candidate within replicate and case.  Deterministic per-block
    # rotation prevents any candidate from always occupying the same row.
    for rep in range(1, replicates + 1):
        for ci, case in enumerate(selected_cases):
            rotated = candidates[(ci + rep - 1) % len(candidates):] + candidates[:(ci + rep - 1) % len(candidates)]
            for cand in rotated:
                entries.append({"well_type": "test", "candidate_id": cand["candidate_id"],
                                "case_id": case["case_id"], "case_role": case["role"],
                                "replicate": rep, "expected": case["modeled_states"][cand["candidate_id"]]["state"]})
    ntc_reps = max(0, min(int(controls.get("ntc_replicates", 1)), 12))
    for rep in range(1, ntc_reps + 1):
        for cand in candidates:
            entries.append({"well_type": "no_template_control", "candidate_id": cand["candidate_id"],
                            "case_id": "NTC", "case_role": "control", "replicate": rep,
                            "expected": "not_amplified"})
    for control in controls.get("positive_controls", []) or []:
        if not isinstance(control, dict):
            continue
        for cand in candidates:
            entries.append({"well_type": "positive_control", "candidate_id": cand["candidate_id"],
                            "case_id": str(control.get("control_id") or "positive-control"),
                            "case_role": "control", "replicate": 1, "expected": "amplified"})
    for control in controls.get("extraction_controls", []) or []:
        if not isinstance(control, dict):
            continue
        entries.append({"well_type": "extraction_control", "candidate_id": str(control.get("candidate_id") or "panel"),
                        "case_id": str(control.get("control_id") or "extraction-control"),
                        "case_role": "control", "replicate": 1,
                        "expected": str(control.get("expected") or "amplified")})
    wells = _well_names(plate_format)
    if not use_edge_wells:
        wells = [w for w in wells if not _is_edge(w, plate_format)]
    if len(entries) > len(wells):
        raise ValueError("planned wells (%d) exceed usable %d-well capacity (%d)" %
                         (len(entries), plate_format, len(wells)))
    rng = random.Random(int(seed))
    # Shuffle case-block ordering while retaining candidate interleaving inside
    # each block.  The recorded seed makes the result exactly reproducible.
    block = len(candidates)
    test_n = len(selected_cases) * replicates * block
    blocks = [entries[i:i + block] for i in range(0, test_n, block)]
    rng.shuffle(blocks)
    ordered = [x for b in blocks for x in b] + entries[test_n:]
    placed = []
    for well, entry in zip(wells, ordered):
        row = dict(entry)
        row["well"] = well
        row["edge_well"] = _is_edge(well, plate_format)
        placed.append(row)
    edge_tests = sum(1 for x in placed if x["edge_well"] and x["well_type"] == "test")
    warnings = []
    if edge_tests:
        warnings.append("%d test wells are on a plate edge; inspect edge effects before interpretation" % edge_tests)
    return {
        "schema_version": "oligoforge-plate-layout/v1", "plate_format": plate_format,
        "replicates": replicates, "randomization_seed": int(seed),
        "candidate_interleaved": True, "use_edge_wells": bool(use_edge_wells),
        "n_wells": len(placed), "wells": placed, "warnings": warnings,
    }


def plate_csv(layout):
    fields = ["well", "well_type", "candidate_id", "case_id", "case_role", "replicate",
              "expected", "observed", "cq", "efficiency", "melt_abnormal", "unexpected_product",
              "probe_signal", "notes"]
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in layout.get("wells", []):
        writer.writerow({k: _safe_cell(row.get(k)) for k in fields})
    return buf.getvalue()


def _canonical_plan_payload(plan):
    return {k: v for k, v in plan.items() if k not in {"plan_sha256"}}


def create_plan(candidates, cases, objective="candidate comparison", reaction_conditions=None,
                plate_format=96, replicates=3, controls=None, acceptance_criteria=None,
                model=None, max_cases=12, seed=0, use_edge_wells=True, existing_evidence=None):
    cands = normalize_candidates(candidates)
    recs = normalize_cases(cases)
    rows, declared_model = evaluate_disagreements(cands, recs, model)
    selected = select_informative_cases(rows, max_cases=max_cases)
    if not selected:
        raise ValueError("no supplied case distinguishes the candidate assays under the declared model")
    layout = build_plate_layout(cands, selected, plate_format=plate_format, replicates=replicates,
                                controls=controls, seed=seed, use_edge_wells=use_edge_wells)
    plan = {
        "schema_version": SCHEMA_VERSION, "validation_studio_version": VALIDATION_STUDIO_VERSION,
        "objective": str(objective or "candidate comparison"),
        "scope_statement": "Bounded candidate-disagreement experiment; laboratory confirmation required.",
        "candidates": cands,
        "case_inventory": [{k: v for k, v in x.items() if k != "sequence"} for x in recs],
        "evaluations": rows, "selected_cases": selected,
        "selection_status": {"method": "deterministic greedy pair-coverage with group diversity",
                             "globally_optimal": False, "n_evaluated": len(rows),
                             "n_informative": sum(1 for x in rows if x["informative"]),
                             "n_selected": len(selected)},
        "reaction_conditions": dict(reaction_conditions or {}),
        "declared_model": declared_model,
        "acceptance_criteria": dict(acceptance_criteria or {}),
        "existing_evidence": list(existing_evidence or []),
        "plate_layout": layout,
        "limitations": [
            "Modeled products do not establish amplification efficiency, LOD, fluorescence, inhibition, or matrix performance.",
            "Case selection is bounded and diversity-aware, not a proof of global experimental optimality.",
            "Candidate comparison is valid only for the declared cases, controls, conditions, and acceptance criteria.",
        ],
    }
    plan["plan_sha256"] = sha256_value(_canonical_plan_payload(plan))
    return plan


def parse_results_csv(text, plan):
    if not isinstance(text, str) or len(text.encode("utf-8")) > 2_000_000:
        raise ValueError("results CSV is missing or too large")
    expected = {x["well"]: x for x in (plan.get("plate_layout") or {}).get("wells", [])}
    reader = csv.DictReader(io.StringIO(text))
    rows, seen, errors = [], set(), []
    allowed = {"amplified", "not_amplified", "missing", "invalid"}
    for lineno, raw in enumerate(reader, 2):
        well = str(raw.get("well") or "").strip().upper()
        if well not in expected:
            errors.append("line %d: unknown well %s" % (lineno, well or "(blank)")); continue
        if well in seen:
            errors.append("line %d: duplicate well %s" % (lineno, well)); continue
        seen.add(well)
        observed = str(raw.get("observed") or "missing").strip().lower().replace(" ", "_")
        if observed not in allowed:
            errors.append("line %d: observed must be amplified, not_amplified, missing, or invalid" % lineno); continue
        def num(name, lo=None, hi=None):
            val = str(raw.get(name) or "").strip()
            if not val:
                return None
            try:
                out = float(val)
            except ValueError:
                raise ValueError("line %d: %s is not numeric" % (lineno, name))
            if (lo is not None and out < lo) or (hi is not None and out > hi):
                raise ValueError("line %d: %s is outside the permitted range" % (lineno, name))
            return out
        try:
            cq = num("cq", 0, 80); efficiency = num("efficiency", 0, 300)
        except ValueError as exc:
            errors.append(str(exc)); continue
        base = expected[well]
        rows.append({**base, "observed": observed, "cq": cq, "efficiency": efficiency,
                     "melt_abnormal": str(raw.get("melt_abnormal") or "").lower() in {"1", "true", "yes"},
                     "unexpected_product": str(raw.get("unexpected_product") or "").lower() in {"1", "true", "yes"},
                     "probe_signal": str(raw.get("probe_signal") or "").strip() or None,
                     "notes": str(raw.get("notes") or "")[:2000]})
    if errors:
        raise ValueError("; ".join(errors[:20]))
    return rows


def interpret_results(plan, result_rows):
    controls, invalid_reasons = [], []
    for row in result_rows:
        if row["well_type"] == "no_template_control" and row["observed"] == "amplified":
            invalid_reasons.append("no-template control amplified in %s" % row["well"])
        if row["well_type"] == "positive_control" and row["observed"] != "amplified":
            invalid_reasons.append("positive control failed in %s" % row["well"])
        if row["well_type"] != "test":
            controls.append(row)
    by_key = defaultdict(list)
    for row in result_rows:
        if row["well_type"] == "test":
            by_key[(row["candidate_id"], row["case_id"])].append(row)
    aggregates, candidate_scores = [], defaultdict(lambda: {"supported": 0, "contradicted": 0, "missing": 0})
    for (candidate_id, case_id), reps in sorted(by_key.items()):
        observed = [x["observed"] for x in reps if x["observed"] not in {"missing", "invalid"}]
        cqs = [x["cq"] for x in reps if x["cq"] is not None]
        predicted = reps[0].get("expected")
        pred_amp = predicted in {"product", "signal_product"}
        n_amp = sum(1 for x in observed if x == "amplified")
        n_not = sum(1 for x in observed if x == "not_amplified")
        if not observed:
            agreement = "missing"
        elif (pred_amp and n_amp == len(observed)) or ((not pred_amp) and n_not == len(observed)):
            agreement = "supported"
        elif n_amp and n_not:
            agreement = "mixed"
        else:
            agreement = "contradicted"
        candidate_scores[candidate_id]["supported" if agreement == "supported" else
                                       "missing" if agreement in {"missing", "mixed"} else "contradicted"] += 1
        mean = sum(cqs) / len(cqs) if cqs else None
        sd = math.sqrt(sum((x - mean) ** 2 for x in cqs) / (len(cqs) - 1)) if len(cqs) > 1 else None
        aggregates.append({"candidate_id": candidate_id, "case_id": case_id, "n_replicates": len(reps),
                           "n_observed": len(observed), "n_amplified": n_amp, "n_not_amplified": n_not,
                           "mean_cq": round(mean, 3) if mean is not None else None,
                           "cq_sd": round(sd, 3) if sd is not None else None,
                           "prediction": predicted, "prediction_agreement": agreement,
                           "melt_abnormal": any(x["melt_abnormal"] for x in reps),
                           "unexpected_product": any(x["unexpected_product"] for x in reps)})
    valid = not invalid_reasons
    score_rows = [{"candidate_id": k, **v} for k, v in sorted(candidate_scores.items())]
    ranked = sorted(score_rows, key=lambda x: (x["contradicted"], -x["supported"], x["missing"], x["candidate_id"]))
    if not valid:
        conclusion = "invalid controls; candidate comparison is not interpretable"
        strength = "invalid"
    elif not ranked or all((x["supported"] + x["contradicted"]) == 0 for x in ranked):
        conclusion = "inconclusive; observations are insufficient"
        strength = "inconclusive"
    elif len(ranked) > 1 and (ranked[0]["contradicted"], -ranked[0]["supported"]) == (ranked[1]["contradicted"], -ranked[1]["supported"]):
        conclusion = "inconclusive; leading candidates are not separated by this experiment"
        strength = "inconclusive"
    else:
        lead = ranked[0]
        margin = (ranked[1]["contradicted"] - lead["contradicted"]) + (lead["supported"] - ranked[1]["supported"])
        strength = "strong for this declared experiment" if margin >= 3 and lead["contradicted"] == 0 else "moderate for this declared experiment"
        conclusion = "%s favors %s; this is not a global assay-performance claim" % (strength.capitalize(), lead["candidate_id"])
    report = {
        "schema_version": "oligoforge-validation-interpretation/v1",
        "validation_studio_version": VALIDATION_STUDIO_VERSION,
        "plan_sha256": plan.get("plan_sha256"), "controls_valid": valid,
        "invalid_control_reasons": invalid_reasons, "control_rows": controls,
        "candidate_case_results": aggregates, "candidate_summary": score_rows,
        "conclusion_strength": strength, "conclusion": conclusion,
        "predictions_supported": [x for x in aggregates if x["prediction_agreement"] == "supported"],
        "predictions_contradicted": [x for x in aggregates if x["prediction_agreement"] == "contradicted"],
        "uncertainties_remaining": [x for x in aggregates if x["prediction_agreement"] in {"mixed", "missing"}],
        "limitations": ["Interpretation applies only to this declared experiment and acceptance criteria.",
                        "Experimental evidence is reported separately and does not silently retrain the deterministic ranker."],
    }
    report["report_sha256"] = sha256_value(report)
    return report


def plan_json(plan):
    return json.dumps(plan, sort_keys=True, indent=2, ensure_ascii=False) + "\n"

