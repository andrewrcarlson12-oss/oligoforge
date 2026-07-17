#!/usr/bin/env python3
"""Regenerate the frozen biological ranking trace from local fixtures.

This is intentionally separate from the fast synthetic evidence-order benchmark.
It executes the real bounded design pipeline and records a compact candidate trace,
the complete attrition ledger, and the self-hashing scientific run manifest.  It
uses no network services and makes no wet-lab optimality claim.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from oligoforge import __version__
from oligoforge import autodesign as AD
from oligoforge import profiles as P

FIXTURES = ROOT / "tests" / "fixtures"
OUT = ROOT / "tests" / "benchmark" / "plasmodium_ranking_trace.json.gz"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate(row):
    e = row.get("evidence") or {}
    off = e.get("offtarget") or {}
    target = e.get("target") or {}
    robust = e.get("condition_robustness") or {}
    return {
        "rank": row.get("rank"),
        "finalist_category": row.get("finalist_category"),
        "display_score": row.get("display_score"),
        "equivalence_group": row.get("equivalence_group"),
        "assay": row.get("assay"),
        "evidence_summary": {
            "hard_valid": e.get("hard_valid"),
            "hard_failures": e.get("hard_failures"),
            "target_coverage": e.get("target_coverage"),
            "target_product_subjects": target.get("product_subjects"),
            "target_signal_subjects": target.get("signal_subjects"),
            "worst_isolate_3prime": e.get("worst_isolate_3prime"),
            "probe_mean_identity": e.get("probe_mean_identity"),
            "offtarget_product_subjects": off.get("product_subjects"),
            "offtarget_signal_subjects": off.get("signal_subjects"),
            "condition_valid_fraction": robust.get("valid_fraction"),
            "condition_failure_reasons": sorted({
                reason
                for scenario in robust.get("scenarios") or []
                for reason in scenario.get("failure_reasons") or []
            }),
            "triplet_penalty": e.get("triplet_penalty"),
            "practical_penalty": e.get("practical_penalty"),
            "degeneracy_fold": e.get("degeneracy_fold"),
            "panel_risk": e.get("panel_risk"),
            "evaluations": e.get("evaluations"),
        },
        "ranking_evidence_vector": row.get("ranking_evidence_vector"),
        "rank_trace": row.get("rank_trace"),
        "rank_explanation": row.get("rank_explanation"),
    }


def main() -> int:
    plas_path = FIXTURES / "plasmodium_cytb.json"
    haem_path = FIXTURES / "haemoproteus_cytb.json"
    expected_path = FIXTURES / "autodesign_expected_v1.37.json"
    plas = json.loads(plas_path.read_text())["sequences"]
    haem = json.loads(haem_path.read_text())["sequences"]
    expected = json.loads(expected_path.read_text())
    started = time.perf_counter()
    result = AD.design_from_sequences(
        plas, P.PROFILES["parasite_mtdna"], offs=haem,
        n_candidates=5, objective="discrimination"
    )
    elapsed = round(time.perf_counter() - started, 3)
    error = result.get("error")
    ids = [
        [c["assay"]["forward"], c["assay"]["reverse"], c["assay"].get("probe")]
        for c in result.get("candidates") or []
    ]
    output = {
        "schema": "oligoforge-biological-ranking-trace/v3",
        "generated_for": "OligoForge %s" % __version__,
        "purpose": (
            "Versioned computational trace on the frozen Plasmodium/Haemoproteus fixtures; "
            "not wet-lab optimality evidence."
        ),
        "fixture_sha256": {
            plas_path.name: _sha(plas_path),
            haem_path.name: _sha(haem_path),
            expected_path.name: _sha(expected_path),
        },
        "runtime_seconds": elapsed,
        "error": error,
        "expected_ordering_match": ids == expected.get("discrimination_ids"),
        "n_targets": len(plas),
        "n_offs": len(haem),
        "n_candidates_screened": result.get("n_candidates_screened"),
        "n_candidates": result.get("n_candidates"),
        "search_status": result.get("search_status"),
        "objective_profile": result.get("objective_profile"),
        "ranker_manifest": result.get("ranker_manifest"),
        "candidate_attrition": result.get("candidate_attrition"),
        "ranking_statement": result.get("ranking_statement"),
        "candidates": [_candidate(c) for c in result.get("candidates") or []],
    }
    rendered = (json.dumps(output, indent=2, sort_keys=False) + "\n").encode("utf-8")
    # A blank embedded filename and epoch mtime make the compressed evidence
    # byte-reproducible while retaining the complete machine trace.
    with OUT.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=9,
                           fileobj=raw, mtime=0) as compressed:
            compressed.write(rendered)
    print(json.dumps({
        "output": str(OUT.relative_to(ROOT)),
        "runtime_seconds": elapsed,
        "expected_ordering_match": output["expected_ordering_match"],
        "run_id": (output.get("ranker_manifest") or {}).get("run_id"),
        "candidates": len(output["candidates"]),
    }, indent=2))
    return 0 if not error and output["expected_ordering_match"] else 1


if __name__ == "__main__":
    # primer3-py can stall during interpreter finalization after exhaustive work
    # on CPython 3.13. Flush the result and terminate this isolated process.
    code = main()
    sys.stdout.flush(); sys.stderr.flush()
    os._exit(code)
