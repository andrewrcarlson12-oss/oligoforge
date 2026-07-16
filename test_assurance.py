"""End-to-end Assurance record -> snapshots -> delta -> scan -> OFVR -> package."""
import copy
import json
import os
import random
import subprocess
import sys
import tempfile
sys.path.insert(0, ".")

from oligoforge import validation_studio as VS
from oligoforge.assurance import *
from oligoforge.assurance import snapshots as SNAPSHOTS
from oligoforge.assurance.evidence_package import verify_evidence_package
from oligoforge.specificity import _rc_iupac


F = "ACGTACGTACGTACGTACGT"
R = "TTTTGGGGCCCCAAAATTTT"
P = "GATTACAGATTACAGATTAC"
random.seed(808)
bg = lambda n: "".join(random.choice("ACGT") for _ in range(n))
clean = bg(25) + F + bg(45) + P + bg(45) + _rc_iupac(R) + bg(25)
dropout = bg(25) + F[:-1] + ("A" if F[-1] != "A" else "C") + bg(45) + P + bg(45) + _rc_iupac(R) + bg(25)
off_clean = bg(220)
off_signal = bg(25) + F + bg(45) + P + bg(45) + _rc_iupac(R) + bg(25)

legacy = {"name": "Synthetic example", "version": "1", "assays": [{"assay_id": "assay-A",
          "chemistry": "idt_taqman", "forward": F, "reverse": R,
          "probe": "/56-FAM/" + P + "/3IABkFQ/", "intended_target_groups": ["target"]}]}
sbom = build_assaysbom(legacy)
assert sbom["schema_version"] == ASSAYSBOM_SCHEMA and sbom["assaysbom_id"].startswith("ofsbom_")
assert build_assaysbom(legacy)["assaysbom_id"] == sbom["assaysbom_id"]
probe = next(x for x in sbom["assays"][0]["components"] if x["role"] == "probe")
assert probe["sequence"] == P and probe["order_sequence"].startswith("/56-FAM/")
assert "<script" not in assaysbom_html({**legacy, "name": "<script>alert(1)</script>"}).lower()

base_t = build_snapshot(">t1 baseline\n" + clean + "\n", name="target baseline", role="target",
                        metadata="record_id,group\nt1,lineage-1\n")
current_t = build_snapshot(">t1 unchanged\n" + clean + "\n>t2 new\n" + dropout + "\n",
                           name="target followup", role="target", baseline_snapshot_id=base_t["snapshot_id"],
                           metadata="record_id,group\nt1,lineage-1\nt2,lineage-2\n")
base_o = build_snapshot(">o1\n" + off_clean, name="off baseline", role="off_target")
current_o = build_snapshot(">o1\n" + off_clean + "\n>o2\n" + off_signal,
                           name="off followup", role="off_target", baseline_snapshot_id=base_o["snapshot_id"])
assert validate_snapshot(base_t)["valid"]
delta = snapshot_delta(base_t, current_t)
assert delta["counts"] == {"added": 1, "removed": 0, "unchanged": 1}

# Exact duplicates remain in the raw ledger but are deduplicated for haplotype metrics.
dedup = build_snapshot(">a\n" + clean + "\n>b\n" + clean, role="target")
assert dedup["metrics"]["raw_record_count"] == 2 and dedup["metrics"]["unique_sequence_count"] == 1
assert dedup["accepted_records"][1]["disposition"] == "accepted_duplicate"

scan = scan_drift(sbom, base_t, current_t, baseline_offtarget=base_o, current_offtarget=current_o,
                  model={"max_mm": 2})
assert scan["state"] == "Possible signal-generating off-target", scan["state"]
codes = {x["code"] for x in scan["reason_records"]}
assert "new_target_lost_coherent_product" in codes and "new_signal_generating_off_target" in codes
ofvrs = generate_ofvrs(scan, issuance_year=2026)
assert len(ofvrs) == 2 and all(x["ofvr_id"].startswith("OFVR-2026-") for x in ofvrs)
assert generate_ofvrs(scan, issuance_year=2026) == ofvrs

plan = VS.create_plan([
    {"candidate_id": "A", "assay": {"forward": F, "reverse": R, "probe": P}},
    {"candidate_id": "B", "assay": {"forward": F[:-1] + "A", "reverse": R, "probe": P}},
], [{"case_id": "case-1", "role": "target", "sequence": clean, "source_type": "synthetic"}],
    model={"max_mm": 2}, replicates=1, controls={"ntc_replicates": 1})
package = build_evidence_package(assaysbom=sbom, snapshots=[base_t, current_t, base_o, current_o],
                                 deltas=[delta], drift_scans=[scan], vulnerabilities=ofvrs,
                                 validation_plans=[plan])
assert verify_evidence_package(package)["valid"]
tampered = copy.deepcopy(package); tampered["artifacts"]["assaysbom"]["portfolio_name"] = "changed"
assert not verify_evidence_package(tampered)["valid"]
assert "<script" not in evidence_package_html(package).lower()

# Parser bounds/rejections and CSV metadata formula-like strings remain inert data.
try:
    build_snapshot(">bad\nACGTX", role="target")
except ValueError:
    pass
else:
    raise AssertionError("all-invalid snapshot should fail")

old_cap = SNAPSHOTS.MAX_RECORDS
try:
    SNAPSHOTS.MAX_RECORDS = 2
    try:
        build_snapshot(">r1\nA\n>r2\nC\n>r3\nG\n", role="target")
    except ValueError as exc:
        assert "exceeds 2 records" in str(exc)
    else:
        raise AssertionError("snapshot record cap allowed one extra final record")
finally:
    SNAPSHOTS.MAX_RECORDS = old_cap

# CLI performs an offline build and emits deterministic JSON.
with tempfile.TemporaryDirectory() as td:
    inp = os.path.join(td, "input.json"); out = os.path.join(td, "sbom.json")
    with open(inp, "w", encoding="utf-8") as h: json.dump(legacy, h)
    p = subprocess.run([sys.executable, "-m", "oligoforge.assurance_cli", "build-assaysbom", inp, out],
                       capture_output=True, text=True, check=True)
    assert json.load(open(out, encoding="utf-8"))["assaysbom_id"] == sbom["assaysbom_id"]

print("ASSURANCE WORKFLOW TESTS PASS")
