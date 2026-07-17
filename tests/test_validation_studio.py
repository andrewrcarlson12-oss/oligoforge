"""Validation Studio selection, layout, export and conservative interpretation."""
import csv
import io
import random
import sys
sys.path.insert(0, ".")

from oligoforge import validation_studio as VS
from oligoforge.specificity import _rc_iupac


random.seed(707)
bg = lambda n: "".join(random.choice("ACGT") for _ in range(n))
F1 = "ACGTACGTACGTACGTACGT"
F2 = F1[:-1] + "A"  # fails on the F1 target at the extension-critical 3' end
R = "TTTTGGGGCCCCAAAATTTT"
P = "GATTACAGATTACAGATTAC"
target = bg(30) + F1 + bg(40) + P + bg(40) + _rc_iupac(R) + bg(30)
target2 = bg(30) + F2 + bg(40) + P + bg(40) + _rc_iupac(R) + bg(30)

candidates = [
    {"candidate_id": "A", "chemistry": "probe", "assay": {"forward": F1, "reverse": R, "probe": P}},
    {"candidate_id": "B", "chemistry": "probe", "assay": {"forward": F2, "reverse": R, "probe": P}},
]
cases = [
    {"case_id": "t1", "role": "target", "group": "lineage-1", "sequence": target, "source_type": "synthetic"},
    {"case_id": "t2", "role": "target", "group": "lineage-2", "sequence": target2, "source_type": "synthetic"},
]

plan = VS.create_plan(candidates, cases, plate_format=96, replicates=2,
                      controls={"ntc_replicates": 1, "positive_controls": [{"control_id": "PC"}]},
                      model={"max_mm": 2}, seed=42, use_edge_wells=True)
assert plan["selection_status"]["n_selected"] == 2, plan["selection_status"]
assert plan["selection_status"]["globally_optimal"] is False
assert plan["plate_layout"]["candidate_interleaved"] is True
assert plan["plate_layout"]["randomization_seed"] == 42
assert plan["plan_sha256"] == VS.create_plan(candidates, cases, plate_format=96, replicates=2,
                      controls={"ntc_replicates": 1, "positive_controls": [{"control_id": "PC"}]},
                      model={"max_mm": 2}, seed=42, use_edge_wells=True)["plan_sha256"]

# Candidate chemistry participates in identity; exact duplicates are suppressed.
norm = VS.normalize_candidates(candidates + [{"candidate_id": "dup", **candidates[0]}])
assert len(norm) == 2
try:
    VS.normalize_candidates([candidates[0], dict(candidates[1], candidate_id="A")])
except ValueError:
    pass
else:
    raise AssertionError("duplicate candidate IDs should fail when molecular identities differ")

csv_template = VS.plate_csv(plan["plate_layout"])
parsed = list(csv.DictReader(io.StringIO(csv_template)))
assert len(parsed) == plan["plate_layout"]["n_wells"]
assert {x["candidate_id"] for x in parsed if x["well_type"] == "test"} == {"A", "B"}

# Fill a result set that agrees with predictions and keeps all controls valid.
out = io.StringIO(newline="")
fields = list(parsed[0].keys())
w = csv.DictWriter(out, fieldnames=fields, lineterminator="\n"); w.writeheader()
for row in parsed:
    if row["well_type"] == "no_template_control":
        row["observed"] = "not_amplified"
    elif row["well_type"] == "positive_control":
        row["observed"] = "amplified"
    else:
        row["observed"] = "amplified" if row["expected"] in {"product", "signal_product"} else "not_amplified"
        row["cq"] = "25.0" if row["observed"] == "amplified" else ""
    w.writerow(row)
results = VS.parse_results_csv(out.getvalue(), plan)
report = VS.interpret_results(plan, results)
assert report["controls_valid"] is True
assert report["predictions_supported"] and not report["predictions_contradicted"]

# An amplified NTC invalidates interpretation rather than merely warning.
bad = [dict(x) for x in results]
next(x for x in bad if x["well_type"] == "no_template_control")["observed"] = "amplified"
invalid = VS.interpret_results(plan, bad)
assert invalid["controls_valid"] is False and invalid["conclusion_strength"] == "invalid"

# CSV formula injection is neutralized.
layout = {"wells": [{"well": "A1", "well_type": "test", "candidate_id": "=CMD()", "case_id": "x"}]}
assert "'=CMD()" in VS.plate_csv(layout)

print("VALIDATION STUDIO TESTS PASS")
