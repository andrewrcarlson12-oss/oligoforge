"""Offline API integration checks for long jobs, Validation Studio, and Assurance."""
import csv
import io
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from fastapi.testclient import TestClient

import app as server


class StubJobs:
    def __init__(self):
        self.payload = None

    def stats(self):
        return {"backend": "in_memory_single_worker", "queue_capacity": 2,
                "primary_timeout_seconds": 240, "blast_timeout_seconds": 360,
                "terminal_ttl_seconds": 1800}

    def submit(self, payload, idempotency_key=None):
        self.payload = payload
        return {"job_id": "job123", "status": "queued", "stages": [],
                "result": None, "idempotency_seen": bool(idempotency_key)}

    def get(self, job_id):
        if job_id != "job123":
            return None
        return {"job_id": job_id, "status": "succeeded", "stages": [],
                "result": {"candidates": []}}

    def cancel(self, job_id):
        if job_id != "job123":
            return None
        return {"job_id": job_id, "status": "cancelled", "stages": [],
                "result": None}

    def retry_blast(self, job_id, idempotency_key=None, **values):
        if job_id != "job123":
            return None
        return {"job_id": "retry123", "status": "queued", "stages": [],
                "result": {"candidates": []}}


def must(label, condition):
    if not condition:
        raise AssertionError(label)
    print("PASS " + label)


def completed_csv(csv_text):
    rows = list(csv.DictReader(io.StringIO(csv_text)))
    for row in rows:
        row["observed"] = ("amplified" if row["expected"] in
                           {"product", "signal_product", "amplified"}
                           else "not_amplified")
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=rows[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue()


def main():
    original = server.DESIGN_JOBS
    jobs = StubJobs()
    server.DESIGN_JOBS = jobs
    client = TestClient(server.app)
    try:
        r = client.get("/api/autodesign/limits")
        must("automatic-design limits are visible", r.status_code == 200 and
             r.json()["queue_capacity"] == 2 and r.json()["durable"] is False)
        payload = {"target_query": "Synthetic marker", "n_fetch": 3}
        r = client.post("/api/autodesign/jobs", json=payload,
                        headers={"Idempotency-Key": "browser-action"})
        must("automatic design returns a job promptly", r.status_code == 202 and
             r.json()["job_id"] == "job123" and r.json()["status_url"].endswith("job123"))
        must("validated payload reaches the job manager", jobs.payload["target_query"] == "Synthetic marker")
        must("job capability can be polled", client.get("/api/autodesign/jobs/job123").json()["status"] == "succeeded")
        must("unknown capability is a non-reflective 404", client.get("/api/autodesign/jobs/missing").status_code == 404)
        must("job capability can be cancelled", client.delete("/api/autodesign/jobs/job123").json()["status"] == "cancelled")
        must("specificity-only retry gets a new capability",
             client.post("/api/autodesign/jobs/job123/retry-blast", json={}).status_code == 202)
        previous = server.HOSTED_MODE
        server.HOSTED_MODE = True
        try:
            denied = client.post("/api/autodesign/jobs", json={**payload, "run_blast": True,
                                                                  "blast_mode": "local",
                                                                  "blast_db_path": "/private/db"})
            must("hosted automatic design rejects local server paths", denied.status_code == 403 and
                 "/private/db" not in denied.text)
        finally:
            server.HOSTED_MODE = previous

        fwd = "ACGTACGTACGTACGTACGT"
        rev = "TGCATGCATGCATGCATGCA"
        target = fwd + "A" * 50 + rev
        candidates = [
            {"candidate_id": "supported", "forward": fwd, "reverse": rev,
             "chemistry": "sybr"},
            {"candidate_id": "unsupported", "forward": "G" * 20,
             "reverse": "C" * 20, "chemistry": "sybr"},
        ]
        cases = [{"case_id": "target-1", "role": "target", "sequence": target,
                  "source_type": "synthetic", "group": "fixture"}]
        r = client.post("/api/validation-studio/plan", json={"candidates": candidates,
                                                               "cases": cases,
                                                               "replicates": 1,
                                                               "controls": {"ntc_replicates": 1},
                                                               "seed": 17})
        plan_body = r.json()
        must("Validation Studio returns a hashed plan and fillable CSV", r.status_code == 200 and
             plan_body["plan"]["plan_sha256"] and "well,well_type" in plan_body["plate_csv"])
        r = client.post("/api/validation-studio/interpret",
                        json={"plan": plan_body["plan"],
                              "results_csv": completed_csv(plan_body["plate_csv"])})
        must("Validation Studio conservatively interprets completed results",
             r.status_code == 200 and r.json()["controls_valid"] is True and
             r.json()["conclusion_strength"] in {"inconclusive", "moderate for this declared experiment",
                                                  "strong for this declared experiment"} and
             any("declared experiment" in x for x in r.json()["limitations"]))

        assay = {"portfolio_name": "Fixture", "assays": [{"assay_id": "A1",
                 "chemistry": "sybr", "forward": fwd, "reverse": rev}]}
        sbom_r = client.post("/api/assurance/assaysbom", json={"assay": assay})
        sbom = sbom_r.json()
        must("AssaySBOM endpoint returns a deterministic identifier",
             sbom_r.status_code == 200 and sbom["assaysbom_id"].startswith("ofsbom_"))
        fasta = ">target-1\n" + target + "\n"
        base = client.post("/api/assurance/snapshots", json={"fasta": fasta,
                                                               "name": "baseline",
                                                               "role": "target",
                                                               "source": {"adapter": "offline_fixture"},
                                                               "metadata": "record_id,group\ntarget-1,fixture\n"}).json()
        follow = client.post("/api/assurance/snapshots", json={"fasta": fasta,
                                                                 "name": "follow-up",
                                                                 "role": "target"}).json()
        must("sequence snapshots are immutable and distinct by declared metadata",
             base["snapshot_id"].startswith("ofsnap_") and base["snapshot_id"] != follow["snapshot_id"])
        must("snapshot HTTP metadata and source remain structured",
             base["source"]["adapter"] == "offline_fixture" and
             base["unique_records"][0]["metadata"]["group"] == "fixture")
        delta = client.post("/api/assurance/snapshots/delta",
                            json={"baseline": base, "followup": follow}).json()
        must("snapshot delta reports exact unchanged haplotypes", delta["counts"]["unchanged"] == 1)
        scan_r = client.post("/api/assurance/drift-scan", json={"assaysbom": sbom,
                                                                 "baseline_target": base,
                                                                 "current_target": follow})
        scan = scan_r.json()
        must("DriftGuard reconstructs complete products", scan_r.status_code == 200 and
             scan["state"] == "Stable" and
             scan["assay_results"][0]["current_target"]["coherent_products"] == 1)
        records = client.post("/api/assurance/ofvr", json={"drift_scan": scan,
                                                             "issuance_year": 2026}).json()["records"]
        must("stable scan produces no cosmetic vulnerability record", records == [])
        package_r = client.post("/api/assurance/package", json={"assaysbom": sbom,
                                                                   "snapshots": [base, follow],
                                                                   "deltas": [delta],
                                                                   "drift_scans": [scan]})
        must("evidence package self-verifies and has escaped HTML", package_r.status_code == 200 and
             package_r.json()["verification"]["valid"] is True and
             "<script" not in package_r.json()["html"].lower())
    finally:
        server.DESIGN_JOBS = original
        original.shutdown(wait=False)


if __name__ == "__main__":
    main()
    print("NEW API INTEGRATION TESTS PASS")
