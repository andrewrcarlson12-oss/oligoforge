"""Focused offline checks for the asynchronous automatic-design backend.

Run from the repository root with:
    /tmp/oligoforge-venv/bin/python tests/integration/test_autodesign_jobs.py

The repository's historical runner executes root-level standalone scripts only,
so this file keeps the same no-pytest convention while living with integration
tests as the job/backend surface grows.
"""

from copy import deepcopy
import json
import os
import sys
import threading
import time
import traceback

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from oligoforge import autodesign as AD
from oligoforge.jobs import (DesignJobManager, IdempotencyConflict,
                             QueueFull, TERMINAL_STATUSES)


FAILURES = []


def check(name, condition, detail=""):
    ok = bool(condition)
    print(("  PASS " if ok else "  FAIL ") + name +
          (("  -> " + str(detail)) if detail and not ok else ""))
    if not ok:
        FAILURES.append(name)


def stage_map(stages):
    return {
        "resolve_fetch": stages.resolve,
        "design": stages.design,
        "enrich": stages.enrich,
        "blast": stages.blast,
    }


class FakeStages:
    """Deterministic stage functions with optional blocking/failure controls."""

    def __init__(self):
        self.calls = []
        self.resolve_started = threading.Event()
        self.resolve_gate = None
        self.blast_started = threading.Event()
        self.blast_gate = None
        self.fail_blast_count = 0
        self.extra_result = {}

    def resolve(self, target_query, off_query, n_fetch):
        self.calls.append("resolve_fetch")
        self.resolve_started.set()
        if self.resolve_gate is not None:
            self.resolve_gate.wait(3.0)
        raw = "ACGT" * 40
        return {
            "target_query": target_query,
            "off_query": off_query,
            "targets": [raw],
            "offs": [],
            "target_pairs": [("ACC1 synthetic", raw)],
            "off_pairs": [],
            "organism_name": "Synthetic",
            "gene_name": "marker",
            "resolved": None,
        }

    def design(self, context, profile, min_ident, objective):
        self.calls.append("design")
        return {
            "profile_used": profile,
            "candidates": [{"assay": {
                "forward": "ACGTACGTACGTACGTACGT",
                "reverse": "TGCATGCATGCATGCATGCA",
                "probe": "CGTACGTACGTACGTACGTA",
            }}],
            # Deliberately faulty leakage; the manager must remove it.
            "targets": list(context["targets"]),
        }

    def enrich(self, out, context, *, min_ident, prefer_junction, nested, objective):
        self.calls.append("enrich")
        result = deepcopy(out)
        result["primary_complete"] = True
        result.update(deepcopy(self.extra_result))
        return result

    def blast(self, out, *, blast_mode, blast_db, blast_db_path, organism,
              suppress_errors):
        self.calls.append("blast")
        self.blast_started.set()
        if self.blast_gate is not None:
            self.blast_gate.wait(3.0)
        if self.fail_blast_count:
            self.fail_blast_count -= 1
            raise RuntimeError("BLAST failed with PRIVATE-EXCEPTION-DETAIL")
        result = deepcopy(out)
        result["specificity"] = {"mode": blast_mode, "hits": 0}
        return result


def make_manager(stages, **kwargs):
    defaults = dict(queue_capacity=4, ttl_s=5.0, primary_timeout_s=2.0,
                    blast_timeout_s=2.0, poll_interval_s=0.01,
                    stage_functions=stage_map(stages))
    defaults.update(kwargs)
    return DesignJobManager(**defaults)


def wait_terminal(manager, job_id, timeout=1.5):
    snapshot = manager.wait(job_id, timeout)
    return snapshot if snapshot and snapshot["status"] in TERMINAL_STATUSES else snapshot


def test_prompt_submission_and_stage_records():
    stages = FakeStages()
    stages.resolve_gate = threading.Event()
    manager = make_manager(stages)
    try:
        before = time.perf_counter()
        submitted = manager.submit({"target_query": "Synthetic marker", "profile": "auto"})
        elapsed = time.perf_counter() - before
        check("submission returns promptly while retrieval is blocked", elapsed < 0.25, elapsed)
        check("capability ID has high entropy", len(submitted["id"]) >= 40, submitted["id"])
        check("snapshot exposes both stable ID field names",
              submitted["id"] == submitted["job_id"])
        check("initial snapshot contains no result", submitted["result"] is None)
        check("stage records declare required versus optional",
              [s["required"] for s in submitted["stages"]] == [True, True, True, False])
        check("retrieval stage really started", stages.resolve_started.wait(1.0))
        running = manager.get(submitted["id"])
        check("polling reports a real running stage",
              running["status"] == "running" and running["stages"][0]["status"] == "running",
              running)
        stages.resolve_gate.set()
        final = wait_terminal(manager, submitted["id"])
        check("required stages succeed", final["status"] == "succeeded", final)
        check("required stage order is explicit",
              stages.calls == ["resolve_fetch", "design", "enrich"], stages.calls)
        check("unrequested BLAST is skipped",
              final["stages"][3]["status"] == "skipped", final["stages"][3])
        check("primary result is returned", final["result"].get("primary_complete") is True)
        check("internal fetched corpus is stripped from result", "targets" not in final["result"])
        final["result"]["primary_complete"] = False
        check("callers cannot mutate retained job state through a snapshot",
              manager.get(submitted["id"])["result"]["primary_complete"] is True)
    finally:
        stages.resolve_gate.set()
        manager.shutdown()


def test_idempotency_and_conflict():
    stages = FakeStages()
    manager = make_manager(stages)
    try:
        payload = {"target_query": "Synthetic marker", "profile": "auto"}
        first = manager.submit(payload, idempotency_key="one-client-action")
        second = manager.submit(payload, idempotency_key="one-client-action")
        check("same idempotency key and input return the same job",
              first["id"] == second["id"], (first["id"], second["id"]))
        conflicted = False
        try:
            manager.submit({"target_query": "Different marker"},
                           idempotency_key="one-client-action")
        except IdempotencyConflict:
            conflicted = True
        check("idempotency key reuse with different input is rejected", conflicted)
        wait_terminal(manager, first["id"])
        check("idempotent replay executes retrieval only once",
              stages.calls.count("resolve_fetch") == 1, stages.calls)
    finally:
        manager.shutdown()


def test_queue_cap_and_queued_cancel():
    stages = FakeStages()
    stages.resolve_gate = threading.Event()
    manager = make_manager(stages, queue_capacity=1)
    try:
        first = manager.submit({"target_query": "first"})
        check("first job occupies worker", stages.resolve_started.wait(1.0))
        second = manager.submit({"target_query": "second"})
        full = False
        try:
            manager.submit({"target_query": "third"})
        except QueueFull:
            full = True
        check("bounded waiting queue rejects excess work", full)
        cancelled = manager.cancel(second["id"])
        check("queued cancellation is immediately terminal",
              cancelled["status"] == "cancelled", cancelled)
        check("queued cancellation has no partial result", cancelled["result"] is None)
        stages.resolve_gate.set()
        check("active job can still finish after queued cancellation",
              wait_terminal(manager, first["id"])["status"] == "succeeded")
    finally:
        stages.resolve_gate.set()
        manager.shutdown()


def test_running_cancel_is_prompt_and_has_no_partial_result():
    stages = FakeStages()
    stages.resolve_gate = threading.Event()
    manager = make_manager(stages)
    try:
        job = manager.submit({"target_query": "cancel me"})
        check("cancellable stage started", stages.resolve_started.wait(1.0))
        requested = manager.cancel(job["id"])
        check("cancel request is recorded", requested["cancel_requested"] is True)
        final = wait_terminal(manager, job["id"], 0.75)
        check("running cancellation becomes visible promptly",
              final and final["status"] == "cancelled", final)
        check("required-stage cancellation never exposes a finalist",
              final and final["result"] is None)
        check("active stage is marked cancelled",
              final and final["stages"][0]["status"] == "cancelled", final)
    finally:
        stages.resolve_gate.set()  # lets the sole worker drain the legacy call
        manager.shutdown()


def test_required_timeout_is_terminal_without_result():
    stages = FakeStages()
    stages.resolve_gate = threading.Event()
    manager = make_manager(stages, primary_timeout_s=0.06)
    try:
        job = manager.submit({"target_query": "slow target"})
        check("slow required stage started", stages.resolve_started.wait(1.0))
        final = wait_terminal(manager, job["id"], 0.75)
        check("primary deadline reports timed_out promptly",
              final and final["status"] == "timed_out", final)
        check("primary timeout has no result", final and final["result"] is None)
        check("timed-out stage is explicit",
              final and final["stages"][0]["status"] == "timed_out", final)
    finally:
        stages.resolve_gate.set()
        manager.shutdown()


def test_optional_blast_failure_retains_primary_and_retry_reuses_it():
    stages = FakeStages()
    stages.fail_blast_count = 1
    manager = make_manager(stages)
    try:
        job = manager.submit({"target_query": "target", "run_blast": True,
                              "blast_mode": "remote"})
        failed_blast = wait_terminal(manager, job["id"])
        check("optional BLAST failure is a warning success",
              failed_blast["status"] == "succeeded_with_warnings", failed_blast)
        check("primary result survives optional failure",
              failed_blast["result"].get("primary_complete") is True,
              failed_blast["result"])
        check("optional stage is marked failed",
              failed_blast["stages"][3]["status"] == "failed",
              failed_blast["stages"][3])
        encoded = json.dumps(failed_blast)
        check("native/network exception detail is never serialized",
              "PRIVATE-EXCEPTION-DETAIL" not in encoded)

        retried = manager.retry_blast(job["id"], idempotency_key="blast-retry-once")
        replay = manager.retry_blast(job["id"], idempotency_key="blast-retry-once")
        check("BLAST retry is idempotent", retried["id"] == replay["id"])
        retry_final = wait_terminal(manager, retried["id"])
        check("BLAST-only retry can succeed", retry_final["status"] == "succeeded",
              retry_final)
        check("retry stages prove primary work was reused",
              [s["status"] for s in retry_final["stages"][:3]] ==
              ["skipped", "skipped", "skipped"], retry_final["stages"])
        check("retry never refetches or redesigns",
              stages.calls.count("resolve_fetch") == 1 and
              stages.calls.count("design") == 1 and
              stages.calls.count("enrich") == 1 and
              stages.calls.count("blast") == 2, stages.calls)
        check("successful retry adds specificity to retained primary",
              retry_final["result"].get("specificity", {}).get("hits") == 0,
              retry_final["result"])
    finally:
        manager.shutdown()


def test_optional_blast_timeout_retains_primary():
    stages = FakeStages()
    stages.blast_gate = threading.Event()
    manager = make_manager(stages, blast_timeout_s=0.06)
    try:
        job = manager.submit({"target_query": "target", "run_blast": True})
        check("optional BLAST stage started", stages.blast_started.wait(1.0))
        final = wait_terminal(manager, job["id"], 0.75)
        check("optional timeout yields succeeded_with_warnings",
              final and final["status"] == "succeeded_with_warnings", final)
        check("optional timeout retains primary result",
              final and final["result"].get("primary_complete") is True, final)
        check("optional timeout is explicit in stage status",
              final and final["stages"][3]["status"] == "timed_out", final)
    finally:
        stages.blast_gate.set()
        manager.shutdown()


def test_snapshot_sanitization():
    stages = FakeStages()
    email = "private.person@example.test"
    key = "NCBI-PRIVATE-KEY-123"
    db_path = "/private/genomes/person/db"
    raw_sequence = "ACGTRYSWKMBDHVN" * 6
    stages.extra_result = {
        "email": email,
        "ncbi_key": key,
        "blast_db_path": db_path,
        "targets": [raw_sequence],
        "note": " | ".join([email, key, db_path, raw_sequence]),
    }
    manager = make_manager(stages)
    try:
        job = manager.submit({
            "target_query": raw_sequence,
            "email": email,
            "ncbi_key": key,
            "blast_db_path": db_path,
        })
        final = wait_terminal(manager, job["id"])
        serialized = json.dumps(final, sort_keys=True)
        check("snapshot omits contact email", email not in serialized)
        check("snapshot omits NCBI key", key not in serialized)
        check("snapshot omits local database path", db_path not in serialized)
        check("snapshot omits raw private target sequence", raw_sequence not in serialized)
        check("snapshot remains JSON serializable", json.loads(serialized)["id"] == job["id"])
        check("manager intentionally provides no job-list capability",
              not hasattr(manager, "list"))
    finally:
        manager.shutdown()


def test_terminal_ttl_and_idempotency_expiry():
    stages = FakeStages()
    manager = make_manager(stages, ttl_s=0.06)
    try:
        payload = {"target_query": "short lived"}
        first = manager.submit(payload, idempotency_key="expires-too")
        check("short-lived job completes",
              wait_terminal(manager, first["id"])["status"] == "succeeded")
        time.sleep(0.09)
        check("terminal capability is forgotten after TTL", manager.get(first["id"]) is None)
        second = manager.submit(payload, idempotency_key="expires-too")
        check("expired idempotency entry can create fresh work",
              second["id"] != first["id"], (first["id"], second["id"]))
        wait_terminal(manager, second["id"])
    finally:
        manager.shutdown()


def test_legacy_entrypoint_composes_reusable_stages():
    originals = (AD.resolve_and_fetch_query, AD.design_query_corpus,
                 AD.enrich_query_design, AD.blast_winner)
    calls = []

    def resolve(target_query, off_query, n_fetch):
        calls.append("resolve_fetch")
        return {"target_query": target_query, "targets": ["A" * 100]}

    def design(context, profile, min_ident, objective):
        calls.append("design")
        return {"candidates": [{"assay": {"forward": "A", "reverse": "T"}}]}

    def enrich(out, context, **kwargs):
        calls.append("enrich")
        out["enriched"] = True
        return out

    def blast(out, **kwargs):
        calls.append("blast")
        check("legacy entrypoint requests BLAST error suppression",
              kwargs.get("suppress_errors") is True, kwargs)
        out["specificity"] = {"hits": 0}
        return out

    try:
        AD.resolve_and_fetch_query = resolve
        AD.design_query_corpus = design
        AD.enrich_query_design = enrich
        AD.blast_winner = blast
        out = AD.design_from_query("Synthetic marker", run_blast=True)
        check("legacy entrypoint composes all four stages in order",
              calls == ["resolve_fetch", "design", "enrich", "blast"], calls)
        check("legacy composed result is returned",
              out.get("enriched") and out.get("specificity", {}).get("hits") == 0, out)
    finally:
        (AD.resolve_and_fetch_query, AD.design_query_corpus,
         AD.enrich_query_design, AD.blast_winner) = originals


TESTS = [
    test_prompt_submission_and_stage_records,
    test_idempotency_and_conflict,
    test_queue_cap_and_queued_cancel,
    test_running_cancel_is_prompt_and_has_no_partial_result,
    test_required_timeout_is_terminal_without_result,
    test_optional_blast_failure_retains_primary_and_retry_reuses_it,
    test_optional_blast_timeout_retains_primary,
    test_snapshot_sanitization,
    test_terminal_ttl_and_idempotency_expiry,
    test_legacy_entrypoint_composes_reusable_stages,
]


if __name__ == "__main__":
    for test in TESTS:
        print("\n== " + test.__name__ + " ==")
        try:
            test()
        except BaseException:
            FAILURES.append(test.__name__ + " (unexpected exception)")
            traceback.print_exc()
    if FAILURES:
        print("\nFAIL: " + ", ".join(FAILURES))
        sys.exit(1)
    print("\nALL AUTODESIGN JOB ASSERTS PASS")
