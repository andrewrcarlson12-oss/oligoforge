"""Bounded, capability-addressed automatic-design jobs.

This module is intentionally a small single-process backend for OligoForge's
current one-worker deployment.  It provides prompt HTTP-friendly submission and
status snapshots without pretending to be a durable distributed queue:

* jobs disappear on process restart and after their terminal TTL;
* exactly one scientific stage is allowed to run at a time;
* capability-style random job IDs are the only lookup mechanism (there is no
  public job-list operation);
* fetched corpora, credentials, local database paths, and exception details are
  retained only in private in-memory job state and never included in snapshots.

Cancellation and deadlines are observed promptly by the manager, between the
cooperative stage boundaries introduced in :mod:`oligoforge.autodesign`.  Python
cannot safely kill a thread that is inside Primer3, NCBI, or legacy BLAST code.
When such a call is cancelled or times out, the public job becomes terminal at
once, but the sole worker waits for that call to return before starting another
scientific stage.  This preserves the one-active-design safety invariant and is
strictly more honest than reporting a partial candidate as a completed design.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import queue
import secrets
import threading
import time
from typing import Any, Callable, Mapping, Optional

from . import autodesign as AD


TERMINAL_STATUSES = frozenset({
    "succeeded", "succeeded_with_warnings", "failed", "timed_out", "cancelled"
})

_STAGE_DEFS = (
    ("resolve_fetch", True),
    ("design", True),
    ("enrich", True),
    ("blast", False),
)

_DROP_KEYS = frozenset({
    "email", "ncbi_key", "ncbikey", "api_key", "apikey", "token",
    "password", "secret", "blast_db_path", "db_path",
    # Internal stage-corpus names.  Normal autodesign results do not contain
    # these, but dropping them makes snapshots safe even with a faulty stage.
    "targets", "offs", "off_targets", "target_pairs", "off_pairs",
    "raw_sequence", "raw_sequences", "sequences", "template",
})


class QueueFull(RuntimeError):
    """The bounded waiting queue has no free slot."""


class IdempotencyConflict(RuntimeError):
    """An idempotency key was reused for different input."""


class RetryNotAvailable(RuntimeError):
    """A BLAST-only retry was requested without a completed primary result."""


@dataclass
class _Job:
    job_id: str
    payload: dict[str, Any]
    fingerprint: str
    idempotency_hash: Optional[str]
    action: str
    created_wall: float
    created_mono: float
    stages: dict[str, dict[str, Any]]
    status: str = "queued"
    updated_wall: float = 0.0
    started_wall: Optional[float] = None
    finished_wall: Optional[float] = None
    finished_mono: Optional[float] = None
    expires_wall: Optional[float] = None
    warnings: list[str] = field(default_factory=list)
    error: Optional[dict[str, str]] = None
    result: Optional[dict[str, Any]] = None
    primary_result: Optional[dict[str, Any]] = None
    cancel_requested: bool = False
    cancel_event: threading.Event = field(default_factory=threading.Event, repr=False)
    terminal_event: threading.Event = field(default_factory=threading.Event, repr=False)


def _utc(ts: Optional[float]) -> Optional[str]:
    if ts is None:
        return None
    return datetime.fromtimestamp(ts, timezone.utc).isoformat().replace("+00:00", "Z")


def _key_name(value: Any) -> str:
    return "".join(ch for ch in str(value).lower() if ch.isalnum() or ch == "_")


def _looks_like_sequence(value: str) -> bool:
    compact = "".join(value.split()).upper()
    return len(compact) >= 40 and all(ch in "ACGTRYSWKMBDHVN.-" for ch in compact)


class DesignJobManager:
    """One-worker, bounded in-memory job manager for automatic design.

    ``stage_functions`` is an optional mapping used by focused tests and local
    embedders.  Its keys are ``resolve_fetch``, ``design``, ``enrich``, and
    ``blast`` and its functions must have the signatures of the corresponding
    public helpers in :mod:`oligoforge.autodesign`.
    """

    def __init__(self, queue_capacity: int = 8, ttl_s: float = 1800.0,
                 primary_timeout_s: float = 240.0, blast_timeout_s: float = 360.0,
                 poll_interval_s: float = 0.025,
                 stage_functions: Optional[Mapping[str, Callable[..., Any]]] = None):
        if int(queue_capacity) < 1:
            raise ValueError("queue_capacity must be at least 1")
        if float(ttl_s) <= 0:
            raise ValueError("ttl_s must be positive")
        if float(primary_timeout_s) <= 0 or float(blast_timeout_s) <= 0:
            raise ValueError("stage timeouts must be positive")
        self.queue_capacity = int(queue_capacity)
        self.ttl_s = float(ttl_s)
        self.primary_timeout_s = float(primary_timeout_s)
        self.blast_timeout_s = float(blast_timeout_s)
        self.poll_interval_s = max(0.005, float(poll_interval_s))
        self._stage_functions = dict(stage_functions or {})
        unknown = set(self._stage_functions) - {x[0] for x in _STAGE_DEFS}
        if unknown:
            raise ValueError("unknown stage function(s): " + ", ".join(sorted(unknown)))

        self._lock = threading.RLock()
        self._jobs: dict[str, _Job] = {}
        self._idempotency: dict[str, tuple[str, str]] = {}
        self._queue: queue.Queue[_Job] = queue.Queue(maxsize=self.queue_capacity)
        self._stop = threading.Event()
        self._worker = threading.Thread(target=self._worker_loop,
                                        name="oligoforge-design-worker", daemon=True)
        self._worker.start()

    # ---- public lifecycle and lookup API ---------------------------------

    def submit(self, payload: Mapping[str, Any],
               idempotency_key: Optional[str] = None) -> dict[str, Any]:
        """Queue a design and return immediately with a sanitized snapshot."""
        if self._stop.is_set():
            raise RuntimeError("job manager is shut down")
        if not isinstance(payload, Mapping):
            raise TypeError("payload must be a mapping")
        private_payload = deepcopy(dict(payload))
        fingerprint = self._fingerprint("design", private_payload)
        idem_hash = self._idempotency_hash(idempotency_key)
        now_wall, now_mono = time.time(), time.monotonic()
        with self._lock:
            self._cleanup_locked(now_mono)
            existing = self._idempotent_existing_locked(idem_hash, fingerprint)
            if existing is not None:
                return self._snapshot_locked(existing)
            job = self._new_job(private_payload, fingerprint, idem_hash, "design",
                                now_wall, now_mono)
            self._jobs[job.job_id] = job
            if idem_hash:
                self._idempotency[idem_hash] = (fingerprint, job.job_id)
            try:
                self._queue.put_nowait(job)
            except queue.Full:
                self._remove_locked(job)
                raise QueueFull("automatic-design queue is full")
            return self._snapshot_locked(job)

    def retry_blast(self, job_id: str, idempotency_key: Optional[str] = None,
                    **overrides: Any) -> Optional[dict[str, Any]]:
        """Queue BLAST alone from a retained primary result.

        A new capability ID is returned.  Resolve/fetch/design/enrich are marked
        skipped in that job, making it mechanically impossible for this retry to
        rerun them.
        """
        if self._stop.is_set():
            raise RuntimeError("job manager is shut down")
        allowed = {"blast_mode", "blast_db", "blast_db_path", "organism"}
        unknown = set(overrides) - allowed
        if unknown:
            raise TypeError("unknown BLAST override(s): " + ", ".join(sorted(unknown)))
        now_wall, now_mono = time.time(), time.monotonic()
        with self._lock:
            self._cleanup_locked(now_mono)
            source = self._jobs.get(job_id) if isinstance(job_id, str) else None
            if source is None:
                return None
            if source.primary_result is None:
                raise RetryNotAvailable("primary design is not complete")
            payload = {
                "blast_mode": source.payload.get("blast_mode", "remote"),
                "blast_db": source.payload.get("blast_db", "nt"),
                "blast_db_path": source.payload.get("blast_db_path"),
                "organism": source.payload.get("organism"),
                "email": source.payload.get("email"),
                "ncbi_key": source.payload.get("ncbi_key"),
            }
            payload.update(overrides)
            # The source capability is hashed into the fingerprint, not exposed.
            fp_payload = dict(payload, source_job=job_id)
            fingerprint = self._fingerprint("blast_retry", fp_payload)
            idem_hash = self._idempotency_hash(idempotency_key)
            existing = self._idempotent_existing_locked(idem_hash, fingerprint)
            if existing is not None:
                return self._snapshot_locked(existing)
            job = self._new_job(payload, fingerprint, idem_hash, "blast_retry",
                                now_wall, now_mono)
            for name in ("resolve_fetch", "design", "enrich"):
                job.stages[name]["status"] = "skipped"
                job.stages[name]["message"] = "reused completed primary result"
            job.primary_result = deepcopy(source.primary_result)
            job.result = deepcopy(source.primary_result)
            self._jobs[job.job_id] = job
            if idem_hash:
                self._idempotency[idem_hash] = (fingerprint, job.job_id)
            try:
                self._queue.put_nowait(job)
            except queue.Full:
                self._remove_locked(job)
                raise QueueFull("automatic-design queue is full")
            return self._snapshot_locked(job)

    def get(self, job_id: str) -> Optional[dict[str, Any]]:
        """Return a defensive public snapshot, or ``None`` for unknown/expired."""
        now = time.monotonic()
        with self._lock:
            self._cleanup_locked(now)
            job = self._jobs.get(job_id) if isinstance(job_id, str) else None
            return self._snapshot_locked(job) if job is not None else None

    def wait(self, job_id: str, timeout: Optional[float] = None) -> Optional[dict[str, Any]]:
        """Testing/embedding convenience: wait for terminal state, then snapshot."""
        with self._lock:
            self._cleanup_locked(time.monotonic())
            job = self._jobs.get(job_id) if isinstance(job_id, str) else None
            event = job.terminal_event if job is not None else None
        if event is None:
            return None
        event.wait(timeout)
        return self.get(job_id)

    def cancel(self, job_id: str) -> Optional[dict[str, Any]]:
        """Request cancellation; queued jobs become terminal synchronously."""
        now_wall, now_mono = time.time(), time.monotonic()
        with self._lock:
            self._cleanup_locked(now_mono)
            job = self._jobs.get(job_id) if isinstance(job_id, str) else None
            if job is None:
                return None
            if job.status in TERMINAL_STATUSES:
                return self._snapshot_locked(job)
            job.cancel_requested = True
            job.cancel_event.set()
            self._touch_locked(job, now_wall)
            if job.status == "queued":
                first_pending = next((name for name, stage in job.stages.items()
                                      if stage["status"] == "pending"), None)
                if first_pending is not None:
                    self._set_stage_locked(job, first_pending, "cancelled",
                                           "cancelled before stage start")
                self._skip_pending_locked(job)
                self._finish_locked(job, "cancelled", now_wall, now_mono,
                                    error={"code": "cancelled",
                                           "message": "automatic design was cancelled"})
            return self._snapshot_locked(job)

    def cleanup(self) -> int:
        """Forget terminal jobs older than their configured TTL."""
        with self._lock:
            return self._cleanup_locked(time.monotonic())

    def stats(self) -> dict[str, Any]:
        """Non-identifying operational counters suitable for a health endpoint."""
        with self._lock:
            self._cleanup_locked(time.monotonic())
            counts: dict[str, int] = {}
            for job in self._jobs.values():
                counts[job.status] = counts.get(job.status, 0) + 1
            return {
                "backend": "memory_single_process",
                "worker_count": 1,
                "queue_capacity": self.queue_capacity,
                "queue_depth": self._queue.qsize(),
                "jobs_by_status": counts,
                "terminal_ttl_seconds": self.ttl_s,
                "primary_timeout_seconds": self.primary_timeout_s,
                "blast_timeout_seconds": self.blast_timeout_s,
                "durable": False,
            }

    def shutdown(self, wait: bool = True, timeout: float = 2.0) -> None:
        """Stop accepting work and cancel queued/running jobs.

        Threads are daemons, so a non-cooperative native/network call never holds
        process shutdown hostage.
        """
        self._stop.set()
        with self._lock:
            jobs = list(self._jobs.values())
        for job in jobs:
            if job.status not in TERMINAL_STATUSES:
                self.cancel(job.job_id)
        if wait and self._worker is not threading.current_thread():
            self._worker.join(max(0.0, float(timeout)))

    close = shutdown

    # ---- worker and stage execution --------------------------------------

    def _worker_loop(self) -> None:
        while not self._stop.is_set():
            try:
                job = self._queue.get(timeout=0.1)
            except queue.Empty:
                continue
            try:
                with self._lock:
                    live = self._jobs.get(job.job_id) is job
                    runnable = live and job.status == "queued" and not job.cancel_event.is_set()
                if not runnable:
                    continue
                credentials = []
                try:
                    credentials = self._apply_job_credentials(job)
                    if job.action == "blast_retry":
                        self._run_blast_retry(job)
                    else:
                        self._run_design(job)
                except BaseException:
                    # Direct embedders can bypass Pydantic validation.  Bad
                    # input or an orchestration defect must not kill the sole
                    # worker or expose its exception detail.
                    with self._lock:
                        if job.status not in TERMINAL_STATUSES:
                            job.result = None
                            job.primary_result = None
                            self._skip_pending_locked(job)
                            self._finish_locked(job, "failed", error={
                                "code": "job_failed",
                                "message": "automatic design could not complete",
                            })
                finally:
                    self._restore_job_credentials(credentials)
            finally:
                self._queue.task_done()

    @staticmethod
    def _apply_job_credentials(job: _Job):
        """Scope process-global Biopython credentials to the sole active job.

        Biopython's Entrez interface stores contact details globally.  The
        manager deliberately has one worker, so saving/restoring them here
        prevents a later queued browser submission from changing the active
        job's contact/API-key context.
        """
        holders = (AD.N.Entrez, AD.SP.Entrez)
        old = [(holder, getattr(holder, "email", None), getattr(holder, "api_key", None))
               for holder in holders]
        email, key = job.payload.get("email"), job.payload.get("ncbi_key")
        for holder in holders:
            if email:
                holder.email = email
            if key:
                holder.api_key = key
        return old

    @staticmethod
    def _restore_job_credentials(old):
        for holder, email, key in old or []:
            holder.email = email
            holder.api_key = key

    def _run_design(self, job: _Job) -> None:
        self._start_job(job)
        deadline = time.monotonic() + self.primary_timeout_s
        p = job.payload

        outcome, context, lingering = self._invoke_stage(
            job, "resolve_fetch", self._stage_function("resolve_fetch"), deadline,
            p.get("target_query"), p.get("off_query"), min(int(p.get("n_fetch", 20)), 30))
        if outcome != "complete":
            self._abort_required(job, outcome, "resolve_fetch", lingering)
            return
        if self._domain_error(job, "resolve_fetch", context):
            return

        outcome, out, lingering = self._invoke_stage(
            job, "design", self._stage_function("design"), deadline,
            context, p.get("profile", "auto"), float(p.get("min_ident", 0.6)),
            p.get("objective", "balanced"))
        if outcome != "complete":
            self._abort_required(job, outcome, "design", lingering)
            return
        if self._domain_error(job, "design", out):
            return

        outcome, out, lingering = self._invoke_stage(
            job, "enrich", self._stage_function("enrich"), deadline,
            out, context, min_ident=float(p.get("min_ident", 0.6)),
            prefer_junction=bool(p.get("prefer_junction", False)),
            nested=bool(p.get("nested", False)),
            objective=p.get("objective", "balanced"))
        if outcome != "complete":
            self._abort_required(job, outcome, "enrich", lingering)
            return
        if self._domain_error(job, "enrich", out):
            return

        with self._lock:
            job.primary_result = deepcopy(out)
            job.result = deepcopy(out)
            self._touch_locked(job)

        if not bool(p.get("run_blast", False)):
            with self._lock:
                self._skip_stage_locked(job, "blast", "not requested")
                self._finish_locked(job, "succeeded")
            return
        self._run_optional_blast(job, p)

    def _run_blast_retry(self, job: _Job) -> None:
        self._start_job(job)
        self._run_optional_blast(job, job.payload)

    def _run_optional_blast(self, job: _Job, params: Mapping[str, Any]) -> None:
        with self._lock:
            primary = deepcopy(job.primary_result)
        if primary is None:
            with self._lock:
                self._finish_locked(job, "failed", error={
                    "code": "primary_result_missing",
                    "message": "primary design result is unavailable",
                })
            return
        deadline = time.monotonic() + self.blast_timeout_s
        outcome, result, lingering = self._invoke_stage(
            job, "blast", self._stage_function("blast"), deadline, primary,
            blast_mode=params.get("blast_mode", "remote"),
            blast_db=params.get("blast_db", "nt"),
            blast_db_path=params.get("blast_db_path"),
            organism=params.get("organism"), suppress_errors=False)
        # A custom/legacy stage may encode failure instead of raising.  Treat it
        # as optional failure while retaining the untouched primary result.
        encoded_error = (isinstance(result, Mapping) and
                         isinstance(result.get("specificity"), Mapping) and
                         bool(result["specificity"].get("error")))
        if outcome == "complete" and not encoded_error:
            with self._lock:
                job.result = deepcopy(result)
                self._finish_locked(job, "succeeded")
            return

        if outcome == "cancelled":
            warning = "Specificity analysis was cancelled; the primary design is retained."
        elif outcome == "timed_out":
            warning = "Specificity analysis timed out; the primary design is retained."
        else:
            warning = "Specificity analysis was unavailable; the primary design is retained."
            if encoded_error:
                with self._lock:
                    stage = job.stages["blast"]
                    stage["status"] = "failed"
                    stage["message"] = "optional specificity analysis was unavailable"
                    stage["finished_at"] = _utc(time.time())
                    self._touch_locked(job)
        with self._lock:
            job.result = deepcopy(job.primary_result)
            job.warnings.append(warning)
            self._finish_locked(job, "succeeded_with_warnings")
        self._drain_lingering(lingering)

    def _invoke_stage(self, job: _Job, stage_name: str, func: Callable[..., Any],
                      deadline: float, *args: Any, **kwargs: Any):
        with self._lock:
            if job.cancel_event.is_set():
                self._set_stage_locked(job, stage_name, "cancelled",
                                       "cancelled before stage start")
                return "cancelled", None, None
            self._set_stage_locked(job, stage_name, "running")

        done = threading.Event()
        box: dict[str, Any] = {}

        def call() -> None:
            try:
                box["value"] = func(*args, **kwargs)
            except BaseException as exc:  # never serialize/log native or network details
                box["exception"] = exc
            finally:
                box["ended_mono"] = time.monotonic()
                done.set()

        stage_thread = threading.Thread(target=call,
                                        name="oligoforge-stage-" + stage_name,
                                        daemon=True)
        stage_thread.start()
        outcome = "complete"
        while not done.wait(self.poll_interval_s):
            if job.cancel_event.is_set():
                outcome = "cancelled"
                break
            if time.monotonic() >= deadline:
                outcome = "timed_out"
                break
        if outcome == "complete" and box.get("ended_mono", time.monotonic()) > deadline:
            outcome = "timed_out"
        if outcome == "complete" and job.cancel_event.is_set():
            outcome = "cancelled"

        if outcome == "cancelled":
            with self._lock:
                self._set_stage_locked(job, stage_name, "cancelled",
                                       "cancellation requested")
            return outcome, None, stage_thread if stage_thread.is_alive() else None
        if outcome == "timed_out":
            with self._lock:
                self._set_stage_locked(job, stage_name, "timed_out",
                                       "stage exceeded its deadline")
            return outcome, None, stage_thread if stage_thread.is_alive() else None
        if "exception" in box:
            box.pop("exception", None)  # discard detail before any snapshot/log path
            with self._lock:
                self._set_stage_locked(job, stage_name, "failed",
                                       "stage failed")
            return "failed", None, None
        value = box.pop("value", None)
        with self._lock:
            self._set_stage_locked(job, stage_name, "complete")
        return "complete", value, None

    # ---- state helpers ----------------------------------------------------

    def _stage_function(self, name: str) -> Callable[..., Any]:
        if name in self._stage_functions:
            return self._stage_functions[name]
        return {
            "resolve_fetch": AD.resolve_and_fetch_query,
            "design": AD.design_query_corpus,
            "enrich": AD.enrich_query_design,
            "blast": AD.blast_winner,
        }[name]

    def _start_job(self, job: _Job) -> None:
        with self._lock:
            if job.status != "queued":
                return
            now = time.time()
            job.status = "running"
            job.started_wall = now
            self._touch_locked(job, now)

    def _abort_required(self, job: _Job, outcome: str, stage_name: str,
                        lingering: Optional[threading.Thread]) -> None:
        status = outcome if outcome in {"cancelled", "timed_out"} else "failed"
        messages = {
            "cancelled": ("cancelled", "automatic design was cancelled"),
            "timed_out": ("primary_timeout", "automatic design exceeded its primary deadline"),
            "failed": ("stage_failed", "automatic design could not complete"),
        }
        code, message = messages[status]
        with self._lock:
            job.result = None
            job.primary_result = None
            self._skip_pending_locked(job)
            self._finish_locked(job, status,
                                error={"code": code, "message": message})
        self._drain_lingering(lingering)

    def _domain_error(self, job: _Job, stage_name: str, value: Any) -> bool:
        if isinstance(value, Mapping) and value.get("error"):
            message = self._safe_domain_message(str(value.get("error")), job)
            with self._lock:
                self._set_stage_locked(job, stage_name, "failed",
                                       "stage could not produce a result")
                self._skip_pending_locked(job)
                job.result = None
                job.primary_result = None
                self._finish_locked(job, "failed", error={
                    "code": stage_name + "_failed",
                    "message": message,
                })
            return True
        if not isinstance(value, Mapping):
            with self._lock:
                self._set_stage_locked(job, stage_name, "failed", "invalid stage result")
                self._skip_pending_locked(job)
                self._finish_locked(job, "failed", error={
                    "code": "invalid_stage_result",
                    "message": "automatic design returned an invalid stage result",
                })
            return True
        return False

    def _new_job(self, payload: dict[str, Any], fingerprint: str,
                 idem_hash: Optional[str], action: str,
                 now_wall: float, now_mono: float) -> _Job:
        stages = {
            name: {
                "name": name,
                "required": required,
                "status": "pending",
                "started_at": None,
                "finished_at": None,
                "message": None,
            }
            for name, required in _STAGE_DEFS
        }
        job_id = None
        for _attempt in range(8):
            candidate = secrets.token_urlsafe(32)
            if candidate not in self._jobs:
                job_id = candidate
                break
        if job_id is None:  # cryptographically implausible unless the RNG is broken
            raise RuntimeError("could not allocate a unique job capability")
        return _Job(job_id=job_id, payload=payload,
                    fingerprint=fingerprint, idempotency_hash=idem_hash,
                    action=action, created_wall=now_wall, created_mono=now_mono,
                    updated_wall=now_wall, stages=stages)

    def _set_stage_locked(self, job: _Job, name: str, status: str,
                          message: Optional[str] = None) -> None:
        stage = job.stages[name]
        now = time.time()
        if status == "running" and stage["started_at"] is None:
            stage["started_at"] = _utc(now)
        if status in {"complete", "failed", "timed_out", "cancelled", "skipped"}:
            stage["finished_at"] = _utc(now)
        stage["status"] = status
        stage["message"] = message
        self._touch_locked(job, now)

    def _skip_stage_locked(self, job: _Job, name: str, message: str) -> None:
        if job.stages[name]["status"] == "pending":
            self._set_stage_locked(job, name, "skipped", message)

    def _skip_pending_locked(self, job: _Job) -> None:
        for name, stage in job.stages.items():
            if stage["status"] == "pending":
                self._skip_stage_locked(job, name, "not reached")

    def _finish_locked(self, job: _Job, status: str,
                       now_wall: Optional[float] = None,
                       now_mono: Optional[float] = None,
                       error: Optional[dict[str, str]] = None) -> None:
        if job.status in TERMINAL_STATUSES:
            return
        now_wall = time.time() if now_wall is None else now_wall
        now_mono = time.monotonic() if now_mono is None else now_mono
        job.status = status
        job.error = error
        job.finished_wall = now_wall
        job.finished_mono = now_mono
        job.expires_wall = now_wall + self.ttl_s
        # Sanitize retained terminal values *before* dropping credentials from
        # the private request.  Otherwise a faulty stage that copied a secret
        # into an innocently named result field could become visible on the
        # next poll after the redaction source had been removed.
        redactions = self._private_redactions(job.payload)
        if job.result is not None:
            job.result = self._sanitize(job.result, redactions)
        if job.primary_result is not None:
            job.primary_result = self._sanitize(job.primary_result, redactions)
        if job.error is not None:
            job.error = self._sanitize(job.error, redactions)
        # Credentials were needed only while the stage ran.  Keep the input
        # fingerprint for idempotency, but do not retain secret values for the
        # remainder of the terminal job TTL.
        job.payload.pop("email", None)
        job.payload.pop("ncbi_key", None)
        self._touch_locked(job, now_wall)
        job.terminal_event.set()

    def _touch_locked(self, job: _Job, now_wall: Optional[float] = None) -> None:
        job.updated_wall = time.time() if now_wall is None else now_wall

    def _drain_lingering(self, stage_thread: Optional[threading.Thread]) -> None:
        # Keep the worker occupied until the non-cooperative call really exits.
        # The job is already terminal and visible to pollers while this waits.
        if stage_thread is not None:
            stage_thread.join()

    # ---- idempotency, expiry, and sanitization ---------------------------

    @staticmethod
    def _fingerprint(action: str, payload: Mapping[str, Any]) -> str:
        raw = json.dumps({"action": action, "payload": payload}, sort_keys=True,
                         separators=(",", ":"), default=str)
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    @staticmethod
    def _idempotency_hash(key: Optional[str]) -> Optional[str]:
        if key is None or str(key) == "":
            return None
        value = str(key)
        if len(value) > 512:
            raise ValueError("idempotency key is too long")
        return hashlib.sha256(("oligoforge-job:" + value).encode("utf-8")).hexdigest()

    def _idempotent_existing_locked(self, idem_hash: Optional[str],
                                    fingerprint: str) -> Optional[_Job]:
        if not idem_hash:
            return None
        entry = self._idempotency.get(idem_hash)
        if entry is None:
            return None
        old_fingerprint, job_id = entry
        job = self._jobs.get(job_id)
        if job is None:
            self._idempotency.pop(idem_hash, None)
            return None
        if old_fingerprint != fingerprint:
            raise IdempotencyConflict("idempotency key belongs to different input")
        return job

    def _cleanup_locked(self, now_mono: float) -> int:
        expired = [job for job in self._jobs.values()
                   if job.finished_mono is not None and
                   now_mono - job.finished_mono >= self.ttl_s]
        for job in expired:
            self._remove_locked(job)
        return len(expired)

    def _remove_locked(self, job: _Job) -> None:
        self._jobs.pop(job.job_id, None)
        if job.idempotency_hash:
            entry = self._idempotency.get(job.idempotency_hash)
            if entry and entry[1] == job.job_id:
                self._idempotency.pop(job.idempotency_hash, None)

    def _safe_domain_message(self, message: str, job: _Job) -> str:
        text = " ".join(message.split())[:500]
        for secret in self._private_redactions(job.payload):
            text = text.replace(secret, "[redacted]")
        return text or "automatic design could not produce a result"

    def _private_redactions(self, payload: Mapping[str, Any]) -> tuple[str, ...]:
        values: list[str] = []
        for key, value in payload.items():
            name = _key_name(key)
            if value is None:
                continue
            if name in _DROP_KEYS or "password" in name or "secret" in name or "token" in name:
                if isinstance(value, str) and value:
                    values.append(value)
            elif isinstance(value, str) and _looks_like_sequence(value):
                values.append(value)
        return tuple(sorted(set(values), key=len, reverse=True))

    def _sanitize(self, value: Any, redactions: tuple[str, ...]) -> Any:
        if isinstance(value, Mapping):
            out = {}
            for key, child in value.items():
                name = _key_name(key)
                if (name in _DROP_KEYS or "password" in name or "secret" in name or
                        "token" in name or name.endswith("apikey")):
                    continue
                out[str(key)] = self._sanitize(child, redactions)
            return out
        if isinstance(value, (list, tuple)):
            return [self._sanitize(x, redactions) for x in value]
        if isinstance(value, str):
            text = value
            for secret in redactions:
                text = text.replace(secret, "[redacted]")
            return text
        if value is None or isinstance(value, (bool, int, float)):
            return value
        # Public snapshots must remain JSON serializable and must not call repr()
        # on arbitrary scientific/native objects that could embed private data.
        return "[unavailable]"

    def _snapshot_locked(self, job: _Job) -> dict[str, Any]:
        redactions = self._private_redactions(job.payload)
        snapshot = {
            "id": job.job_id,
            "job_id": job.job_id,
            "status": job.status,
            "created_at": _utc(job.created_wall),
            "updated_at": _utc(job.updated_wall),
            "started_at": _utc(job.started_wall),
            "finished_at": _utc(job.finished_wall),
            "expires_at": _utc(job.expires_wall),
            "cancel_requested": bool(job.cancel_requested),
            "stages": [deepcopy(job.stages[name]) for name, _required in _STAGE_DEFS],
            "warnings": list(job.warnings),
            "error": deepcopy(job.error),
            "result": deepcopy(job.result),
        }
        return self._sanitize(snapshot, redactions)


__all__ = [
    "DesignJobManager", "QueueFull", "IdempotencyConflict",
    "RetryNotAvailable", "TERMINAL_STATUSES",
]
