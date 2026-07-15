"""Deterministic scientific provenance for design and ranking runs.

The manifest deliberately excludes wall-clock time, host names and local paths so
identical scientific inputs under the same software/model configuration create
the same run identifier on another machine. Runtime/library versions are part of
the identity because they can change native numerical behaviour.  External database dates/versions are recorded only when the caller
actually knows them; absence is represented explicitly rather than guessed.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from importlib import metadata
from typing import Any, Dict, Iterable, Mapping, Optional

from . import __version__ as APP_VERSION
from . import thermo as T

PROVENANCE_SCHEMA_VERSION = "1.0.0"
THERMODYNAMIC_MODEL_VERSION = "oligoforge-nn-primer3-v1.33"
SPECIFICITY_MODEL_VERSION = "paired-all-site-ambiguity-conservative-v1.33"


def _package_version(name: str) -> Optional[str]:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None
    except Exception:
        return None


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def sha256_value(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def sequence_hashes(sequences: Iterable[str], prefix: str = "sequence") -> Dict[str, Any]:
    rows = [str(x or "").upper() for x in sequences]
    return {
        "count": len(rows),
        "individual_sha256": [hashlib.sha256(x.encode("utf-8")).hexdigest() for x in rows],
        "corpus_sha256": hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest(),
        "label": prefix,
    }


def software_versions() -> Dict[str, Any]:
    return {
        "oligoforge": APP_VERSION,
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "primer3_py": _package_version("primer3-py"),
        "viennarna": _package_version("ViennaRNA"),
        "biopython": _package_version("biopython"),
        "fastapi": _package_version("fastapi"),
        "pydantic": _package_version("pydantic"),
    }


def build_manifest(*, ranker_version: str, ranking_schema: str,
                   profile_version: str, objective: Optional[str],
                   candidate_limits: Optional[Mapping[str, Any]] = None,
                   input_hashes: Optional[Mapping[str, Any]] = None,
                   external_databases: Optional[Mapping[str, Any]] = None,
                   warnings: Optional[Iterable[str]] = None,
                   fallbacks: Optional[Iterable[str]] = None,
                   constraints: Optional[Mapping[str, Any]] = None,
                   search_version: Optional[str] = None,
                   retention_version: Optional[str] = None,
                   manual_design_version: Optional[str] = None,
                   rescue_version: Optional[str] = None,
                   random_seed: Optional[int] = None) -> Dict[str, Any]:
    """Build a deterministic, self-hashing run manifest.

    ``external_databases`` should contain caller-supplied accession/database
    versions or retrieval dates.  The function never substitutes the current date.
    """
    known_db = dict(external_databases or {})
    db_state = "declared" if known_db else "not_recorded"
    body = {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "application_version": APP_VERSION,
        "ranker_version": ranker_version,
        "ranking_schema": ranking_schema,
        "scoring_profile_version": profile_version,
        "search_version": search_version,
        "retention_version": retention_version,
        "manual_design_version": manual_design_version,
        "assay_rescue_version": rescue_version,
        "objective": objective,
        "software_versions": software_versions(),
        "scientific_models": {
            "thermodynamics": THERMODYNAMIC_MODEL_VERSION,
            "specificity": SPECIFICITY_MODEL_VERSION,
            "reaction_condition_snapshot": {
                "mv_conc_mM": T._snapshot()[0],
                "dv_conc_mM": T._snapshot()[1],
                "dntp_conc_mM": T._snapshot()[2],
                "total_oligo_conc_nM": T._snapshot()[3],
                "anneal_c": T._snapshot()[4],
            },
        },
        "candidate_limits": dict(candidate_limits or {}),
        "constraints": dict(constraints or {}),
        "input_hashes": dict(input_hashes or {}),
        "external_databases": known_db,
        "external_database_state": db_state,
        "warnings": sorted({str(x) for x in (warnings or []) if str(x).strip()}),
        "fallbacks": sorted({str(x) for x in (fallbacks or []) if str(x).strip()}),
        "deterministic": random_seed is None,
        "random_seed": random_seed,
        "timestamp_policy": "excluded_from_run_id",
        "host_policy": "host names, local paths and wall-clock timestamps excluded",
    }
    body["run_id"] = "ofrun_" + sha256_value(body)[:24]
    body["manifest_sha256"] = sha256_value(body)
    return body


def verify_manifest(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    """Verify both the deterministic run identifier and manifest digest.

    Returns a diagnostic record instead of raising so imported reports can be
    audited without trusting their contents.  Unknown extra fields are covered by
    the digest and therefore invalidate a modified manifest.
    """
    supplied = dict(manifest or {})
    stored_digest = supplied.pop("manifest_sha256", None)
    calculated_digest = sha256_value(supplied)
    stored_run_id = supplied.pop("run_id", None)
    calculated_run_id = "ofrun_" + sha256_value(supplied)[:24]
    digest_valid = bool(stored_digest) and stored_digest == calculated_digest
    run_id_valid = bool(stored_run_id) and stored_run_id == calculated_run_id
    return {
        "valid": bool(digest_valid and run_id_valid),
        "digest_valid": digest_valid,
        "run_id_valid": run_id_valid,
        "stored_manifest_sha256": stored_digest,
        "calculated_manifest_sha256": calculated_digest,
        "stored_run_id": stored_run_id,
        "calculated_run_id": calculated_run_id,
    }
