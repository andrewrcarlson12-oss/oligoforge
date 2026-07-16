"""Shared deterministic helpers for the offline Assurance records.

Assurance objects are content addressed.  Wall-clock time, host names, local
paths, and mutable server state are deliberately excluded from their identity.
"""
from __future__ import annotations

import copy
import json
from typing import Any, Iterable, Mapping

from .. import provenance as PROV


ASSURANCE_MODEL_VERSION = "oligoforge-assurance-offline-v1"


def canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      default=str)


def sha256_value(value: Any) -> str:
    return PROV.sha256_value(value)


def content_address(body: Mapping[str, Any], prefix: str, hash_field: str) -> dict:
    """Return a defensive copy with deterministic identifier and self-hash.

    ``body`` must not already contain the derived identifier or hash.  The ID is
    derived before the final digest, matching the existing run-manifest policy.
    """
    out = copy.deepcopy(dict(body))
    ident_field = {
        "ofsbom": "sbom_id", "ofsnap": "snapshot_id", "ofdelta": "delta_id",
        "ofscan": "scan_id", "ofpkg": "package_id",
    }.get(prefix, prefix + "_id")
    out.pop(ident_field, None)
    out.pop(hash_field, None)
    out[ident_field] = prefix + "_" + sha256_value(out)[:24]
    out[hash_field] = sha256_value(out)
    return out


def verify_content_address(record: Mapping[str, Any], prefix: str, hash_field: str) -> dict:
    supplied = copy.deepcopy(dict(record or {}))
    ident_field = {
        "ofsbom": "sbom_id", "ofsnap": "snapshot_id", "ofdelta": "delta_id",
        "ofscan": "scan_id", "ofpkg": "package_id",
    }.get(prefix, prefix + "_id")
    stored_hash = supplied.pop(hash_field, None)
    stored_id = supplied.pop(ident_field, None)
    expected_id = prefix + "_" + sha256_value(supplied)[:24]
    supplied[ident_field] = stored_id
    expected_hash = sha256_value(supplied)
    return {
        "valid": bool(stored_id == expected_id and stored_hash == expected_hash),
        "id_valid": stored_id == expected_id,
        "hash_valid": stored_hash == expected_hash,
        "stored_id": stored_id,
        "expected_id": expected_id,
        "stored_sha256": stored_hash,
        "expected_sha256": expected_hash,
    }


def require_keys(value: Mapping[str, Any], keys: Iterable[str], area: str) -> None:
    missing = [key for key in keys if value.get(key) in (None, "", [])]
    if missing:
        raise ValueError("%s requires %s" % (area, ", ".join(missing)))

