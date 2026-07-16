"""Bounded, offline, immutable sequence snapshots and deterministic deltas."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re

from ..provenance import sha256_value


SNAPSHOT_VERSION = "1.0.0"
SNAPSHOT_SCHEMA = "oligoforge-sequence-snapshot/v1"
MAX_COMPRESSED_BYTES = 20_000_000
MAX_UNCOMPRESSED_BYTES = 100_000_000
MAX_RECORDS = 10_000
MAX_RECORD_LENGTH = 500_000
MAX_TOTAL_BASES = 10_000_000
_ALPHABET = set("ACGTRYSWKMBDHVN")


def _decode_fasta(payload):
    if isinstance(payload, bytes):
        if len(payload) > MAX_COMPRESSED_BYTES:
            raise ValueError("snapshot input exceeds the compressed-byte limit")
        if payload[:2] == b"\x1f\x8b":
            try:
                data = gzip.GzipFile(fileobj=io.BytesIO(payload)).read(MAX_UNCOMPRESSED_BYTES + 1)
            except (OSError, EOFError) as exc:
                raise ValueError("invalid FASTA.GZ input") from exc
            if len(data) > MAX_UNCOMPRESSED_BYTES:
                raise ValueError("decompressed FASTA exceeds the safety limit")
            return data.decode("utf-8", errors="replace")
        return payload.decode("utf-8", errors="replace")
    text = str(payload or "")
    if len(text.encode("utf-8")) > MAX_UNCOMPRESSED_BYTES:
        raise ValueError("FASTA input exceeds the safety limit")
    return text


def _metadata_table(text):
    if not text:
        return {}
    sample = str(text)[:4096]
    delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
    reader = csv.DictReader(io.StringIO(str(text)), delimiter=delimiter)
    out = {}
    for row in reader:
        rid = str(row.get("record_id") or row.get("id") or row.get("accession") or "").strip()
        if rid:
            out[rid] = {str(k): str(v or "") for k, v in row.items() if k is not None}
    return out


def _parse_records(text):
    records, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
                if len(records) > MAX_RECORDS:
                    raise ValueError("snapshot exceeds %d records" % MAX_RECORDS)
            header, seq = line[1:].strip(), []
        elif header is None:
            # Plain single-sequence input is supported with a deterministic ID.
            header, seq = "record-1", [line]
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
        if len(records) > MAX_RECORDS:
            raise ValueError("snapshot exceeds %d records" % MAX_RECORDS)
    return records


def build_snapshot(fasta, *, name="Sequence snapshot", source=None, metadata=None,
                   role="target", baseline_snapshot_id=None, retrieval=None):
    role = str(role).lower().replace("-", "_")
    if role not in {"target", "off_target", "background"}:
        raise ValueError("snapshot role must be target, off_target, or background")
    text = _decode_fasta(fasta)
    raw_hash = hashlib.sha256(text.encode("utf-8")).hexdigest()
    meta = _metadata_table(metadata)
    accepted, rejected, unique, seen_hash = [], [], [], {}
    total = 0
    for index, (header, raw_seq) in enumerate(_parse_records(text), 1):
        rid = (header.split(None, 1)[0] if header else "record-%d" % index)[:200]
        sequence = re.sub(r"\s+", "", raw_seq.upper().replace("U", "T"))
        reasons = []
        if not sequence:
            reasons.append("empty_sequence")
        invalid = sorted(set(sequence) - _ALPHABET)
        if invalid:
            reasons.append("invalid_alphabet")
        if len(sequence) > MAX_RECORD_LENGTH:
            reasons.append("record_too_long")
        total += len(sequence)
        if total > MAX_TOTAL_BASES:
            raise ValueError("snapshot exceeds %d total bases" % MAX_TOTAL_BASES)
        if reasons:
            rejected.append({"record_id": rid, "input_index": index, "reasons": reasons,
                             "length": len(sequence)})
            continue
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        row = {"record_id": rid, "header": header[:1000], "sequence": sequence,
               "sequence_sha256": digest, "length": len(sequence), "metadata": meta.get(rid, {})}
        if digest in seen_hash:
            row["disposition"] = "accepted_duplicate"
            row["duplicate_of"] = seen_hash[digest]
        else:
            row["disposition"] = "accepted_unique"
            seen_hash[digest] = rid
            unique.append(row)
        accepted.append(row)
    if not accepted:
        raise ValueError("snapshot contains no accepted sequence records")
    group_counts = {}
    for row in unique:
        group = str(row["metadata"].get("group") or row["metadata"].get("lineage") or "ungrouped")
        group_counts[group] = group_counts.get(group, 0) + 1
    body = {
        "schema_version": SNAPSHOT_SCHEMA, "snapshot_version": SNAPSHOT_VERSION,
        "name": str(name), "role": role, "source": dict(source or {}),
        "retrieval": dict(retrieval or {}), "input_sha256": raw_hash,
        "baseline_snapshot_id": baseline_snapshot_id,
        "accepted_records": accepted, "rejected_records": rejected,
        "unique_records": unique,
        "metrics": {"raw_record_count": len(accepted), "unique_sequence_count": len(unique),
                    "rejected_record_count": len(rejected), "raw_total_bases": sum(x["length"] for x in accepted),
                    "unique_total_bases": sum(x["length"] for x in unique), "exact_duplicate_count": len(accepted) - len(unique),
                    "group_counts": group_counts},
        "weighting_note": "Raw, unique-haplotype and metadata-group counts are reported separately; no sampling representativeness is assumed.",
        "immutable": True,
    }
    body["content_sha256"] = sha256_value(body)
    body["snapshot_id"] = "ofsnap_" + body["content_sha256"][:24]
    return body


def validate_snapshot(value):
    if not isinstance(value, dict) or value.get("schema_version") != SNAPSHOT_SCHEMA:
        return {"valid": False, "errors": ["unsupported snapshot schema"]}
    body = {k: v for k, v in value.items() if k not in {"snapshot_id", "content_sha256"}}
    digest = sha256_value(body)
    valid = value.get("content_sha256") == digest and value.get("snapshot_id") == "ofsnap_" + digest[:24]
    return {"valid": valid, "errors": [] if valid else ["snapshot content hash or identifier does not verify"],
            "calculated_content_sha256": digest}


def snapshot_delta(baseline, followup):
    vb, vf = validate_snapshot(baseline), validate_snapshot(followup)
    if not vb["valid"] or not vf["valid"]:
        raise ValueError("both snapshots must have valid immutable hashes")
    if baseline.get("role") != followup.get("role"):
        raise ValueError("snapshot roles differ")
    b = {x["sequence_sha256"]: x for x in baseline.get("unique_records", [])}
    f = {x["sequence_sha256"]: x for x in followup.get("unique_records", [])}
    added, removed, unchanged = sorted(set(f) - set(b)), sorted(set(b) - set(f)), sorted(set(b) & set(f))
    body = {
        "schema_version": "oligoforge-snapshot-delta/v1",
        "baseline_snapshot_id": baseline["snapshot_id"], "followup_snapshot_id": followup["snapshot_id"],
        "role": baseline["role"],
        "added": [{"sequence_sha256": h, "record_id": f[h]["record_id"], "length": f[h]["length"]} for h in added],
        "removed": [{"sequence_sha256": h, "record_id": b[h]["record_id"], "length": b[h]["length"]} for h in removed],
        "unchanged": [{"sequence_sha256": h, "baseline_record_id": b[h]["record_id"],
                       "followup_record_id": f[h]["record_id"]} for h in unchanged],
        "counts": {"added": len(added), "removed": len(removed), "unchanged": len(unchanged)},
        "incremental_equivalent_to_full_set_difference": True,
    }
    body["delta_sha256"] = sha256_value(body)
    return body
