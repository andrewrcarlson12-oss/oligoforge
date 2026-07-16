"""Offline command-line workflow for OligoForge Assurance.

No command performs network retrieval.  Input and output paths are explicit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .assurance import (build_assaysbom, build_snapshot, snapshot_delta, scan_drift,
                        generate_ofvrs, build_evidence_package)


def _read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path, value):
    text = json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def build_parser():
    p = argparse.ArgumentParser(prog="oligoforge-assurance",
                                description="Offline, deterministic assay-evidence lifecycle tools")
    sub = p.add_subparsers(dest="command", required=True)
    sb = sub.add_parser("build-assaysbom"); sb.add_argument("input"); sb.add_argument("output")
    sn = sub.add_parser("build-snapshot"); sn.add_argument("fasta"); sn.add_argument("output")
    sn.add_argument("--name", default="Sequence snapshot"); sn.add_argument("--role", choices=["target", "off_target", "background"], default="target")
    sn.add_argument("--metadata"); sn.add_argument("--source-label")
    de = sub.add_parser("delta"); de.add_argument("baseline"); de.add_argument("followup"); de.add_argument("output")
    ds = sub.add_parser("drift-scan"); ds.add_argument("assaysbom"); ds.add_argument("baseline_target")
    ds.add_argument("current_target"); ds.add_argument("output"); ds.add_argument("--baseline-offtarget"); ds.add_argument("--current-offtarget")
    vr = sub.add_parser("ofvr"); vr.add_argument("scan"); vr.add_argument("output"); vr.add_argument("--year", type=int)
    ep = sub.add_parser("package"); ep.add_argument("assaysbom"); ep.add_argument("output")
    ep.add_argument("--snapshots", nargs="*", default=[]); ep.add_argument("--deltas", nargs="*", default=[])
    ep.add_argument("--scans", nargs="*", default=[]); ep.add_argument("--vulnerabilities", nargs="*", default=[])
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "build-assaysbom":
        value = build_assaysbom(_read_json(args.input))
    elif args.command == "build-snapshot":
        fasta = Path(args.fasta).read_bytes()
        metadata = Path(args.metadata).read_text(encoding="utf-8") if args.metadata else None
        value = build_snapshot(fasta, name=args.name, role=args.role, metadata=metadata,
                               source={"label": args.source_label} if args.source_label else {})
    elif args.command == "delta":
        value = snapshot_delta(_read_json(args.baseline), _read_json(args.followup))
    elif args.command == "drift-scan":
        bo = _read_json(args.baseline_offtarget) if args.baseline_offtarget else None
        co = _read_json(args.current_offtarget) if args.current_offtarget else None
        value = scan_drift(_read_json(args.assaysbom), _read_json(args.baseline_target),
                           _read_json(args.current_target), baseline_offtarget=bo, current_offtarget=co)
    elif args.command == "ofvr":
        value = generate_ofvrs(_read_json(args.scan), issuance_year=args.year)
    elif args.command == "package":
        load_many = lambda paths: [_read_json(x) for x in paths]
        value = build_evidence_package(assaysbom=_read_json(args.assaysbom),
                                       snapshots=load_many(args.snapshots), deltas=load_many(args.deltas),
                                       drift_scans=load_many(args.scans), vulnerabilities=load_many(args.vulnerabilities))
    else:  # pragma: no cover
        raise AssertionError(args.command)
    _write_json(args.output, value)
    print(json.dumps({"ok": True, "command": args.command, "output": str(args.output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
