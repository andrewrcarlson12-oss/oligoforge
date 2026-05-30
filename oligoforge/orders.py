"""Build IDT-compatible order files.

Oligo bulk format (verified against IDT requirements): CSV columns
Name, Sequence, Scale, Purification. gBlocks are a separate product, emitted as
FASTA. ZEN probes are best ordered through IDT's PrimeTime qPCR-probe entry
(double-quench config + BOGO); the CSV row carries the modified sequence string
for reference using IDT modification codes.
"""
import csv, io

SCALE = {"primer": "25nm", "probe_zen": "100nm", "probe_lna": "100nm"}
PURIF = {"primer": "STD", "probe_zen": "HPLC", "probe_lna": "HPLC"}


def _wrap(kind, seq):
    if kind == "probe_zen":
        # 5' 6-FAM, 3' Iowa Black FQ. Internal ZEN is added by the PrimeTime
        # double-quench config at order time.
        return f"/56-FAM/{seq}/3IABkFQ/"
    return seq


def oligo_csv(entries):
    """entries: [{name, seq, kind}] -> CSV text (Name,Sequence,Scale,Purification)."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["Name", "Sequence", "Scale", "Purification"])
    for e in entries:
        k = e.get("kind", "primer")
        w.writerow([e["name"], _wrap(k, e["seq"].upper()),
                    SCALE.get(k, "25nm"), PURIF.get(k, "STD")])
    return buf.getvalue()


def gblock_fasta(blocks):
    """blocks: [{name, seq}] -> FASTA text for the gBlocks/eBlocks entry."""
    return "\n".join(f">{b['name']}\n{b['seq'].upper()}" for b in blocks)
