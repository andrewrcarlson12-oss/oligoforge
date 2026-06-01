"""Thermodynamics for oligo QC at qPCR salt conditions.

Wraps primer3-py (SantaLucia 1998 nearest-neighbor + Owczarzy salt correction)
so every value matches what IDT OligoAnalyzer reports to within ~1-2 C absolute,
and tracks it exactly for relative comparisons. Running locally means these are
the REAL primer3 numbers, not a browser approximation.
"""
import primer3
import itertools
import threading
from functools import lru_cache

_P3 = threading.Lock()        # primer3's C core is not safe under the sync-endpoint threadpool
_P3_MAXLEN = 60               # primer3 thermo alignment refuses (raises) on sequences > 60 nt

# Typical TaqMan/SYBR master-mix conditions. Edit here to match your kit.
COND = dict(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=200.0)
# full IUPAC complement so a degenerate primer (e.g. the W in Plas_F) never crashes
_COMP = {"A":"T","T":"A","G":"C","C":"G","N":"N","U":"A",
         "R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K",
         "B":"V","V":"B","D":"H","H":"D"}
# a representative concrete base per ambiguity code: primer3 only accepts ACGT
_RESOLVE = {"R":"A","Y":"C","S":"G","W":"A","K":"G","M":"A",
            "B":"C","D":"A","H":"A","V":"A","N":"A","U":"T"}


_IUPAC = set("ACGTRYSWKMBDHVN")


def clean_seq(raw):
    """Normalize a user-pasted nucleotide sequence for safe analysis.

    Returns (clean, notes, error). IUPAC ambiguity codes are preserved. A leading
    FASTA header is dropped; whitespace and numbering from a pasted formatted block
    are removed; RNA (U) is converted to DNA (T). Any other non-IUPAC character is
    reported as an error rather than silently dropped, so a corrupt paste can never
    poison a result."""
    notes = []
    if not raw or not raw.strip():
        return "", notes, "empty sequence"
    raw_lines = raw.splitlines()
    body = [ln for ln in raw_lines if not ln.lstrip().startswith(">")]
    if len(body) != len(raw_lines):
        notes.append("removed a FASTA header line")
    joined = "".join(body)
    no_fmt = "".join(ch for ch in joined if not (ch.isspace() or ch.isdigit()))
    if len(no_fmt) != len(joined):
        notes.append("removed spacing/numbering from a formatted paste")
    up = no_fmt.upper()
    bad = sorted({ch for ch in up if ch not in _IUPAC and ch != "U"})
    if bad:
        return "", notes, "invalid character(s): " + " ".join(bad) + " (expected DNA / IUPAC codes)"
    if "U" in up:
        up = up.replace("U", "T")
        notes.append("converted RNA (U) to DNA (T)")
    if not up:
        return "", notes, "no sequence left after cleaning"
    return up, notes, None


def revcomp(seq):
    return "".join(_COMP.get(b, "N") for b in reversed(seq.upper()))


def _resolve(seq):
    return "".join(_RESOLVE.get(b, b) for b in seq.upper())


def has_degenerate(seq):
    return any(b not in "ACGTN" for b in seq.upper())


def gc_percent(seq):
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def amplicon_tm(seq):
    """Estimated melt Tm of a PCR product. This is the empirical salt-adjusted
    formula Tm = 81.5 + 16.6*log10[Na+] + 0.41*%GC - 600/length, NOT a
    nearest-neighbor calc — a rough predictor for the SYBR melt-curve peak, to be
    confirmed against the observed melt. Tracks the reaction monovalent salt set
    via set_conditions()."""
    import math
    seq = seq.upper().replace("U", "T")
    n = len(seq)
    if n < 1:
        return 0.0
    na = max(COND.get("mv_conc", 50.0), 1.0) / 1000.0
    return 81.5 + 16.6 * math.log10(na) + 0.41 * gc_percent(seq) - 600.0 / n


@lru_cache(maxsize=8192)
def tm(seq):
    with _P3:
        return primer3.calc_tm(_resolve(seq), tm_method="santalucia",
                               salt_corrections_method="owczarzy", **COND)


_IUPAC_SETS = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC",
               "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
               "V": "ACG", "N": "ACGT"}


def set_conditions(mv_conc=None, dv_conc=None, dntp_conc=None, dna_conc=None):
    """Update, in place, the reaction salt used by every Tm and structure calc.

    Only the arguments you pass are changed; the rest keep their current value.
    Lets the Tm match your actual master mix (free Mg2+, monovalent, dNTP, oligo nM)
    instead of a generic default. Session-global for this local single-user server.
    """
    if mv_conc is not None:
        COND["mv_conc"] = float(mv_conc)
    if dv_conc is not None:
        COND["dv_conc"] = float(dv_conc)
    if dntp_conc is not None:
        COND["dntp_conc"] = float(dntp_conc)
    if dna_conc is not None:
        COND["dna_conc"] = float(dna_conc)
    tm.cache_clear(); hairpin.cache_clear(); self_dimer.cache_clear(); hetero_dimer.cache_clear()
    return dict(COND)


def tm_range(seq, cap=64):
    """Tm across every resolution of a degenerate oligo.

    A degenerate base (e.g. W = A/T) has two different Tm values depending on which
    base is present in a given template copy. Reporting one number hides that spread;
    this enumerates the resolutions (capped) and returns the min/max so the real range
    is visible. Non-degenerate oligos return min == max.
    """
    s = seq.upper().replace("U", "T")
    degen = [(i, _IUPAC_SETS.get(b, b)) for i, b in enumerate(s) if len(_IUPAC_SETS.get(b, b)) > 1]
    if not degen:
        t = round(tm(s), 1)
        return dict(min=t, max=t, n=1, degenerate=False, capped=False)
    ncomb = 1
    for _, opts in degen:
        ncomb *= len(opts)
    if ncomb > cap:
        t = round(tm(s), 1)
        return dict(min=t, max=t, n=ncomb, degenerate=True, capped=True)
    idxs = [i for i, _ in degen]
    base = list(s)
    tms = []
    for combo in itertools.product(*[opts for _, opts in degen]):
        for j, b in zip(idxs, combo):
            base[j] = b
        tms.append(tm("".join(base)))
    return dict(min=round(min(tms), 1), max=round(max(tms), 1), n=ncomb, degenerate=True, capped=False)


@lru_cache(maxsize=8192)
def hairpin(seq):
    s = _resolve(seq)
    if len(s) > _P3_MAXLEN:              # not an oligo; primer3 would raise. Neutral (no hairpin).
        return 0.0, 0.0
    with _P3:
        r = primer3.calc_hairpin(s, **COND)
    return r.dg / 1000.0, r.tm           # (dG kcal/mol, melting Tm C)


@lru_cache(maxsize=8192)
def self_dimer(seq):
    s = _resolve(seq)
    if len(s) > _P3_MAXLEN:
        return 0.0
    with _P3:
        return primer3.calc_homodimer(s, **COND).dg / 1000.0


@lru_cache(maxsize=8192)
def hetero_dimer(a, b):
    sa, sb = _resolve(a), _resolve(b)
    if len(sa) > _P3_MAXLEN or len(sb) > _P3_MAXLEN:
        return 0.0
    with _P3:
        return primer3.calc_heterodimer(sa, sb, **COND).dg / 1000.0


def max_run(seq, base=None):
    seq = seq.upper()
    if base:
        best = cur = 0
        for b in seq:
            cur = cur + 1 if b == base else 0
            best = max(best, cur)
        return best
    return max(max_run(seq, b) for b in "ATGC")


def last5_gc(seq):
    tail = seq.upper()[-5:]
    return tail.count("G") + tail.count("C")


# ---- LNA-aware Tm (Affinity Plus / "+N" notation) ----
def strip_lna(seq):
    """'+A' LNA notation -> (dna_backbone_seq, n_lna, [1-based LNA positions])."""
    s = seq.upper().strip()
    dna, pos, idx, i = [], [], 0, 0
    while i < len(s):
        if s[i] == "+" and i + 1 < len(s):
            dna.append(s[i + 1]); idx += 1; pos.append(idx); i += 2
        else:
            dna.append(s[i]); idx += 1; i += 1
    return "".join(dna), len(pos), pos


def suggest_lna(seq, snp_pos=None, max_lna=None, spacing=3):
    """Heuristic LNA (+) placement on a DNA probe -> '+N' notation, positions, count, Tm range.

    Two modes:
      - Tm boost (default): space LNAs ~every `spacing` bases, nudged onto G/C for the biggest
        gain, clear of the ends, never 3 in a row, capped.
      - Genotyping (snp_pos, 0-based): LNA centred on the polymorphic base (+ its neighbours),
        which maximises the match-vs-mismatch Tm gap an allele-specific probe relies on.

    This is a STARTING point only. Optimal LNA placement is strongly context-dependent and IDT's
    Affinity Plus / OligoAnalyzer model is authoritative; always confirm there, and for allele
    discrimination validate the match-vs-mismatch dTm empirically.
    """
    s = "".join(c for c in seq.upper().replace("U", "T") if c in "ACGT")
    n = len(s)
    if n < 6:
        return dict(error="probe too short for LNA placement (need >= 6 nt)")
    pos = set()
    if snp_pos is not None and 0 <= snp_pos < n:
        for p in (snp_pos - 1, snp_pos, snp_pos + 1):
            if 0 <= p < n:
                pos.add(p)
        mode = "genotyping (LNA on the discriminating base)"
    else:
        cap = max_lna if max_lna else min(8, max(3, n // 4))
        i = 2
        while i < n - 2 and len(pos) < cap:
            cand = i
            for off in (0, 1, -1):
                j = i + off
                if 2 <= j < n - 2 and s[j] in "GC":
                    cand = j; break
            pos.add(cand)
            i = cand + spacing
        mode = "Tm boost (LNA spaced on G/C)"
    pos = sorted(pos)
    notation = "".join((("+" + b) if i in pos else b) for i, b in enumerate(s))
    out = dict(dna=s, lna_pos=[p + 1 for p in pos], n_lna=len(pos), notation=notation, mode=mode)
    out.update(tm_lna(notation))
    out["note"] = ("Heuristic starting placement — optimal LNA positions are context-dependent. "
                   "Confirm in IDT OligoAnalyzer (Affinity Plus). For allele/SNP discrimination, "
                   "place an LNA on the polymorphic base and validate match-vs-mismatch dTm.")
    return out


def tm_lna(seq, n_lna=None, per_lna_low=2.0, per_lna_high=8.0):
    """DNA-backbone Tm plus an honest effective-Tm RANGE for LNA oligos.

    LNA Tm is not modeled exactly here. IDT OligoAnalyzer's calibrated LNA model
    is authoritative; this just keeps the tool from reporting the bare-DNA Tm as
    if it were the probe Tm. The 2-8 C per-LNA span reflects context dependence
    (McTigue 2004)."""
    dna, n, pos = strip_lna(seq)
    if n_lna is not None:
        n = n_lna
    base = tm(dna)
    out = dict(dna_backbone_tm=round(base, 1), n_lna=n, lna_pos=pos)
    if n:
        out.update(est_tm_low=round(base + per_lna_low * n, 1),
                   est_tm_high=round(base + per_lna_high * n, 1),
                   note="LNA strongly stabilizes; effective Tm is a range, not a point. "
                        "Confirm in IDT OligoAnalyzer (calibrated LNA model).")
    return out
