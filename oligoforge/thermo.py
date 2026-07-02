"""Thermodynamics for oligo QC at qPCR salt conditions.

Wraps primer3-py (SantaLucia 1998 nearest-neighbor + Owczarzy salt correction)
so every value matches what IDT OligoAnalyzer reports to within ~1-2 C absolute,
and tracks it exactly for relative comparisons. Running locally means these are
the REAL primer3 numbers, not a browser approximation.
"""
import primer3
from Bio.SeqUtils import MeltingTemp as _mt
import itertools
import threading
import re
from functools import lru_cache

_P3 = threading.Lock()        # primer3's C core is not safe under the sync-endpoint threadpool
_P3_MAXLEN = 60               # primer3 thermo alignment refuses (raises) on sequences > 60 nt

# Typical TaqMan/SYBR master-mix conditions. Edit here to match your kit.
COND = dict(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=200.0)
# qPCR annealing/extension temperature. Secondary structure is evaluated HERE (where priming
# actually happens) in addition to 37 C: a hairpin or dimer with a strong 37 C dG can be fully
# melted by ~60 C, so the 37 C value alone over-calls structure. 37 C is kept because it is the
# IDT OligoAnalyzer-comparable convention. Not part of COND (which is splatted into primer3 Tm
# calls that do not accept a temp_c kwarg); set via set_conditions(anneal_c=...).
ANNEAL_C = 60.0
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
    no_fmt = "".join(ch for ch in joined if not (ch.isspace() or ch.isdigit() or ch in "-."))
    if len(no_fmt) != len(joined):
        notes.append("removed spacing / numbering / alignment gaps from a formatted paste")
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


_MOD_BLOCK = re.compile(r"/[^/]*/")   # IDT 5'/3'/internal modification codes: /56-FAM/, /3IABkFQ/, /ZEN/, ...


def strip_mods(seq):
    """Recover the bare nucleotide sequence from an oligo written in modification / LNA notation.
    e.g. an IDT order string '/56-FAM/CTTA+CA+A+GATAT+CC+ACCACA/3IABkFQ/' -> 'CTTACAAGATATCCACCACA'.
    Removes /.../ modification blocks (FAM, ZEN, IABkFQ, internal mods), the LNA '+' prefix (one '+'
    marks the FOLLOWING base as a locked nucleic acid; the base itself stays), the phosphorothioate
    '*' linkage marker, and whitespace. IUPAC letters are left intact and case is preserved; callers
    do their own upper/U->T. A bare sequence passes through unchanged."""
    if not seq:
        return seq
    s = _MOD_BLOCK.sub("", seq)
    return s.replace("+", "").replace("*", "").replace(" ", "").replace("\t", "")


def _resolve(seq):
    # primer3 accepts only ACGT: map each IUPAC code / U to a concrete base and drop anything
    # else (whitespace, alignment gaps, stray punctuation) so Tm never crashes on a dirty oligo.
    return "".join(_RESOLVE.get(b, b) for b in seq.upper() if b in "ACGT" or b in _RESOLVE)


def has_degenerate(seq):
    return any(b not in "ACGTN" for b in seq.upper())


def gc_percent(seq):
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def amplicon_tm(seq):
    """Estimated melt Tm of a PCR product. Empirical salt-adjusted formula
    Tm = 81.5 + 16.6*log10[Na+]eq + 0.41*%GC - 600/length, NOT a nearest-neighbor
    calc -- a rough predictor for the SYBR melt-curve peak, to be confirmed against
    the observed melt. Monovalent salt is taken as the von Ahsen Mg2+-equivalent
    [Na+]eq = [Na+] + 120*sqrt([Mg2+]_free), [Mg2+]_free = [Mg2+] - [dNTP] (1:1
    chelation), so the divalent contribution that dominates SYBR product stability
    is no longer dropped. Tracks the reaction salt set via set_conditions()."""
    import math
    seq = seq.upper().replace("U", "T")
    n = len(seq)
    if n < 1:
        return 0.0
    na_mM = max(COND.get("mv_conc", 50.0), 0.0)
    mg_free_mM = max(COND.get("dv_conc", 0.0) - COND.get("dntp_conc", 0.0), 0.0)
    na_eq = max(na_mM + 120.0 * math.sqrt(mg_free_mM), 1.0) / 1000.0   # von Ahsen, to M
    return 81.5 + 16.6 * math.log10(na_eq) + 0.41 * gc_percent(seq) - 600.0 / n


# ---- Two Tm scales, on purpose ----
#  * tm()      -> primer3 SantaLucia + Owczarzy-2004. The DESIGN/SELECTION Tm that every candidate
#                 ranking, acceptance window, and golden test is calibrated to. Kept stable so the
#                 validated designs (locked panel, autodesign winners) never silently drift.
#  * tm_acc()  -> SantaLucia 1998 NN + Owczarzy 2008 (free Mg2+ via a Mg:dNTP dissociation constant,
#                 R=[Mg2+]/[Na+] regime switch) + Ct/4 -- i.e. exactly IDT OligoAnalyzer's method
#                 (Biopython Tm_NN, nn_table=DNA_NN3, saltcorr=7; dnac1==dnac2 gives Ct/4). This is
#                 the ACCURATE Tm shown to the user in QC / pair / viewer / report. It reads ~1.6 C
#                 above the primer3 path, which was the systematic gap users saw vs OligoAnalyzer
#                 (benchmarks put this method within ~0.5 C of OligoAnalyzer). DNA_NN3 = Allawi &
#                 SantaLucia 1997, the set underlying the 1998 unified parameters.
# Selection and display are separated deliberately: improving the displayed accuracy must never
# change which primers the tool designs.
#
# SALT-MODEL AUDIT (v1.27.0): the qPCR-Mg2+ salt-correction question was reviewed and the model
# left unchanged, because it is already the correct one where it matters. Every Tm the user READS
# and REPORTS (QC / pair / viewer / report / MIQE, all via tm_acc + nn.params) uses the
# divalent-aware Owczarzy-2008 correction with free Mg2+ = [Mg2+] - [dNTP] (von Ahsen 1:1
# chelation) -- the quantity that actually sets duplex stability in a PCR master mix. This was
# verified to within 0.03 C against an independent from-scratch Owczarzy-2008 reimplementation
# across all 14 locked-panel oligos, and it responds correctly to Mg2+ (e.g. Mg 3->6 mM raises a
# 22-mer ~1.3 C) and to dNTP chelation (0.8 mM dNTP lowers it ~0.5 C). primer3's selection Tm is
# near-insensitive to Mg2+ (Owczarzy-2004 monovalent-equivalent path; Mg 1.5->6 mM moves it <0.3 C
# and dNTP chelation is ignored) -- acceptable BECAUSE it is used ONLY to rank/gate candidates
# against a fixed window, never shown as an accurate number. Switching the displayed Tm's salt
# model would change nothing (it is already Owczarzy-2008); switching the SELECTION Tm's model
# would shift the hand-validated locked panel and autodesign goldens for no accuracy gain, so it
# was deliberately not done. Mg/dNTP sensitivity of the reported Tm is pinned in test_regression.
TM_METHOD = "santalucia"      # primer3 NN set for the selection Tm
SALT_METHOD = "owczarzy"      # primer3 monovalent salt model (Owczarzy 2004) for the selection Tm
_NN_TABLE = _mt.DNA_NN3       # Biopython NN set for the accurate (display) Tm


def _calc_tm(seq, salt_method=None):
    """Selection Tm via primer3 under an explicit salt model (uncached; the calibration harness
    sweeps salt models without cache clashes)."""
    r = _resolve(seq)
    if not r:
        return 0.0
    with _P3:
        return primer3.calc_tm(r, tm_method=TM_METHOD,
                               salt_corrections_method=salt_method or SALT_METHOD, **COND)


@lru_cache(maxsize=8192)
def tm(seq):
    return _calc_tm(seq, SALT_METHOD)


@lru_cache(maxsize=8192)
def tm_acc(seq):
    """Accurate, IDT-method Tm (deg C): SantaLucia 1998 NN + Owczarzy 2008 + Ct/4. Shown to the
    user wherever a Tm is displayed. Tracks IDT OligoAnalyzer to ~0.5 C for plain DNA."""
    r = _resolve(seq)
    if len(r) < 2:
        return 0.0
    try:
        return float(_mt.Tm_NN(r, nn_table=_NN_TABLE,
                               dnac1=COND["dna_conc"], dnac2=COND["dna_conc"], selfcomp=False,
                               Na=COND["mv_conc"], K=0.0, Tris=0.0,
                               Mg=COND["dv_conc"], dNTPs=COND["dntp_conc"], saltcorr=7))
    except Exception:
        return 0.0


_IUPAC_SETS = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC",
               "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
               "V": "ACG", "N": "ACGT"}


def set_conditions(mv_conc=None, dv_conc=None, dntp_conc=None, dna_conc=None, anneal_c=None):
    """Update, in place, the reaction salt used by every Tm and structure calc.

    Only the arguments you pass are changed; the rest keep their current value.
    Lets the Tm match your actual master mix (free Mg2+, monovalent, dNTP, oligo nM)
    instead of a generic default. anneal_c sets the qPCR annealing temperature used for
    the anneal-temperature structure profiles. Session-global for this local single-user server.
    """
    global ANNEAL_C
    # Validate physical ranges BEFORE mutating: this salt is session-global and feeds every Tm
    # and structure calc, so a bad value (negative, non-finite, zero oligo, no salt, or absurdly
    # high) must not be allowed to corrupt it and silently poison downstream Tm.
    LIMITS = {"mv_conc": (0.0, 2000.0), "dv_conc": (0.0, 200.0),
              "dntp_conc": (0.0, 100.0), "dna_conc": (1e-6, 1e6)}
    given = {"mv_conc": mv_conc, "dv_conc": dv_conc, "dntp_conc": dntp_conc, "dna_conc": dna_conc}
    for _name, _val in given.items():
        if _val is None:
            continue
        try:
            _v = float(_val)
        except (TypeError, ValueError):
            return dict(error="%s must be a number" % _name)
        if _v != _v or _v in (float("inf"), float("-inf")):
            return dict(error="%s must be a finite number" % _name)
        _lo, _hi = LIMITS[_name]
        if _v < _lo or _v > _hi:
            return dict(error="%s out of range (%g..%g)" % (_name, _lo, _hi))
    if anneal_c is not None:
        try:
            _a = float(anneal_c)
        except (TypeError, ValueError):
            return dict(error="anneal_c must be a number")
        if _a != _a or _a in (float("inf"), float("-inf")):
            return dict(error="anneal_c must be a finite number")
        if _a < 30.0 or _a > 85.0:
            return dict(error="anneal_c out of range (30..85)")
    _eff_mv = float(mv_conc) if mv_conc is not None else COND["mv_conc"]
    _eff_dv = float(dv_conc) if dv_conc is not None else COND["dv_conc"]
    if _eff_mv + _eff_dv <= 0:
        return dict(error="need some salt: monovalent + divalent must be > 0 mM")
    if mv_conc is not None:
        COND["mv_conc"] = float(mv_conc)
    if dv_conc is not None:
        COND["dv_conc"] = float(dv_conc)
    if dntp_conc is not None:
        COND["dntp_conc"] = float(dntp_conc)
    if dna_conc is not None:
        COND["dna_conc"] = float(dna_conc)
    if anneal_c is not None:
        ANNEAL_C = float(anneal_c)
    tm.cache_clear(); tm_acc.cache_clear(); hairpin.cache_clear(); self_dimer.cache_clear(); hetero_dimer.cache_clear()
    _hairpin_full_at.cache_clear(); _self_dimer_full_at.cache_clear()
    _hetero_dimer_full_at.cache_clear(); end_stability.cache_clear()
    return dict(COND, anneal_c=ANNEAL_C)


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
        t = round(tm_acc(s), 1)
        return dict(min=t, max=t, n=1, degenerate=False, capped=False)
    ncomb = 1
    for _, opts in degen:
        ncomb *= len(opts)
    if ncomb > cap:
        t = round(tm_acc(s), 1)
        return dict(min=t, max=t, n=ncomb, degenerate=True, capped=True)
    idxs = [i for i, _ in degen]
    base = list(s)
    tms = []
    for combo in itertools.product(*[opts for _, opts in degen]):
        for j, b in zip(idxs, combo):
            base[j] = b
        tms.append(tm_acc("".join(base)))
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


# ---- anneal-temperature structure profiles (additive QC; the design gate still uses the
# 37 C functions above so validated rankings are unchanged) ----
# Each returns (dG@37 C, dG@anneal C, structure melting Tm C). dG@37 is byte-identical to the
# 37 C function above; dG@anneal is the temperature-correct value; the melting Tm answers
# "does this structure even exist at the anneal temperature" directly (Tm < anneal => melted).
# The anneal temperature is an EXPLICIT argument (part of the cache key), not the module global,
# so an assay is gated at ITS OWN anneal temp: 60 C for the host TaqMan/SYBR panel, 54 C for the
# AT-rich parasite mtDNA assays. Judging a 54 C assay's structure at 60 C melts ~6 C more than
# reality and under-rejects (a self-dimer that survives 54 C is scored as gone) -- the mirror of the
# 37 C over-rejection this release fixed. Design gates pass the profile's anneal_c (design.py); QC/
# display callers pass none and fall back to the session global via the public wrappers below.
@lru_cache(maxsize=8192)
def _hairpin_full_at(seq, anneal_c):
    s = _resolve(seq)
    if len(s) > _P3_MAXLEN:
        return 0.0, 0.0, 0.0
    with _P3:
        r37 = primer3.calc_hairpin(s, temp_c=37.0, **COND)
        ran = primer3.calc_hairpin(s, temp_c=anneal_c, **COND)
    return r37.dg / 1000.0, ran.dg / 1000.0, r37.tm


def hairpin_full(seq, anneal_c=None):
    return _hairpin_full_at(seq, ANNEAL_C if anneal_c is None else float(anneal_c))


@lru_cache(maxsize=8192)
def _self_dimer_full_at(seq, anneal_c):
    s = _resolve(seq)
    if len(s) > _P3_MAXLEN:
        return 0.0, 0.0, 0.0
    with _P3:
        r37 = primer3.calc_homodimer(s, temp_c=37.0, **COND)
        ran = primer3.calc_homodimer(s, temp_c=anneal_c, **COND)
    return r37.dg / 1000.0, ran.dg / 1000.0, r37.tm


def self_dimer_full(seq, anneal_c=None):
    return _self_dimer_full_at(seq, ANNEAL_C if anneal_c is None else float(anneal_c))


@lru_cache(maxsize=8192)
def _hetero_dimer_full_at(a, b, anneal_c):
    sa, sb = _resolve(a), _resolve(b)
    if len(sa) > _P3_MAXLEN or len(sb) > _P3_MAXLEN:
        return 0.0, 0.0, 0.0
    with _P3:
        r37 = primer3.calc_heterodimer(sa, sb, temp_c=37.0, **COND)
        ran = primer3.calc_heterodimer(sa, sb, temp_c=anneal_c, **COND)
    return r37.dg / 1000.0, ran.dg / 1000.0, r37.tm


def hetero_dimer_full(a, b, anneal_c=None):
    return _hetero_dimer_full_at(a, b, ANNEAL_C if anneal_c is None else float(anneal_c))


@lru_cache(maxsize=8192)
def end_stability(a, b):
    """dG (kcal/mol) of the most stable 3'-end-anchored duplex of primer a against b
    (primer3 calc_end_stability). The 3'-anchored case is the dangerous one: a 3'-engaged
    primer-dimer can be EXTENDED into an artifact that competes with the target, whereas an
    internal dimer of equal global dG only titrates primer away. Used to flag which flagged
    cross-dimers are 3'-engaged. More negative = stronger 3' anchoring."""
    sa, sb = _resolve(a), _resolve(b)
    if not sa or not sb or len(sa) > _P3_MAXLEN or len(sb) > _P3_MAXLEN:
        return 0.0
    with _P3:
        return primer3.calc_end_stability(sa, sb, **COND).dg / 1000.0


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


# McTigue 2004 per-LNA Tm increments (deg C), averaged over neighbor context. The full McTigue
# model is 32 nearest-neighbor terms (ddH/ddS for 5'-MX(L) and 5'-X(L)N) giving ~2 C accuracy;
# these base-averaged increments are a point estimate that sits inside McTigue's spread. LNA
# pyrimidines and A stabilize most; purine neighbors add stability. IDT OligoAnalyzer (Affinity
# Plus) uses the full McTigue model and is authoritative -- always confirm an LNA Tm there.
_LNA_DTM = {"A": 4.0, "T": 3.0, "C": 4.0, "G": 3.0}


def tm_lna(seq, n_lna=None):
    """Effective Tm of an LNA oligo: the accurate all-DNA backbone Tm (SantaLucia + Owczarzy
    2008) plus a McTigue-informed per-LNA increment. Returns a point estimate plus a +/-2 C
    band (McTigue's stated accuracy). Consecutive LNAs are slightly less than additive, so the
    estimate is mildly conservative for runs of LNA."""
    dna, n, pos = strip_lna(seq)
    if n_lna is not None:
        n = n_lna
    base = tm_acc(dna)
    out = dict(dna_backbone_tm=round(base, 1), n_lna=n, lna_pos=pos)
    if n:
        bases = [dna[p - 1] for p in pos] if (pos and all(0 < p <= len(dna) for p in pos)) else []
        inc = sum(_LNA_DTM.get(b, 3.5) for b in bases) if bases else 3.5 * n
        est = base + inc
        out.update(est_tm=round(est, 1), est_tm_low=round(est - 2.0, 1), est_tm_high=round(est + 2.0, 1),
                   note="LNA Tm = accurate DNA-backbone Tm + McTigue-informed per-LNA increment "
                        "(point estimate, ~+/-2 C). Consecutive LNAs are less than additive. IDT "
                        "OligoAnalyzer (Affinity Plus) uses the full McTigue model -- confirm there.")
    return out
