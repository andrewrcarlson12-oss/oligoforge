"""Unified nearest-neighbor thermodynamics for DNA/DNA duplexes — the physical-chemistry core.

Returns the quantities the rest of the app (and the upcoming mismatch/LNA work) actually need and
that a bare Tm number hides: duplex enthalpy ΔH°, entropy ΔS°, free energy ΔG°37, a salt-corrected
melting temperature at the reaction conditions, and the fraction of oligo hybridized at the
annealing temperature.

Parameters: SantaLucia (1998) unified nearest-neighbor set for Watson-Crick pairs (ΔH° kcal/mol,
ΔS° cal/(mol·K) at 1 M Na+). These are *identical* to the table Biopython's MeltingTemp.Tm_NN uses
(asserted in tests/test_nn.py), so this engine and the display Tm (thermo.tm_acc) share one source
of truth — nothing drifts.

Salt: monovalent via Owczarzy 2004; magnesium via Owczarzy 2008, with FREE magnesium computed as
[Mg2+]_total - [dNTP]_total (von Ahsen 2001 — each dNTP chelates ~1 Mg2+ 1:1), which is the
quantity that actually sets duplex stability in a PCR master mix. The monovalent/divalent regime is
selected by the ratio R = sqrt([Mg2+]_free)/[Mon+] exactly as in Owczarzy 2008.

This module is additive: it does not touch the selection Tm (thermo.tm, primer3) or the display Tm
(thermo.tm_acc), so every locked-panel / golden number is unchanged. Validated against Biopython in
tests/test_nn.py.
"""
import math
from . import thermo as T

R = 1.987  # gas constant, cal/(mol·K)

# SantaLucia (1998) unified NN parameters. Key = 5'-dinucleotide / its base-complement (NOT reversed),
# e.g. 5'-CA-3' over 3'-GT-5' is "CA/GT". (ΔH° kcal/mol, ΔS° cal/(mol·K)).
NN = {
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}
INIT = {"init": (0.0, 0.0), "A/T": (2.3, 4.1), "G/C": (0.1, -2.8), "sym": (0.0, -1.4)}
_C = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _rc(s):
    return "".join(_C[b] for b in reversed(s))


def duplex_nn(seq):
    """(ΔH° kcal/mol, ΔS° cal/(mol·K)) for the perfectly complementary duplex of `seq` at 1 M Na+.
    Returns None for anything that isn't a >=2 nt pure-ACGT sequence (mismatch / degenerate / modified
    oligos are out of scope for this first engine and handled elsewhere)."""
    s = (seq or "").upper()
    if len(s) < 2 or any(b not in "ACGT" for b in s):
        return None
    dh, ds = INIT["init"]
    for end in (s[0], s[-1]):                      # terminal initiation, per end
        h, sd = INIT["A/T"] if end in "AT" else INIT["G/C"]
        dh += h; ds += sd
    if s == _rc(s):                                # self-complementary symmetry term
        dh += INIT["sym"][0]; ds += INIT["sym"][1]
    for i in range(len(s) - 1):
        di = s[i:i + 2]
        pair = NN.get(di + "/" + _C[di[0]] + _C[di[1]])
        if pair is None:                           # try the reverse orientation of the key
            k = (di + "/" + _C[di[0]] + _C[di[1]])[::-1]
            pair = NN.get(k)
        dh += pair[0]; ds += pair[1]
    return dh, ds


def _tm_1M(dh, ds, dnac_M, selfcomp):
    """Melting temperature (°C) at 1 M Na+ from ΔH/ΔS and oligo strand concentration dnac_M (mol/L).
    Total-strand term follows the SantaLucia convention (= Biopython Tm_NN with dnac1=dnac2): the
    effective concentration is dnac for a self-complementary strand and dnac/2 for two complementary
    strands at equal concentration."""
    k = dnac_M if selfcomp else dnac_M / 2.0
    return (1000.0 * dh) / (ds + R * math.log(k)) - 273.15


def _salt_tm(tm1M_c, fgc, nbp, mono_M, mg_free_M):
    """Correct a 1 M-Na+ Tm (°C) to the actual monovalent + free-Mg2+ conditions.
    Owczarzy 2004 (monovalent) / Owczarzy 2008 (magnesium), regime by R = sqrt([Mg])/[Mon]."""
    inv = 1.0 / (tm1M_c + 273.15)

    def _mono(inv):
        ln = math.log(mono_M)
        return inv + (4.29 * fgc - 3.95) * 1e-5 * ln + 9.40e-6 * ln * ln

    if mg_free_M <= 1e-9:
        return (1.0 / _mono(inv) - 273.15) if mono_M > 0 else tm1M_c
    a, b, c, d, e, f, g = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
    if mono_M > 0:
        ratio = math.sqrt(mg_free_M) / mono_M
        if ratio < 0.22:                            # monovalent dominates
            return 1.0 / _mono(inv) - 273.15
        if ratio < 6.0:                             # both ions contribute: modify a, d, g
            lm = math.log(mono_M)
            a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mono_M) * lm)
            d = 1.42e-5 * (1.279 - 4.03e-3 * lm - 8.03e-3 * lm * lm)
            g = 8.31e-5 * (0.486 - 0.258 * lm + 5.25e-3 * lm ** 3)
    lnMg = math.log(mg_free_M)
    inv = (inv + a + b * lnMg + fgc * (c + d * lnMg)
           + (1.0 / (2.0 * (nbp - 1))) * (e + f * lnMg + g * lnMg * lnMg))
    return 1.0 / inv - 273.15


def _conds(cond):
    c = dict(T.COND)
    if cond:
        c.update(cond)
    return c


def gibbs(seq, temp_c, cond=None):
    """ΔG (kcal/mol) of duplex formation at temperature `temp_c`, ΔG(T) = ΔH - T·ΔS. 1 M Na+
    reference (salt shifts Tm, captured separately in `params`); negative = favorable."""
    dn = duplex_nn(seq)
    if dn is None:
        return None
    dh, ds = dn
    return round(dh - (temp_c + 273.15) * ds / 1000.0, 2)


def params(seq, cond=None, anneal_c=None):
    """Full NN readout at the reaction conditions. Returns dict(dh, ds, dg37, tm, frac_bound,
    anneal_c, mono_mM, mg_free_mM, gc) or None if `seq` is not a clean ACGT oligo (>=2 nt)."""
    dn = duplex_nn(seq)
    if dn is None:
        return None
    dh, ds = dn
    c = _conds(cond)
    ac = T.ANNEAL_C if anneal_c is None else float(anneal_c)
    nbp = len(seq)
    fgc = T.gc_percent(seq) / 100.0
    selfc = seq.upper() == _rc(seq.upper())
    tm1 = _tm_1M(dh, ds, max(c["dna_conc"], 1e-3) * 1e-9, selfc)
    mono = max(c["mv_conc"], 0.0) / 1000.0
    mg_free = max(c["dv_conc"] - c["dntp_conc"], 0.0) / 1000.0
    tm = _salt_tm(tm1, fgc, nbp, mono, mg_free)
    dg37 = dh - 310.15 * ds / 1000.0
    # fraction of oligo hybridized at the annealing temperature (primer-excess two-state, anchored at
    # Tm where it is 0.5): θ(T) = 1 / (1 + exp[(ΔH/R)(1/T - 1/Tm)]).
    tmK, taK = tm + 273.15, ac + 273.15
    frac = 1.0 / (1.0 + math.exp((dh * 1000.0 / R) * (1.0 / taK - 1.0 / tmK))) if tmK > 0 else None
    return dict(dh=round(dh, 1), ds=round(ds, 1), dg37=round(dg37, 2), tm=round(tm, 1),
                frac_bound=(round(frac, 3) if frac is not None else None), anneal_c=ac,
                mono_mM=round(mono * 1000, 1), mg_free_mM=round(mg_free * 1000, 2),
                gc=round(fgc * 100, 1))


# ---------------------------------------------------------------------------
# LNA (locked nucleic acid) layer — McTigue, Peterson & Kahn (2004) Biochemistry
# 43:5388-5405. Per-substitution nearest-neighbor INCREMENTS (ΔΔH cal/mol, ΔΔS
# cal/mol·K) added on top of the SantaLucia DNA sum: ΔH = ΔH_DNA + ΣΔΔH_LNA,
# ΔS = ΔS_DNA + ΣΔΔS_LNA (the additive model; see IDT US2012/0029891). Values
# transcribed from the McTigue 2004 set as distributed with EBI MELTING; this
# engine's DNA ΔH/ΔS were verified byte-identical to MELTING's, and the combined
# model reproduces McTigue's 100 published duplex Tm's at ~1.7 °C RMSE
# (tests/test_nn.py). Key: "XLY/.." = LNA on the 5' base of the step,
# "XYL/.." = LNA on the 3' base; bottom strand is the base-complement.
_LNA_INC = {
    "ALA/TT": (707.0, 2.5),  "ALC/TG": (1131.0, 4.1), "ALG/TC": (264.0, 2.6),  "ALT/TA": (2282.0, 7.5),
    "AAL/TT": (992.0, 4.1),  "ACL/TG": (2890.0, 10.6),"AGL/TC": (-1200.0, -1.8),"ATL/TA": (1816.0, 6.9),
    "CLA/GT": (1049.0, 4.3), "CLC/GG": (2096.0, 8.0), "CLG/GC": (785.0, 3.7),  "CLT/GA": (708.0, 4.2),
    "CAL/GT": (1358.0, 4.4), "CCL/GG": (2063.0, 7.6), "CGL/GC": (-276.0, -0.7),"CTL/GA": (-1671.0, -4.1),
    "GLA/CT": (3162.0, 10.5),"GLC/CG": (-360.0, -0.3),"GLG/CC": (-2844.0, -6.7),"GLT/CA": (-212.0, 0.1),
    "GAL/CT": (444.0, 2.9),  "GCL/CG": (-925.0, -1.1),"GGL/CC": (-943.0, -0.9),"GTL/CA": (-635.0, -0.3),
    "TLA/AT": (-46.0, 1.6),  "TLC/AG": (1893.0, 6.7), "TLG/AC": (-1540.0, -3.0),"TLT/AA": (1528.0, 5.3),
    "TAL/AT": (1591.0, 5.3), "TCL/AG": (609.0, 3.2),  "TGL/AC": (2165.0, 7.2), "TTL/AA": (2326.0, 8.1),
}


def parse_lna(seq):
    """Parse an oligo written with IDT '+' LNA notation (and other mod blocks) into a base list and a
    parallel list of LNA flags. '/5.../' mod blocks, phosphorothioate '*', and whitespace are removed;
    '+X' marks base X as LNA. Returns (bases, lna_flags)."""
    import re
    s = re.sub(r"/[^/]*/", "", seq or "").replace("*", "").replace(" ", "").upper()
    bases, lna, pend = [], [], False
    for ch in s:
        if ch == "+":
            pend = True
        elif ch in "ACGT":
            bases.append(ch); lna.append(pend); pend = False
    return bases, lna


def _duplex_nn_lna(bases, lna):
    """(ΔH kcal/mol, ΔS cal/mol·K, warnings) for an LNA-substituted oligo = SantaLucia DNA core +
    McTigue increments. `warnings` flags conditions outside McTigue's single-internal-substitution
    parameterization (adjacent LNAs, terminal LNA, high density)."""
    core = duplex_nn("".join(bases))
    if core is None:
        return None
    dh, ds = core
    ih = isum = 0.0
    n = len(bases)
    for i in range(n - 1):
        b0, b1 = bases[i], bases[i + 1]
        comp = _C[b0] + _C[b1]
        keys = []
        if lna[i]:
            keys.append(b0 + "L" + b1 + "/" + comp)      # LNA is 5' base of this step
        if lna[i + 1]:
            keys.append(b0 + b1 + "L" + "/" + comp)       # LNA is 3' base of this step
        for k in keys:
            p = _LNA_INC.get(k)
            if p:
                ih += p[0]; isum += p[1]
    warn = []
    if any(lna[i] and lna[i + 1] for i in range(n - 1)):
        warn.append("adjacent LNAs")
    if lna and (lna[0] or lna[-1]):
        warn.append("terminal LNA")
    if sum(lna) > n / 2.0:
        warn.append("high LNA density")
    return dh + ih / 1000.0, ds + isum, warn


def params_lna(seq, cond=None, anneal_c=None):
    """NN thermodynamics for an LNA-substituted oligo at the reaction conditions, computed (not
    estimated) via McTigue 2004 increments. Returns dict(dh, ds, dg37, tm, frac_bound, lna_n, gc,
    anneal_c, mono_mM, mg_free_mM, beyond_param, note) or None if there is no '+' LNA, fewer than two
    bases, or any non-ACGT base (degenerate / unmodelled). Uses the same salt and concentration
    convention as params(), so the LNA Tm is directly comparable to the DNA-backbone NN Tm."""
    bases, lna = parse_lna(seq)
    import re as _re
    _clean = _re.sub(r"/[^/]*/", "", seq or "").replace("*", "").replace(" ", "").replace("+", "").upper()
    if any(ch not in "ACGT" for ch in _clean):   # degenerate/unknown base present -> not modelled here
        return None
    if len(bases) < 2 or not any(lna) or any(b not in "ACGT" for b in bases):
        return None
    res = _duplex_nn_lna(bases, lna)
    if res is None:
        return None
    dh, ds, warn = res
    s = "".join(bases)
    c = _conds(cond)
    ac = T.ANNEAL_C if anneal_c is None else float(anneal_c)
    fgc = T.gc_percent(s) / 100.0
    selfc = s == _rc(s)
    tm1 = _tm_1M(dh, ds, max(c["dna_conc"], 1e-3) * 1e-9, selfc)
    mono = max(c["mv_conc"], 0.0) / 1000.0
    mg_free = max(c["dv_conc"] - c["dntp_conc"], 0.0) / 1000.0
    tm = _salt_tm(tm1, fgc, len(s), mono, mg_free)
    dg37 = dh - 310.15 * ds / 1000.0
    tmK, taK = tm + 273.15, ac + 273.15
    frac = 1.0 / (1.0 + math.exp((dh * 1000.0 / R) * (1.0 / taK - 1.0 / tmK))) if tmK > 0 else None
    note = ("McTigue 2004 parameters cover single internal LNA substitutions; this oligo has "
            + " + ".join(warn) + ", so the value is an approximation — confirm in IDT OligoAnalyzer.") if warn else ""
    return dict(dh=round(dh, 1), ds=round(ds, 1), dg37=round(dg37, 2), tm=round(tm, 1),
                frac_bound=(round(frac, 3) if frac is not None else None), lna_n=sum(lna),
                anneal_c=ac, mono_mM=round(mono * 1000, 1), mg_free_mM=round(mg_free * 1000, 2),
                gc=round(fgc * 100, 1), beyond_param=bool(warn), note=note)


# ---------------------------------------------------------------------------
# Mismatch layer — Allawi & SantaLucia (1997), Allawi & SantaLucia (1998),
# Peyret et al. (1999): nearest-neighbor parameters for internal and terminal
# single mismatches in DNA. Used to predict the stability of an oligo bound to a
# non-perfect target (an off-by-N isolate), turning "N mismatches" into a
# predicted duplex Tm, ΔTm vs the perfect match, and fraction bound at the
# annealing temperature. Parameters are taken from Biopython's DNA_IMM1
# (internal) and DNA_TMM1 (terminal) tables (same Allawi/SantaLucia/Peyret
# provenance as our Watson-Crick set); the walk replicates Biopython's Tm_NN
# term logic exactly and is validated against it in tests/test_nn.py.
_COMPL = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _complement(s):
    return "".join(_COMPL.get(b, "N") for b in s.upper())


def _duplex_general(seq, c_seq):
    """(ΔH kcal/mol, ΔS cal/mol·K) for seq (5'->3') paired with c_seq (3'->5', base-aligned, same
    length), allowing internal and terminal mismatches. Replicates Biopython Tm_NN's term logic with
    SantaLucia WC + Allawi/SantaLucia/Peyret mismatch parameters. None if lengths differ / <2 / unknown
    step. (For the DNA_NN3 set, init reduces to terminal A/T & G/C terms — the others are zero.)"""
    from Bio.SeqUtils import MeltingTemp as _mt
    IMM, TMM = _mt.DNA_IMM1, _mt.DNA_TMM1
    s, cs = seq.upper(), c_seq.upper()
    if len(s) != len(cs) or len(s) < 2:
        return None
    dh = ds = 0.0
    ts, tcs = s, cs
    left = tcs[:2][::-1] + "/" + ts[:2][::-1]          # terminal mismatch, 5' end
    if left in TMM:
        dh += TMM[left][0]; ds += TMM[left][1]; ts = ts[1:]; tcs = tcs[1:]
    right = ts[-2:] + "/" + tcs[-2:]                    # terminal mismatch, 3' end
    if right in TMM:
        dh += TMM[right][0]; ds += TMM[right][1]; ts = ts[:-1]; tcs = tcs[:-1]
    ends = s[0] + s[-1]                                  # terminal initiation (Biopython uses original seq)
    at = ends.count("A") + ends.count("T"); gc = ends.count("G") + ends.count("C")
    dh += INIT["A/T"][0] * at + INIT["G/C"][0] * gc
    ds += INIT["A/T"][1] * at + INIT["G/C"][1] * gc
    for i in range(len(ts) - 1):                        # zipping: mismatch table first, then WC
        nb = ts[i:i + 2] + "/" + tcs[i:i + 2]
        if nb in IMM:
            dh += IMM[nb][0]; ds += IMM[nb][1]
        elif nb[::-1] in IMM:
            dh += IMM[nb[::-1]][0]; ds += IMM[nb[::-1]][1]
        elif nb in NN:
            dh += NN[nb][0]; ds += NN[nb][1]
        elif nb[::-1] in NN:
            dh += NN[nb[::-1]][0]; ds += NN[nb[::-1]][1]
        else:
            return None
    return dh, ds


# IUPAC degeneracy -> the set of ACGT bases each code can represent.
_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
          "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
          "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _parse_oligo_sets(seq):
    """Parse an oligo (IDT '+' LNA notation, '/.../' mod blocks, '*', IUPAC degeneracy) into per-position
    ACGT option tuples and parallel LNA flags. An unrecognised base letter yields an empty tuple (the
    caller treats the oligo as unscoreable)."""
    import re
    s = re.sub(r"/[^/]*/", "", seq or "").replace("*", "").replace(" ", "").upper()
    sets, lna, pend = [], [], False
    for ch in s:
        if ch == "+":
            pend = True
            continue
        opts = _IUPAC.get(ch)
        if opts is None:
            if ch.isalpha():
                sets.append(()); lna.append(pend); pend = False
            continue
        sets.append(tuple(opts)); lna.append(pend); pend = False
    return sets, lna


def _variant_dh_ds(variant, lna, t):
    """ΔH/ΔS for a concrete ACGT `variant` (5'->3') vs target `t` (same length): mismatch-aware NN walk
    plus McTigue LNA increments added only on FULLY MATCHED steps (within McTigue's matched-context
    parameterization). Returns (dh, ds, lna_at_mismatch) or None."""
    core = _duplex_general(variant, _complement(t))
    if core is None:
        return None
    dh, ds = core
    ih = isum = 0.0
    n = len(variant)
    for i in range(n - 1):
        if not (variant[i] == t[i] and variant[i + 1] == t[i + 1]):
            continue                                   # LNA increment only where the whole step is WC
        b0, b1 = variant[i], variant[i + 1]
        comp = _COMPL[b0] + _COMPL[b1]
        if lna[i]:
            p = _LNA_INC.get(b0 + "L" + b1 + "/" + comp)
            if p:
                ih += p[0]; isum += p[1]
        if lna[i + 1]:
            p = _LNA_INC.get(b0 + b1 + "L" + "/" + comp)
            if p:
                ih += p[0]; isum += p[1]
    lna_mm = any(lna[i] and variant[i] != t[i] for i in range(n))
    return dh + ih / 1000.0, ds + isum, lna_mm


def _best_variant(sets, lna, t):
    """The concrete ACGT variant of a (possibly degenerate) oligo that best matches target `t`, plus the
    total number of synthesised variants. Covered positions resolve to the target base (the matching
    variant in the synthesised mix binds); uncovered degenerate positions are enumerated (bounded) to
    maximise duplex stability. Returns (variant_str, n_var) or (None, n_var)."""
    n = len(t)
    n_var = 1
    for S in sets:
        n_var *= max(len(S), 1)
    if any(len(S) == 0 for S in sets):
        return None, n_var
    fixed = [None] * n
    free = []
    for i in range(n):
        if t[i] in sets[i]:
            fixed[i] = t[i]
        elif len(sets[i]) == 1:
            fixed[i] = sets[i][0]
        else:
            free.append(i)
    if not free:
        return "".join(fixed), n_var
    import itertools
    opts = [sets[i] for i in free]
    total = 1
    for o in opts:
        total *= len(o)
    combos = itertools.product(*opts) if total <= 4096 else [tuple(o[0] for o in opts)]
    best, best_tm = None, -1e18
    for combo in combos:
        v = fixed[:]
        for idx, b in zip(free, combo):
            v[idx] = b
        vs = "".join(v)
        r = _variant_dh_ds(vs, lna, t)
        if r is None:
            continue
        tm = _tm_1M(r[0], r[1], 1e-7, vs == _rc(vs))   # 1 M Tm, fixed conc, for ranking only
        if tm > best_tm:
            best_tm, best = tm, vs
    return best, n_var


def mismatch_params(oligo, target, cond=None, anneal_c=None):
    """Predict the stability of `oligo` bound to `target` — the genomic region it anneals to, in the
    oligo's own 5'->3' frame (target == the oligo's matching variant for a perfect match), same length.
    Handles plain DNA, IDT '+' LNA, and IUPAC-degenerate oligos: a degenerate oligo is resolved to its
    best-matching synthesised variant against the target, and McTigue LNA increments are added on matched
    steps. Returns dict(tm, tm_perfect, dtm, frac_bound, frac_bound_perfect, n_mm, mm_pos, mm_from_3p, dh,
    ds, dg37, anneal_c, lna_n, n_var, degenerate, lna_mismatch, note) or None if the target is not ACGT,
    lengths differ, the oligo is < 2 nt, or it contains an unrecognised base. mm_from_3p is the nearest
    mismatch's distance from the 3' end (1 = terminal) — the number that decides primer extension."""
    t = (target or "").upper()
    sets, lna = _parse_oligo_sets(oligo)
    n = len(sets)
    if n < 2 or len(t) != n or any(b not in "ACGT" for b in t):
        return None
    variant, n_var = _best_variant(sets, lna, t)
    if variant is None:
        return None
    vd = _variant_dh_ds(variant, lna, t)
    if vd is None:
        return None
    dh, ds, lna_mm = vd
    c = _conds(cond)
    ac = T.ANNEAL_C if anneal_c is None else float(anneal_c)
    dnac = max(c["dna_conc"], 1e-3) * 1e-9
    mono = max(c["mv_conc"], 0.0) / 1000.0
    mg_free = max(c["dv_conc"] - c["dntp_conc"], 0.0) / 1000.0
    fgc = T.gc_percent(variant) / 100.0
    selfc = variant == _rc(variant)
    tm = _salt_tm(_tm_1M(dh, ds, dnac, selfc), fgc, n, mono, mg_free)
    pm = _variant_dh_ds(variant, lna, variant)         # perfect: best variant vs its own complement
    tm_pm = _salt_tm(_tm_1M(pm[0], pm[1], dnac, selfc), fgc, n, mono, mg_free) if pm else None
    taK = ac + 273.15

    def _frac(dHkcal, tmc):
        if tmc is None:
            return None
        tmK = tmc + 273.15
        return (1.0 / (1.0 + math.exp((dHkcal * 1000.0 / R) * (1.0 / taK - 1.0 / tmK)))) if tmK > 0 else None
    mm_pos = [i for i in range(n) if variant[i] != t[i]]
    from3p = (n - max(mm_pos)) if mm_pos else None
    degen = any(len(S) > 1 for S in sets)
    bits = []
    if degen:
        bits.append("degenerate oligo resolved to its best-matching variant (mix of %d)" % n_var)
    if lna_mm:
        bits.append("an LNA base sits on a mismatch (beyond McTigue parameterization)")
    return dict(tm=round(tm, 1), tm_perfect=(round(tm_pm, 1) if tm_pm is not None else None),
                dtm=(round(tm - tm_pm, 1) if tm_pm is not None else None),
                frac_bound=(round(_frac(dh, tm), 3) if _frac(dh, tm) is not None else None),
                frac_bound_perfect=(round(_frac(pm[0], tm_pm), 3) if pm and tm_pm is not None else None),
                n_mm=len(mm_pos), mm_pos=mm_pos, mm_from_3p=from3p,
                dh=round(dh, 1), ds=round(ds, 1), dg37=round(dh - 310.15 * ds / 1000.0, 2),
                anneal_c=ac, lna_n=sum(lna), n_var=n_var, degenerate=degen, lna_mismatch=lna_mm,
                note=("; ".join(bits) if bits else ""))
