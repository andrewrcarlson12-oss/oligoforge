"""Candidate enumeration + scoring, mirroring the manual QC workflow:
enumerate windows, gate each on Tm/GC/3'-end/structure, pair on amplicon size
and Tm match, then search the interior for a clean hydrolysis probe.
"""
from . import thermo as T


def _ok_primer(seq, c):
    if not (c["len_min"] <= len(seq) <= c["len_max"]): return False
    if not (c["gc_min"] <= T.gc_percent(seq) <= c["gc_max"]): return False
    if not (c["tm_min"] <= T.tm(seq) <= c["tm_max"]): return False
    if c.get("no_three_prime_T") and seq[-1] == "T": return False
    if T.last5_gc(seq) > c["max_3prime_gc"]: return False
    if T.max_run(seq, "G") >= c["max_g_run"]: return False
    if T.max_run(seq) >= c["max_any_run"]: return False
    if T.hairpin(seq)[0] <= c["hairpin_min"]: return False
    if T.self_dimer(seq) <= c["self_dimer_min"]: return False
    return True


def enumerate_primers(template, c):
    template = template.upper()
    fwd, rev = [], []
    n = len(template)
    for i in range(n):
        for k in range(c["len_min"], c["len_max"] + 1):
            if i + k > n: break
            w = template[i:i + k]
            if _ok_primer(w, c): fwd.append((i, i + k, w))
            rv = T.revcomp(w)
            if _ok_primer(rv, c): rev.append((i, i + k, rv))
    return fwd, rev


_PAIR_CAP = 300   # rank pairs by Tm-fit first, then dimer-test down the list only this far.
                  # Callers consume the best pair(s); on AT-rich templates hundreds of primers
                  # pass _ok_primer, so the naive O(F x R) hetero-dimer pass cost ~20 s/window.

def pair_primers(fwd, rev, c):
    """Best primer pairs, ranked by Tm-fit. Reverse-primer Tms are computed ONCE (were
    recomputed for every f x r pair, twice each); size/Tm-gap-valid pairs are scored cheaply,
    then the expensive hetero-dimer check is applied walking the ranked list until _PAIR_CAP
    dimer-clean pairs are collected. Identical output to the prior filter-then-sort for the
    top pairs every caller actually uses; just orders of magnitude fewer dimer calls."""
    amin, amax = c["amp_min"], c["amp_max"]
    gmin, gapmax, topt = c["min_probe_gap"], c["pair_tm_gap_max"], c["tm_opt"]
    rtm = [(rs, re, r, T.tm(r)) for (rs, re, r) in rev]
    prelim = []
    for (fs, fe, f) in fwd:
        ftm = T.tm(f)
        for (rs, re, r, rt) in rtm:
            amp = re - fs
            if not (amin <= amp <= amax): continue
            if rs - fe < gmin: continue
            gap = abs(ftm - rt)
            if gap > gapmax: continue
            score = abs((ftm + rt) / 2 - topt) + gap
            prelim.append((score, fs, fe, f, rs, re, r, amp, gap))
    prelim.sort(key=lambda x: x[0])
    pairs = []
    for (score, fs, fe, f, rs, re, r, amp, gap) in prelim:
        if T.hetero_dimer(f, r) <= c["pair_dimer_min"]: continue
        pairs.append(dict(score=score, fstart=fs, fend=fe, f=f,
                          rstart=rs, rend=re, r=r, amp=amp, gap=gap))
        if len(pairs) >= _PAIR_CAP: break
    return pairs


def find_probe(template, fend, rstart, fseq, rseq, c):
    template = template.upper()
    inner = template[fend:rstart]
    maxp = max(T.tm(fseq), T.tm(rseq))
    best = None
    for a in range(len(inner)):
        for k in range(c["probe_len_min"], c["probe_len_max"] + 1):
            b = a + k
            if b > len(inner): break
            for strand, sub in (("+", inner[a:b]), ("-", T.revcomp(inner[a:b]))):
                if "G" in sub[:3]: continue                       # 5'-FAM quench guard
                if sub.count("C") < sub.count("G"): continue      # more C than G
                if T.max_run(sub, "G") >= c["max_g_run"]: continue
                if T.max_run(sub) >= c["max_any_run"]: continue
                t = T.tm(sub)
                if not (maxp + c["probe_offset_min"] <= t <= maxp + c["probe_offset_max"]): continue
                hdg, htm = T.hairpin(sub)
                if hdg <= c["probe_hairpin_min"]: continue
                if T.self_dimer(sub) <= c["self_dimer_min"]: continue
                if T.hetero_dimer(sub, fseq) <= c["pair_dimer_min"]: continue
                if T.hetero_dimer(sub, rseq) <= c["pair_dimer_min"]: continue
                mid = (c["probe_offset_min"] + c["probe_offset_max"]) / 2
                key = (-hdg, abs((t - maxp) - mid))
                if best is None or key < best[0]:
                    best = (key, dict(probe=sub, strand=strand, tm=t, offset=t - maxp,
                                      hairpin_dg=hdg, hairpin_tm=htm,
                                      self_dimer=T.self_dimer(sub),
                                      dimer_f=T.hetero_dimer(sub, fseq),
                                      dimer_r=T.hetero_dimer(sub, rseq)))
    return best[1] if best else None


def design_assay(template, c):
    """Full pipeline: best primer pair + its probe + gBlock, against a template.
    For no-probe (SYBR) profiles, returns the best primer pair with probe=None."""
    fwd, rev = enumerate_primers(template, c)
    pairs = pair_primers(fwd, rev, c)
    if c.get("no_probe"):
        if not pairs:
            return None
        p = pairs[0]
        gb, gs, ge = build_gblock(template, p["fstart"], p["rend"])
        return dict(forward=p["f"], reverse=p["r"], probe=None, amplicon=p["amp"],
                    amplicon_tm=round(T.amplicon_tm(template.upper()[p["fstart"]:p["rend"]]), 1),
                    pair_tm_gap=p["gap"], f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]),
                    probe_info=None, gblock=gb, gblock_span=(gs, ge),
                    f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
    for p in pairs[:10]:
        probe = find_probe(template, p["fend"], p["rstart"], p["f"], p["r"], c)
        if probe:
            gb, gs, ge = build_gblock(template, p["fstart"], p["rend"])
            return dict(forward=p["f"], reverse=p["r"], probe=probe["probe"],
                        amplicon=p["amp"], pair_tm_gap=p["gap"],
                        amplicon_tm=round(T.amplicon_tm(template.upper()[p["fstart"]:p["rend"]]), 1),
                        f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]), probe_info=probe,
                        gblock=gb, gblock_span=(gs, ge),
                        f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
    return None


_GBLOCK_FILLER = ("GCATAGCTGACTGATCAGCTAGTCGATCAGTACGATCGTAGCATCGTAGCTAGCT"
                  "GACGTACGATCGATCAGCTAGCATCGATGCATCGATCGATCGTACGTAGCATGCA")  # GC~50%, no long runs


def build_gblock(template, fstart, rend, flank=40, min_len=125):
    """Amplicon + flanking sequence as a synthesizable standard, at least IDT's
    gBlock minimum (125 bp). Flanks first widen within the template; only if the
    template itself can't reach the minimum is a balanced, low-binding filler
    appended (IDT rejects gene-fragment orders under 125 bp)."""
    s = max(0, fstart - flank)
    e = min(len(template), rend + flank)
    while (e - s) < min_len and (s > 0 or e < len(template)):
        if s > 0:
            s -= 1
        if e < len(template):
            e += 1
    block = template[s:e]
    if len(block) < min_len:
        need = min_len - len(block)
        block += (_GBLOCK_FILLER * (need // len(_GBLOCK_FILLER) + 1))[:need]
    return block, s, e
