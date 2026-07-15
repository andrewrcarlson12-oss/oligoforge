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
    # Structure REJECTION is judged at THIS ASSAY'S annealing temperature (c["anneal_c"]: 60 C host,
    # 54 C parasite), not primer3's 37 C default and not a single session global: a hairpin/self-dimer
    # fully melted at Ta does not exist during priming and must not disqualify an otherwise-good primer,
    # and a 54 C assay must be judged at 54 C, not 60. hairpin_full/self_dimer_full return
    # (dG@37, dG@anneal, structure_Tm); we gate on the anneal-temperature dG. Thresholds are unchanged.
    # (Ranking/tie-breaks elsewhere stay on the 37 C basis so validated picks are stable.)
    ac = c.get("anneal_c", T.ANNEAL_C)
    if T.hairpin_full(seq, ac)[1] <= c["hairpin_min"]: return False
    if T.self_dimer_full(seq, ac)[1] <= c["self_dimer_min"]: return False
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
    ac = c.get("anneal_c", T.ANNEAL_C)
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
    _sd = {}
    def _self(s):
        v = _sd.get(s)
        if v is None:
            v = T.self_dimer(s); _sd[s] = v
        return v
    for (score, fs, fe, f, rs, re, r, amp, gap) in prelim:
        # REJECT cross-dimers that persist at the annealing temperature; a cross-dimer melted at Ta
        # is not a real interference. The reported/ranked `worst` dimer stays on the 37 C basis so
        # the dimer-floor tie-break (and thus every validated pick) is unchanged.
        if T.hetero_dimer_full(f, r, ac)[1] <= c["pair_dimer_min"]: continue
        worst = min(_self(f), _self(r), T.hetero_dimer(f, r))
        pairs.append(dict(score=score, fstart=fs, fend=fe, f=f,
                          rstart=rs, rend=re, r=r, amp=amp, gap=gap, dimer=round(worst, 2)))
        if len(pairs) >= _PAIR_CAP: break
    # Final pick is dimer-aware, not Tm-fit alone: a penalty that bites only below a safe dimer
    # floor (-5.5 kcal/mol) lets a clean pair overtake a marginally-better-Tm pair that self- or
    # hetero-dimerizes. Pairs at/above the floor keep their Tm-fit order, so existing designs whose
    # worst dimer is milder than -5.5 are unchanged.
    pairs.sort(key=lambda p: p["score"] + 2.0 * max(0.0, -5.5 - p["dimer"]))
    return pairs


def find_probe(template, fend, rstart, fseq, rseq, c):
    template = template.upper()
    inner = template[fend:rstart]
    maxp = max(T.tm(fseq), T.tm(rseq))
    ac = c.get("anneal_c", T.ANNEAL_C)
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
                # REJECT probe structure at the annealing temperature (a probe hairpin/dimer melted
                # at Ta will not block hybridization); hdg (dG@37) is retained for reporting and for
                # the weakest-hairpin tie-break below so the validated probe pick is stable.
                hdg, hdg_an, htm = T.hairpin_full(sub, ac)
                if hdg_an <= c["probe_hairpin_min"]: continue
                if T.self_dimer_full(sub, ac)[1] <= c["self_dimer_min"]: continue
                if T.hetero_dimer_full(sub, fseq, ac)[1] <= c["pair_dimer_min"]: continue
                if T.hetero_dimer_full(sub, rseq, ac)[1] <= c["pair_dimer_min"]: continue
                # Prefer a probe ~9 C over the hotter primer -- the standard 8-10 C TaqMan placement --
                # whenever one is reachable, clamped to the profile's allowed window (AT-rich profiles
                # whose window caps below 9 fall back to their ceiling). Offset distance is bucketed to
                # whole degrees so near-equal offsets defer to the weakest-hairpin probe.
                target = min(max(c["probe_offset_min"], 9.0), c["probe_offset_max"])
                key = (round(abs((t - maxp) - target)), -hdg)
                if best is None or key < best[0]:
                    best = (key, dict(probe=sub, strand=strand, tm=t, offset=t - maxp,
                                      hairpin_dg=hdg, hairpin_tm=htm,
                                      self_dimer=T.self_dimer(sub),
                                      dimer_f=T.hetero_dimer(sub, fseq),
                                      dimer_r=T.hetero_dimer(sub, rseq)))
    return best[1] if best else None


def design_assay(template, c):
    """Full pipeline: best primer pair + its probe + gBlock, against a template.
    For no-probe (SYBR) profiles, returns the best primer pair with probe=None.

    CONTRACT: returns a dict on success, or **None** when no assay clears the profile
    (empty / too-short / all-N template, or no pair+probe passes the gates). Callers MUST
    guard for None before subscripting or .get() -- e.g. `a = design_assay(...); if not a:`.
    Every in-tree caller (batch_design, autodesign._design_one, design_candidates) already
    does, and /api/design length-checks the template first and returns a clean {error}.
    Do NOT "fix" this to return {"error": ...}: an error-dict is truthy and would silently
    break the existing `if not a` guards. Pinned by test_fuzz.py (no design endpoint 500s)."""
    fwd, rev = enumerate_primers(template, c)
    pairs = pair_primers(fwd, rev, c)
    if c.get("no_probe"):
        if not pairs:
            return None
        p = pairs[0]
        gb, gs, ge = build_gblock(template, p["fstart"], p["rend"])
        return dict(forward=p["f"], reverse=p["r"], probe=None, amplicon=p["amp"],
                    amplicon_tm=round(T.amplicon_tm(template.upper()[p["fstart"]:p["rend"]]), 1),
                    pair_tm_gap=p["gap"], pair_dimer=p.get("dimer"),
                    f_tm=T.tm(p["f"]), r_tm=T.tm(p["r"]),
                    probe_info=None, gblock=gb, gblock_span=(gs, ge),
                    f_xy=(p["fstart"], p["fend"]), r_xy=(p["rstart"], p["rend"]))
    for p in pairs[:10]:
        probe = find_probe(template, p["fend"], p["rstart"], p["f"], p["r"], c)
        if probe:
            gb, gs, ge = build_gblock(template, p["fstart"], p["rend"])
            return dict(forward=p["f"], reverse=p["r"], probe=probe["probe"],
                        amplicon=p["amp"], pair_tm_gap=p["gap"], pair_dimer=p.get("dimer"),
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


def build_offtarget_gblock(amplicon, off_seqs, flank=50, min_len=125):
    """Worst-case off-target discrimination control. Finds the off-target sequence whose
    region homologous to the assay's amplicon is the CLOSEST match -- the off-target the
    assay is most likely to fail to reject -- and returns it plus flanks as a synthesizable
    gene fragment. Ordering this next to the positive-control gBlock lets the bench test
    whether off-target rejection holds against the MOST similar off-target, not an average
    one (an average off-target makes a leaky block look clean).
    Returns {seq, off_index, amplicon_identity, span} or None."""
    if not amplicon or not off_seqs:
        return None
    L = len(amplicon)
    if L < 20:
        return None
    best = None  # (identity, off_index, oriented_seq, pos)
    for oi, raw in enumerate(off_seqs):
        for s in (raw, T.revcomp(raw)):
            if len(s) < L:
                continue
            for i in range(len(s) - L + 1):
                cutoff = L if best is None else int(L * (1 - best[0] / 100.0))
                mm = 0
                for a, b in zip(amplicon, s[i:i + L]):
                    if a != b:
                        mm += 1
                        if mm > cutoff:
                            break
                else:
                    ident = 100.0 * (L - mm) / L
                    if best is None or ident > best[0]:
                        best = (ident, oi, s, i)
    if best is None:
        return None
    ident, oi, s, pos = best
    a = max(0, pos - flank)
    e = min(len(s), pos + L + flank)
    while (e - a) < min_len and (a > 0 or e < len(s)):
        if a > 0:
            a -= 1
        if e < len(s):
            e += 1
    block = s[a:e]
    if len(block) < min_len:
        need = min_len - len(block)
        block += (_GBLOCK_FILLER * (need // len(_GBLOCK_FILLER) + 1))[:need]
    return dict(seq=block, off_index=oi, amplicon_identity=round(ident, 1), span=(a, e))


def probe_span(template, assay):
    """Exact base coordinates of the probe within the template -- the literal location the probe
    was taken from in design_assay (find_probe scans template[fend:rstart] on both strands).
    Returns [start, end] (so template[start:end] == probe for a +-strand probe, or == revcomp(probe)
    for a --strand probe), or None. No estimation: it is a substring search for the real oligo."""
    pi = assay.get("probe_info"); probe = assay.get("probe")
    if not pi or not probe:
        return None
    template = template.upper()
    needle = probe if pi.get("strand", "+") == "+" else T.revcomp(probe)
    fx, rx = assay["f_xy"][1], assay["r_xy"][0]
    i = template[fx:rx].find(needle)
    if i >= 0:
        return [fx + i, fx + i + len(needle)]
    j = template.find(needle)
    return [j, j + len(needle)] if j >= 0 else None


def _spread_starts(starts):
    """Return coordinates in an endpoint/midpoint order so a time budget samples the
    whole template before filling local gaps."""
    starts = sorted(set(int(x) for x in starts))
    out, seen = [], set()

    def visit(lo, hi):
        if lo > hi:
            return
        mid = (lo + hi) // 2
        for idx in (lo, hi, mid):
            if idx not in seen:
                seen.add(idx)
                out.append(starts[idx])
        visit(lo + 1, mid - 1)
        visit(mid + 1, hi - 1)

    if starts:
        visit(0, len(starts) - 1)
    return out


def _candidate_rank(a, c):
    """Lower is better. Rank complete assays globally rather than returning the
    first acceptable design window."""
    ft = float(a.get("f_tm", 0.0)); rt = float(a.get("r_tm", 0.0))
    score = abs((ft + rt) / 2.0 - float(c.get("tm_opt", (ft + rt) / 2.0)))
    score += 1.5 * abs(ft - rt)
    amp = int(a.get("amplicon") or 0)
    if amp > 150:
        score += 0.20 * (amp - 150)
    pair_dimer = a.get("pair_dimer")
    if isinstance(pair_dimer, (int, float)):
        score += 2.0 * max(0.0, -5.5 - float(pair_dimer))
    pi = a.get("probe_info")
    if pi:
        target = min(max(float(c.get("probe_offset_min", 0.0)), 9.0),
                     float(c.get("probe_offset_max", 10.0)))
        score += 0.35 * abs(float(pi.get("offset", target)) - target)
        score += 1.5 * max(0.0, -3.0 - float(pi.get("hairpin_dg", 0.0)))
    return score


def design_candidates(template, c, n=5, window=400, step=120, budget=12.0):
    """Return globally ranked, coordinate-exact assays across the full template.

    The old implementation stopped after the first ``n`` successful 5'-to-3'
    windows.  This created positional bias and could hide stronger downstream
    assays.  The revised search samples the entire target in a spread order,
    collects a bounded pool, remaps every coordinate to the full template, then
    ranks complete assays before trimming to ``n``.
    """
    import time
    template = template.upper()
    L = len(template)
    n = max(1, int(n))
    window = max(int(window), min(int(c.get("amp_max", 150)) + 120, 2200))
    if L <= window:
        starts = [0]
    else:
        span = L - window
        max_windows = 30
        even = [round(i * span / (max_windows - 1)) for i in range(max_windows)]
        grid = list(range(0, span + 1, max(1, int(step)))) + [span]
        merged = sorted(set(even + grid))
        if len(merged) > max_windows:
            merged = [merged[round(i * (len(merged) - 1) / (max_windows - 1))]
                      for i in range(max_windows)]
        starts = _spread_starts(merged)

    out, seen = [], set()
    t0 = time.monotonic()
    for idx, start in enumerate(starts):
        if idx and time.monotonic() - t0 >= float(budget):
            break
        sub = template[start:start + window]
        try:
            a = design_assay(sub, c)
        except Exception:
            a = None
        if not a:
            continue
        key = (a.get("forward"), a.get("reverse"), a.get("probe"))
        if key in seen:
            continue
        seen.add(key)
        pxy = probe_span(sub, a)
        a["search_window_start"] = start
        a["f_xy"] = [a["f_xy"][0] + start, a["f_xy"][1] + start]
        a["r_xy"] = [a["r_xy"][0] + start, a["r_xy"][1] + start]
        a["probe_xy"] = [pxy[0] + start, pxy[1] + start] if pxy else None
        a["amplicon_xy"] = [a["f_xy"][0], a["r_xy"][1]]
        if a.get("gblock_span"):
            a["gblock_span"] = [a["gblock_span"][0] + start,
                                  a["gblock_span"][1] + start]
        a["candidate_rank"] = round(_candidate_rank(a, c), 4)
        out.append(a)

    out.sort(key=lambda a: (a.get("candidate_rank", float("inf")),
                            a.get("amplicon", 10**9),
                            a.get("search_window_start", 0)))
    return out[:n]
