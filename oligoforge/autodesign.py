"""One-call assay design from an NCBI target.

Flow: fetch sequences for a target query -> design a primer/probe set on the most
complete sequence -> score every candidate for conservation across the WHOLE fetched
set and discrimination against an off-target set -> rank. This turns the manual
fetch/design/validate loop into one call.

It is a SCREEN, not a final answer. It designs on a representative sequence and
validates the result across the set; it does not replace looking at an alignment,
confirming finals in OligoAnalyzer, or bench validation. A bare organism name also
under-specifies the problem: it needs the gene (so the fetched sequences share a
locus) and, for a detection assay, what to discriminate against.
"""
from . import design as D, conservation as C, ncbi as N, thermo as T, profiles as PROF, specificity as SP, structure as STR, refmarkers as RM
from . import candidate_search as CSEARCH
from . import candidate_retention as CRET
from . import ranking as RANK
from . import ranking_explain as REXPLAIN
import hashlib
import json
import threading
from collections import OrderedDict
from copy import deepcopy

# Mitochondrial / ribosomal barcode words: their GenBank records (e.g. MalAvi cytb
# barcodes) are deposited as partial sequences not linked to a Gene record, so a precise
# [Gene Name] fetch would MISS them — keep these as free-text. Protein-coding genes
# (IFNG, IL4, RPL13...) instead get the precise fetch so "interferon gamma" lands on IFNG,
# not IFNGR1 (interferon-gamma receptor 1).
# Candidate generation makes many native Primer3 calls.  Bounded defensive-copy
# caches prevent repeated identical requests from re-entering that native search,
# improve interactive reproducibility, and avoid a primer3-py long-run stall seen
# after back-to-back exhaustive searches.  Keys include reaction conditions and all
# ranking/search inputs; callers never receive the cached mutable object itself.
_CACHE_LOCK = threading.RLock()
_SEARCH_CACHE = OrderedDict()
_DESIGN_CACHE = OrderedDict()
_SEARCH_CACHE_MAX = 8
_DESIGN_CACHE_MAX = 8


def _stable_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(raw.encode()).hexdigest()


def _cache_get(cache, key):
    with _CACHE_LOCK:
        if key not in cache:
            return None
        value = cache.pop(key)
        cache[key] = value
        return deepcopy(value)


def _cache_put(cache, key, value, limit):
    with _CACHE_LOCK:
        cache.pop(key, None)
        cache[key] = deepcopy(value)
        while len(cache) > limit:
            cache.popitem(last=False)


def clear_design_caches():
    """Clear bounded sequence-design memoization (tests/admin diagnostics)."""
    with _CACHE_LOCK:
        _SEARCH_CACHE.clear()
        _DESIGN_CACHE.clear()


def design_cache_info():
    with _CACHE_LOCK:
        return {"search_entries": len(_SEARCH_CACHE), "design_entries": len(_DESIGN_CACHE),
                "search_limit": _SEARCH_CACHE_MAX, "design_limit": _DESIGN_CACHE_MAX}


_MARKER_WORDS = {"cytb", "cob", "cytochrome", "oxidase", "nadh", "cox1", "cox2", "cox3", "coi", "co1", "coii", "coiii",
                 "nad1", "nad5", "nad4", "atp6", "18s", "28s", "16s", "23s", "12s",
                 "its", "its1", "its2", "rrna", "ribosomal", "rbcl", "matk", "trnl",
                 "trnh", "psba", "rdrp", "ssu", "lsu", "d-loop",
                 "spacer", "intergenic", "transcribed", "minicircle", "kinetoplast"}


def _is_marker(gene):
    g = (gene or "").lower()
    toks = set(g.replace("-", " ").split())
    return bool(toks & _MARKER_WORDS) or "control region" in g or "ribosomal" in g


def _split_query(q):
    """Best-effort split of a free-text target into (organism, gene): a leading capitalised
    genus (+ lowercase species epithet) is the organism, the remainder is the gene. Used to
    resolve the official gene symbol and to fill the workbench organism/gene fields separately."""
    toks = (q or "").split()
    if not toks:
        return "", ""
    org, i = [], 0
    if toks[0][:1].isupper() and toks[0].isalpha():
        org.append(toks[0]); i = 1
        if i < len(toks) and toks[i].islower() and toks[i].isalpha() and toks[i] not in _MARKER_WORDS:
            org.append(toks[i]); i += 1
    return " ".join(org), " ".join(toks[i:])


def _reference(seqs):
    """Most complete sequence at the transcript scale: longest at or below ~12 kb,
    which admits real mRNAs (immune-gene transcripts run well past 2 kb) while still
    excluding whole mitogenomes (~16 kb+) that would otherwise get sliced in the
    wrong gene. If everything is larger, fall back to the shortest."""
    capped = [s for s in seqs if len(s) <= 12000]
    return max(capped, key=len) if capped else min(seqs, key=len)


def _spread_order(starts):
    """Visit a sorted coordinate list in a target-spanning order.

    The sequence 5' -> 3' is a poor order when a runtime budget may expire: it
    recreates the exact 5'-bias the sliding-window search is intended to remove.
    Recursive endpoint/midpoint sampling touches the full target early, then fills
    the gaps.
    """
    starts = sorted(set(starts))
    out, seen = [], set()
    def add(lo, hi):
        if lo > hi:
            return
        for idx in (lo, hi, (lo + hi) // 2):
            if idx not in seen:
                seen.add(idx); out.append(starts[idx])
        mid = (lo + hi) // 2
        add(lo + 1, mid - 1)
        add(mid + 1, hi - 1)
    if starts:
        add(0, len(starts) - 1)
    return out


def _candidates_with_ledger(reference, profile, n=3, window=420, step=140, budget_s=35.0):
    """Target-wide complete-triplet search plus machine-readable attrition ledger."""
    # ``n`` is the requested number of displayed finalists, not the number allowed
    # to survive preliminary search.  A broad pool must reach full annotation.
    retained_limit = max(48, min(120, int(n) * 12))
    key = _stable_hash({"reference": reference, "profile": profile, "n": int(n),
                        "window": int(window), "step": int(step), "budget_s": float(budget_s),
                        "retained_limit": retained_limit, "conditions": T._snapshot(),
                        "search_version": getattr(CSEARCH, "SEARCH_VERSION", "unknown")})
    cached = _cache_get(_SEARCH_CACHE, key)
    if cached is not None:
        return cached
    result = CSEARCH.search(reference, profile, window=window, step=step,
                            budget_s=budget_s, pair_limit=12, probes_per_pair=3,
                            triplets_per_window=30, retained_limit=retained_limit,
                            max_windows=18)
    _cache_put(_SEARCH_CACHE, key, result, _SEARCH_CACHE_MAX)
    return deepcopy(result)


def _candidates(reference, profile, n=3, window=420, step=140, budget_s=35.0):
    """Backward-compatible candidate-list API; use _candidates_with_ledger for audit data."""
    rows, _ledger = _candidates_with_ledger(reference, profile, n, window, step, budget_s)
    return rows

def _score(assay, targets, offs, min_ident):
    oligos = {"F": assay["forward"], "R": assay["reverse"]}
    if assay.get("probe"):
        oligos["P"] = assay["probe"]
    cons = {k: C.conservation(s, targets, min_ident) for k, s in oligos.items()}
    disc = {k: C.discrimination(s, offs) for k, s in oligos.items()} if offs else None
    s = (cons["F"]["worst_3prime"] + cons["R"]["worst_3prime"]) / 2
    if "P" in cons:
        s += cons["P"]["mean_ident"]
    if disc:
        # A primer whose 3' end mismatches the off-target won't extend on it, so a single such
        # primer stops the off-target amplicon -- reward the STRONGEST 3'-end blocker, and
        # separately penalize broad off-target similarity across all oligos.
        ids, blocks = [], []
        for k in ("F", "R"):
            dk = disc.get(k)
            if dk and dk.get("n"):
                ids.append(dk.get("max_ident", 0))
                blocks.append(min(dk.get("min_3prime_mismatch", 0), 3))
        dp = disc.get("P")
        if dp and dp.get("n"):
            ids.append(dp.get("max_ident", 0))
        if ids:
            s -= sum(ids) / len(ids)
        if blocks:
            # 9*sum + 5*max equals the old 14*max for a SINGLE blocker (unchanged), and adds 9 per
            # extra 3'-block: a pair that fails the off-target at BOTH 3' ends is far more robust --
            # both terminal mismatches must independently read through to make a false product.
            s += 9.0 * sum(blocks) + 5.0 * max(blocks)
    # Cross-candidate dimer preference: among windows, demote a candidate whose worst self/hetero
    # dimer is worse than a safe floor (-5.5 kcal/mol) so a comparable dimer-clean candidate wins.
    # Candidates milder than the floor get zero penalty, so existing picks are unchanged.
    _wd = min(T.self_dimer(assay["forward"]), T.self_dimer(assay["reverse"]),
              T.hetero_dimer(assay["forward"], assay["reverse"]))
    s -= 4.0 * max(0.0, -5.5 - _wd)
    return round(s, 1), cons, disc


# IUPAC degenerate lookup (reverse of thermo._IUPAC_SETS): frozenset of bases -> code letter
_DEG_CODE = {frozenset(v): k for k, v in T._IUPAC_SETS.items()}


def _seq_quality(seq):
    """Sequence red flags that primer3's per-oligo filters still let through: long
    mononucleotide runs, dinucleotide repeats (slippage / homodimer), extreme GC.
    Returns (flag_strings, score_penalty)."""
    s = (seq or "").upper(); n = len(s)
    if not n:
        return [], 0.0
    flags, pen = [], 0.0
    gc = 100.0 * sum(c in "GC" for c in s) / n
    run = mx = 1
    for i in range(1, n):
        run = run + 1 if s[i] == s[i - 1] else 1
        mx = max(mx, run)
    if mx >= 6:
        flags.append("%d-base run" % mx); pen += 2.0
    elif mx == 5:
        flags.append("5-base run"); pen += 0.5
    best_di = 0
    for i in range(n - 3):
        a, b = s[i], s[i + 1]
        if a == b:
            continue
        j, units = i, 0
        while j + 1 < n and s[j] == a and s[j + 1] == b:
            units += 1; j += 2
        best_di = max(best_di, units)
    if best_di >= 4:
        flags.append("%dx dinucleotide repeat" % best_di); pen += 2.0
    elif best_di == 3:
        flags.append("3x dinucleotide repeat"); pen += 0.8
    if gc > 72:
        flags.append("GC %.0f%% high" % gc); pen += 1.5
    elif gc < 28:
        flags.append("GC %.0f%% low" % gc); pen += 1.0
    return flags, pen


def _best_window(oligo, seq):
    """Lowest-mismatch ungapped placement of oligo in seq; returns (window, mismatches)."""
    L = len(oligo); best, bestmm = None, L + 1
    for i in range(0, len(seq) - L + 1):
        mm = 0
        for a, b in zip(oligo, seq[i:i + L]):
            if a != b:
                mm += 1
                if mm >= bestmm:
                    break
        if mm < bestmm:
            bestmm, best = mm, seq[i:i + L]
            if mm == 0:
                break
    return best, bestmm


def _degenerate(oligo, targets, max_mm_frac=0.25, min_minor=0.18, min_count=2, cap=40):
    """Place the oligo against every target (both orientations), tally the best-matching
    window per target, and emit IUPAC codes where targets disagree -- this is what lets a
    genus design carry a W/R/Y for pan-genus coverage rather than a single consensus base.
    Returns (degenerate_seq, n_degenerate_positions, n_targets_used)."""
    o = (oligo or "").upper(); L = len(o)
    if L < 8:
        return o, 0, 0
    cols = [{} for _ in range(L)]; used = 0
    for t in targets[:cap]:
        t = t.upper()
        w1, m1 = _best_window(o, t)
        try:
            w2, m2 = _best_window(o, T.revcomp(t))
        except Exception:
            w2, m2 = None, L + 1
        w, mm = (w1, m1) if m1 <= m2 else (w2, m2)
        if w is None or mm > L * max_mm_frac:
            continue
        used += 1
        for i, ch in enumerate(w):
            if ch in "ACGT":
                cols[i][ch] = cols[i].get(ch, 0) + 1
    if used < 2:
        return o, 0, used
    out, ndeg = [], 0
    for i, ch0 in enumerate(o):
        c = cols[i]; tot = sum(c.values())
        if not tot:
            out.append(ch0); continue
        maj = max(c, key=c.get)
        keep = {b for b, k in c.items() if k >= min_count and k / tot >= min_minor}
        keep.add(maj)
        out.append(_DEG_CODE.get(frozenset(keep), maj))
        if len(keep) > 1:
            ndeg += 1
    return "".join(out), ndeg, used


def _disc_candidates(reference, profile, offs, want=6, window=420, screen=80, max_windows=12):
    """When an off-target set is given, surface primer pairs that DISCRIMINATE it.
    The normal path keeps only the single Tm-optimal pair per window, so a pair whose
    primer 3' end mismatches the off-target -- and would block its amplicon -- is thrown
    away before scoring even when it sits right there in the enumeration. Here we enumerate
    per window, keep every pair whose forward OR reverse primer has a 3'-end mismatch to the
    off-target, rank by 3'-block strength / low off-target identity / Tm-fit, attach a probe,
    and hand them to the scorer next to the Tm-optimal candidates. Best-effort and bounded;
    only runs when offs are provided, so the no-off-target path is unchanged."""
    ref = (reference or "").upper(); L = len(ref)
    offs = [o for o in (offs or []) if o][:12]
    if not offs or L < 40:
        return []
    if L <= window:
        starts = [0]
    else:
        span = L - window
        starts = sorted(set(round(i * span / max(1, max_windows - 1)) for i in range(max_windows)))
        starts = _spread_order(starts)
    dcache = {}
    def blk3(seq):
        v = dcache.get(seq)
        if v is None:
            if len(dcache) >= 500:            # bound total alignment work on long, many-window references
                return (0, 100)
            d = C.discrimination(seq, offs)
            v = ((d.get("min_3prime_mismatch", 0), d.get("max_ident", 100)) if d.get("n") else (0, 100))
            dcache[seq] = v
        return v
    ranked = []
    for s in starts:
        win = ref[s:s + window]
        try:
            fwd, rev = D.enumerate_primers(win, profile)
            pairs = D.pair_primers(fwd, rev, profile)
        except Exception:
            continue
        for p in pairs[:screen]:
            fb, fid = blk3(p["f"]); rb, rid = blk3(p["r"])
            blk = max(fb, rb)
            if blk < 1:                       # neither primer 3'-blocks the off-target -> not a discriminator
                continue
            # rank by TOTAL 3'-block first (a pair that fails the off-target at BOTH ends beats one that
            # leans on a single primer), then strongest single block, then low off-target identity, Tm-fit.
            ranked.append((-(fb + rb), -blk, round((fid + rid) / 2.0, 1), p["score"], s, win, p))
    ranked.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    out, seen = [], set()
    no_probe = profile.get("no_probe")
    # Preserve primer-pair diversity before spending additional slots on probe
    # alternatives.  A naive pair-major loop let the first two pairs consume the
    # complete specialist budget with three probes each, discarding a stronger
    # later pair before full annotation.  Select one triplet per unique pair first,
    # then add second/third probes round-robin.
    pair_groups = []
    pair_seen = set()
    pair_budget = max(8, min(int(want), (3 * int(want) + 3) // 4))
    for _sum, _b, _id, _sc, start, win, pair in ranked:
        pident = (pair["f"], pair["r"])
        if pident in pair_seen:
            continue
        if no_probe:
            probes = [None]
        else:
            probes = D.enumerate_probe_candidates(
                win, pair["fend"], pair["rstart"], pair["f"], pair["r"], profile, limit=3)
            if not probes:
                continue
        pair_seen.add(pident)
        pair_groups.append((start, win, pair, probes))
        if len(pair_groups) >= pair_budget:
            break

    max_depth = max((len(x[3]) for x in pair_groups), default=0)
    for depth in range(max_depth):
        for start, win, pair, probes in pair_groups:
            if depth >= len(probes):
                continue
            probe = probes[depth]
            ident = (pair["f"], pair["r"], probe.get("probe") if probe else None)
            if ident in seen:
                continue
            assay = D.assay_from_pair(win, pair, probe=probe)
            assay["search_window_start"] = start
            assay["f_xy"] = [assay["f_xy"][0] + start, assay["f_xy"][1] + start]
            assay["r_xy"] = [assay["r_xy"][0] + start, assay["r_xy"][1] + start]
            if assay.get("probe_xy"):
                assay["probe_xy"] = [assay["probe_xy"][0] + start, assay["probe_xy"][1] + start]
            assay["amplicon_xy"] = [assay["f_xy"][0], assay["r_xy"][1]]
            if assay.get("gblock_span"):
                assay["gblock_span"] = [assay["gblock_span"][0] + start,
                                          assay["gblock_span"][1] + start]
            seen.add(ident); out.append(assay)
            if len(out) >= want:
                return out
    return out


def epcr_offline(forward, reverse, sequences, probe=None, min_product=40, max_product=3000, min_ident=0.75):
    """Deterministic in-silico PCR of a primer pair over caller-supplied sequences -- NO BLAST, no network.
    Places each primer on each sequence (ungapped, IUPAC- and mismatch-aware via conservation.best_placement),
    converts the placement to + -strand coordinates, and reuses the shipped epcr() convergence/size logic to
    predict products. For each product, also reports whether the probe binds inside it (the difference between a
    harmless mis-priming and a false-positive signal). Returns a list of dict(subject, size, span, probe_binds).
    The 3'-end of each primer must match its binding window (q3) for the primer to be extension-competent."""
    hits = []
    for idx, seq in enumerate(sequences or []):
        su = (seq or "").upper().replace("U", "T")
        if not su:
            continue
        for primer, nm in ((forward, "F"), (reverse, "R")):
            if not primer:
                continue
            p = C.best_placement(primer, su)
            if not p or p["ident"] < min_ident:
                continue
            L = len(primer); i = p["offset"]; w = p.get("window") or ""
            q3 = bool(w) and C._match(primer[-1], w[-1])      # extension needs a matched 3' terminal base
            if p["strand"] == "+":
                lo, hi, strand = i + 1, i + L, "+"            # 1-based inclusive on the + strand
            else:                                             # matched on rc(su) -> map back to + coords
                lo, hi, strand = len(su) - (i + L) + 1, len(su) - i, "-"
            hits.append(dict(primer=nm, subject=idx, lo=lo, hi=hi, strand=strand, q3=q3,
                             ident=round(100 * p["ident"], 1)))
    products = SP.epcr(hits, min_product=min_product, max_product=max_product, require_3prime=True)
    out = []
    for pr in products:
        binds = None
        if probe:
            sub = (sequences[pr["subject"]] or "").upper().replace("U", "T")
            seg = sub[max(0, pr["span"][0] - 1): pr["span"][1]]
            pp = C.best_placement(probe, seg)
            binds = bool(pp and pp["ident"] >= 0.80)
        out.append(dict(subject=pr["subject"], size=pr["size"], span=pr["span"], probe_binds=binds))
    return out


def _junction_record(assay, junctions):
    if not junctions:
        return None
    def crosses(sp):
        return [j for j in junctions if sp and sp[0] < j < sp[1]]
    f = crosses(assay.get("f_xy")); r = crosses(assay.get("r_xy")); p = crosses(assay.get("probe_xy"))
    amp = crosses(assay.get("amplicon_xy"))
    level = "strong" if (f or r or p) else ("size" if amp else "none")
    return dict(level=level, forward=bool(f), reverse=bool(r), probe=bool(p),
                amplicon=bool(amp), junctions=sorted(set(f + r + p + amp)))



def _augment_objective_probes(candidates, reference, profile, targets, offs, objective,
                              min_ident=0.6, pair_limit=10, probe_scan_limit=24,
                              additions_per_pair=2):
    """Recover target-aware probe alternatives that a cheap local beam can miss.

    Probe enumeration initially has only the reference sequence and thermodynamic
    evidence.  In multi-isolate or discrimination designs, a probe with slightly
    worse local Tm preference can be decisively better after conservation or
    off-target annotation.  This bounded second-stage scan revisits a diverse set
    of retained primer pairs, ranks additional probes using cheap corpus evidence,
    and sends only a few alternatives per pair to the authoritative full ranker.
    """
    if profile.get("no_probe") or not candidates or (len(targets or []) < 2 and not offs):
        return [], dict(stage="objective_probe_augmentation", entered=0, retained=0,
                        rejected=0, hard_gate=False, reversible=True,
                        reason="not applicable")
    # One representative per primer pair, ordered by the existing preliminary
    # pair/triplet evidence.  Pair diversity was already preserved upstream.
    pair_rows = {}
    for assay in candidates:
        key = (assay.get("forward"), assay.get("reverse"),
               tuple(assay.get("f_xy") or ()), tuple(assay.get("r_xy") or ()))
        if not key[0] or not key[1] or len(key[2]) != 2 or len(key[3]) != 2:
            continue
        old = pair_rows.get(key)
        if old is None or float(assay.get("candidate_rank", 1e12)) < float(old.get("candidate_rank", 1e12)):
            pair_rows[key] = assay
    ordered_pairs = sorted(pair_rows.values(), key=lambda a:(float(a.get("candidate_rank", 1e12)),
                                                              (a.get("f_xy") or [0])[0],
                                                              a.get("forward"), a.get("reverse")))[:max(1,int(pair_limit))]
    existing = {(a.get("forward"), a.get("reverse"), a.get("probe")) for a in candidates}
    added=[]; scanned=0; accepted=0; decisions=[]
    for assay in ordered_pairs:
        fxy, rxy = assay["f_xy"], assay["r_xy"]
        pair = dict(f=assay["forward"], r=assay["reverse"],
                    fstart=int(fxy[0]), fend=int(fxy[1]),
                    rstart=int(rxy[0]), rend=int(rxy[1]),
                    amp=int(assay.get("amplicon") or (int(rxy[1])-int(fxy[0]))),
                    gap=float(assay.get("pair_tm_gap") or abs(T.tm(assay["forward"])-T.tm(assay["reverse"]))),
                    dimer=float(assay.get("pair_dimer") or T.hetero_dimer(assay["forward"],assay["reverse"])),
                    score=float(assay.get("preliminary_pair_score") or 0.0))
        try:
            probes = D.enumerate_probe_candidates(reference, pair["fend"], pair["rstart"],
                                                   pair["f"], pair["r"], profile,
                                                   limit=max(4,int(probe_scan_limit)))
        except Exception as exc:
            decisions.append(dict(pair=[pair["f"],pair["r"]],decision="error",reason=type(exc).__name__))
            continue
        scored=[]
        for probe in probes:
            scanned += 1
            ident=(pair["f"],pair["r"],probe.get("probe"))
            if ident in existing:
                decisions.append(dict(pair=[pair["f"],pair["r"]],probe=probe.get("probe"),
                                      decision="rejected",reason="existing_triplet_duplicate"))
                continue
            tc=C.conservation(probe["probe"],targets,min_ident)
            dc=C.discrimination(probe["probe"],offs) if offs else None
            target_mean=float(tc.get("mean_ident") or 0.0)
            target_min=float(tc.get("min_ident") or 0.0)
            off_max=float((dc or {}).get("max_ident") or 0.0)
            if objective in {"discrimination","confirmatory"}:
                corpus_key=(-target_min,-target_mean,off_max,float(probe.get("preliminary_penalty",1e9)))
            else:
                corpus_key=(-target_min,-target_mean,float(probe.get("preliminary_penalty",1e9)),off_max)
            scored.append((corpus_key,probe,tc,dc))
        ordered_scored=sorted(scored,key=lambda x:(x[0],x[1].get("start",0),x[1].get("probe","")))
        keep_n=max(1,int(additions_per_pair))
        for idx,(corpus_key,probe,tc,dc) in enumerate(ordered_scored):
            if idx >= keep_n:
                decisions.append(dict(pair=[pair["f"],pair["r"]],probe=probe.get("probe"),
                                      decision="rejected",reason="objective_probe_budget",
                                      target_mean_identity=tc.get("mean_ident"),target_min_identity=tc.get("min_ident"),
                                      off_target_max_identity=(dc or {}).get("max_ident"),
                                      selection_key=list(corpus_key)))
                continue
            new=D.assay_from_pair(reference,pair,probe=probe)
            new["search_window_start"]=0
            new["amplicon_xy"]=[new["f_xy"][0],new["r_xy"][1]]
            new["objective_probe_augmentation"]={"target_conservation":tc,"offtarget_discrimination":dc,
                                                  "selection_key":list(corpus_key),
                                                  "status":"cheap_corpus_screen_only"}
            ident=(new.get("forward"),new.get("reverse"),new.get("probe"))
            if ident not in existing:
                existing.add(ident);added.append(new);accepted += 1
                decisions.append(dict(pair=[pair["f"],pair["r"]],probe=probe.get("probe"),
                                      decision="retained",reason="objective_corpus_probe_alternative",
                                      target_mean_identity=tc.get("mean_ident"),target_min_identity=tc.get("min_ident"),
                                      off_target_max_identity=(dc or {}).get("max_ident"),
                                      selection_key=list(corpus_key)))
    reason_counts={}
    for decision in decisions:
        if decision.get("decision") in {"retained","rejected"}:
            reason=decision.get("reason") or "unspecified"
            reason_counts[reason]=reason_counts.get(reason,0)+1
    ledger=dict(stage="objective_probe_augmentation",unit="probe_candidates",
                entered=scanned,retained=accepted,rejected=max(0,scanned-accepted),
                hard_gate=False,reversible=True,pairs_scanned=len(ordered_pairs),
                pair_limit=int(pair_limit),probe_scan_limit=int(probe_scan_limit),
                additions_per_pair=int(additions_per_pair),candidate_decisions=decisions,
                rejection_reasons={k:v for k,v in reason_counts.items() if k != "objective_corpus_probe_alternative"},
                reason="target/off-target corpus evidence can rescue probes outside the local thermodynamic beam")
    return added,ledger

def design_from_sequences(targets, profile, offs=None, min_ident=0.6, n_candidates=5,
                          objective="balanced", junctions=None, panel=None,
                          search_budget_s=35.0):
    """Design and rank complete assays under an explicit objective profile.

    Every retained triplet receives the same target coverage, supplied off-target,
    conservation, and condition-robustness evaluations before final sorting.
    """
    _clean = []
    for t in (targets or []):
        c, _n, err = T.clean_seq(t) if isinstance(t, str) and t.strip() else ("", None, "empty")
        if c and not err:
            _clean.append(c)
    targets = [t for t in _clean if len(t) > 60]
    if not targets:
        return dict(error="no usable target sequences (need >=1, ideally several)")
    _off_clean = []
    for t in (offs or []):
        c, _n, err = T.clean_seq(t) if isinstance(t, str) and t.strip() else ("", None, "empty")
        if c and not err and len(c) > 20:
            _off_clean.append(c)
    offs = _off_clean
    search_budget_s = max(3.0, min(float(search_budget_s), 120.0))
    design_key = _stable_hash({"targets": targets, "offs": offs, "profile": profile,
                               "min_ident": float(min_ident), "n_candidates": int(n_candidates),
                               "objective": objective, "junctions": junctions or [],
                               "panel": panel or [], "conditions": T._snapshot(),
                               "search_budget_s": float(search_budget_s),
                               "ranker_version": getattr(RANK, "RANKING_SCHEMA_VERSION", "unknown")})
    cached_design = _cache_get(_DESIGN_CACHE, design_key)
    if cached_design is not None:
        return cached_design
    ref = _reference(targets)
    cands, attrition = _candidates_with_ledger(ref, profile, n_candidates,
                                               budget_s=search_budget_s)
    attrition.setdefault("candidate_limits", {})["search_budget_seconds"] = search_budget_s
    if offs:
        # Preserve discrimination-specialist pairs that a generic preliminary
        # thermodynamic beam may not surface.  They still undergo the same final
        # full annotation and cannot bypass hard constraints.
        try:
            _dc = _disc_candidates(ref, profile, offs, want=max(16, int(n_candidates) * 3))
            _have = {(c.get("forward"), c.get("reverse"), c.get("probe")) for c in cands}
            _before = len(_have)
            for _a in _dc:
                _key = (_a.get("forward"), _a.get("reverse"), _a.get("probe"))
                if _key not in _have:
                    cands.append(_a); _have.add(_key)
            attrition["discrimination_specialists_added"] = len(_have) - _before
        except Exception as exc:
            attrition["discrimination_specialist_error"] = type(exc).__name__
    augmented_probes, augmentation_ledger = _augment_objective_probes(
        cands, ref, profile, targets, offs, objective, min_ident=min_ident,
        pair_limit=max(8, int(n_candidates) * 2), probe_scan_limit=24,
        additions_per_pair=2)
    if augmented_probes:
        cands.extend(augmented_probes)
    attrition.setdefault("stages", []).append(augmentation_ledger)
    attrition["objective_probe_alternatives_added"] = len(augmented_probes)
    if not cands:
        return dict(error="no primer/probe set met this chemistry's Tm window on the reference. "
                          "AT-rich targets like parasite mtDNA can't reach a high-Tm probe — use the "
                          "Auto setting, or pick a low-Tm / MGB profile. (Very short fetched sequences "
                          "can also lack a full amplicon window.)", candidate_attrition=attrition)
    multi = len(targets) >= 2
    prelim = []
    for a in cands:
        sc, cons, disc = _score(a, targets, offs, min_ident)
        qf, pen = _seq_quality(a["forward"])
        rqf, rpen = _seq_quality(a["reverse"]); qf = list(qf) + list(rqf); pen += rpen
        if a.get("probe"):
            pqf, ppen = _seq_quality(a["probe"]); qf += pqf; pen += ppen
        a["quality_flags"] = qf
        _amp = a.get("amplicon")
        amp_pen = 0.25 * (_amp - 150) if isinstance(_amp, int) and _amp > 150 else 0.0
        if multi:
            df, ndf, nu = _degenerate(a["forward"], targets)
            dr, ndr, _ = _degenerate(a["reverse"], targets)
            dp, ndp = (None, 0)
            if a.get("probe"):
                dp, ndp, _ = _degenerate(a["probe"], targets)
            if nu >= 2 and (ndf + ndr + ndp) > 0:
                a["forward_deg"], a["reverse_deg"] = df, dr
                if dp is not None:
                    a["probe_deg"] = dp
                a["n_degenerate"] = ndf + ndr + ndp
                a["deg_targets"] = nu
        prelim.append(dict(score=round(sc - pen - amp_pen, 1), score_raw=sc,
                           quality_penalty=round(pen, 1), amplicon_penalty=round(amp_pen, 1),
                           assay=a, conservation=cons, discrimination=disc))

    # Primer3's native structure engine can become unstable after very large
    # uninterrupted batches, and fully annotating every cheap-screen survivor is
    # not necessary for an interactive bounded search.  Retain a broad,
    # objective-aware and regionally diverse pool before expensive all-site PCR
    # and condition-envelope structure analysis.  This stage is fully auditable;
    # discarded candidates remain visible in the attrition ledger.
    full_annotation_limit = max(20, min(28, int(n_candidates) * 5))
    for row in prelim:
        # retain_diverse is a lower-is-better beam; the legacy preliminary score
        # is higher-is-better and already includes objective-relevant cheap
        # conservation/discrimination evidence.
        row["assay"]["candidate_rank"] = -float(row.get("score", 0.0))
    retained_assays, annotation_ledger = CRET.retain_diverse(
        [row["assay"] for row in prelim], limit=full_annotation_limit,
        region_size=max(100, int(profile.get("amp_max", 150))),
        per_region=max(4, full_annotation_limit // 4), per_near=3)
    rows_by_identity = {}
    for row in sorted(prelim, key=lambda x: (-float(x.get("score", 0.0)),
                                             CRET.identity_key(x["assay"]))):
        rows_by_identity.setdefault(CRET.identity_key(row["assay"]), row)
    prelim = [rows_by_identity[CRET.identity_key(a)] for a in retained_assays]
    annotation_ledger = dict(annotation_ledger)
    annotation_ledger.update(stage="full_annotation_diversity_retention",
                             unit="complete_assays",
                             reason="objective-aware regional beam before all-site PCR and condition robustness")
    attrition.setdefault("candidate_limits", {})["full_annotation_pool"] = full_annotation_limit
    attrition.setdefault("stages", []).append(annotation_ledger)

    jmap = {}
    if junctions:
        for row in prelim:
            a = row["assay"]
            jmap[(a.get("forward"), a.get("reverse"), a.get("probe"))] = _junction_record(a, junctions)
    ranked, objective_profile = RANK.rank_candidates(prelim, targets, offs, profile,
                                                      objective_name=objective,
                                                      junction_by_identity=jmap,
                                                      panel=panel)
    finalists = RANK.select_finalists(ranked, n=max(1, int(n_candidates)))
    for idx, row in enumerate(finalists):
        comp = ranked[row["rank"]] if row.get("rank", 0) < len(ranked) else None
        row["rank_explanation"] = REXPLAIN.explain(row, comp)
        # Compatibility: ``score`` remains visible, but is explicitly subordinate.
        row["score"] = row["display_score"]
    attrition.setdefault("stages", []).append({
        "stage": "full_annotation", "entered": len(prelim), "retained": len(ranked),
        "rejected": 0, "hard_gate": False, "reversible": False,
        "evaluations": ["target_epcr", "offtarget_epcr", "conservation", "condition_robustness"]})
    attrition["stages"].append({
        "stage": "finalist_selection", "entered": len(ranked), "retained": len(finalists),
        "rejected": len(ranked) - len(finalists), "hard_gate": False, "reversible": True,
        "reason": "diverse finalist display budget"})
    ih = hashlib.sha256(("\n".join(targets) + "\n--OFF--\n" + "\n".join(offs)).encode()).hexdigest()
    out = dict(n_targets=len(targets), n_offs=len(offs), reference_len=len(ref),
               n_candidates=len(finalists), n_candidates_screened=len(ranked),
               n_requested=n_candidates, candidates=finalists,
               objective_profile=objective_profile, candidate_attrition=attrition,
               ranker_manifest=RANK.manifest(objective_profile, attrition.get("candidate_limits", {}),
                                             input_hashes={"sequence_corpus_sha256": ih}),
               search_status="heuristic_bounded",
               cache_policy=("bounded defensive-copy memoization keyed by sequence corpus, chemistry, "
                             "reaction conditions, objective, and constraints"),
               ranking_statement=("Strongest computational support among the fully evaluated retained pool "
                                  "under the declared objective; not a universal wet-lab optimum."))
    if len(finalists) < n_candidates:
        out["constraint_note"] = ("only %d finalist set(s) survived the %s search and ranking constraints. "
            "Try a different objective/profile or a longer target region; do not relax true off-target or "
            "geometry failures merely to fill the list." % (len(finalists), profile.get("name", "selected")))
    _cache_put(_DESIGN_CACHE, design_key, out, _DESIGN_CACHE_MAX)
    return deepcopy(out)


AUTO_ORDER = ["idt_taqman", "parasite_mtdna", "gc_rich", "parasite_sybr"]   # normal-composition default
_PRETTY = {"idt_taqman": "IDT PrimeTime (ZEN double-quenched probe)",
           "idt_affinity": "IDT Affinity Plus (LNA probe)",
           "thermo_taqman": "Thermo TaqMan (MGB)",
           "parasite_mtdna": "low-Tm TaqMan \u2014 order the probe as IDT Affinity Plus (LNA)",
           "parasite_sybr": "low-Tm SYBR (primers only, no probe)",
           "gc_rich": "GC-rich high-Tm TaqMan (short primers; consider DMSO / 7-deaza-dGTP)",
           "gc_rich_sybr": "GC-rich high-Tm SYBR (no probe)",
           "biorad_probe": "Bio-Rad PrimePCR", "sybr_generic": "SYBR generic"}


def _auto_order(ref):
    """Composition-aware chemistry order for Auto. AT-rich and GC-rich targets each get their
    tuned low-/high-Tm profile tried FIRST (the same first-class treatment), instead of spending
    the first attempt on a generic profile that cannot meet the Tm/GC window \u2014 and, on a GC-rich
    target, failing every profile. 40-62% GC keeps the standard IDT-first order.
    Returns (ordered_profile_keys, gc_percent)."""
    s = (ref or "").upper()
    gc = 100.0 * sum(c in "GC" for c in s) / max(1, len(s))
    if gc < 40.0:
        return ["parasite_mtdna", "idt_taqman", "gc_rich", "parasite_sybr"], gc
    if gc > 62.0:
        return ["gc_rich", "idt_taqman", "parasite_mtdna", "sybr_generic"], gc
    return list(AUTO_ORDER), gc


def _amplicon_on(seq, fwd, amplen):
    """Amplicon start = where the forward primer sits; length is the DESIGNED amplicon
    length. Deliberately does NOT search for the reverse primer: on a long, AT-rich
    whole-mitogenome reference the reverse primer can match far away and yield a multi-kb
    span, and folding that is O(n^3) (minutes). Returns (start, end) or None."""
    if not seq or not amplen:
        return None
    fs = SP._locate(fwd, seq)
    if fs is None:
        return None
    return (fs, min(len(seq), fs + int(amplen)))


def _annotate(out, ref, prefer_junction):
    """Add template-structure (always, offline) and, when prefer_junction is set,
    exon-junction-spanning to each candidate. Structure is folded on the actual designed
    amplicon (from the reference). For the junction call the exon table and the amplicon
    must share a coordinate frame, so junctions are read from whatever transcript the
    gene_table returns and the amplicon is re-located on THAT transcript -- not assumed to
    match the reference (isoform variants differ in exon structure)."""
    cands = out.get("candidates") or []
    _prof = PROF.PROFILES.get(out.get("profile_used")) or {}   # anneal temp for the fold; falls back to the global default
    junctions = mrna = None
    if prefer_junction and out.get("source_accession"):
        try:
            junctions, jinfo, used = SP.exon_junctions_mrna("", "", out["source_accession"])
        except Exception:
            junctions, jinfo, used = None, "couldn't read exon annotation for the source accession", None
        out["junction_info"] = jinfo
        if junctions and used:
            try:
                recs = N.fetch_accessions([used], "fasta")
                mrna = str(recs[0].seq) if recs else None
            except Exception:
                mrna = None
    for c in cands:
        a = c["assay"]
        rspan = _amplicon_on(ref, a["forward"], a.get("amplicon")) if ref else None
        c["amp_span"] = list(rspan) if rspan else None
        if STR.available() and ref and rspan:
            amp = ref[rspan[0]:rspan[1]]
            f = STR.fold_ensemble(amp, anneal_c=_prof.get("anneal_c", T.ANNEAL_C))
            if f:
                fl, rl = len(a["forward"]), len(a["reverse"])
                pp = SP._locate(a["probe"], amp) if a.get("probe") else None
                if a.get("probe") and pp is None:
                    pp = SP._locate(SP._rc_iupac(a["probe"]), amp)
                _pr, _pa = f.get("paired_prob"), f.get("paired_anneal")
                _has_p = bool(a.get("probe")) and pp is not None
                c["structure"] = dict(
                    mfe=f["mfe"], mfe_per_nt=f["mfe_per_nt"], dna=f["dna_params"],
                    f_paired=STR.site_paired_fraction(f["paired"], 0, fl),
                    r_paired=STR.site_paired_fraction(f["paired"], len(amp) - rl, len(amp)),
                    p_paired=(STR.site_paired_fraction(f["paired"], pp, pp + len(a["probe"]))
                              if _has_p else None),
                    # ensemble (partition-function) paired probability -- more honest than MFE alone
                    f_paired_prob=STR.site_paired_prob(_pr, 0, fl),
                    r_paired_prob=STR.site_paired_prob(_pr, len(amp) - rl, len(amp)),
                    p_paired_prob=(STR.site_paired_prob(_pr, pp, pp + len(a["probe"])) if _has_p else None),
                    # structure that survives at the annealing temperature (what priming actually sees)
                    anneal_c=f.get("anneal_c"), mfe_anneal=f.get("mfe_anneal"),
                    f_paired_anneal=STR.site_paired_fraction(_pa, 0, fl),
                    r_paired_anneal=STR.site_paired_fraction(_pa, len(amp) - rl, len(amp)),
                    p_paired_anneal=(STR.site_paired_fraction(_pa, pp, pp + len(a["probe"])) if _has_p else None))
        if prefer_junction:
            jspan = _amplicon_on(mrna, a["forward"], a.get("amplicon")) if (mrna and junctions) else None
            c["spans_junction"] = (any(jspan[0] < j < jspan[1] for j in junctions)
                                   if (junctions and jspan) else None)
    if prefer_junction and junctions:
        cands.sort(key=lambda c: (0 if c.get("spans_junction") else 1, -c["score"]))
    return out


def design_nested(reference, profile, inner_assay, outer_flank_max=600, min_gap=8,
                  outer_amp_max=1400, topk=40):
    """Fully-nested OUTER pair that flanks a given inner (diagnostic) assay on the same
    reference: outer forward upstream of the inner forward, outer reverse downstream of the
    inner reverse, so a second-round reaction only amplifies a correct first-round product
    (the basis of the haemosporidian MalAvi screen). Inner-first, so the diagnostic assay
    keeps its full scoring; outer primer quality is gated by the SAME profile. Returns
    {"outer": {...}} or None when there isn't enough flanking sequence."""
    if not reference or not inner_assay:
        return None
    span = _amplicon_on(reference, inner_assay["forward"], inner_assay.get("amplicon"))
    if not span:
        return None
    I_fs, I_re = span
    L = len(reference)
    up = reference[max(0, I_fs - outer_flank_max):I_fs]
    dn = reference[I_re:min(L, I_re + outer_flank_max)]
    U0 = max(0, I_fs - outer_flank_max)
    if len(up) < profile["len_min"] or len(dn) < profile["len_min"]:
        return None
    fwd, _ = D.enumerate_primers(up, profile)
    _, rev = D.enumerate_primers(dn, profile)
    topt = profile["tm_opt"]
    F = sorted(((U0 + i, U0 + j, w, T.tm(w)) for (i, j, w) in fwd if (U0 + j) <= I_fs - min_gap),
               key=lambda x: abs(x[3] - topt))[:topk]
    R = sorted(((I_re + a, I_re + b, w, T.tm(w)) for (a, b, w) in rev if (I_re + a) >= I_re + min_gap),
               key=lambda x: abs(x[3] - topt))[:topk]
    inner_amp = I_re - I_fs
    best = None
    for (fs, fe, f, ftm) in F:
        for (rs, re, r, rtm) in R:
            amp = re - fs
            if amp < inner_amp + 2 * min_gap or amp > outer_amp_max:
                continue
            gap = abs(ftm - rtm)
            if gap > profile["pair_tm_gap_max"]:
                continue
            if T.hetero_dimer(f, r) <= profile["pair_dimer_min"]:
                continue
            score = abs((ftm + rtm) / 2 - topt) + gap
            if best is None or score < best["score"]:
                best = dict(score=round(score, 2), forward=f, reverse=r, amplicon=amp,
                            f_tm=round(ftm, 1), r_tm=round(rtm, 1), pair_tm_gap=round(gap, 1),
                            f_outside=I_fs - fe, r_outside=rs - I_re)
    return dict(outer=best) if best else None


def resolve_and_fetch_query(target_query, off_query=None, n_fetch=20):
    """Resolve a query and fetch the target/off-target sequence corpora.

    The returned mapping is an *internal stage value*: it intentionally contains
    the fetched sequences needed by the later design stage and must not be used as
    a public job-status payload.  Splitting this stage out lets the asynchronous
    runner record retrieval separately and, importantly, avoids fetching again
    when only optional specificity analysis is retried.
    """
    target_query = (target_query or "").strip()
    if not target_query:
        return dict(error='enter a target: an organism plus a gene or marker (e.g. "Plasmodium cytochrome b"), or paste sequences in the Design tab')
    _org, _gene = _split_query(target_query)
    _resolved, _fetch_q = None, target_query
    if _gene and _org and not _is_marker(_gene) and not RM.gold_standard_query(_org, _gene):
        try:
            _resolved = N.resolve_gene(_gene, _org)
        except Exception:
            _resolved = None
        if _resolved and _resolved.get("found"):
            if _resolved.get("in_requested_organism") and _resolved.get("clean"):
                _fetch_q = f'{_resolved["symbol"]}[Gene Name] AND {_org}[Organism]'
            elif _resolved.get("symbol") and not _resolved.get("in_requested_organism"):
                return dict(error=f"{_resolved['symbol']} ({_resolved.get('description','')}) isn't "
                                  f"annotated in {_org} — NCBI has it in "
                                  f"{_resolved.get('organism','another species')}. Design there, or "
                                  f"paste a {_gene} sequence in the Design tab. (Designing on the loosely "
                                  f"matching records in {_org} would risk a paralog like a receptor.)")
    _pairs = N.search_fetch_fasta(_fetch_q, n_fetch)
    if not _pairs and _fetch_q != target_query:
        _pairs = N.search_fetch_fasta(target_query, n_fetch)      # precise found nothing -> free text
    tg = [seq for _, seq in _pairs]
    if not tg:
        msg = f"NCBI returned nothing for: {target_query}"
        if _resolved and _resolved.get("found") and not _resolved.get("in_requested_organism"):
            msg += (f". {_resolved['symbol']} ({_resolved.get('description','')}) is annotated in "
                    f"{_resolved.get('organism','another species')}, not in {_org}.")
        return dict(error=msg)
    _off_pairs = N.search_fetch_fasta(off_query, max(8, n_fetch // 2)) if off_query else []
    off = [seq for _, seq in _off_pairs] if off_query else None

    return dict(target_query=target_query, off_query=off_query,
                organism_name=_org, gene_name=_gene, resolved=_resolved,
                fetch_query=_fetch_q, target_pairs=_pairs, off_pairs=_off_pairs,
                targets=tg, offs=off)


def design_query_corpus(query_context, profile_key="auto", min_ident=0.6,
                        objective="balanced"):
    """Run profile selection and fully rank a previously fetched query corpus."""
    if query_context.get("error"):
        return dict(error=query_context["error"])
    target_query = query_context["target_query"]
    off_query = query_context.get("off_query")
    tg = query_context.get("targets") or []
    off = query_context.get("offs")

    if profile_key == "auto":
        out, tried = None, []
        _refseq = _reference([t for t in tg if t and len(t) > 60] or tg)
        order, _autogc = _auto_order(_refseq)
        for pk in order:
            tried.append(pk)
            r = design_from_sequences(tg, PROF.PROFILES[pk], off, min_ident, objective=objective)
            if r.get("candidates"):
                out = r
                out["profile_used"] = pk
                break
        if out is None:
            return dict(error="Auto could not place a clean assay in any chemistry for this target "
                              "(tried %s at %.0f%% GC). The region may be too variable, or the fetched "
                              "sequences too short for a full amplicon; try a different gene/region, a "
                              "wider amplicon, or pick a profile manually." % (", ".join(tried), _autogc),
                              tried=tried, target_query=target_query, off_query=off_query)
        out["auto_gc"] = round(_autogc, 1)
    else:
        prof = PROF.PROFILES.get(profile_key)
        if not prof:
            return dict(error=f"unknown profile: {profile_key}")
        out = design_from_sequences(tg, prof, off, min_ident, objective=objective)
        if out.get("error"):
            out["target_query"] = target_query
            out["off_query"] = off_query
            return out
        out["profile_used"] = profile_key

    return out


def enrich_query_design(out, query_context, min_ident=0.6,
                        prefer_junction=False, nested=False,
                        objective="balanced"):
    """Attach query provenance, structure/junction data, controls, and nesting.

    This is the final required stage of automatic design.  Remote/local BLAST is
    deliberately excluded so a valid primary assay remains independently usable
    when that optional external stage is unavailable.
    """
    if out.get("error"):
        return out
    target_query = query_context["target_query"]
    off_query = query_context.get("off_query")
    _org = query_context.get("organism_name") or ""
    _gene = query_context.get("gene_name") or ""
    _resolved = query_context.get("resolved")
    _pairs = query_context.get("target_pairs") or []
    _off_pairs = query_context.get("off_pairs") or []
    tg = query_context.get("targets") or []
    off = query_context.get("offs")

    out["profile_pretty"] = _PRETTY.get(out.get("profile_used"), out.get("profile_used"))
    out["target_query"] = target_query
    out["resolved_gene"] = (_resolved.get("symbol") if (_resolved and _resolved.get("found")) else _gene) or ""
    out["resolved_organism"] = _org or (_resolved.get("organism") if _resolved else "") or ""
    out["off_query"] = off_query
    def _subjects(pairs):
        rows = []
        for desc, _sq in (pairs or [])[:40]:
            parts = (desc or "").split(None, 1)
            rows.append({"acc": parts[0] if parts else "", "title": parts[1] if len(parts) > 1 else "", "len": len(_sq)})
        return rows
    out["target_subjects"] = _subjects(_pairs)   # actual accessions/species behind the conservation call
    out["off_subjects"] = _subjects(_off_pairs)   # ...and the off-target set behind the discrimination call
    # Rebuild the deterministic manifest with the concrete accession set.  NCBI
    # retrieval dates are not invented; accession identifiers and query strings
    # are still enough to expose exactly what external evidence was requested.
    if out.get("ranker_manifest"):
        _old_manifest = out.get("ranker_manifest") or {}
        out["ranker_manifest"] = RANK.manifest(
            out.get("objective_profile") or RANK.get_profile(objective),
            (out.get("candidate_attrition") or {}).get("candidate_limits", {}),
            input_hashes=_old_manifest.get("input_hashes") or {},
            external_databases={
                "ncbi_nucleotide": {
                    "target_query": target_query,
                    "target_accessions": [x.get("acc") for x in out["target_subjects"] if x.get("acc")],
                    "offtarget_query": off_query,
                    "offtarget_accessions": [x.get("acc") for x in out["off_subjects"] if x.get("acc")],
                    "retrieval_date": None,
                    "database_snapshot": "live service; exact snapshot not reported by NCBI response",
                }
            },
            warnings=(["no supplied off-target corpus; exclusivity remains unresolved"] if not off else []),
            fallbacks=(out.get("fallbacks") or []),
            constraints={"workflow": "automatic_design", "profile": out.get("profile_used"),
                         "min_identity": min_ident, "prefer_junction": bool(prefer_junction),
                         "nested_requested": bool(nested)},
        )

    if out.get("candidates"):
        try:
            _ref = _reference(tg)
        except Exception:
            _ref = None
        _desc = next((a for a, sq in _pairs if sq == _ref), None) if (_ref and _pairs) else None
        out["source_accession"] = _desc.split()[0] if _desc else None
        _annotate(out, _ref, prefer_junction)
        if out.get("profile_used") in ("idt_affinity", "parasite_mtdna"):
            for c in out["candidates"]:                       # LNA chemistry -> suggest + positions
                pb = (c.get("assay") or {}).get("probe")
                if pb:
                    try:
                        c["assay"]["probe_lna"] = T.suggest_lna(pb)
                    except Exception:
                        pass
        # discrimination-control gBlocks for the lead candidate: positive control (target
        # amplicon + flanks) and the worst-case off-target, so the block can be tested on the bench.
        try:
            _a0 = out["candidates"][0].get("assay") or {}
            _fi = _ref.find(_a0.get("forward", "")) if _ref else -1
            _rcr = T.revcomp(_a0["reverse"]) if _a0.get("reverse") else ""
            _ri = _ref.find(_rcr) if (_ref and _rcr) else -1
            if _ref and _fi >= 0 and _ri >= 0:
                _amp = _ref[_fi:_ri + len(_a0["reverse"])]
                if not _a0.get("gblock"):
                    _gb, _gs, _ge = D.build_gblock(_ref, _fi, _ri + len(_a0["reverse"]))
                    _a0["gblock"] = _gb
                if off:
                    _og = D.build_offtarget_gblock(_amp, off)
                    if _og:
                        _a0["offtarget_gblock"] = _og["seq"]
                        _a0["offtarget_gblock_identity"] = _og["amplicon_identity"]
        except Exception:
            pass
        if nested:
            _prof = PROF.PROFILES.get(out.get("profile_used")) or {}
            try:
                nz = design_nested(_ref, _prof, out["candidates"][0]["assay"]) if (_ref and _prof) else None
            except Exception:
                nz = None
            if nz and nz.get("outer"):
                _inner = out["candidates"][0]["assay"]
                out["nested"] = dict(outer=nz["outer"], inner_amplicon=_inner["amplicon"],
                                     inner_forward=_inner["forward"], inner_reverse=_inner["reverse"])
            else:
                out["nested"] = dict(note="couldn't place a flanking outer pair — the reference needs "
                                          "roughly 150+ bp of usable sequence beyond the inner amplicon "
                                          "on each side. The inner assay above runs fine on its own.")

    return out


def blast_winner(out, blast_mode="remote", blast_db="nt", blast_db_path=None,
                 organism=None, suppress_errors=False):
    """Run optional in-silico PCR for the lead assay without redoing design.

    ``suppress_errors`` exists solely for the legacy synchronous API contract,
    which historically encoded BLAST failure in ``specificity.error``.  The job
    runner leaves it false so it can mark the optional stage failed/timed out,
    retain the primary result, and expose an explicit warning.
    """
    if out.get("candidates"):
        a = out["candidates"][0]["assay"]
        try:
            out["specificity"] = SP.in_silico_pcr(a["forward"], a["reverse"], mode=blast_mode,
                                                  db=blast_db, db_path=blast_db_path, organism=organism)
        except Exception as e:
            if not suppress_errors:
                raise
            out["specificity"] = dict(error=f"in-silico PCR could not run: {e}")
    return out


def design_from_query(target_query, profile_key="auto", off_query=None, n_fetch=20,
                      min_ident=0.6, run_blast=False, blast_mode="remote",
                      blast_db="nt", blast_db_path=None, organism=None, prefer_junction=False,
                      nested=False, objective="balanced"):

    """Fetch -> design -> enrich -> (optional) in-silico-PCR the winning pair.

    profile_key may be a specific profile, or "auto": Auto tries the IDT-orderable
    chemistries in order and returns the first that yields a clean assay, so an
    AT-rich parasite target lands on the low-Tm TaqMan instead of just failing.

    The synchronous compatibility entry point now composes the same reusable
    stages as the job runner; its result and error behavior remain unchanged.
    """
    context = resolve_and_fetch_query(target_query, off_query, n_fetch)
    if context.get("error"):
        return context
    out = design_query_corpus(context, profile_key, min_ident, objective)
    if out.get("error"):
        return out
    out = enrich_query_design(out, context, min_ident=min_ident,
                              prefer_junction=prefer_junction, nested=nested,
                              objective=objective)
    if run_blast:
        out = blast_winner(out, blast_mode=blast_mode, blast_db=blast_db,
                           blast_db_path=blast_db_path, organism=organism,
                           suppress_errors=True)
    return out
